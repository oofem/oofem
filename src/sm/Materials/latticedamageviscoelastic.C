/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "latticedamageviscoelastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "Elements/latticestructuralelement.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "math.h"


namespace oofem {
REGISTER_Material(LatticeDamageViscoelastic);

LatticeDamageViscoelastic :: LatticeDamageViscoelastic(int n, Domain *d) : LatticeDamage(n, d)
{
    slaveMat = 0;
}


LatticeDamageViscoelastic :: ~LatticeDamageViscoelastic()
//
// destructor
//
{}

int
LatticeDamageViscoelastic :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( ( mode == _2dLattice ) || ( mode == _3dLattice ) ) {
        return 1;
    }

    return 0;
}


IRResultType
LatticeDamageViscoelastic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro


    LatticeDamage :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, slaveMat, _IFT_LatticeDamageViscoelastic_slaveMat); // number of slave material

    RheoChainMaterial *rChM = giveViscoelasticMaterial();

    // override elastic modulus by the one given by the compliance function
    // provided values of the strain at peak stress must be changed

    double E28 = 1. / rChM->computeCreepFunction(28.01, 28.,NULL,NULL); // modulus of elasticity evaluated at 28 days, duration of loading 15 min
    this->e0Mean *= ( this->eNormalMean / E28 ); // transform strain at peak stress accordingly
    this->eNormalMean = E28; // swap elastic modulus/stiffness


    // must match with the viscoelastic material
    if ( rChM->giveAlphaOne() != this->alphaOne ) {
        OOFEM_ERROR("a1 must be set to the same value in both master and viscoelastic slave materials");
    }

    // must match with the viscoelastic material
    if ( rChM->giveAlphaTwo() != this->alphaTwo ) {
        OOFEM_ERROR("a2 must be set to the same value in both master and viscoelastic slave materials");
    }

    return IRRT_OK;
}


MaterialStatus *
LatticeDamageViscoelastic :: CreateStatus(GaussPoint *gp) const
{
    LatticeDamageViscoelasticStatus *answer = new LatticeDamageViscoelasticStatus(1, LatticeDamageViscoelastic :: domain, gp, slaveMat);
    return answer;
}



void
LatticeDamageViscoelastic :: giveRealStressVector(FloatArray &answer,
                                                  GaussPoint *gp,
                                                  const FloatArray &totalStrain,
                                                  TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    LatticeDamageViscoelasticStatus *status = static_cast< LatticeDamageViscoelasticStatus * >( this->giveStatus(gp) );

    FloatArray reducedStrain;
    FloatArray quasiElasticStrain;


    RheoChainMaterial *rChM = giveViscoelasticMaterial();

    MaterialMode mMode = gp->giveMaterialMode();
    StressVector tempEffectiveStress(mMode);

    GaussPoint *rChGP = status->giveViscoelasticGaussPoint();

    // subtract stress independent part = temperature
    // slave (viscoelastic) material uses incremental formulation
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain, tStep, VM_Total);

    rChM->giveRealStressVector(tempEffectiveStress, rChGP, reducedStrain, tStep);

    const double e0 = this->give(e0_ID, gp) * this->e0Mean;
    status->setE0(e0);

    double f, equivStrain, tempKappa, omega = 0.;


    // transform effective stress computed by the viscoelastic model into strain by
    // dividing it by a stiffness matrix (which is diagonal)
    FloatMatrix elasticStiffnessMatrix;

    if ( mMode == _2dLattice ) {
        LatticeLinearElastic :: give2dLatticeStiffMtrx(elasticStiffnessMatrix, ElasticStiffness, gp, tStep);
    } else { // 3d
        LatticeLinearElastic :: give3dLatticeStiffMtrx(elasticStiffnessMatrix, ElasticStiffness, gp, tStep);
    }

    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );

    quasiElasticStrain.resize(rsize);
    quasiElasticStrain.zero();
    for ( int i = 1; i <= rsize; i++ ) { // only diagonal terms matter
        quasiElasticStrain.at(i) = tempEffectiveStress.at(i) / elasticStiffnessMatrix.at(i, i);
    }

    // compute equivalent strain from "quasielastic" strain
    this->computeEquivalentStrain(equivStrain, quasiElasticStrain, gp, tStep);


    // compute value of loading function if strainLevel crit apply
    f = equivStrain - status->giveKappa();

    if ( f <= 0.0 ) {
        // damage does not grow
        tempKappa = status->giveKappa();
        omega = status->giveDamage();
        if ( status->giveCrackFlag() != 0 ) {
            status->setTempCrackFlag(2);
        } else {
            status->setTempCrackFlag(0);
        }
    } else {
        // damage grows
        tempKappa = equivStrain;

        // evaluate damage parameter

        this->computeDamageParam(omega, tempKappa, gp);
        if ( omega > 0 ) {
            status->setTempCrackFlag(1);
        }
    }

    answer.resize(rsize);
    answer.zero();
    answer.add(tempEffectiveStress);
    answer.times(1. - omega);


    //Compute crack width
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();

    double crackWidth;

    // use total formulation to get mechanical strain necessary to get approximate value of the crack width
    // how to consider creep? it is also needed to subtract total shrinkage strain
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain, tStep, VM_Total);

    FloatArray shrinkageStrain;
    rChM->giveShrinkageStrainVector(shrinkageStrain, rChGP, tStep, VM_Total);
    reducedStrain.subtract(shrinkageStrain);

    // crack width is only approximate!!!
    if ( gp->giveMaterialMode() == _2dLattice ) {
        crackWidth = omega * sqrt( pow(reducedStrain.at(1), 2.) + pow(reducedStrain.at(2), 2.) ) * length;
    } else {
        crackWidth = omega * sqrt( pow(reducedStrain.at(1), 2.) + pow(reducedStrain.at(2), 2.) + pow(reducedStrain.at(3), 2.) ) * length;
    }


    // compute dissipation and store it. How does it work for a viscoelastic material with damage?

    // reduced strain = total strain - temperature effects - shrinkage

    status->setTempEquivalentStrain(equivStrain);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempReducedStrainBe(reducedStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);

    status->setTempNormalStress( answer.at(1) );
    status->setTempCrackWidth(crackWidth);


    return;
}



RheoChainMaterial *
LatticeDamageViscoelastic :: giveViscoelasticMaterial() {
    Material *mat;
    RheoChainMaterial *rChMat;
    mat = domain->giveMaterial(slaveMat);

    rChMat = dynamic_cast< RheoChainMaterial * >(mat);

    return rChMat;
}


void
LatticeDamageViscoelastic :: give2dLatticeStiffMtrx(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
    /* Returns elastic moduli in reduced stress-strain space*/

    LatticeDamageViscoelasticStatus *status = static_cast< LatticeDamageViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;
    RheoChainMaterial *rChMat;
    double Eincr;

    slaveGp = status->giveViscoelasticGaussPoint();
    rChMat = giveViscoelasticMaterial();

    LatticeDamage :: give2dLatticeStiffMtrx(answer, rmode, gp, tStep);

    Eincr = rChMat->giveEModulus(slaveGp, tStep);

    // lattice damage model uses normal stiffness different from unity, it is necessary to rescale it appropriately
    answer.times(Eincr / this->eNormalMean);
}

void
LatticeDamageViscoelastic :: give3dLatticeStiffMtrx(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
    /* Returns elastic moduli in reduced stress-strain space*/

    LatticeDamageViscoelasticStatus *status = static_cast< LatticeDamageViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;
    RheoChainMaterial *rChMat;
    double Eincr;

    slaveGp = status->giveViscoelasticGaussPoint();
    rChMat = giveViscoelasticMaterial();

    LatticeDamage :: give3dLatticeStiffMtrx(answer, rmode, gp, tStep);

    Eincr = rChMat->giveEModulus(slaveGp, tStep);

    answer.times(Eincr / this->eNormalMean);
}

int
LatticeDamageViscoelastic :: giveIPValue(FloatArray &answer,
                                         GaussPoint *gp,
                                         InternalStateType type,
                                         TimeStep *tStep)
{
    return LatticeDamage :: giveIPValue(answer, gp, type, tStep);
}


LatticeDamageViscoelasticStatus :: LatticeDamageViscoelasticStatus(int n, Domain *d, GaussPoint *g, int s) :
    LatticeDamageStatus(g), slaveMat(s)
{
    viscoelasticGP = new GaussPoint( g->giveIntegrationRule(), g->giveNumber(), g->giveNaturalCoordinates(), g->giveWeight(), g->giveMaterialMode() );
}

void
LatticeDamageViscoelasticStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticeDamageStatus :: initTempStatus();
}

void
LatticeDamageViscoelasticStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    MaterialStatus *mS = this->giveViscoelasticMatStatus();

    LatticeDamageStatus :: printOutputAt(file, tStep);
    fprintf(file, "\nViscoelastic material:");

    mS->printOutputAt(file, tStep);

    fprintf(file, "\n");
}

MaterialStatus *
LatticeDamageViscoelasticStatus :: giveViscoelasticMatStatus() {
    Material *mat;
    RheoChainMaterial *rChMat;
    GaussPoint *rChGP;

    //    mat = domain->giveMaterial(slaveMat);
    //    rChMat = dynamic_cast< RheoChainMaterial * >(mat);

    rChGP = this->giveViscoelasticGaussPoint();

    MaterialStatus *mS = rChMat->giveStatus(rChGP);

    return mS;
}


void
LatticeDamageViscoelasticStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    MaterialStatus *mS = this->giveViscoelasticMatStatus();
    mS->updateYourself(tStep);

    LatticeDamageStatus :: updateYourself(tStep);
}

void
LatticeDamageViscoelasticStatus :: saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{

    // save parent class status
  LatticeDamageStatus :: saveContext(stream, mode);

  this->giveViscoelasticMatStatus()->saveContext(stream, mode);
  
}

void
LatticeDamageViscoelasticStatus :: restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    // read parent class status
  LatticeDamageStatus :: restoreContext(stream, mode);

  this->giveViscoelasticMatStatus()->restoreContext(stream, mode);
}
  
}     // end namespace oofem
