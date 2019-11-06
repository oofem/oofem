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

#include "latticeplasticitydamageviscoelastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(LatticePlasticityDamageViscoelastic);

LatticePlasticityDamageViscoelastic :: LatticePlasticityDamageViscoelastic(int n, Domain *d) : LatticePlasticityDamage(n, d)
{
    slaveMat = 0;
}


LatticePlasticityDamageViscoelastic :: ~LatticePlasticityDamageViscoelastic()
//
// destructor
//
{}

int
LatticePlasticityDamageViscoelastic :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    if ( mode == _3dLattice ) {
        return 1;
    }

    return 0;
}


void
LatticePlasticityDamageViscoelastic :: initializeFrom(InputRecord &ir)
{

    LatticePlasticityDamage :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, slaveMat, _IFT_LatticePlasticityDamageViscoelastic_slaveMat); // number of slave material

    RheoChainMaterial *rChM = giveViscoelasticMaterial();

    // override elastic modulus by the one given by the compliance function
    // provided values of the strain at peak stress must be changed

    double E28 = 1. / rChM->computeCreepFunction(28.01, 28., NULL, NULL); // modulus of elasticity evaluated at 28 days, duration of loading 15 min
    this->eNormalMean = E28; // swap elastic modulus/stiffness


    // must match with the viscoelastic material
    if ( rChM->giveAlphaOne() != this->alphaOne ) {
        OOFEM_ERROR("a1 must be set to the same value in both master and viscoelastic slave materials");
    }

    // must match with the viscoelastic material
    if ( rChM->giveAlphaTwo() != this->alphaTwo ) {
        OOFEM_ERROR("a2 must be set to the same value in both master and viscoelastic slave materials");
    }

}


MaterialStatus *
LatticePlasticityDamageViscoelastic :: CreateStatus(GaussPoint *gp) const
{
    LatticePlasticityDamageViscoelasticStatus *answer = new LatticePlasticityDamageViscoelasticStatus(1, LatticePlasticityDamageViscoelastic :: domain, gp, slaveMat);
    return answer;
}

void
LatticePlasticityDamageViscoelastic :: giveReducedStrain(FloatArray &answer,
                                                         GaussPoint *gp,
                                                         TimeStep *tStep)
{
    LatticePlasticityDamageStatus *status = ( LatticePlasticityDamageStatus * ) ( this->giveStatus(gp) );

    FloatMatrix elasticStiffnessMatrix;
    LatticeLinearElastic :: give3dLatticeStiffMtrx(elasticStiffnessMatrix, ElasticStiffness, gp, tStep);

    answer = status->giveStressVector();

    double omega = status->giveDamage();
    answer.times(1. / ( 1. - omega ) );


    FloatArray oldPlasticStrain;
    oldPlasticStrain = status->givePlasticStrain();

    for ( int i = 1; i <= 6; i++ ) { // only diagonal terms matter
        answer.at(i) /= elasticStiffnessMatrix.at(i, i);
    }

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i) += oldPlasticStrain.at(i);
    }
}

void
LatticePlasticityDamageViscoelastic :: giveRealStressVector(FloatArray &answer,
                                                            GaussPoint *gp,
                                                            const FloatArray &totalStrain,
                                                            TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    double tol = 1.e-12;

    LatticePlasticityDamageViscoelasticStatus *status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );

    //    this->initGpForNewStep(gp);

    // total deformation - temperature - plastic strain
    FloatArray reducedStrainForViscoMat;

    // VE stress / stiffness
    FloatArray quasiElasticStrain;

    // VE stress / stiffness + plastic strain
    FloatArray reducedStrain;

    RheoChainMaterial *rChM = giveViscoelasticMaterial();

    MaterialMode mMode = gp->giveMaterialMode();
    StressVector tempEffectiveStress(mMode);

    GaussPoint *rChGP = status->giveViscoelasticGaussPoint();

    FloatArray tempPlasticStrain;

    FloatMatrix elasticStiffnessMatrix;
    LatticeLinearElastic :: give3dLatticeStiffMtrx(elasticStiffnessMatrix, ElasticStiffness, gp, tStep);

    FloatArray strain(3);

    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() );

    int itercount = 0;
    double prevIterKappaP = 0.;

    while ( ( fabs(status->giveTempKappaP() - prevIterKappaP) >= tol ) || ( itercount == 0 ) ) {
        itercount++;

        if ( itercount > 100 ) {
            OOFEM_ERROR("Algorithm not converging");
        }

        prevIterKappaP = status->giveTempKappaP();

        // subtract stress independent part = temperature
        // slave (viscoelastic) material uses incremental formulation
        // must subtract not only temperature but also increment of plastic deformation
        this->giveStressDependentPartOfStrainVector(reducedStrainForViscoMat, gp, totalStrain, tStep, VM_Total);

        tempPlasticStrain = status->giveTempPlasticStrain();

        for ( int i = 1; i <= 3; i++ ) {
            reducedStrainForViscoMat.at(i) -= tempPlasticStrain.at(i);
        }

        rChM->giveRealStressVector(tempEffectiveStress, rChGP, reducedStrainForViscoMat, tStep);

        // transform effective stress computed by the viscoelastic model into strain by
        // dividing it by a stiffness matrix (which is diagonal)
        quasiElasticStrain.resize(rsize);
        for ( int i = 1; i <= rsize; i++ ) { // only diagonal terms matter
            quasiElasticStrain.at(i) = tempEffectiveStress.at(i) / elasticStiffnessMatrix.at(i, i);
        }

        // total deformation is needed = visco stress / stiffness + plastic strain
        for ( int i = 1; i <= 3; i++ ) {
            strain.at(i) = quasiElasticStrain.at(i) + tempPlasticStrain.at(i);
        }

        // this (parent) method needs a total, temperature-free strain (with plastic strain)
        // from now on new = updated plastic strain
        this->performPlasticityReturn(answer, gp, strain, tStep);
    } // end while loop

    reducedStrain = quasiElasticStrain;
    for ( int i = 1; i <= 3; i++ ) {
        reducedStrain.at(i) += tempPlasticStrain.at(i);
    }

    // performDamageEvolution input "strain" is used only to compute crack width. Therefore
    // it should correspond to stress-dependent strain deformation = total - shrinkage - temperature
    // but with plastic deformation

    FloatArray shrinkageStrain;
    rChM->giveShrinkageStrainVector(shrinkageStrain, rChGP, tStep, VM_Total);

    strain = reducedStrainForViscoMat;
    strain.subtract(shrinkageStrain);
    strain.resizeWithValues(3);
    tempPlasticStrain = status->giveTempPlasticStrain();
    strain.add(tempPlasticStrain);

    // strain is passed only to compute crack width, nothing else
    this->performDamageEvaluation(gp, strain);

    double omega;
    omega = status->giveTempDamage();

    answer.resizeWithValues(6);
    answer.at(1) = ( 1. - omega ) * answer.at(1);
    answer.at(2) = ( 1. - omega ) * answer.at(2);
    answer.at(3) = ( 1. - omega ) * answer.at(3);
    answer.at(4) = ( 1. - omega ) * quasiElasticStrain.at(4) * this->alphaTwo * this->eNormalMean;
    answer.at(5) = ( 1. - omega ) * quasiElasticStrain.at(5) * this->alphaTwo * this->eNormalMean;
    answer.at(6) = ( 1. - omega ) * quasiElasticStrain.at(6) * this->alphaTwo * this->eNormalMean;

    status->letTempStrainVectorBe(totalStrain);
    status->letTempReducedStrainBe(reducedStrain);
    status->letTempStressVectorBe(answer);

    return;
}

RheoChainMaterial *
LatticePlasticityDamageViscoelastic :: giveViscoelasticMaterial() {
    Material *mat;
    RheoChainMaterial *rChMat;
    mat = domain->giveMaterial(slaveMat);

    rChMat = dynamic_cast< RheoChainMaterial * >( mat );

    return rChMat;
}

void
LatticePlasticityDamageViscoelastic :: give2dLatticeStiffMtrx(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
    /* Returns elastic moduli in reduced stress-strain space*/
    OOFEM_ERROR("2d mode not implemented");
}

void
LatticePlasticityDamageViscoelastic :: give3dLatticeStiffMtrx(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
    /* Returns elastic moduli in reduced stress-strain space*/

    LatticePlasticityDamageViscoelasticStatus *status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;
    RheoChainMaterial *rChMat;
    double Eincr;

    slaveGp = status->giveViscoelasticGaussPoint();
    rChMat = giveViscoelasticMaterial();

    LatticePlasticityDamage :: give3dLatticeStiffMtrx(answer, rmode, gp, tStep);

    Eincr = rChMat->giveEModulus(slaveGp, tStep);

    answer.times(Eincr / this->eNormalMean);
}

int
LatticePlasticityDamageViscoelastic :: giveIPValue(FloatArray &answer,
                                                   GaussPoint *gp,
                                                   InternalStateType type,
                                                   TimeStep *tStep)
{
    return LatticePlasticityDamage :: giveIPValue(answer, gp, type, tStep);
}


LatticePlasticityDamageViscoelasticStatus :: LatticePlasticityDamageViscoelasticStatus(int n, Domain *d, GaussPoint *g, int s) :
    LatticePlasticityDamageStatus(n, d, g), slaveMat(s)
{
    viscoelasticGP = new GaussPoint(g->giveIntegrationRule(), g->giveNumber(), g->giveNaturalCoordinates(), g->giveWeight(), g->giveMaterialMode() );
}

void
LatticePlasticityDamageViscoelasticStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticePlasticityDamageStatus :: initTempStatus();
}

void
LatticePlasticityDamageViscoelasticStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    MaterialStatus *mS = this->giveViscoelasticMatStatus();

    LatticePlasticityDamageStatus :: printOutputAt(file, tStep);
    fprintf(file, "\nViscoelastic material:");

    mS->printOutputAt(file, tStep);

    fprintf(file, "\n");
}

MaterialStatus *
LatticePlasticityDamageViscoelasticStatus :: giveViscoelasticMatStatus() {
    Material *mat;
    RheoChainMaterial *rChMat;
    GaussPoint *rChGP;

    // @todo: Don't now how to fix this!
    //    mat = domain->Givematerial(slaveMat);
    //    rChMat = dynamic_cast< RheoChainMaterial * >(mat);

    rChGP = this->giveViscoelasticGaussPoint();

    MaterialStatus *mS = rChMat->giveStatus(rChGP);

    return mS;
}


void
LatticePlasticityDamageViscoelasticStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    MaterialStatus *mS = this->giveViscoelasticMatStatus();
    mS->updateYourself(tStep);

    LatticePlasticityDamageStatus :: updateYourself(tStep);
}

void
LatticePlasticityDamageViscoelasticStatus :: saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    LatticePlasticityDamageStatus :: saveContext(stream, mode);


    this->giveViscoelasticMatStatus()->saveContext(stream, mode);

    return;
}

void
LatticePlasticityDamageViscoelasticStatus :: restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticePlasticityDamageStatus :: restoreContext(stream, mode);

    this->giveViscoelasticMatStatus()->restoreContext(stream, mode);
}
}     // end namespace oofem
