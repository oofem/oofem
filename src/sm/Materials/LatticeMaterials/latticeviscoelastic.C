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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#include "latticeviscoelastic.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticeViscoelastic);

LatticeViscoelastic :: LatticeViscoelastic(int n, Domain *d) : LatticeLinearElastic(n, d)
{
    slaveMat = 0;
}


bool
LatticeViscoelastic :: hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _2dLattice ) || ( mode == _3dLattice );
}


void
LatticeViscoelastic :: initializeFrom(InputRecord &ir)
{

    LatticeLinearElastic :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, slaveMat, _IFT_LatticeViscoelastic_slaveMat); // number of slave material

    RheoChainMaterial *rChM = giveViscoelasticMaterial();

    if ( rChM->giveAlphaOne() != this->alphaOne ) {
        OOFEM_ERROR("a1 must be set to the same value in both master and viscoelastic slave materials");
    }


    if ( rChM->giveAlphaTwo() != this->alphaTwo ) {
        OOFEM_ERROR("a2 must be set to the same value in both master and viscoelastic slave materials");
    }

}


MaterialStatus *
LatticeViscoelastic :: CreateStatus(GaussPoint *gp) const
{
    LatticeViscoelasticStatus *answer = new LatticeViscoelasticStatus(gp);
    return answer;
}



void
LatticeViscoelastic :: giveRealStressVector(FloatArray &answer,
                                            GaussPoint *gp,
                                            const FloatArray &totalStrain,
                                            TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    LatticeViscoelasticStatus *status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );

    //    this->initGpForNewStep(gp);

    FloatArray reducedStrain;


    // evaluate viscoelastic material:
    // pass: viscoelasticGP, temperature-free increment of strain

    RheoChainMaterial *rChM = giveViscoelasticMaterial();

    MaterialMode mode = gp->giveMaterialMode();
    StressVector tempEffectiveStress(mode);

    GaussPoint *rChGP = status->giveSlaveGaussPointVisco();

    // subtract stress independent part = temperature

    this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain, tStep, VM_Incremental);

    rChM->giveRealStressVector(tempEffectiveStress, rChGP, reducedStrain, tStep);

    answer.zero();
    answer.add(tempEffectiveStress);

    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);

    return;
}



RheoChainMaterial *
LatticeViscoelastic :: giveViscoelasticMaterial() {
    Material *mat;
    RheoChainMaterial *rChMat;
    mat = domain->giveMaterial(slaveMat);

    rChMat = dynamic_cast< RheoChainMaterial * >( mat );

    return rChMat;
}




void
LatticeViscoelastic :: give2dLatticeStiffMtrx(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
    /* Returns elastic moduli in reduced stress-strain space*/

    LatticeViscoelasticStatus *status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;
    RheoChainMaterial *rChMat;
    double Eincr;

    slaveGp = status->giveSlaveGaussPointVisco();
    rChMat = giveViscoelasticMaterial();

    LatticeLinearElastic :: give2dLatticeStiffMtrx(answer, rmode, gp, tStep);

    Eincr = rChMat->giveEModulus(slaveGp, tStep);

    answer.times(Eincr / this->eNormalMean);
}

void
LatticeViscoelastic :: give3dLatticeStiffMtrx(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep)
{
    /* Returns elastic moduli in reduced stress-strain space*/

    LatticeViscoelasticStatus *status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;
    RheoChainMaterial *rChMat;
    double Eincr;

    slaveGp = status->giveSlaveGaussPointVisco();
    rChMat = giveViscoelasticMaterial();

    LatticeLinearElastic :: give3dLatticeStiffMtrx(answer, rmode, gp, tStep);

    Eincr = rChMat->giveEModulus(slaveGp, tStep);

    answer.times(Eincr / this->eNormalMean);
}

int
LatticeViscoelastic :: giveIPValue(FloatArray &answer,
                                   GaussPoint *gp,
                                   InternalStateType type,
                                   TimeStep *tStep)
{
    return LatticeLinearElastic :: giveIPValue(answer, gp, type, tStep);
}


LatticeViscoelasticStatus :: LatticeViscoelasticStatus(GaussPoint *g) :
    LatticeMaterialStatus(g),
    slaveGpVisco(new GaussPoint(gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode() ) )
{}

void
LatticeViscoelasticStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticeMaterialStatus :: initTempStatus();
}

void
LatticeViscoelasticStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    MaterialStatus *mS = this->giveViscoelasticMatStatus();

    LatticeMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "\nViscoelastic material:");

    mS->printOutputAt(file, tStep);

    fprintf(file, "\n");
}

MaterialStatus *
LatticeViscoelasticStatus :: giveViscoelasticMatStatus() const {
  //    Material *mat;
  //    RheoChainMaterial *rChMat;
  //    GaussPoint *rChGP;

    //    mat = domain->giveMaterial(slaveMat);
    //    rChMat = dynamic_cast< RheoChainMaterial * >(mat);

  //    rChGP = this->giveSlaveGaussPointVisco();

    //    MaterialStatus *mS = rChMat->giveStatus(rChGP);

    //    return mS;
    return NULL;
}


void
LatticeViscoelasticStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    MaterialStatus *mS = this->giveViscoelasticMatStatus();
    mS->updateYourself(tStep);

    LatticeMaterialStatus :: updateYourself(tStep);
}

void
LatticeViscoelasticStatus :: saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    // save parent class status
    LatticeMaterialStatus :: saveContext(stream, mode);

    this->giveViscoelasticMatStatus()->saveContext(stream, mode);
}

void
LatticeViscoelasticStatus :: restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticeMaterialStatus :: restoreContext(stream, mode);

    this->giveViscoelasticMatStatus()->restoreContext(stream, mode);
}
}     // end namespace oofem
