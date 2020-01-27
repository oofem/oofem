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
{}

void
LatticePlasticityDamageViscoelastic :: initializeFrom(InputRecord &ir)
{
    LatticePlasticityDamage :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, viscoMat, _IFT_LatticePlasticityDamageViscoelastic_viscoMat); // number of slave material

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    // override elastic modulus by the one given by the compliance function
    double timeFactor = 0.;
    IR_GIVE_FIELD(ir, timeFactor, _IFT_LatticePlasticityDamageViscoelastic_timeFactor); // timeConversion equal to 1 day in current time units of the analysis

    double E28 = 1. / rheoMat->computeCreepFunction(timeFactor * 28.01, timeFactor * 28., NULL, NULL); // modulus of elasticity evaluated at 28 days, duration of loading 15 min

    this->eNormalMean = E28; // swap elastic modulus/stiffness
}


MaterialStatus *
LatticePlasticityDamageViscoelastic :: CreateStatus(GaussPoint *gp) const
{
    return new LatticePlasticityDamageViscoelasticStatus(1, LatticePlasticityDamageViscoelastic :: domain, gp);
}


FloatArrayF< 6 >
LatticePlasticityDamageViscoelastic :: giveLatticeStress3d(const FloatArrayF< 6 > &totalStrain,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep)
{
    double tol = 1.e-12; // error in order of Pascals

    auto status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    GaussPoint *rChGP = status->giveSlaveGaussPointVisco();
    // just a dummy yet essential call to create a status of viscomaterial. Otherwise initTempStatus() would fail.
    rheoMat->giveStatus(rChGP);

    status->initTempStatus();

    FloatArrayF< 6 >reducedStrainForViscoMat;
    FloatArrayF< 6 >reducedStrain;
    FloatArrayF< 6 >quasiReducedStrain;
    FloatArrayF< 6 >tempPlasticStrain;

    FloatArray viscoStress;
    FloatArrayF< 6 >tempStress;
    FloatArrayF< 6 >stress;

    int itercount = 1;

    FloatArray indepStrain;
    indepStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);

    auto elasticStiffnessMatrix = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    double tolerance = 1.;


    do {
        if ( itercount > 100 ) {
            printf("tolerance = %e and itercount = %d\n", tolerance, itercount);
            OOFEM_ERROR("Algorithm not converging");
        }

        reducedStrainForViscoMat = totalStrain;

        if ( indepStrain.giveSize() > 0 ) {
            reducedStrainForViscoMat -= FloatArrayF< 6 >(indepStrain);
        }

        tempPlasticStrain = status->giveTempPlasticLatticeStrain();
        reducedStrainForViscoMat -= tempPlasticStrain;

        FloatArrayF< 6 >tempDamageLatticeStrain = status->giveTempDamageLatticeStrain();
        reducedStrainForViscoMat -= FloatArrayF< 6 >(tempDamageLatticeStrain);

        rheoMat->giveRealStressVector(viscoStress, rChGP, reducedStrainForViscoMat, tStep);
        tempStress = FloatArrayF< 6 >(viscoStress);


        for ( int i = 1; i <= 6; i++ ) {
            quasiReducedStrain.at(i) =  tempStress.at(i) / elasticStiffnessMatrix.at(i, i) + tempPlasticStrain.at(i) + tempDamageLatticeStrain.at(i);
        }

        stress = this->performPlasticityReturn(gp, quasiReducedStrain, tStep);
        this->performDamageEvaluation(gp, quasiReducedStrain);
        double tempDamage = status->giveTempDamage();
        stress *= ( 1. - tempDamage );

        tolerance = norm(stress - tempStress) / this->eNormalMean;
        itercount++;

        //	printf("tolerance = %e, itercount = %d\n", tolerance, itercount);
    } while ( tolerance >= tol );

    status->letTempLatticeStrainBe(totalStrain);
    status->letTempLatticeStressBe(tempStress);
    status->letTempReducedLatticeStrainBe(quasiReducedStrain);

    return tempStress;
}


FloatMatrixF< 6, 6 >
LatticePlasticityDamageViscoelastic :: give3dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
{
    LatticePlasticityDamageViscoelasticStatus *status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;

    // get status of the slave viscoelastic material
    slaveGp = status->giveSlaveGaussPointVisco();

    // get viscoelastic material
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    double Eincr = rheoMat->giveEModulus(slaveGp, tStep);

    auto answer = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    answer *= ( Eincr / this->eNormalMean );

    return answer;
}

int
LatticePlasticityDamageViscoelastic :: giveIPValue(FloatArray &answer,
                                                   GaussPoint *gp,
                                                   InternalStateType type,
                                                   TimeStep *tStep)
{
    if ( ( type == IST_DryingShrinkageTensor ) ||
         ( type == IST_AutogenousShrinkageTensor ) ||
         ( type == IST_TotalShrinkageTensor ) ||
         ( type == IST_CreepStrainTensor ) ||
         ( type == IST_DryingShrinkageTensor ) ||
         ( type == IST_ThermalStrainTensor ) ) {
        LatticePlasticityDamageViscoelasticStatus *status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );

        RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );
        return rheoMat->giveIPValue(answer, status->giveSlaveGaussPointVisco(), type, tStep);
    }

    return LatticePlasticityDamage :: giveIPValue(answer, gp, type, tStep);
}


int LatticePlasticityDamageViscoelastic :: checkConsistency()
{
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    if ( rheoMat->giveAlphaOne() != this->alphaOne ) {
        OOFEM_ERROR("a1 must be set to the same value in both master and viscoelastic slave materials");
    }

    if ( rheoMat->giveAlphaTwo() != this->alphaTwo ) {
        OOFEM_ERROR("a2 must be set to the same value in both master and viscoelastic slave materials");
    }

    GaussPoint *noGP = NULL;
    if ( rheoMat->give(tAlpha, noGP) != 0. ) {
        OOFEM_ERROR("tAlpha must be set to 0. in slave viscoelastic material");
    }


    return FEMComponent :: checkConsistency();
}


LatticePlasticityDamageViscoelasticStatus :: LatticePlasticityDamageViscoelasticStatus(int n, Domain *d, GaussPoint *gp) :
    LatticePlasticityDamageStatus(n, d, gp), slaveGpVisco(std :: make_unique< GaussPoint >(gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode() ) )
{}

void
LatticePlasticityDamageViscoelasticStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
{
    LatticePlasticityDamageStatus :: initTempStatus();

    // at this point rheomat :: giveStatus (viscoGP) has to be called first
    RheoChainMaterialStatus *rheoStatus = static_cast< RheoChainMaterialStatus * >( this->giveSlaveGaussPointVisco()->giveMaterialStatus() );

    rheoStatus->initTempStatus();
}

void
LatticePlasticityDamageViscoelasticStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticePlasticityDamageStatus :: printOutputAt(file, tStep);
    fprintf(file, "\nViscoelastic material:");

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->printOutputAt(file, tStep);

    fprintf(file, "\n");
}


void
LatticePlasticityDamageViscoelasticStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->updateYourself(tStep);

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

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->saveContext(stream, mode);
}

void
LatticePlasticityDamageViscoelasticStatus :: restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticePlasticityDamageStatus :: restoreContext(stream, mode);

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->restoreContext(stream, mode);
}
}     // end namespace oofem
