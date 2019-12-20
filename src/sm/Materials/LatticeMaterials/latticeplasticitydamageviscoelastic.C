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
    double E28 = 1. / rheoMat->computeCreepFunction(28.01, 28., NULL, NULL); // modulus of elasticity evaluated at 28 days, duration of loading 15 min
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
    status->initTempStatus();

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    GaussPoint *rChGP = status->giveSlaveGaussPointVisco();

    // total deformation - temperature - plastic strain
    FloatArrayF< 6 >reducedStrainForViscoMat;

    // VE stress / stiffness + plastic strain
    FloatArrayF< 6 >reducedStrain;

    // components 1-3 of the reducedStrain
    FloatArrayF< 6 >quasiTotalStrain;

    // temporary value of the plastic strain
    FloatArrayF< 6 >tempPlasticStrain;

    // return stress from viscoelastic material - effective stress
    FloatArray viscoStress;
    // answer - viscoelastic stress multiplied by 1-omega - nominal stress
    FloatArrayF< 6 >tempStress;

    // answer from the plasticity stress-return algorithm
    FloatArrayF< 6 >stress;

    auto elasticStiffnessMatrix = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    int itercount = 1;

    auto indepStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);

    double tolerance = 1.;


    do {
        if ( itercount > 100 ) {
            OOFEM_ERROR("Algorithm not converging");
        }

        reducedStrainForViscoMat = totalStrain;

        if ( indepStrain.giveSize() > 0 ) {
            reducedStrainForViscoMat -= FloatArrayF< 6 >(indepStrain);
        }

        tempPlasticStrain = status->giveTempPlasticLatticeStrain();
        for ( int i = 1; i <= 3; i++ ) {
            reducedStrainForViscoMat.at(i) -= tempPlasticStrain.at(i);
        }
	rheoMat->giveRealStressVector(viscoStress, rChGP, reducedStrainForViscoMat, tStep);
	tempStress = FloatArrayF< 6 >(viscoStress);
        for ( int i = 1; i <= 6; i++ ) {
            quasiTotalStrain.at(i) =  tempStress.at(i) / elasticStiffnessMatrix.at(i, i) + tempPlasticStrain.at(i) + indepStrain.at(i);
        }

        stress = LatticePlasticityDamage :: giveLatticeStress3d(quasiTotalStrain, gp, tStep);

        tolerance = norm(stress - tempStress) / this->eNormalMean;

        itercount++;
    } while ( tolerance >= tol );

    status->letTempLatticeStrainBe(totalStrain);

    return tempStress;
}


FloatMatrixF< 6, 6 >
LatticePlasticityDamageViscoelastic :: give3dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
{
    LatticePlasticityDamageViscoelasticStatus *status = static_cast< LatticePlasticityDamageViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;


    auto answer = LatticePlasticityDamage :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    // get status of the slave viscoelastic material
    slaveGp = status->giveSlaveGaussPointVisco();
    // get viscoelastic material
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    double Eincr = rheoMat->giveEModulus(slaveGp, tStep);

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


    return FEMComponent :: checkConsistency();
}


LatticePlasticityDamageViscoelasticStatus :: LatticePlasticityDamageViscoelasticStatus(int n, Domain *d, GaussPoint *gp) :
    LatticePlasticityDamageStatus(n, d, gp), slaveGpVisco(std :: make_unique< GaussPoint >(gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode() ) )
{}

void
LatticePlasticityDamageViscoelasticStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticePlasticityDamageStatus :: initTempStatus();

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
    //    this->slaveGpVisco->giveMaterialStatus()->updateYourself(tStep);
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
