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
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "math.h"


namespace oofem {
REGISTER_Material(LatticeDamageViscoelastic);

LatticeDamageViscoelastic :: LatticeDamageViscoelastic(int n, Domain *d) : LatticeDamage(n, d)
{}


void
LatticeDamageViscoelastic :: initializeFrom(InputRecord &ir)
{
    LatticeDamage :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, viscoMat, _IFT_LatticeDamageViscoelastic_viscoMat); // number of slave material

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    // override elastic modulus by the one given by the compliance function

    // fix to work at arbitrary time units!
    double E28 = 1. / rheoMat->computeCreepFunction(28.01, 28., NULL, NULL); // modulus of elasticity evaluated at 28 days, duration of loading 15 min
    this->e0Mean *= ( this->eNormalMean / E28 ); // transform strain at peak stress accordingly
    this->eNormalMean = E28; // swap elastic modulus/stiffness
}


MaterialStatus *
LatticeDamageViscoelastic :: CreateStatus(GaussPoint *gp) const
{
    return new LatticeDamageViscoelasticStatus(1, LatticeDamageViscoelastic :: domain, gp);
}


FloatArrayF< 6 >
LatticeDamageViscoelastic :: giveLatticeStress3d(const FloatArrayF< 6 > &totalStrain,
                                                 GaussPoint *gp,
                                                 TimeStep *tStep)
{
    LatticeDamageViscoelasticStatus *status = static_cast< LatticeDamageViscoelasticStatus * >( this->giveStatus(gp) );

    // the value from the status seems to be unused except for printout
    status->setE0(this->give(e0_ID, gp) * this->e0Mean);
    status->initTempStatus();

    FloatArrayF< 6 >reducedStrain;
    FloatArrayF< 6 >quasiElasticStrain;

    FloatArray viscoStress;
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );
    GaussPoint *rChGP = status->giveSlaveGaussPointVisco();

    reducedStrain = totalStrain;
    // subtract stress independent part = temperature (total temperature)
    auto indepStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( indepStrain.giveSize() > 0 ) {
        reducedStrain -= FloatArrayF< 6 >(indepStrain);
    }

    // size of viscoStress will correspond to the correct material mode: need to correct 2D case!
    FloatArrayF< 6 >tempStress;
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dLattice:
        OOFEM_ERROR("mode (%s) not implemented", __MaterialModeToString(mMode) );
        break;
    case _2dLattice:
        rheoMat->giveRealStressVector(viscoStress, rChGP, reducedStrain [ { 0, 1, 5 } ], tStep);
        //      tempStress.assemble<6> ( viscoStress, { 0, 1, 5 });
        // HELP!!!
        tempStress.at(1) = viscoStress.at(1);
        tempStress.at(2) = viscoStress.at(2);
        tempStress.at(6) = viscoStress.at(3);
        break;
    case _3dLattice:
        rheoMat->giveRealStressVector(viscoStress, rChGP, reducedStrain, tStep);
        tempStress = FloatArrayF< 6 >(viscoStress);
        break;
    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }

    // transform effective stress computed by the viscoelastic model into strain by
    // dividing it by a stiffness matrix (which is diagonal)
    auto elasticStiffnessMatrix = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    // loop over only first 3 components?
    for ( int i = 1; i <= 6; i++ ) { // only diagonal terms matter
        quasiElasticStrain.at(i) = tempStress.at(i) / elasticStiffnessMatrix.at(i, i);
    }

    // compute equivalent strain from "quasielastic" strain
    auto equivStrain =  this->computeEquivalentStrain(quasiElasticStrain, gp, tStep);

    // compute value of loading function if strainLevel crit apply
    auto f = equivStrain - status->giveKappa();

    double tempKappa, omega = 0.;
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
        omega = this->computeDamageParam(tempKappa, gp);
        if ( omega > 0. ) {
            status->setTempCrackFlag(1);
        }
    }

    tempStress *= ( 1. - omega );

    //Compute crack width - here it is assumed that also creep strain contributes to the crack width
    // another option would be to use quasiElasticStrain computed earlier here
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    double crackWidth = 0.;

    FloatArray shrinkageStrain;
    rheoMat->giveShrinkageStrainVector(shrinkageStrain, rChGP, tStep, VM_Total);

    switch ( mMode ) {
    case _1dLattice:
        OOFEM_ERROR("mode (%s) not implemented", __MaterialModeToString(mMode) );
        break;
    case _2dLattice:
        // HELP!!!
        reducedStrain.at(1) -= shrinkageStrain.at(1);
        reducedStrain.at(2) -= shrinkageStrain.at(2);
        reducedStrain.at(6) -= shrinkageStrain.at(3);
        crackWidth = omega * sqrt(pow(reducedStrain.at(1), 2.) + pow(reducedStrain.at(2), 2.) ) * length;
        break;
    case _3dLattice:
        reducedStrain -= FloatArrayF< 6 >(shrinkageStrain);
        crackWidth = omega * sqrt(pow(reducedStrain.at(1), 2.) + pow(reducedStrain.at(2), 2.) + pow(reducedStrain.at(3), 2.) ) * length;
        break;
    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }

    // compute dissipation and store it. How does it work for a viscoelastic material with damage?

    status->setTempEquivalentStrain(equivStrain);
    status->letTempLatticeStrainBe(totalStrain);
    status->letTempReducedLatticeStrainBe(reducedStrain);
    status->letTempLatticeStressBe(tempStress);
    status->setTempKappa(tempKappa);
    status->setTempDamage(omega);

    status->setTempNormalLatticeStress(tempStress.at(1) );
    status->setTempCrackWidth(crackWidth);

    return tempStress;
}


FloatMatrixF< 6, 6 >
LatticeDamageViscoelastic :: give3dLatticeStiffnessMatrix(MatResponseMode rmode,
                                                          GaussPoint *gp,
                                                          TimeStep *tStep) const
{
    LatticeDamageViscoelasticStatus *status = static_cast< LatticeDamageViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;


    auto answer = LatticeDamage :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    // get status of the slave viscoelastic material
    slaveGp = status->giveSlaveGaussPointVisco();
    // get viscoelastic material
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    double Eincr = rheoMat->giveEModulus(slaveGp, tStep);

    answer *= ( Eincr / this->eNormalMean );

    return answer;
}

FloatMatrixF< 3, 3 >
LatticeDamageViscoelastic :: give2dLatticeStiffnessMatrix(MatResponseMode rmode,
                                                          GaussPoint *gp,
                                                          TimeStep *tStep) const
{
    LatticeDamageViscoelasticStatus *status = static_cast< LatticeDamageViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;


    auto answer = LatticeDamage :: give2dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    // get status of the slave viscoelastic material
    slaveGp = status->giveSlaveGaussPointVisco();
    // get viscoelastic material
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    double Eincr = rheoMat->giveEModulus(slaveGp, tStep);

    answer *= ( Eincr / this->eNormalMean );

    return answer;
}



int
LatticeDamageViscoelastic :: giveIPValue(FloatArray &answer,
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
        LatticeDamageViscoelasticStatus *status = static_cast< LatticeDamageViscoelasticStatus * >( this->giveStatus(gp) );

        RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );
        return rheoMat->giveIPValue(answer, status->giveSlaveGaussPointVisco(), type, tStep);
    }

    return LatticeDamage :: giveIPValue(answer, gp, type, tStep);
}


int LatticeDamageViscoelastic :: checkConsistency()
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


LatticeDamageViscoelasticStatus :: LatticeDamageViscoelasticStatus(int n, Domain *d, GaussPoint *g) :
    LatticeDamageStatus(g), slaveGpVisco(std :: make_unique< GaussPoint >(gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode() ) )
{}


void
LatticeDamageViscoelasticStatus :: initTempStatus()
{
    LatticeDamageStatus :: initTempStatus();

    RheoChainMaterialStatus *rheoStatus = static_cast< RheoChainMaterialStatus * >( this->giveSlaveGaussPointVisco()->giveMaterialStatus() );
    rheoStatus->initTempStatus();
}


void
LatticeDamageViscoelasticStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeDamageStatus :: printOutputAt(file, tStep);
    fprintf(file, "\nViscoelastic material:");

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->printOutputAt(file, tStep);

    fprintf(file, "\n");
}


void
LatticeDamageViscoelasticStatus :: updateYourself(TimeStep *tStep)
{
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->updateYourself(tStep);

    LatticeDamageStatus :: updateYourself(tStep);
}

void
LatticeDamageViscoelasticStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    // save parent class status
    LatticeDamageStatus :: saveContext(stream, mode);

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->saveContext(stream, mode);
}

void
LatticeDamageViscoelasticStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    // read parent class status
    LatticeDamageStatus :: restoreContext(stream, mode);

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->restoreContext(stream, mode);
}
}     // end namespace oofem
