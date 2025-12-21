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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

LatticeViscoelastic::LatticeViscoelastic(int n, Domain *d) : LatticeLinearElastic(n, d)
{}


void
LatticeViscoelastic::initializeFrom(InputRecord &ir)
{
    LatticeLinearElastic::initializeFrom(ir);

    IR_GIVE_FIELD(ir, viscoMat, _IFT_LatticeViscoelastic_viscoMat); // number of slave material

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    // fix to work at arbitrary time units!
    double E28 = 1. / rheoMat->computeCreepFunction(28.01, 28., NULL, NULL); // modulus of elasticity evaluated at 28 days, duration of loading 15 min
    this->eNormalMean = E28;  // swap elastic modulus/stiffness
}


std::unique_ptr<MaterialStatus> 
LatticeViscoelastic::CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<LatticeViscoelasticStatus>(gp);
}


FloatArrayF< 6 >
LatticeViscoelastic::giveLatticeStress3d(const FloatArrayF< 6 > &totalStrain,
                                         GaussPoint *gp,
                                         TimeStep *tStep)
{
    auto status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    FloatArray viscoStress;
    FloatArray partialStrain;

    GaussPoint *rChGP = status->giveSlaveGaussPointVisco();
    double Eincr = rheoMat->giveEModulus(rChGP, tStep);
    this->eNormalMean = Eincr;

    FloatArrayF< 6 >reducedStrainForViscoMat;
    FloatArrayF< 6 >quasiTotalStrain;

    FloatArray indepStrain;

    FloatArrayF< 6 >tempStress;
    reducedStrainForViscoMat = totalStrain;

    indepStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);

    if ( indepStrain.giveSize() > 0 ) {
        reducedStrainForViscoMat -= FloatArrayF< 6 >(indepStrain);
    }


    rheoMat->giveRealStressVector(viscoStress, rChGP, reducedStrainForViscoMat, tStep);
    tempStress = FloatArrayF< 6 >(viscoStress);


    status->letTempLatticeStrainBe(totalStrain);
    status->letTempLatticeStressBe(tempStress);

    return tempStress;
}

FloatMatrixF< 6, 6 >
LatticeViscoelastic::give3dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const

{
    LatticeViscoelasticStatus *status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;

    // get status of the slave viscoelastic material
    slaveGp = status->giveSlaveGaussPointVisco();
    // get viscoelastic material
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    auto answer = LatticeLinearElastic::give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    double Eincr = rheoMat->giveEModulus(slaveGp, tStep);

    answer *= ( Eincr / this->eNormalMean );

    return answer;
}

FloatMatrixF< 3, 3 >
LatticeViscoelastic::give2dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
{
    LatticeViscoelasticStatus *status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;

    // get status of the slave viscoelastic material
    slaveGp = status->giveSlaveGaussPointVisco();
    // get viscoelastic material
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    auto answer = LatticeLinearElastic::give2dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    double Eincr = rheoMat->giveEModulus(slaveGp, tStep);
    answer *= ( Eincr / this->eNormalMean );

    return answer;
}


int
LatticeViscoelastic::giveIPValue(FloatArray &answer,
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
        LatticeViscoelasticStatus *status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );
        RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );
        return rheoMat->giveIPValue(answer, status->giveSlaveGaussPointVisco(), type, tStep);
    }

    return LatticeLinearElastic::giveIPValue(answer, gp, type, tStep);
}



int LatticeViscoelastic::checkConsistency()
{
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );


    if ( rheoMat->giveAlphaOne() != this->alphaOne ) {
        OOFEM_ERROR("a1 must be set to the same value in both master and viscoelastic slave materials");
    }

    if ( rheoMat->giveAlphaTwo() != this->alphaTwo ) {
        OOFEM_ERROR("a2 must be set to the same value in both master and viscoelastic slave materials");
    }


    return FEMComponent::checkConsistency();
}





LatticeViscoelasticStatus::LatticeViscoelasticStatus(GaussPoint *gp) :
    LatticeMaterialStatus(gp),
    slaveGpVisco( std::make_unique< GaussPoint >( gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode() ) )
    //    slaveGpVisco(std::make_unique<GaussPoint>( gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode()) )


{}

void
LatticeViscoelasticStatus::initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticeMaterialStatus::initTempStatus();

    RheoChainMaterialStatus *rheoStatus = static_cast< RheoChainMaterialStatus * >( this->giveSlaveGaussPointVisco()->giveMaterialStatus() );
    rheoStatus->initTempStatus();
}

void
LatticeViscoelasticStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    //    MaterialStatus *mS = this->giveViscoelasticMatStatus();

    LatticeMaterialStatus::printOutputAt(file, tStep);
    fprintf(file, "\nViscoelastic material:");

    //this->slaveGpVisco->giveMaterialStatus()->printOutputAt(file, tStep);
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->printOutputAt(file, tStep);

    fprintf(file, "\n");
}

void
LatticeViscoelasticStatus::updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    //  this->slaveGpVisco->giveMaterialStatus()->updateYourself(tStep);
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->updateYourself(tStep);

    LatticeMaterialStatus::updateYourself(tStep);
}

void
LatticeViscoelasticStatus::saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    // save parent class status
    LatticeMaterialStatus::saveContext(stream, mode);

    //    this->slaveGpVisco->giveMaterialStatus()->saveContext(stream, mode);
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->saveContext(stream, mode);
}

void
LatticeViscoelasticStatus::restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticeMaterialStatus::restoreContext(stream, mode);

    //  this->slaveGpVisco->giveMaterialStatus()->restoreContext(stream, mode);
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->restoreContext(stream, mode);
}
}     // end namespace oofem
