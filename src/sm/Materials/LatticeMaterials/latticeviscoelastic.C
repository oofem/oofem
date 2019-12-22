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
{}


void
LatticeViscoelastic :: initializeFrom(InputRecord &ir)
{
    LatticeLinearElastic :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, viscoMat, _IFT_LatticeViscoelastic_viscoMat); // number of slave material
}


MaterialStatus *
LatticeViscoelastic :: CreateStatus(GaussPoint *gp) const
{
    LatticeViscoelasticStatus *answer = new LatticeViscoelasticStatus(gp);
    return answer;
}


FloatArrayF< 6 >
LatticeViscoelastic :: giveLatticeStress3d(const FloatArrayF< 6 > &totalStrain,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
{
  double tol = 1.e-12; //
  
    auto status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    FloatArray viscoStress;
    FloatArray partialStrain;

    GaussPoint *rChGP = status->giveSlaveGaussPointVisco();

    FloatArrayF< 6 > reducedStrainForViscoMat;
    FloatArrayF< 6 >quasiTotalStrain;
    
    FloatArray indepStrain;

    FloatArray thermalStrain;
    
    FloatArrayF< 6> stress;
    FloatArrayF< 6 >tempStress;
    
    int itercount = 1;
    double tolerance = 1.;

    do {
        if ( itercount > 100 ) {
            OOFEM_ERROR("Algorithm not converging");
        }

	reducedStrainForViscoMat = totalStrain;
	
	indepStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
	
	if ( indepStrain.giveSize() > 0 ) {
	  reducedStrainForViscoMat -= FloatArrayF< 6 >(indepStrain);
	}

		
	rheoMat->giveRealStressVector(viscoStress, rChGP, reducedStrainForViscoMat, tStep);
	tempStress = FloatArrayF< 6 >(viscoStress);

	auto elasticStiffnessMatrix = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

	
    for ( int i = 1; i <= 6; i++ ) { // only diagonal terms matter
      quasiTotalStrain.at(i) = tempStress.at(i) / elasticStiffnessMatrix.at(i, i);
    }

    thermalStrain = LatticeLinearElastic :: computeStressIndependentStrainVector(gp, tStep, VM_Total);
    
    if ( thermalStrain.giveSize() > 0 ) {
        quasiTotalStrain += FloatArrayF< 6 >(thermalStrain);
    }
    
    stress = LatticeLinearElastic :: giveLatticeStress3d(quasiTotalStrain, gp, tStep);

    tolerance = norm(stress - tempStress) / this->eNormalMean;
    printf("tolerance = %e\n", tolerance);
    
     itercount++;
    } while ( tolerance >= tol );

	
    status->letTempLatticeStrainBe(totalStrain);
    status->letTempLatticeStressBe(tempStress);

    return tempStress;
}

FloatMatrixF< 6, 6 >
LatticeViscoelastic :: give3dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
{
    LatticeViscoelasticStatus *status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;


    auto answer = LatticeLinearElastic :: give3dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    // get status of the slave viscoelastic material
    slaveGp = status->giveSlaveGaussPointVisco();
    // get viscoelastic material
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    double Eincr = rheoMat->giveEModulus(slaveGp, tStep);

    answer *= ( Eincr / this->eNormalMean );

    return answer;
}

FloatMatrixF< 3, 3 >
LatticeViscoelastic :: give2dLatticeStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep) const
{
    LatticeViscoelasticStatus *status = static_cast< LatticeViscoelasticStatus * >( this->giveStatus(gp) );
    GaussPoint *slaveGp;


    auto answer = LatticeLinearElastic :: give2dLatticeStiffnessMatrix(ElasticStiffness, gp, tStep);

    // get status of the slave viscoelastic material
    slaveGp = status->giveSlaveGaussPointVisco();
    // get viscoelastic material
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    double Eincr = rheoMat->giveEModulus(slaveGp, tStep);

    // check if it works!
    answer *= ( Eincr / this->eNormalMean );

    return answer;
}


int
LatticeViscoelastic :: giveIPValue(FloatArray &answer,
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

    return LatticeLinearElastic :: giveIPValue(answer, gp, type, tStep);
}



int LatticeViscoelastic :: checkConsistency()
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





LatticeViscoelasticStatus :: LatticeViscoelasticStatus(GaussPoint *gp) :
    LatticeMaterialStatus(gp),
    slaveGpVisco(std :: make_unique< GaussPoint >(gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode() ) )
    //    slaveGpVisco(std::make_unique<GaussPoint>( gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode()) )


{}

void
LatticeViscoelasticStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticeMaterialStatus :: initTempStatus();

    RheoChainMaterialStatus *rheoStatus = static_cast< RheoChainMaterialStatus * >( this->giveSlaveGaussPointVisco()->giveMaterialStatus() );
    rheoStatus->initTempStatus();
}

void
LatticeViscoelasticStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    //    MaterialStatus *mS = this->giveViscoelasticMatStatus();

    LatticeMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "\nViscoelastic material:");

    //this->slaveGpVisco->giveMaterialStatus()->printOutputAt(file, tStep);
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->printOutputAt(file, tStep);

    fprintf(file, "\n");
}

void
LatticeViscoelasticStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    //  this->slaveGpVisco->giveMaterialStatus()->updateYourself(tStep);
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->updateYourself(tStep);

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

    //    this->slaveGpVisco->giveMaterialStatus()->saveContext(stream, mode);
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->saveContext(stream, mode);
}

void
LatticeViscoelasticStatus :: restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    LatticeMaterialStatus :: restoreContext(stream, mode);

    //  this->slaveGpVisco->giveMaterialStatus()->restoreContext(stream, mode);
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->restoreContext(stream, mode);
}
}     // end namespace oofem
