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
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "../sm/Elements/Bars/truss3dnl.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Materials/structuralms.h"
#include "fei3dlinelin.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(Truss3dnl);

Truss3dnl :: Truss3dnl(int n, Domain *aDomain) : Truss3d(n, aDomain)
{
}


IRResultType
Truss3dnl :: initializeFrom(InputRecord *ir)
{
  IRResultType result = Truss3d :: initializeFrom(ir);
  if ( result != IRRT_OK ) {
    return result;
  }
  initialStretch = 1;
  IR_GIVE_OPTIONAL_FIELD(ir, initialStretch, _IFT_Truss3dnl_initialStretch);
  
  return IRRT_OK;      
}



  
void
Truss3dnl :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  FloatMatrix B, Be;
  FloatArray vStress, vStrain, u;
  
  // This function can be quite costly to do inside the loops when one has many slave dofs.
  this->computeVectorOf(VM_Total, tStep, u);
  // subtract initial displacements, if defined
  if ( initialDisplacements ) {
    u.subtract(* initialDisplacements);
  }
  
  // zero answer will resize accordingly when adding first contribution
  answer.clear();
  
  for ( auto &gp: *this->giveDefaultIntegrationRulePtr() ) {
    StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
    this->computeBmatrixAt(gp, B, tStep, true);
    this->computeBmatrixAt(gp, Be, tStep);
    if ( useUpdatedGpRecord == 1 ) {
      vStress = matStat->giveStressVector();
    } else {
      ///@todo Is this really what we should do for inactive elements?
      if ( !this->isActivated(tStep) ) {
	vStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
	vStrain.zero();
      }
      vStrain.beProductOf(Be, u);
      // add influence of initial stress/stretch
      double l2 = initialStretch*initialStretch;
      vStrain.times(l2);
      FloatArray E0(1);
      E0.at(1) = (l2-1.)/2.;
      vStrain.add(E0);
      //
      this->computeStressVector(vStress, vStrain, gp, tStep);
    }
    
    if ( vStress.giveSize() == 0 ) { /// @todo is this check really necessary?
      break;
    }
    
    // Compute nodal internal forces at nodes as f = B^T*Stress dV
    double dV  = this->computeVolumeAround(gp);
    
    if ( vStress.giveSize() == 6 ) {
      // It may happen that e.g. plane strain is computed
      // using the default 3D implementation. If so,
      // the stress needs to be reduced.
      // (Note that no reduction will take place if
      //  the simulation is actually 3D.)
      FloatArray stressTemp;
      StructuralMaterial :: giveReducedSymVectorForm( stressTemp, vStress, gp->giveMaterialMode() );
      answer.plusProduct(B, stressTemp, dV);
    } else   {
      answer.plusProduct(B, vStress, dV);
    }
    
    
    // If inactive: update fields but do not give any contribution to the internal forces
    if ( !this->isActivated(tStep) ) {
      answer.zero();
      return;
    }
  }
}
  
  
  
  
void
Truss3dnl :: computeStiffnessMatrix(FloatMatrix &answer,
				    MatResponseMode rMode, TimeStep *tStep)
{
  StructuralCrossSection *cs = this->giveStructuralCrossSection();
  bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);
  
  answer.clear();
  
  if ( !this->isActivated(tStep) ) {
    return;
  }
  
  // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
  if ( integrationRulesArray.size() == 1 ) {
    FloatMatrix B, D, DB, Ksigma;
    for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
      this->computeBmatrixAt(gp, B, tStep, true);
      this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
      double dV = this->computeVolumeAround(gp);
      DB.beProductOf(D, B);
      if ( matStiffSymmFlag ) {
	answer.plusProductSymmUpper(B, DB, dV);
      } else {
	answer.plusProductUnsym(B, DB, dV);
      }
      this->computeInitialStressStiffness(Ksigma, gp, tStep);
      Ksigma.times(dV);
      answer.add(Ksigma);
      
    }
    
    if ( matStiffSymmFlag ) {
      answer.symmetrized();
    }
  }
}
  
  
void
Truss3dnl :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool lin)
{
  FloatMatrix Bl, Bnl;
  this->computeBlMatrixAt(gp, Bl);
  this->computeBnlMatrixAt(gp, Bnl, tStep, lin);
  answer = Bl;
  answer.add(Bnl);
}



void
Truss3dnl :: computeBlMatrixAt(GaussPoint *gp, FloatMatrix &answer)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
  Truss3d::computeBmatrixAt(gp, answer);
}



void
Truss3dnl :: computeBnlMatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool lin)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
  FloatArray d;
  this->computeVectorOf(VM_Total, tStep, d);
    
  FloatMatrix Bnl, A(6,6);
  A.at(1,1) = A.at(2,2) = A.at(3,3) = A.at(4,4) = A.at(5,5) = A.at(6,6) =  1.0;
  A.at(1,4) = A.at(2,5) = A.at(3,6) = A.at(4,1) = A.at(5,2) = A.at(6,3) = -1.0;
  double l0 = this->computeLength();
  double factor = 1/l0/l0;
  if(!lin) {
    factor /= 2;
  } 
  Bnl.beProductOf(A,d);
  Bnl.times(factor);
  answer.beTranspositionOf(Bnl);
  
}
  
  
void
Truss3dnl :: computeInitialStressStiffness(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
  
  answer.resize(6,6);
  answer.at(1,1) = answer.at(2,2) = answer.at(3,3) = answer.at(4,4) = answer.at(5,5) = answer.at(6,6) =  1.0;
  answer.at(1,4) = answer.at(2,5) = answer.at(3,6) = answer.at(4,1) = answer.at(5,2) = answer.at(6,3) = -1.0;
  
  FloatArray d, stress, strain;
  FloatMatrix B;
  this->computeVectorOf(VM_Total, tStep, d);
  this->computeBmatrixAt(gp, B, tStep);	  
  strain.beProductOf(B, d);
  // add influence of initial stress/stretch
  double l2 = initialStretch*initialStretch;
  strain.times(l2);
  FloatArray E0(1);
  E0.at(1) = (l2-1.)/2;
  strain.add(E0);
  /////////////////////////////////////////////////////////////////////////////////////////
  this->giveStructuralCrossSection()->giveRealStress_1d(stress, gp, strain, tStep);
  double l0 = this->computeLength();	
  double factor = 1/l0/l0;
  // prevent zero initial stress stiffness
  if(stress.at(1) == 0) {
    stress.at(1) = 1;
  }
  answer.times(stress.at(1));
  answer.times(factor);
  
}
    

} // end namespace oofem

