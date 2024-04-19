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

#include "../sm/Elements/Bars/truss3dnl2.h"
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
REGISTER_Element(Truss3dnl2);

Truss3dnl2 :: Truss3dnl2(int n, Domain *aDomain) : Truss3d(n, aDomain)
{
  cellGeometryWrapper = NULL;
}


void
Truss3dnl2 :: initializeFrom(InputRecord &ir)
{
  Truss3d :: initializeFrom(ir);
  X = this-> giveCellGeometryWrapper()->giveVertexCoordinates( 1 );
  X.append(this-> giveCellGeometryWrapper()->giveVertexCoordinates( 2 ));
  L = this->computeLength();
}

  
void
Truss3dnl2 :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
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
    this->computeBmatrixAt(gp, B, tStep, u);
    if ( useUpdatedGpRecord == 1 ) {
      vStress = matStat->givePVector();
    } else {
      ///@todo Is this really what we should do for inactive elements?
      if ( !this->isActivated(tStep) ) {
	vStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
	vStrain.zero();
      }
      // compute strain tensor, i.e., Biot strain
      auto vStrain = this->computeStrainVector(gp, u);
      // compute stress tensor, i.e., firt Piola-Kirchhoff
      vStress = this->giveStructuralCrossSection()->giveFirstPKStresses(vStrain, gp, tStep);
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


double
Truss3dnl2 :: computeLength()
{
  FloatArray X12;
  X12.beProductOf(this->givePmatrix(), X);
  return X12.computeNorm();
}
  
double
Truss3dnl2 :: computeDeformedLength(const FloatArray &d)
{
  FloatArray x12, x(X);
  x.add(d);
  x12.beProductOf(this->givePmatrix(), x);
  return x12.computeNorm();  
}
    



  
FloatArray
Truss3dnl2 :: computeStrainVector(GaussPoint *gp, const FloatArray &d)
{
  FloatArray answer(1);
  auto l = this->computeDeformedLength(d);
  answer.at(1) = l/this->computeLength();
  return answer;
}
  

void
Truss3dnl2 :: computeStiffnessMatrix(FloatMatrix &answer,
				    MatResponseMode rMode, TimeStep *tStep)
{
  StructuralCrossSection *cs = this->giveStructuralCrossSection();
  bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);
  
  answer.clear();
  
  if ( !this->isActivated(tStep) ) {
    return;
  }
  FloatArray u;
  this->computeVectorOf(VM_Total, tStep, u);
 
  // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
  if ( integrationRulesArray.size() == 1 ) {
    FloatMatrix B, D, DB, Ksigma;
    for ( auto &gp : *this->giveDefaultIntegrationRulePtr() ) {
      this->computeBmatrixAt(gp, B, tStep, u);
      this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
      double dV = this->computeVolumeAround(gp);
      DB.beProductOf(D, B);
      if ( matStiffSymmFlag ) {
	answer.plusProductSymmUpper(B, DB, dV);
      } else {
	answer.plusProductUnsym(B, DB, dV);
      }
      
      this->computeInitialStressStiffness(Ksigma, rMode, gp, tStep,B,u);
      Ksigma.times(dV);
      answer.add(Ksigma);   
    }
    
    if ( matStiffSymmFlag ) {
      answer.symmetrized();
    }
  }
}


void
Truss3dnl2 :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_dPdF_1d(rMode, gp, tStep);
}


  
void
Truss3dnl2 :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, const FloatArray &d)
{
  double L = computeLength();
  double l = computeDeformedLength(d);
  FloatMatrixF<6,6> A = this->giveAmatrix();
  FloatArray x(X);
  x.add(d); 
  answer.beTProductOf(x,A);
  answer.times(1./l/L);
}


FloatMatrixF<6,6> 
Truss3dnl2 :: giveAmatrix()
{
  FloatMatrix A;
  A = {{1, 0, 0, -1, 0, 0}, {0, 1, 0, 0, -1, 0}, {0, 0, 1, 0, 0, -1}, {-1, 0, 0, 1, 0, 0}, {0, -1, 0, 0, 1, 0}, {0, 0, -1, 0, 0, 1}};
  return A;
}

FloatMatrixF<3,6> 
Truss3dnl2 :: givePmatrix()
{
  FloatMatrix P;
  P = {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1}};
  return P;
}




  
void
Truss3dnl2 :: computeInitialStressStiffness(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep, const FloatMatrix &B, const FloatArray &d)
//
{
    
  FloatMatrix BB, A = this->giveAmatrix();
  double L = computeLength();
  double l = computeDeformedLength(d);
  FloatArray x(X);
  FloatMatrix xx, Axx, AxxA;
  x.add(d); 

  xx.beProductTOf(x,x);
  Axx.beProductOf(A,xx);
  AxxA.beProductOf(Axx,A);
  AxxA.times(1./l/l);
  
  answer = A;
  answer.subtract(AxxA);
  answer.times(1./l/L);
  auto stress = this->giveStructuralCrossSection()->giveFirstPKStresses(this->computeStrainVector(gp, d), gp, tStep);

   // prevent zero initial stress stiffness
   if ( stress.at(1) == 0 ) {
     FloatMatrix D;
     this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
     stress.at(1) = D.at(1,1) * initStressFactor;
     
   }
   answer.times(stress.at(1));
}
  


void
Truss3dnl2 :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 2);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


FEICellGeometry *
Truss3dnl2 :: giveCellGeometryWrapper()
{
    if ( !cellGeometryWrapper ) {
        cellGeometryWrapper = new FEIElementGeometryWrapper(this);
    }

    return cellGeometryWrapper;
}  

} // end namespace oofem

