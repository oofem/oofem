/* $Header: /home/cvs/bp/oofem/oofemlib/src/ppdeelement.C,v 1.14.4.1 2004/04/05 15:19:43 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/


//   file PPDEELEMENT.C

#include "ppdeelement.h"
#include "domain.h"
#include "timestep.h"
#include "node.h"
#include "load.h"
#include "dof.h"
#include "material.h"
#include "heatcrosssection.h"
#include "gausspnt.h"
#include "integrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "skyline.h"
#include "debug.h"
#include "verbose.h"
#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <stdio.h>
#endif

PPdeElement :: PPdeElement (int n, Domain* aDomain)
  : Element (n, aDomain)
   // Constructor. Creates an element with number n, belonging to aDomain.
{
  rotationMatrix     = NULL ;
  rotationMatrixDefined = 0 ;
  // unknownType = TemperatureVector
}


PPdeElement :: ~PPdeElement ()
   // Destructor.
{

 delete rotationMatrix ;
}


void
PPdeElement :: computeConstitutiveMatrixAt (FloatMatrix& answer, 
                      MatResponseMode rMode, GaussPoint* gp, TimeStep* tStep)
   // Returns the  material matrix {E} of the receiver.
 // type of matrix is determined by this->giveMaterialMode()
 // rMode parameter determines type of stiffness matrix to be requested
 // (tangent, secant, ...)
{
 ((HeatCrossSection*)this->giveCrossSection())
   ->giveCharMaterialConductivityMatrix(answer, rMode,gp,tStep);
 return;
}

void
PPdeElement :: computeBcRhsVectorAt (FloatArray& answer, TimeStep* stepN, ValueModeType mode)
   // Computes the load vector due to the boundary conditions acting on the
   // receiver's nodes, at stepN. Returns NULL if this array contains only
   // zeroes.
{
 FloatArray  d, dp;
 FloatMatrix s;
/*
 this -> computeVectorOfPrescribed(FluxVector,TotalMode,stepN, d) ;
 if ((stepN->giveLoadResponseMode()==IncrementOfLoad) && (!stepN->isTheFirstStep())) {
  this -> computeVectorOfPrescribed(FluxVector,TotalMode,stepN->givePreviousStep(), dp);
  d.substract (dp);
  // delete dp;
 }
*/
 this -> computeVectorOfPrescribed(EID_ConservationEquation, mode,stepN, d) ;

 if (d.containsOnlyZeroes())
   answer.resize (0);
 else {
   this -> computeConductivityMatrix(s, TangentStiffness, stepN);
   answer.beProductOf (s, d); answer.negated() ;
 }
 
 // delete d ;
 return  ;
}


void
PPdeElement :: computeGeneratorRhsVectorAt (FloatArray& answer, Load* forLoad, TimeStep* stepN, ValueModeType mode)
// Computes numerically the generator Rhs vector of the receiver due to the generator
//  at stepN.
// // load is firrst transformed to local cs.
// // load vector is then transformed to coordinate system in each node.
// // (should be global coordinate system, but there may be defined 
// //  different coordinate system in each node)
{
  int         i ;
  double      dV ;
  GaussPoint* gp ;
  FloatArray  force, ntf ;
  FloatMatrix n, nt ;
 IntegrationRule* iRule = integrationRulesArray[giveDefaultIntegrationRule()];
  
  forLoad -> computeComponentArrayAt(force, stepN, mode);
 force.times(this->giveMaterial() -> give('d'));
  // ask for transformation from global cs to local cs.
  // R = this-> GiveGtoLRotationMatrix ();
  // if (R) { force -> rotatedWith(R,'n') ; delete R;}
  
  answer.resize (0);

  if (force.giveSize()) {

  for (i=0 ; i<iRule->getNumberOfIntegrationPoints() ; i++) {
   gp  = iRule->getIntegrationPoint(i) ;
   this -> computeNmatrixAt(gp, n) ;
   dV  = this -> computeVolumeAround(gp) ;
   nt.beTranspositionOf (n);
   ntf.beProductOf (nt, force);
   ntf.times(dV) ;
   // ntf = nt   -> Times(force) -> times(dV) ;
   answer.add(ntf) ;
  }
  //delete force ;
  // R = this -> giveRotationMatrix() ;
  // if (R) answer -> rotatedWith(R,'t'); 
  return   ;
 }
}



void 
PPdeElement :: computeEdgeRhsVectorAt (FloatArray& answer, Load*,TimeStep*, ValueModeType mode) 
{
  _error ("computeEdgeLoadVectorAt : not implemented for this element type");
}

void
PPdeElement :: computeSurfaceRhsVectorAt (FloatArray& answer, Load*, TimeStep*, ValueModeType mode)
{
  _error ("computeSurfaceLoadVectorAt : not implemented for this element type");
}


void
PPdeElement :: computeCapacityMatrix (FloatMatrix& answer, TimeStep* tStep)
   // Computes numerically the consistent (full) capacity matrix of the receiver.
{
   int         i ;
   double      density,dV ;
   FloatMatrix n ;
   GaussPoint  *gp ;
  IntegrationRule* iRule = integrationRulesArray[giveDefaultIntegrationRule()];

   answer.resize (0,0);
   density = this -> giveMaterial() -> give('d') ;
   //specHeat= this -> giveMaterial() -> give('c') ;
   for (i=0 ; i<iRule->getNumberOfIntegrationPoints() ; i++) {
      gp      = iRule->getIntegrationPoint(i) ;
      this -> computeNmatrixAt(gp, n) ;
      dV      = this -> computeVolumeAround(gp) ;
      answer.plusProduct(n, n,density*dV) ;
  }

   answer.symmetrized() ;
}


void
PPdeElement :: computeRhsVectorAt (FloatArray& answer, TimeStep* stepN, ValueModeType mode)
// Computes the load vector of the receiver, at stepN.
{
  int         i,n,nLoads ;
  bcGeomType    ltype;
  Load*       load;

  FloatArray  helpLoadVector ;
  // FloatArray* temperatureLoadVector= NULL ;
  
  answer.resize (0) ;
  
  nLoads    = this -> giveBodyLoadArray()->giveSize() ;
  for (i=1 ; i<=nLoads ; i++) {
  n     = bodyLoadArray.at(i) ;
  load  = (Load*) domain->giveLoad(n) ;
  ltype = load->giveBCGeoType();
  if (ltype == BodyLoadBGT) {
  this -> computeGeneratorRhsVectorAt(helpLoadVector, load,stepN, mode) ;
  if (helpLoadVector.giveSize()) {
    answer.add(helpLoadVector) ;
  }
 } else {
  _error ("computeRhsVectorAt : unsupported load type class");
  exit (1);
 }
 }
 
  this -> computeBcRhsVectorAt(helpLoadVector, stepN, mode) ;
  if (helpLoadVector.giveSize()) {
  answer.add(helpLoadVector) ;
 }
  
 return  ;
}



void
PPdeElement ::  computeBcLhsDueConvection (FloatMatrix& answer, TimeStep*)
//
// computes Lhs contribution if convection bc is imposed on boundary.
// = \int N^T N h ds
//
{
  answer.resize (0,0);
}


void
PPdeElement :: computeConductivityMatrix (FloatMatrix& answer, MatResponseMode rMode, 
                     TimeStep* tStep)
   // Computes numerically the stiffness matrix of the receiver.
{
   int         i;
   double      dV ;
   FloatMatrix b, db, d ;
   GaussPoint  *gp ;
  IntegrationRule* iRule = integrationRulesArray[giveDefaultIntegrationRule()];
   
   answer.resize (0,0);
   for (i=0 ; i < iRule->getNumberOfIntegrationPoints() ; i++) {
      gp = iRule-> getIntegrationPoint(i) ;
      this -> computeBmatrixAt(gp, b) ;
      //      d  = this -> giveConstitutiveMatrix() ;
      this -> computeConstitutiveMatrixAt(d, rMode, gp, tStep);
      dV = this -> computeVolumeAround(gp) ;
      db.beProductOf (d,b) ;
      answer.plusProduct(b,db,dV) ;

      //delete b ;
      //delete db ;
      //delete d ;
    }
  answer.symmetrized() ;
  if (this->updateRotationMatrix()) answer.rotatedWith(*this->rotationMatrix) ;
   return  ;
}


void
PPdeElement ::  giveCharacteristicMatrix (FloatMatrix& answer, CharType mtrx, TimeStep *tStep) 
// 
// returns characteristics matrix of receiver accordind to mtrx
//
{
  if (mtrx == ConductivityMatrix) 
  this -> computeConductivityMatrix(answer, TangentStiffness, tStep); 
  else if (mtrx == CapacityMatrix) 
  this -> computeCapacityMatrix(answer, tStep);
  else if (mtrx ==  BcLhsDueToConvection) 
  this -> computeBcLhsDueConvection (answer, tStep);
  else _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");

  return ;
}



void
PPdeElement ::  giveCharacteristicVector (FloatArray& answer, CharType mtrx, ValueModeType mode, TimeStep *tStep) 
// 
// returns characteristics vector of receiver accordind to mtrx
//
{
  if (mtrx == ElementPPDELoadVector) this -> computeRhsVectorAt (answer, tStep, mode);
  else _error("giveCharacteristicVector: Unknown Type or Mode of characteristic mtrx.");

  return ;
}


void
PPdeElement :: printOutputAt (FILE * file, TimeStep* stepN)
   // Performs end-of-step operations.
{
}


void  
PPdeElement :: updateYourself (TimeStep* stepN)
   // Updates the receiver at end of step.
{
  this->Element :: updateYourself(stepN);
  return;
}


int
PPdeElement :: checkConsistency ()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Heat" versions
//
{
  int result =1;

  if (!this->giveMaterial()->testMaterialExtension(Material_HeatCapability)) {
    _warning("checkConsistency : material without heat support");
  result =0;
  }
  if (!this->giveCrossSection()->testCrossSectionExtension(CS_HeatCapability)) {
    _warning("checkConsistency : cross-section without heat support");
  result =0;
  }
  return result;
}




int
PPdeElement :: updateRotationMatrix () 
{
/* returns a tranformation matrix between local coordinate system
 and global coordinate system, taking into account possible local 
 coordinate system in nodes.
 if no transformation necessary - returns NULL
*/
 int isT_GtoL, isT_NtoG;
 FloatMatrix T_GtoL, T_NtoG;

 if (rotationMatrixDefined) return (this->rotationMatrix != NULL);
 rotationMatrixDefined = 1;
 isT_GtoL = this->computeGtoLRotationMatrix (T_GtoL);
 isT_NtoG = this->computeGNDofRotationMatrix (T_NtoG, _toNodalCS);

#ifdef DEBUG
 if (isT_GtoL) 
  if ((!T_GtoL.isSquare ()) || 
    (T_GtoL.giveNumberOfRows () != this->computeNumberOfDofs (EID_ConservationEquation)))
   _error ("PPdeElement :: updateRotationMatrix - T_GtoL transformation matrix size mismatch");
 if (isT_NtoG) 
  if ((!T_NtoG.isSquare ()) || 
    (T_NtoG.giveNumberOfRows () != this->computeNumberOfDofs (EID_ConservationEquation)))
   _error ("PPdeElement :: updateRotationMatrix - T_NtoG transformation matrix size mismatch");
#endif

 if (isT_GtoL && T_NtoG.isNotEmpty()) rotationMatrix = T_GtoL.Times (&T_NtoG);
 else if (isT_GtoL) rotationMatrix = T_GtoL.GiveCopy();
 else if (T_NtoG.isNotEmpty()) rotationMatrix = T_NtoG.GiveCopy();
 else rotationMatrix = NULL;
 
 //delete T_GtoL;
 //delete T_GtoNTransp;
 return (this->rotationMatrix != NULL);

}

int  
PPdeElement::computeGNDofRotationMatrix (FloatMatrix& answer, DofManTrasfType mode) 
{
 int i,j,k,lastPos = 0, flag = 0;

 // test if transformation is necessary
 for (i=1; i<= numberOfDofMans; i++)
  flag += this->giveDofManager(i)->requiresTransformation ();
 if (flag == 0) {answer.beEmptyMtrx(); return 0;}

 // initialize answer
 IntArray loc;
 giveLocationArray (loc, EID_ConservationEquation);
 answer.resize (this->computeNumberOfDofs(EID_ConservationEquation),loc.giveSize());
 answer.zero();

 FloatMatrix dofManT;
 IntArray dofIDmask;
 this->giveDofManDofIDMask (1,EID_ConservationEquation,dofIDmask);
 int nr = dofIDmask.giveSize();
 int nc;
 // loop over nodes
 for (i=1; i<= numberOfDofMans; i++) {
  this->giveDofManager(i)->computeDofTransformation (dofManT, &dofIDmask, mode);
  nc = dofManT.giveNumberOfColumns();
  for (j=1; j<=nr; j++)
   for (k=1; k<= nc; k++)
    // localize node contributions 
    answer.at((i-1)*nr+j, lastPos+k) = dofManT.at(j,k);
  lastPos += nc;
 }

/*
 // loop over sides
 for (i=1; i<= numberOfSides; i++) {
  this->giveSide(i)->computeTransformation (dofManT, &dofIDmask);
  nc = dofManT.giveNumberOfColumns();
  for (j=1; j<=nr; j++)
   for (k=1; k<= nc; k++)
    // localize node contributions 
    answer.at((i-1)*nr+j, lastPos+k) = dofManT.at(j,k);
  lastPos += nc;
 }
*/
 return 1;

}


#ifdef __OOFEG
void
PPdeElement ::drawYourself (oofegGraphicContext& gc)
{
}

#endif
 
//





