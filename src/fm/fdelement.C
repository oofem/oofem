/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.C,v 1.3.4.1 2004/04/05 15:19:53 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



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

#include "fdelement.h"
#include "domain.h"
#include "timestep.h"
#include "node.h"
#include "dof.h"
#include "load.h"
#include "boundaryload.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "debug.h"
#include "verbose.h"
#include "cltypes.h"
#include "elementside.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <stdio.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif

FluidDynamicElement :: FluidDynamicElement (int n, Domain* aDomain, ElementMode em)
: Element (n, aDomain)
   // Constructor. Creates an element with number n, belonging to aDomain.
{
  emode = em;
}


FluidDynamicElement :: ~FluidDynamicElement ()
   // Destructor.
{}

void
FluidDynamicElement ::  giveCharacteristicMatrix (FloatMatrix& answer, 
                                                  CharType mtrx, TimeStep *tStep) 
// 
// returns characteristics matrix of receiver accordind to mtrx
//
{
  if (mtrx == MassMatrix) 
    this -> computeMassMatrix(answer, tStep); 
  else if (mtrx == DiffusionMatrix) 
    this -> computeDiffusionMatrix(answer, tStep);
  else if (mtrx == AdvectionMatrix)
    this -> computeAdvectionMatrix (answer, tStep);
  else if (mtrx == CouplingUPMatrix)
    this -> computeCouplingUPMatrix (answer, tStep);

  else _error2("giveCharacteristicMatrix: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx));
  
  return ;
}

void
FluidDynamicElement ::  giveCharacteristicVector (FloatArray& answer, CharType mtrx, ValueModeType mode,
                                                  TimeStep *tStep) 
// 
// returns characteristics vector of receiver according to requested type
//
{
  if (mtrx == ElementBCTransportVector) this->computeBCVectorAt (answer, tStep, mode); 
  else _error2("giveCharacteristicVector: Unknown Type of characteristic mtrx (%s)",
	       __CharTypeToString(mtrx));
  
  return ;
}

int
FluidDynamicElement :: checkConsistency ()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
  int result =1;
  /*
  if (!this->giveMaterial()->testMaterialExtension(Material_TransportCapability)) {
    _warning("checkConsistency : material without support for transport problems");
    result =0;
  }
  */
  return result;
}

void
FluidDynamicElement :: printOutputAt (FILE * file, TimeStep* stepN)
  // Performs end-of-step operations.
{
}

/*
void
FluidDynamicElement::computeMassSubMatrix (FloatMatrix& answer, int ip, int iri, TimeStep* tStep)
{
  int         i ;
  double      dV;
  FloatMatrix n ;
  GaussPoint  *gp ;
  IntegrationRule* iRule = integrationRulesArray[iri];
  
  answer.resize (0,0);
  answer.zero();
  for (i=0 ; i<iRule->getNumberOfIntegrationPoints() ; i++) {
    gp      = iRule->getIntegrationPoint(i) ;
    this -> computeNSubMatrixAt(n, ip, gp->giveCoordinates()) ;
    dV      = this -> computeVolumeAround(gp) ;
    answer.plusProduct(n, n, dV * c) ;
  }
  
  answer.symmetrized() ;
}

void
FluidDynamicElement::computeAdvectionSubMatrix (FloatMatrix& answer, int nsd, int iri, TimeStep* tStep)
{
  int         ip ;
  double      dV;
  FloatMatrix n ;
  GaussPoint  *gp ;
  IntegrationRule* iRule = integrationRulesArray[iri];
  int vsize = this->giveNumberOfVelocityUnknowns();


  answer.resize (0,0);
  answer.zero();
  for (ip=0 ; ip<iRule->getNumberOfIntegrationPoints() ; ip++) {
    gp      = iRule->getIntegrationPoint(ip) ;
    this -> computeNSubMatrixAt(n, 1, gp->giveCoordinates()) ; // 1 - velocity
    this -> computeGradientMatrixAt (b, gp->giveCoordinates()) ;
    this -> computeVelocityAt (v, n);
    for (i=1; i<=vsize; i++) {
      for (j=1; j<=vsize; j++) {
        for (k=1; k<=nsd; k++) {
          answer.at(i,j) += n.at(1,i)*v.at(k)*b.at(1,k)*dV;
        }
      }
    }
  }
}




TransportElement :: updateInternalState (TimeStep* stepN)
   // Updates the receiver at end of step.
{
  int i,j ;
  IntegrationRule* iRule;
  FloatArray f,r;
  FloatMatrix n;
  TransportMaterial* mat = ((TransportMaterial*)this->giveMaterial());
  GaussPoint* gp;
  
  // force updating ip values
  for (i=0 ; i < numberOfIntegrationRules ; i++) {
    iRule = integrationRulesArray[i];
    for (j=0; j < iRule->getNumberOfIntegrationPoints(); j++) {
      gp = iRule->getIntegrationPoint(j);
      this -> computeNmatrixAt(n, gp->giveCoordinates()) ;
      this -> computeVectorOf (FluxVector, VM_Total, stepN, r);
      f.beProductOf(n,r);
      mat->updateInternalState (f, gp, stepN);
    }
  }
}
*/


#ifdef __OOFEG
int 
FluidDynamicElement::giveInternalStateAtNode (FloatArray& answer, InternalStateType type, InternalStateMode mode, 
                                           int node, TimeStep* atTime)
{
  int indx = 1;
  Node* n = this->giveNode(node);

  if (type == IST_Velocity) {
    answer.resize(this->giveSpatialDimension());
    int dofindx;
    if ((dofindx = n->findDofWithDofId (V_u))) answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
    if ((dofindx = n->findDofWithDofId (V_v))) answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
    if ((dofindx = n->findDofWithDofId (V_w))) answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
    return 1;
  } else if (type == IST_Pressure) {
    int dofindx;
    if ((dofindx = n->findDofWithDofId (P_f))) {
      answer.resize(1);
      answer.at(1) = n->giveDof(dofindx)->giveUnknown(EID_ConservationEquation, VM_Total, atTime);
      return 1;
    } else return 0;
  } else return Element::giveInternalStateAtNode (answer, type, mode, node, atTime);
}

#endif

int 
FluidDynamicElement::giveIntVarCompFullIndx (IntArray& answer, InternalStateType type)
{
  if ((type == IST_Velocity) ) {
    IntArray mask;
    int indx = 1;
    answer.resize(3);
    this->giveElementDofIDMask(mask);
    if (mask.findFirstIndexOf (V_u)) answer.at(1) = indx++;
    if (mask.findFirstIndexOf (V_v)) answer.at(2) = indx++;
    if (mask.findFirstIndexOf (V_w)) answer.at(3) = indx++;
    return 1;
  } else if (type == IST_Pressure) {
    answer.resize (1);
    answer.at(1) = 1;
    return 1;
  } else {
    return Element::giveIntVarCompFullIndx (answer, type);
  }
}
