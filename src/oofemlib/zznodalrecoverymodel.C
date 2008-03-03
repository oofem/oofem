/* $Header: /home/cvs/bp/oofem/oofemlib/src/zznodalrecoverymodel.C,v 1.4.4.1 2004/04/05 15:19:44 bp Exp $ */
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

//
// file zznodalrecoverymodel.C
//

#include "zznodalrecoverymodel.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "elementside.h"
#include "integrationrule.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

ZZNodalRecoveryModel :: ZZNodalRecoveryModel (Domain* d) : NodalRecoveryModel (d)
{
}

ZZNodalRecoveryModel :: ~ZZNodalRecoveryModel()
{
}

/***********
int
ZZNodalRecoveryModel :: recoverValues (InternalStateType type, TimeStep* tStep)
{
 int ireg, nregions = this->giveNumberOfRegions ();
 int ielem, nelem = domain->giveNumberOfElements();
 int nnodes = domain->giveNumberOfDofManagers();
 int elemNodes, elemSides;
 int elementValSize, regionValSize;
 int elementNode, elementSide, node, side;
 int regionDofMans;
 int i, neq, eq;
 Element* element;
 ZZNodalRecoveryModelInterface* interface;
 IntArray skipRegionMap (nregions), regionRecSize(nregions);
 IntArray regionNodalNumbers (nnodes);
 IntArray loc;
 FloatArray lhs, rhs, nn, nsig;


 if ((this->valType == type) && (this->stateCounter == tStep->giveSolutionStateCounter())) return 1;

 // clear nodal table
 this->clear();

 // init region table indicating regions to skip
 this->initRegionMap (skipRegionMap, regionRecSize, type);

 // loop over regions 
 for (ireg=1; ireg<= nregions; ireg++) {
  // skip regions 
  if (skipRegionMap.at(ireg)) break;

  // loop over elements and determine local region node numbering and determine and check nodal values size
  if (this->initRegionNodeNumbering (regionNodalNumbers, regionDofMans, ireg )==0) break;

  regionValSize = regionRecSize.at(ireg);
  neq = regionDofMans*regionValSize;
  lhs.resize (neq); lhs.zero();
  rhs.resize (neq); rhs.zero();
  // assemble element contributions
  for (ielem = 1; ielem<= nelem; ielem++) {
   element = domain->giveElement (ielem);
   if ((interface = (ZZNodalRecoveryModelInterface*) element->giveInterface (ZZNodalRecoveryModelInterfaceType)) == NULL) {
    abort ();
   }
   if (giveElementRegion(element) != ireg) continue;

   // ask element contributions
   interface->ZZNodalRecoveryMI_computeNNMatrix (nn, type);
   interface->ZZNodalRecoveryMI_computeNValProduct(nsig, type, tStep);
   // assemble contributions
   elemNodes = element->giveNumberOfNodes ();
   elemSides = element->giveNumberOfSides();

   loc.resize ((elemNodes+elemSides)*regionValSize);
   eq = 1;
   for (elementNode = 1; elementNode<= elemNodes; elementNode++) {
    node = element->giveNode(elementNode)->giveNumber();
    for (i = 1; i<= regionValSize; i++) loc.at(eq++) = (regionNodalNumbers.at(node)-1)*regionValSize + i;
   }
   for (elementSide = 1; elementSide<= elemSides; elementSide++) {
    node = element->giveSide(elementSide)->giveNumber();
    for (i = 1; i<= regionValSize; i++) loc.at(eq++) = (regionNodalNumbers.at(node)-1)*regionValSize + i;
   }
   lhs.assemble (nn, loc);
   rhs.assemble (nsig, loc);
  }  // end assemble element contributions

  // solve for recovered values of active region
  for (i=1; i<= neq; i++) 
   // rhs will be overriden by recovered values
   rhs.at(i) /= lhs.at(i);


  // update recovered values 
  this->updateRegionRecoveredValues (ireg, regionNodalNumbers, regionValSize, rhs);
  
 } // end loop over regions

 this->valType = type;
 this->stateCounter = tStep->giveSolutionStateCounter();
 return 1;
}
***********/

int
ZZNodalRecoveryModel :: recoverValues (InternalStateType type, TimeStep* tStep)
{
 int ireg, nregions = domain->giveNumberOfRegions ();
 int ielem, nelem = domain->giveNumberOfElements();
 int nnodes = domain->giveNumberOfDofManagers();
 int elemNodes;
 int regionValSize;
 int elementNode, node;
 int regionDofMans;
 int i, j, neq, eq;
 Element* element;
 ZZNodalRecoveryModelInterface* interface;
 IntArray skipRegionMap (nregions), regionRecSize(nregions);
 IntArray regionNodalNumbers (nnodes);
// IntArray loc;
 FloatArray lhs, nn, sol;
 FloatMatrix rhs, nsig;


 if ((this->valType == type) && (this->stateCounter == tStep->giveSolutionStateCounter())) return 1;

 // clear nodal table
 this->clear();

 // init region table indicating regions to skip
 this->initRegionMap (skipRegionMap, regionRecSize, type);

 // loop over regions 
 for (ireg=1; ireg<= nregions; ireg++) {
  // skip regions 
  if (skipRegionMap.at(ireg)) continue;

  // loop over elements and determine local region node numbering and determine and check nodal values size
  if (this->initRegionNodeNumbering (regionNodalNumbers, regionDofMans, ireg )==0) break;

  regionValSize = regionRecSize.at(ireg);
  neq = regionDofMans*regionValSize;
  lhs.resize (regionDofMans); lhs.zero();
  rhs.resize (regionDofMans,regionValSize); rhs.zero();
  sol.resize (neq); sol.zero();
  // assemble element contributions
  for (ielem = 1; ielem<= nelem; ielem++) {
   element = domain->giveElement (ielem);
   if (element->giveRegionNumber() != ireg) continue;
   if ((interface = (ZZNodalRecoveryModelInterface*) element->giveInterface (ZZNodalRecoveryModelInterfaceType)) == NULL) {
    abort ();
   }


   // ask element contributions
   interface->ZZNodalRecoveryMI_computeNNMatrix (nn, type);
   interface->ZZNodalRecoveryMI_computeNValProduct(nsig, type, tStep);
   // assemble contributions
   elemNodes = element->giveNumberOfDofManagers ();

   //loc.resize ((elemNodes+elemSides)*regionValSize);
   eq = 1;
   for (elementNode = 1; elementNode<= elemNodes; elementNode++) {
    node = element->giveDofManager(elementNode)->giveNumber();
    lhs.at(regionNodalNumbers.at(node)) += nn.at(eq);
    for (i = 1; i<= regionValSize; i++) 
     rhs.at(regionNodalNumbers.at(node),i) += nsig.at(eq,i);
    eq++;
   }
  }  // end assemble element contributions

  // solve for recovered values of active region
  for (i=1; i<= regionDofMans; i++) { 
   eq = (i-1)*regionValSize;
   for (j=1; j<=regionValSize; j++) {
    // rhs will be overriden by recovered values
    sol.at(eq+j) = rhs.at(i,j)/lhs.at(i);
   }
  }  
  
  // update recovered values 
  this->updateRegionRecoveredValues (ireg, regionNodalNumbers, regionValSize, sol);
  
 } // end loop over regions

 this->valType = type;
 this->stateCounter = tStep->giveSolutionStateCounter();
 return 1;
}


void
ZZNodalRecoveryModel :: initRegionMap (IntArray& regionMap, IntArray& regionValSize, InternalStateType type)
{
 int nregions = domain->giveNumberOfRegions ();
 int ielem, nelem = domain->giveNumberOfElements();
 int i, regionsSkipped = 0;
 Element* element;
 ZZNodalRecoveryModelInterface* interface;

 regionMap.resize (nregions); regionMap.zero();
 regionValSize.resize(nregions); regionValSize.zero();

 // loop over elements and check if implement interface
 for (ielem = 1; ielem<= nelem; ielem++) {
  element = domain->giveElement (ielem);
  if ((interface =  (ZZNodalRecoveryModelInterface*)element->giveInterface (ZZNodalRecoveryModelInterfaceType)) == NULL) {
/*
     printf ("NodalRecoveryModel :: initRegionMap: Element %d does not support required interface", ielem);
   printf ("-skipping region %d\n", element->giveRegionNumber());
*/
   regionsSkipped = 1;
   regionMap.at(element->giveRegionNumber()) = 1;
   continue;
  } else {
   if (regionValSize.at(element->giveRegionNumber())) {
    if (regionValSize.at(element->giveRegionNumber()) != interface->ZZNodalRecoveryMI_giveDofManRecordSize(type)) {
     regionMap.at(element->giveRegionNumber()) = 1;
/*
     printf ("NodalRecoveryModel :: initRegionMap: element %d has incompatible value size, skipping region\n",ielem);
*/
     regionsSkipped = 1;
    }
   } else {
    regionValSize.at(element->giveRegionNumber()) = interface->ZZNodalRecoveryMI_giveDofManRecordSize(type);
    if (regionValSize.at(element->giveRegionNumber()) == 0) {
     regionMap.at(element->giveRegionNumber()) = 1;
     regionsSkipped = 1;
    }
   } 
  }
 }

 if (regionsSkipped) {
   OOFEM_LOG_RELEVANT ("NodalRecoveryModel :: initRegionMap: skipping regions ");
   for ( i =1; i<=nregions; i++)
     if (regionMap.at(i)) OOFEM_LOG_RELEVANT ("%d ", i);
   OOFEM_LOG_RELEVANT ("\n");
 }

 
}






/***************
void
ZZNodalRecoveryModelInterface::ZZNodalRecoveryMI_computeNValProduct (FloatArray& answer, InternalStateType type, 
                                   TimeStep* tStep)
{  // evaluates N^T sigma over element volume
  // N(nsigma, nsigma*nnodes)
  // Definition : sigmaVector = N * nodalSigmaVector
  int i, size;
  double dV;
  FloatArray stressVector, help;
  FloatMatrix n;
 Element* elem  = this->ZZNodalRecoveryMI_giveElement();
  IntegrationRule* iRule = elem->giveDefaultIntegrationRulePtr();
  GaussPoint* gp;

 size = elem->giveNumberOfNodes()*ZZNodalRecoveryMI_giveDofManRecordSize (type);
 if (type == StressVector) answer.resize(size); 
 else return;

 answer.zero();
  for (i=0 ; i < iRule->getNumberOfIntegrationPoints() ; i++) {
  gp  = iRule->getIntegrationPoint(i);
  dV  = elem -> computeVolumeAround(gp) ;
  //this-> computeStressVector(stressVector, gp, stepN);
  elem -> giveIPValue(stressVector, gp, type);

  this-> ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx (n, gp, type);
  help.beTProductOf(n,stressVector);
  answer.add(help.times(dV));
  }
 
  return ;
}



void 
ZZNodalRecoveryModelInterface::ZZNodalRecoveryMI_computeNNMatrix (FloatArray& answer, InternalStateType type)
{
  // 
  // Returns NTN matrix (lumped) for Zienkiewicz-Zhu 
  // The size of N mtrx is (nstresses, nnodes*nstreses)
  // Definition : sigmaVector = N * nodalSigmaVector
  //
  int i,j, size;
  double sum, dV;
  FloatMatrix n, fullAnswer;
 Element* elem  = this->ZZNodalRecoveryMI_giveElement();
  IntegrationRule* iRule = elem->giveDefaultIntegrationRulePtr();
  GaussPoint* gp;

 size = elem->giveNumberOfNodes()*ZZNodalRecoveryMI_giveDofManRecordSize (type);
 fullAnswer.resize(size, size);

  for (i=0 ; i < iRule->getNumberOfIntegrationPoints() ; i++) {
  gp  = iRule->getIntegrationPoint(i);
  dV  = elem -> computeVolumeAround(gp) ;
  this->ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx (n, gp, type);
  fullAnswer.plusProductSymmUpper(n,n,dV);
  }

  fullAnswer.symmetrized();
 answer.resize (size);
 for (i=1; i<=size ; i++) {
  sum = 0.0;
  for (j=1; j<=size; j++)
   sum += fullAnswer.at(i,j);
  answer.at(i) = sum;
 }
  
  return;
}
************/


void
ZZNodalRecoveryModelInterface::ZZNodalRecoveryMI_computeNValProduct (FloatMatrix& answer, InternalStateType type, 
                                   TimeStep* tStep)
{  // evaluates N^T sigma over element volume
  // N(nsigma, nsigma*nnodes)
  // Definition : sigmaVector = N * nodalSigmaVector
  int i, j,k,size;
  double dV;
  FloatArray stressVector, help;
  FloatMatrix n;
 Element* elem  = this->ZZNodalRecoveryMI_giveElement();
  IntegrationRule* iRule = elem->giveDefaultIntegrationRulePtr();
  GaussPoint* gp;

 size = ZZNodalRecoveryMI_giveDofManRecordSize (type);
 answer.resize(elem->giveNumberOfDofManagers(),size); 
 
 answer.zero();
  for (i=0 ; i < iRule->getNumberOfIntegrationPoints() ; i++) {
  gp  = iRule->getIntegrationPoint(i);
  dV  = elem -> computeVolumeAround(gp) ;
  //this-> computeStressVector(stressVector, gp, stepN);
  elem -> giveIPValue(stressVector, gp, type, tStep);

  this-> ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx (n, gp, type);
  for (j=1; j<=elem->giveNumberOfDofManagers(); j++)
   for (k=1; k<=size; k++) {
    answer.at(j,k)+= n.at(1,j)*stressVector.at(k)*dV;
   }
//  help.beTProductOf(n,stressVector);
//  answer.add(help.times(dV));
 }
 
  return ;
}

void 
ZZNodalRecoveryModelInterface::ZZNodalRecoveryMI_computeNNMatrix (FloatArray& answer, InternalStateType type)
{
  // 
  // Returns NTN matrix (lumped) for Zienkiewicz-Zhu 
  // The size of N mtrx is (nstresses, nnodes*nstreses)
  // Definition : sigmaVector = N * nodalSigmaVector
  //
  int i,size;
  double sum, dV, volume=0.0;
  FloatMatrix n, fullAnswer;
 Element* elem  = this->ZZNodalRecoveryMI_giveElement();
  IntegrationRule* iRule = elem->giveDefaultIntegrationRulePtr();
  GaussPoint* gp;

 size = elem->giveNumberOfDofManagers(); //ZZNodalRecoveryMI_giveDofManRecordSize (type);
 fullAnswer.resize(size, size);
 fullAnswer.zero();
 double pok = 0.0;

  for (i=0 ; i < iRule->getNumberOfIntegrationPoints() ; i++) {
  gp  = iRule->getIntegrationPoint(i);
  dV  = elem -> computeVolumeAround(gp) ;
  this->ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx (n, gp, type);
  fullAnswer.plusProductSymmUpper(n,n,dV);
  pok += (n.at(1,1)*dV);
  volume += dV;
  }


  fullAnswer.symmetrized();
  answer.resize (size);
  for (i=1; i<=size ; i++) {
    sum = 0.0;
    for (int j=1; j<=size; j++)
      sum += fullAnswer.at(i,j);
  answer.at(i) = sum;
  }
 //fullAnswer.printYourself();
 //printf ("\nVolume=%lf\n", volume);
/*
 answer.resize (size);
 sum = 0.0;
 for (i=1; i<=size ; i++) {
  sum+=fullAnswer.at(i,i);
  answer.at(i) = fullAnswer.at(i,i);
 }
 answer.times(volume/sum);
 //answer.printYourself();
*/

  return;
}














