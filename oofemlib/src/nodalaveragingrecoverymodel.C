/* $Header: /home/cvs/bp/oofem/oofemlib/src/nodalaveragingrecoverymodel.C,v 1.2.4.1 2004/04/05 15:19:43 bp Exp $ */
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
// file nodalaveragingrecoverymodel.C
//

#include "nodalaveragingrecoverymodel.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "elementside.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

NodalAveragingRecoveryModel :: NodalAveragingRecoveryModel (Domain* d) : NodalRecoveryModel (d)
{
}

NodalAveragingRecoveryModel :: ~NodalAveragingRecoveryModel()
{
}

int
NodalAveragingRecoveryModel :: recoverValues (InternalStateType type, TimeStep* tStep)
{
 int ireg, nregions = domain->giveNumberOfRegions ();
 int ielem, nelem = domain->giveNumberOfElements();
 int inode, nnodes = domain->giveNumberOfDofManagers();
 int elemNodes;
 int regionValSize;
 int elementNode, node;
 int regionDofMans;
 int i, neq, eq;
 Element* element;
 NodalAveragingRecoveryModelInterface* interface;
 IntArray skipRegionMap (nregions), regionRecSize(nregions);
 IntArray regionNodalNumbers (nnodes);
 IntArray regionDofMansConnectivity;
 FloatArray lhs, val;


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
  lhs.resize (neq); lhs.zero();
  regionDofMansConnectivity.resize(regionDofMans); regionDofMansConnectivity.zero();

  // assemble element contributions
  for (ielem = 1; ielem<= nelem; ielem++) {
   element = domain->giveElement (ielem);
   if (element->giveRegionNumber() != ireg) continue;
   if ((interface = (NodalAveragingRecoveryModelInterface*) 
      element->giveInterface (NodalAveragingRecoveryModelInterfaceType)) == NULL) {
    abort ();
   }

   elemNodes = element->giveNumberOfDofManagers ();
   // ask element contributions
   for (elementNode = 1; elementNode<= elemNodes; elementNode++) {
    node = element->giveDofManager(elementNode)->giveNumber();
    interface->NodalAveragingRecoveryMI_computeNodalValue (val, elementNode, type, tStep);
    eq = (regionNodalNumbers.at(node)-1)*regionValSize;
    for (i = 1; i<= regionValSize; i++) lhs.at(eq+i) += val.at(i);
    regionDofMansConnectivity.at(regionNodalNumbers.at(node)) ++;
   }
  }  // end assemble element contributions

  // solve for recovered values of active region
  for (inode=1; inode<= nnodes; inode++) {
   if (regionNodalNumbers.at(inode)) {
    eq = (regionNodalNumbers.at(inode)-1)*regionValSize;
    for (i=1; i<=regionValSize; i++) lhs.at(eq+i) /= regionDofMansConnectivity.at(regionNodalNumbers.at(inode));
   }
  }

  // update recovered values 
  this->updateRegionRecoveredValues (ireg, regionNodalNumbers, regionValSize, lhs);
  
 } // end loop over regions

 this->valType = type;
 this->stateCounter = tStep->giveSolutionStateCounter();
 return 1;
}

void
NodalAveragingRecoveryModel :: initRegionMap (IntArray& regionMap, IntArray& regionValSize, InternalStateType type)
{
 int nregions = domain->giveNumberOfRegions ();
 int ielem, nelem = domain->giveNumberOfElements();
 int i, regionsSkipped = 0;
 Element* element;
 NodalAveragingRecoveryModelInterface* interface;

 regionMap.resize (nregions); regionMap.zero();
 regionValSize.resize(nregions); regionValSize.zero();

 // loop over elements and check if implement interface
 for (ielem = 1; ielem<= nelem; ielem++) {
  element = domain->giveElement (ielem);
  if ((interface =  (NodalAveragingRecoveryModelInterface*)element->
     giveInterface (NodalAveragingRecoveryModelInterfaceType)) == NULL) {
/*
   printf ("NodalAveragingRecoveryModel :: initRegionMap: Element %d does not support required interface", ielem);
   printf ("-skipping region %d\n", element->giveRegionNumber());
*/
   regionsSkipped = 1;
   regionMap.at(element->giveRegionNumber()) = 1;
   continue;
  } else {
   if (regionValSize.at(element->giveRegionNumber())) {
    if (regionValSize.at(element->giveRegionNumber()) != 
      interface->NodalAveragingRecoveryMI_giveDofManRecordSize(type)) {
     regionMap.at(element->giveRegionNumber()) = 1;
/*
         printf ("NodalAveragingRecoveryModel :: initRegionMap: element %d has incompatible value size, skipping region\n",
                 ielem);
*/
     regionsSkipped = 1;
    }
   } else regionValSize.at(element->giveRegionNumber()) = interface->
    NodalAveragingRecoveryMI_giveDofManRecordSize(type);
  }
 }

 if (regionsSkipped) {
   OOFEM_LOG_RELEVANT ("NodalAveragingRecoveryModel :: initRegionMap: skipping regions ");
   for ( i =1; i<=nregions; i++)
     if (regionMap.at(i)) OOFEM_LOG_RELEVANT ("%d ", i);
   OOFEM_LOG_RELEVANT ("\n");
 }
}



