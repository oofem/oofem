/* $Header: /home/cvs/bp/oofem/sm/src/mmashapefunctprojection.C,v 1.7.4.1 2004/04/05 15:19:47 bp Exp $ */
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

#include "mmashapefunctprojection.h"
#include "mathfem.h"
#include "gausspnt.h"
#include "element.h"
#include "node.h"
#include "domain.h"
#include "spatiallocalizer.h"
#include "timestep.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"

MMAShapeFunctProjection::MMAShapeFunctProjection () : MaterialMappingAlgorithm ()
{
 stateCounter = 0;
}

MMAShapeFunctProjection::~MMAShapeFunctProjection() {}

void 
MMAShapeFunctProjection::__init (Domain* dold, Domain* dnew, IntArray& varTypes, FloatArray& coords, int region, TimeStep* tStep)
//(Domain* dold, IntArray& varTypes, GaussPoint* gp, TimeStep* tStep) 
{
 MMAShapeFunctProjectionInterface* interface;
 int ivar, nvar = varTypes.giveSize();
 // check time stemp
 if (stateCounter == tStep->giveSolutionStateCounter()) return;


 // Project Gauss point components to nodes on old mesh
 if (dold->giveSmoother()==NULL) {
  //dold->setSmoother (new ZZNodalRecoveryModel(dold));
  OOFEM_LOG_INFO("MMAShapeFunctProjection: setting NodalAveragingRecoveryModel\n");
  dold->setSmoother (new NodalAveragingRecoveryModel(dold));
 }

 int inode, nnodes = dnew -> giveNumberOfDofManagers();
 int j, nbelemnodes;
 Element *belem;
 MMAShapeFunctProjectionInterface::nodalValContainerType container;
 const FloatArray *nodVal;
 
 this->intVarTypes = varTypes;
 nodalValList.clear();
 nodalValList.growTo (nvar);

 for (ivar = 1; ivar <= nvar; ivar++) {

  nodalValList.put (ivar, new AList<FloatArray> (nnodes));
  for (inode = 1; inode <= nnodes; inode++) nodalValList.at(ivar)->put (inode, new FloatArray);
//  nodalValList.at(ivar)->clear();
//  nodalValList.at(ivar)->growTo (nnodes);
  
  dold -> giveSmoother() -> recoverValues ((InternalStateType)varTypes.at(ivar), tStep);
  
  // Transfer nodal components from the old to a new mesh
  
  for (inode = 1; inode <= nnodes; inode++) {
   // find background element containing new node 
   belem = dold->giveSpatialLocalizer()->giveElementContainingPoint(*dnew->giveNode(inode)->giveCoordinates());
   if (belem == NULL) {
     OOFEM_ERROR2 ("MMAShapeFunctProjection::init node %d not in old mesh",inode);
   }
   // map internal variable to new node
   if ((interface = (MMAShapeFunctProjectionInterface*) 
      belem->giveInterface (MMAShapeFunctProjectionInterfaceType)) == NULL) {
    abort ();
   }
   // set up  belement nodal values
   nbelemnodes = belem->giveNumberOfDofManagers();
   if (container.giveSize() < nbelemnodes) {
    int oldSize = container.giveSize();
    container.growTo(nbelemnodes);
    for (int i = oldSize+1; i<=nbelemnodes; i++) 
     container.put(i, new FloatArray);
   }
   
   for (j=1; j<= nbelemnodes; j++) {
    dold -> giveSmoother() -> giveNodalVector(nodVal, belem->giveDofManager(j)->giveNumber(), 
                                              belem->giveRegionNumber());
    *(container.at(j)) = *nodVal;
   }
   
   interface -> MMAShapeFunctProjectionInterface_interpolateIntVarAt (*nodalValList.at(ivar)->at(inode), 
                                     *dnew->giveNode(inode)->giveCoordinates(),
                                     MMAShapeFunctProjectionInterface::coordType_global,
                                     container, (InternalStateType)varTypes.at(ivar),
                                     tStep);
  }
 }  
 // remember time stemp
 stateCounter = tStep->giveSolutionStateCounter();
}


void
MMAShapeFunctProjection::finish (TimeStep* tStep)
{
 // delete nodalValList
 nodalValList.clear();
}

int
MMAShapeFunctProjection::mapVariable (FloatArray& answer, GaussPoint* gp, InternalStateType type, TimeStep* tStep)
{
 Element *elem = gp->giveElement();
 int inode, nnodes = elem->giveNumberOfDofManagers();
 MMAShapeFunctProjectionInterface::nodalValContainerType container(nnodes);
 MMAShapeFunctProjectionInterface* interface;
 
 if ((interface = (MMAShapeFunctProjectionInterface*) 
    elem->giveInterface (MMAShapeFunctProjectionInterfaceType)) == NULL) {
  abort ();
 }

 int indx = this->intVarTypes.findFirstIndexOf((int)type);
 if (indx) {

  for (inode = 1; inode <= nnodes; inode ++) {
   container.put(inode, new FloatArray);
   *(container.at(inode)) = *(this->nodalValList.at(indx)->at(elem->giveDofManager(inode)->giveNumber()));
  }
  
  interface -> MMAShapeFunctProjectionInterface_interpolateIntVarAt (answer, (*gp->giveCoordinates()),
                                    MMAShapeFunctProjectionInterface::coordType_local,
                                    container, type, tStep);
 } else {
   OOFEM_ERROR ("MMAShapeFunctProjection::mapVariable: var not initialized");
 }
 return 1;
}


int
MMAShapeFunctProjection::__mapVariable (FloatArray& answer, FloatArray& coords, Domain* dnew, 
                                        InternalStateType type, TimeStep* tStep)
{
 Element* elem = dnew -> giveSpatialLocalizer() -> giveElementContainingPoint (coords);
 if (!elem)  {
   OOFEM_ERROR ("MMAShapeFunctProjection::__mapVariable: no suitable source found");
 }

 int inode, nnodes = elem->giveNumberOfDofManagers();
 MMAShapeFunctProjectionInterface::nodalValContainerType container(nnodes);
 MMAShapeFunctProjectionInterface* interface;
 
 if ((interface = (MMAShapeFunctProjectionInterface*) 
      elem->giveInterface (MMAShapeFunctProjectionInterfaceType)) == NULL) {
   abort ();
 }

 int indx = this->intVarTypes.findFirstIndexOf((int)type);
 if (indx) {

  for (inode = 1; inode <= nnodes; inode ++) {
   container.put(inode, new FloatArray);
   *(container.at(inode)) = *(this->nodalValList.at(indx)->at(elem->giveDofManager(inode)->giveNumber()));
  }
  
  interface -> MMAShapeFunctProjectionInterface_interpolateIntVarAt (answer, coords,
                                    MMAShapeFunctProjectionInterface::coordType_global,
                                    container, type, tStep);
 } else {
   OOFEM_ERROR ("MMAShapeFunctProjection::__mapVariable: var not initialized");
 }
 return 1;
}







