/* $Header: /home/cvs/bp/oofem/oofemlib/src/rigidarmnode.C,v 1.15.4.1 2004/04/05 15:19:43 bp Exp $*/
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

#include "rigidarmnode.h"
#include "rigidarmslavedof.h"
#include "nodload.h"
#include "timestep.h"
#include "masterdof.h"

#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "debug.h"
#include "verbose.h"
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif

RigidArmNode :: RigidArmNode (int n, Domain* aDomain)
      : Node (n,aDomain), slaveDofMask()
   // Constructor. Creates a node with number n, belonging to aDomain.
{
 masterDofMngr = 0;
}

RigidArmNode :: ~RigidArmNode()
   // Destructor.
{
}

IRResultType
RigidArmNode :: initializeFrom (InputRecord* ir)
   // Gets from the source line from the data file all the data of the receiver.
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 int j ;
 int hasIc=0,hasBc=0,dofIc=0,dofBc=0;
 IntArray dofIDArry;
 IntArray bc, ic;

 IR_GIVE_FIELD (ir, coordinates, IFT_RigidArmNode_coords, "coords"); // Macro
 
//Node :: instanciateFromString (initString);
 IR_GIVE_FIELD (ir, masterDofMngr, IFT_RigidArmNode_master, "master"); // Macro

 
 loadArray.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, loadArray, IFT_RigidArmNode_load, "load"); // Macro
 if (this->resolveDofIDArray (ir, dofIDArry) != IRRT_OK) IR_IOERR (giveClassName(), __proc, IFT_Unknown, "", ir, result);

 // read bc and ic for primary dofs, for slave dof ignored
 bc.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, bc, IFT_RigidArmNode_bc, "bc"); // Macro

 ic.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, ic, IFT_RigidArmNode_ic, "ic"); // Macro

 // master mask input array stored in slaveDofMask
 // the zero value indicates corresponding dof to be master
 // nonzero indicates corresponding dof to be slave
 slaveDofMask.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, slaveDofMask, IFT_RigidArmNode_mastermask, "mastermask"); // Macro
 if (slaveDofMask.giveSize()) {
  if (slaveDofMask.giveSize() != this->giveNumberOfDofs()) _error ("instanciateFrom: slaveDofMask size mismatch");
  hasIc = !(ic.giveSize() == 0);
  hasBc = !(bc.giveSize() == 0);
  // check sizes
  if (hasBc) if (bc.giveSize() != this->giveNumberOfDofs()) _error ("instanciateFrom: bc size mismatch");
  if (hasIc) if (ic.giveSize() != this->giveNumberOfDofs()) _error ("instanciateFrom: ic size mismatch");
 } else {
  // default all dofs are slaves
  slaveDofMask.resize (this->giveNumberOfDofs());
  for (j=1; j<=this->giveNumberOfDofs(); j++) slaveDofMask.at(j) = 1;
 }

 dofArray = new Dof* [this->giveNumberOfDofs()] ;
 for (j=0 ; j<numberOfDofs ; j++) {
  if (slaveDofMask.at(j+1)) 
   dofArray[j] = new RigidArmSlaveDof(j+1,this,(DofID) dofIDArry.at(j+1)) ;
  else {
   if(hasIc) dofIc = ic.at(j+1) ;
   if(hasBc) dofBc = bc.at(j+1) ;
   dofArray[j] = new MasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
  }
 }
 

// check for own lcs
 if (ir->hasField (IFT_RigidArmNode_lcs, "lcs")) {
  _warning2 ("initializeFrom: lcs on RigidArmNode is not supported", 1);
 }

#ifdef __PARALLEL_MODE
 /* 
    The rigid arm node could not refer to master on other domain.
    So if the partition boundary is between rigid arm and master,
    the master must be duplicated and marked as shared or remote, or whatever to
    preserve consistency.

    In the cuurent implementation, we allow to read globnum, parttion list and parallel mode,
    but we only check, if master exists and that have the same parallel mode as receiver.
    The globnum and partition list are needed on master side only.
 */

 // read globnum, even if it is not needed 
 globalNumber = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, globalNumber, IFT_RigidArmNode_globnum, "globnum"); // Macro

 partitions.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, partitions, IFT_RigidArmNode_partitions, "partitions"); // Macro

 if (ir->hasField (IFT_RigidArmNode_shared, "shared")) parallel_mode = DofManager_shared;
 else if (ir->hasField (IFT_RigidArmNode_remote, "remote")) parallel_mode = DofManager_remote;
 else if (ir->hasField (IFT_RigidArmNode_null, "null")) parallel_mode = DofManager_null;
 else parallel_mode = DofManager_local;

 // in parallel mode,  slaves are allowed, because ((Dr. Rypl promissed)
 // masters have to be in same partition as slaves. They can be again Remote copies.
#endif
 
 return IRRT_OK;
}

int
RigidArmNode::checkConsistency () {
/*
  Checks internal data consistency in node. 
  Current implementation checks (when receiver has slave dofs) if receiver has the same 
  coordinate system as master dofManager of slave dof.
*/
 int result = 1;
 int i;

 result = result && Node::checkConsistency();

 // check if master is RigidArmNode - not supported
 if (giveMasterDofMngr()->giveClassID() == RigidArmNodeClass) {
  _warning2 ("checkConsistency: chaining of RigidArmNodes is not allowed", 1);
  result = 0;
 }

 // check for consistency of receiver and master dofs
 if (this->numberOfDofs != giveMasterDofMngr()->giveNumberOfDofs ()) {
  _warning2 ("checkConsistency: numberOfDofs on RigidArmNode and master differ", 1);
  result = 0;
 }
  
 for (i=1; i<= numberOfDofs; i++)
  if (this->giveDof (i)->giveDofID() != giveMasterDofMngr()->giveDof (i)->giveDofID()) {
   _warning2 ("checkConsistency: dofID mismatch on RigidArmNode and master", 1);
   result = 0;
   break;
 }


#ifdef __PARALLEL_MODE
 // check if master in same mode
 if (parallel_mode != DofManager_local) {
   if (giveMasterDofMngr()->giveParallelMode() != parallel_mode) {
     _warning2 ("checkConsistency: mismatch in parallel mode of RigidArmNode and master", 1);
     result = 0;
   }
 }
#endif

 return result;
}

void
RigidArmNode::computeDofTransformation (FloatMatrix& answer, const IntArray* dofIDArry, DofManTrasfType mode)
{
 // computes trasformation matrix of receiver.
 // transformation should include trasformation from global cs to nodal cs,
 // as well as further necessary transformations (for example in case 
 // rigid arms this must include transformation to master dofs).
 FloatMatrix GNTransf, masterTransf, dofTransf;
 int i, j, indx, ndof;
 RigidArmSlaveDof::RigidArmSlaveDofTransfType ttype;
 Node* master = giveMasterDofMngr();
/*
 // compute dofIDmap in advance
 IntArray dofIDmap (ndof);
 for (i=1; i<= ndof; i++) {
  indx = this->findDofWithDofId (dofIDArry.at(i));
  if (indx) dofIDmap.at(i) = indx; else {
   char buff [80];
   sprintf(buff,"computeTransformation: DofManager %d : uncompatible dof requested",this->giveNumber());
   _error (buff);
  }
 }
*/

/***************** 99
 if (mode == _toGlobalCS) ttype = RigidArmSlaveDof::_toSlave;
 else if (mode == _toNodalCS) ttype = RigidArmSlaveDof::_toMaster;
 else _error ("computeDofTransformation: unknown mode");

 if (dofIDArry) ndof = dofIDArry->giveSize(); else ndof = numberOfDofs;

 // assemble transformation contributions from local dofs
 masterTransf.resize (ndof, numberOfDofs); // master and slave should have the same dofs
 for (i=1; i<= ndof; i++) {
  //uncomment if dofIDmap assembled
  //this->giveDof (dofIDmap.at(i))->computeTransformation (dofTransf);
  if (slaveDofMask->at(i)) { // linked (slave) dof 
   if (dofIDArry == NULL) {
    ((RigidArmSlaveDof*)(this->giveDof (i)))->computeDofTransformation (dofTransf, ttype);
    for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
   } else if (indx = this->findDofWithDofId (dofIDArry->at(i))) {
    ((RigidArmSlaveDof*)(this->giveDof (indx)))->computeDofTransformation (dofTransf, ttype);
    for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
   } else {
    char buff [80];
    sprintf(buff,"computeDofTransformation: DofManager %d : uncompatible dof requested",this->giveNumber());
    _error (buff);
   }
  } else { // primary dof
   masterTransf.at(i,i) = 1.0;
  }
 }

99 *************************/

 if (mode == _toGlobalCS) {
  ttype = RigidArmSlaveDof::_toSlave;
  if (dofIDArry) ndof = dofIDArry->giveSize(); else ndof = numberOfDofs;

  // assemble transformation contributions from local dofs
  masterTransf.resize (ndof, numberOfDofs); // master and slave should have the same dofs
  for (i=1; i<= ndof; i++) {
   //uncomment if dofIDmap assembled
   //this->giveDof (dofIDmap.at(i))->computeTransformation (dofTransf);
   if (slaveDofMask.at(i)) { // linked (slave) dof 
    if (dofIDArry == NULL) {
     ((RigidArmSlaveDof*)(this->giveDof (i)))->computeDofTransformation (dofTransf, ttype);
     for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
    } else if ((indx = this->findDofWithDofId (dofIDArry->at(i)))) {
     ((RigidArmSlaveDof*)(this->giveDof (indx)))->computeDofTransformation (dofTransf, ttype);
     for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
    } else {
      _error2 ("computeDofTransformation: DofManager %d : uncompatible dof requested",this->giveNumber());
    }
   } else { // primary dof
    masterTransf.at(i,i) = 1.0;
   }
  }
 } else if (mode == _toNodalCS) {
  ttype = RigidArmSlaveDof::_toMaster;

  IntArray dofIDmap (numberOfDofs);
  if (dofIDArry) {
   ndof = dofIDArry->giveSize(); 

   // compute dofIDmap in advance
   
   for (i=1; i<= ndof; i++) {
    indx = this->findDofWithDofId (dofIDArry->at(i));
    if (indx) dofIDmap.at(indx) = i; else {
      _error2 ("computeTransformation: DofManager %d : uncompatible dof requested",this->giveNumber());
    }
   }
  } else ndof = numberOfDofs;

  // assemble transformation contributions from local dofs
  masterTransf.resize (numberOfDofs, ndof); // master and slave should have the same dofs
  for (i=1; i<= numberOfDofs; i++) {
   //uncomment if dofIDmap assembled
   //this->giveDof (dofIDmap.at(i))->computeTransformation (dofTransf);
   if (slaveDofMask.at(i)) { // linked (slave) dof 
    ((RigidArmSlaveDof*)(this->giveDof (i)))->computeDofTransformation (dofTransf, ttype);
    if (dofIDArry) {
     for (j=1; j<= numberOfDofs; j++) if ((indx = dofIDmap.at(j))) masterTransf.at(i,j) = dofTransf.at(1,indx);
    } else for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
   } else { // primary dof
    masterTransf.at(i,i) = 1.0;
   }
  }
 } else _error ("computeDofTransformation: unknown mode");
 

 // assemble transformation from global cs to nodal local cs.
 if (master->hasLocalCS()) {
  // the trasformation of all dofs must be assembled
  // because masterTransf depend generally on all master dofs
  master->computeGNDofTransformation (GNTransf, NULL);
  if (mode == _toGlobalCS) {
   FloatMatrix help;
   help.beTranspositionOf (GNTransf);
   GNTransf = help;
  }
 } // end hasLCS
 
 // return result
 if (master->hasLocalCS()) {
  // both transformation apply
  if (mode == _toNodalCS) 
   answer.beProductOf (GNTransf, masterTransf);
  else answer.beProductOf (masterTransf, GNTransf);
 } else {
  answer = masterTransf;
 } 
}


void
RigidArmNode::computeLoadTransformation (FloatMatrix& answer, const IntArray* dofIDArry, DofManTrasfType mode)
{
 // computes trasformation matrix of receiver.
 // transformation should include trasformation from global cs to nodal cs,
 // as well as further necessary transformations (for example in case 
 // rigid arms this must include transformation to master dofs).
 FloatMatrix GNTransf, masterTransf, dofTransf;
 int i, j, indx, ndof;
 RigidArmSlaveDof::RigidArmSlaveDofTransfType ttype;
 Node* master = giveMasterDofMngr();
/*
 // compute dofIDmap in advance
 IntArray dofIDmap (ndof);
 for (i=1; i<= ndof; i++) {
  indx = this->findDofWithDofId (dofIDArry.at(i));
  if (indx) dofIDmap.at(i) = indx; else {
   char buff [80];
   sprintf(buff,"computeTransformation: DofManager %d : uncompatible dof requested",this->giveNumber());
   _error (buff);
  }
 }
*/

/************ 99
 if (mode == _toGlobalCS) ttype = RigidArmSlaveDof::_toSlave;
 else if (mode == _toNodalCS) ttype = RigidArmSlaveDof::_toMaster;
 else _error ("computeLoadTransformation: unknown mode");

 if (dofIDArry) ndof = dofIDArry->giveSize(); else ndof = numberOfDofs;

 // assemble transformation contributions from local dofs
 masterTransf.resize (ndof, numberOfDofs); // master and slave should have the same dofs
 for (i=1; i<= ndof; i++) {
  //uncomment if dofIDmap assembled
  //this->giveDof (dofIDmap.at(i))->computeTransformation (dofTransf);
  if (slaveDofMask->at(i)) { // linked (slave) dof 
   if (dofIDArry == NULL) {
    ((RigidArmSlaveDof*)(this->giveDof (i)))->computeLoadTransformation (dofTransf, ttype);
    for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
   } else if (indx = this->findDofWithDofId (dofIDArry->at(i))) {
    ((RigidArmSlaveDof*)(this->giveDof (indx)))->computeLoadTransformation (dofTransf, ttype);
    for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
   } else {
    char buff [80];
    sprintf(buff,"computeDofTransformation: DofManager %d : uncompatible dof requested",this->giveNumber());
    _error (buff);
   }
  } else { // primary dof
   masterTransf.at(i,i) = 1.0;
  }
 }
99 *************/

 if (mode == _toGlobalCS) {
  ttype = RigidArmSlaveDof::_toSlave;
  if (dofIDArry) ndof = dofIDArry->giveSize(); else ndof = numberOfDofs;

  // assemble transformation contributions from local dofs
  masterTransf.resize (ndof, numberOfDofs); // master and slave should have the same dofs
  for (i=1; i<= ndof; i++) {
   //uncomment if dofIDmap assembled
   //this->giveDof (dofIDmap.at(i))->computeTransformation (dofTransf);
   if (slaveDofMask.at(i)) { // linked (slave) dof 
    if (dofIDArry == NULL) {
     ((RigidArmSlaveDof*)(this->giveDof (i)))->computeLoadTransformation (dofTransf, ttype);
     for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
    } else if ((indx = this->findDofWithDofId (dofIDArry->at(i)))) {
     ((RigidArmSlaveDof*)(this->giveDof (indx)))->computeLoadTransformation (dofTransf, ttype);
     for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
    } else {
      _error2 ("computeLoadTransformation: DofManager %d : uncompatible dof requested",this->giveNumber());
    }
   } else { // primary dof
     if (dofIDArry) {
       if ((indx = this->findDofWithDofId (dofIDArry->at(i)))) masterTransf.at(i,i) = 1.0;
     } else masterTransf.at(i,i) = 1.0;
     //masterTransf.at(i,i) = 1.0;
   }
  }
 } else if (mode == _toNodalCS) {
  ttype = RigidArmSlaveDof::_toMaster;

  IntArray dofIDmap (numberOfDofs);
  if (dofIDArry) {
   ndof = dofIDArry->giveSize(); 

   // compute dofIDmap in advance
   
   for (i=1; i<= ndof; i++) {
    indx = this->findDofWithDofId (dofIDArry->at(i));
    if (indx) dofIDmap.at(indx) = i; else {
      _error2 ("computeTransformation: DofManager %d : uncompatible dof requested",this->giveNumber());
    }
   }
  } else ndof = numberOfDofs;

  // assemble transformation contributions from local dofs
  masterTransf.resize (numberOfDofs, ndof); // master and slave should have the same dofs
  for (i=1; i<= numberOfDofs; i++) {
   //uncomment if dofIDmap assembled
   //this->giveDof (dofIDmap.at(i))->computeTransformation (dofTransf);
   if (slaveDofMask.at(i)) { // linked (slave) dof 
    ((RigidArmSlaveDof*)(this->giveDof (i)))->computeLoadTransformation (dofTransf, ttype);
    if (dofIDArry) {
     for (j=1; j<= numberOfDofs; j++) if ((indx = dofIDmap.at(j))) masterTransf.at(i,j) = dofTransf.at(1,indx);
    } else for (j=1; j<= numberOfDofs; j++) masterTransf.at(i,j) = dofTransf.at(1,j);
   } else { // primary dof
     if (dofIDArry) {
       if ((indx = dofIDmap.at(i))) masterTransf.at(i,i) = 1.0;
     } else masterTransf.at(i,i) = 1.0;
     //masterTransf.at(i,i) = 1.0;
   }
  }
 } else _error ("computeLoadTransformation: unknown mode");
 

 // assemble transformation from global cs to nodal local cs.
 if (master->hasLocalCS()) {
  // the trasformation of all dofs must be assembled
  // because masterTransf depend generally on all master dofs
  master->computeGNDofTransformation (GNTransf, NULL);
  if (mode == _toGlobalCS) {
   FloatMatrix help;
   help.beTranspositionOf (GNTransf);
   GNTransf = help;
  }
 } // end hasLCS
 
 // return result
 if (master->hasLocalCS()) {
  // both transformation apply
  if (mode == _toNodalCS) 
   answer.beProductOf (GNTransf, masterTransf);
  else answer.beProductOf (masterTransf, GNTransf);
 } else {
  answer = masterTransf;
 } 
}



void
RigidArmNode :: computeLoadVectorAt (FloatArray& answer, TimeStep* stepN, ValueModeType mode)
   // Computes the vector of the nodal loads of the receiver.
{
 FloatMatrix masterTransf;
 //FloatArray localAnswer;

 // assemble answer of receiver for "local dofs"
 Node :: computeLoadVectorAt (answer,stepN, mode);

 // transform "local dofs" to master dofs
 if (answer.isNotEmpty ()) {
  // assemble transformation contributions from local dofs
  computeLoadTransformation (masterTransf, NULL, _toNodalCS);
  answer.rotatedWith (masterTransf,'n');
 }
 
}


void  RigidArmNode :: giveCompleteLocationArray (IntArray& locationArray) const
   // Returns the complete location array of the receiver.
   // including all available dofs
{
 Node* master = giveMasterDofMngr();
 int i ;
 // prevents some size problem when connecting different elements with 
 // different number of dofs
 locationArray.resize (numberOfDofs) ;
 for (i=1;i<=numberOfDofs;i++) {
  if (slaveDofMask.at(i)) // linked (slave) dof 
   locationArray.at(i)=master->giveDof(i)->giveEquationNumber();
  else
   locationArray.at(i)=this->giveDof(i)->giveEquationNumber();
 }
 return ;
}

void  RigidArmNode :: giveCompletePrescribedLocationArray (IntArray& locationArray) const
   // Returns the complete location array of prescribed equations of the receiver.
   // including all available dofs
{
 Node* master = giveMasterDofMngr();
 int i ;
 // prevents some size problem when connecting different elements with 
 // different number of dofs
 locationArray.resize (numberOfDofs) ;
 for (i=1;i<=numberOfDofs;i++) {
  if (slaveDofMask.at(i)) // linked (slave) dof 
   locationArray.at(i)=master->giveDof(i)->givePrescribedEquationNumber();
  else
   locationArray.at(i)=this->giveDof(i)->givePrescribedEquationNumber();
 }
 return ;
}


void 
RigidArmNode :: giveLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const
   // Returns the location array of the receiver. Creates this array if it
   // does not exist yet. The location array contains the equation number of
   // every  requested degree of freedom of the receiver.
  // In dofIDArray are stored DofID's of requsted DOFs in receiver.
  // The DofID's are determining the physical meaning of particular DOFs
   //
   // Rigid arm return always full location array
{
 this->giveCompleteLocationArray (locationArray);
}

void 
RigidArmNode :: givePrescribedLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const
   // Returns the location array of prescribed equations of the receiver. Creates this array if it
   // does not exist yet. The location array contains the equation number of
   // every  requested degree of freedom of the receiver.
  // In dofIDArray are stored DofID's of requsted DOFs in receiver.
  // The DofID's are determining the physical meaning of particular DOFs
   //
   // Rigid arm return always full location array
{
 this->giveCompletePrescribedLocationArray (locationArray);
}


void 
RigidArmNode::giveUnknownVector (FloatArray& answer, const IntArray& dofMask, 
                                 EquationID type, ValueModeType mode, TimeStep* stepN)
{
 Node* master = giveMasterDofMngr();
 // IntArray fullDofMask(numberOfDofs);
 int i;

 
 //for (i=1; i<= numberOfDofs; i++) 
 // fullDofMask.at(i) = giveDof(i)->giveDofID();
 answer.resize (numberOfDofs);

 for (i=1; i<= numberOfDofs; i++) {
  if (slaveDofMask.at(i)) { // linked (slave) dof 
   answer.at(i) = master->giveDof(i)->giveUnknown(type, mode, stepN) ;
  } else  { // primary DOF
   answer.at(i) = this->giveDof(i)->giveUnknown(type, mode, stepN) ;
  }
 }
}

void 
RigidArmNode::giveUnknownVector (FloatArray& answer, const IntArray& dofMask, 
                 PrimaryField& field, ValueModeType mode, TimeStep* stepN)
{
 Node* master = giveMasterDofMngr();
 int i;

 answer.resize (numberOfDofs);

 for (i=1; i<= numberOfDofs; i++) {
  if (slaveDofMask.at(i)) { // linked (slave) dof 
   answer.at(i) = master->giveDof(i)->giveUnknown(field, mode, stepN) ;
  } else  { // primary DOF
   answer.at(i) = this->giveDof(i)->giveUnknown(field, mode, stepN) ;
  }
 }
}

void 
RigidArmNode:: givePrescribedUnknownVector (FloatArray& answer, const IntArray& dofMask, 
                      ValueModeType mode, TimeStep* stepN)
{
 Node* master = giveMasterDofMngr();
 // IntArray fullDofMask(numberOfDofs);
 Dof        *dofJ ;
 int i;

 //for (i=1; i<= numberOfDofs; i++) 
 // fullDofMask.at(i) = giveDof(i)->giveDofID();
 answer.resize (numberOfDofs);

 for (i=1; i<= numberOfDofs; i++) {
  if (slaveDofMask.at(i)) { // linked (slave) dof 
   dofJ = master->giveDof(i) ;
   //if (dofJ -> hasBc(stepN) && (dofJ ->giveUnknownType() == type))
   if (dofJ -> hasBc(stepN))
    answer.at(i) = dofJ->giveBcValue(mode,stepN) ;//giveUnknown(u,stepN) ;
     //answer.at(i) = dofJ->giveBcValue(type, mode,stepN) ;//giveUnknown(u,stepN) ;
   else
    answer.at(i) = 0. ;
  } else { // primary DOF
   dofJ = this->giveDof(i) ;
   
   //if (dofJ -> hasBc(stepN) && (dofJ ->giveUnknownType() == type))
   if (dofJ -> hasBc(stepN))
    answer.at(i) = dofJ->giveBcValue(mode,stepN) ;//giveUnknown(u,stepN) ;
    //answer.at(i) = dofJ->giveBcValue(type, mode,stepN) ;//giveUnknown(u,stepN) ;
   else
    answer.at(i) = 0. ;
  }   
 }

 //master->givePrescribedUnknownVector (answer, fullDofMask, type, mode, stepN);
}



double 
RigidArmNode :: giveUpdatedCoordinate (int ic ,TimeStep* tStep, EquationID type, double scale )
//
// returns coordinate + scale * displacement
// displacement is of updMode (UpdateMode) type
//
{
 int i, j;
 FloatMatrix *T;
 
 if ((ic < 1) || (ic > 3)) {
  _error ("giveUpdatedCoordinate: Can't return non-existing coordinate (index not in range 1..3)");
  return 0.;
 }

 if (tStep->isTheCurrentTimeStep ()) {
  double coordinate = this->giveCoordinate (ic);
  if (!this->hasLocalCS ()) {
   // this has no local cs.
   for (i=1 ; i<=numberOfDofs ; i++) {
    j = domain ->  giveCorrespondingCoordinateIndex (i);
    if ((j != 0) && (j == ic))   {
     if (slaveDofMask.at(i)) { // linked (slave) dof 
      coordinate += 
       scale * ((RigidArmSlaveDof*)this->giveDof(i))->giveLocalUnknown (type,VM_Total,tStep);
     } else {
      coordinate += 
       scale * this->giveDof(i)->giveUnknown (type,VM_Total,tStep);
     }
     break;
    }
   }
  } else {
   //
   // this has local cs.
   // We must perform transformation of displacements DOFs
   // in to global c.s and then to add them to global coordinates.
   //
   T = this->giveLocalCoordinateTriplet() ;
   FloatArray displacements (3) ;
   for (i=1 ; i<= 3; i++) displacements.at(i) = 0.;
   for (i=1 ; i<=numberOfDofs; i++) {
    j = domain ->  giveCorrespondingCoordinateIndex (i);
    if (j) // && (this->giveDof(i)->giveUnknownType()==DisplacementVector))
     if (slaveDofMask.at(i)) { // linked (slave) dof 
      displacements.at(j) = scale * ((RigidArmSlaveDof*)this->giveDof(i))->
       giveLocalUnknown (type,VM_Total,tStep);
     } else {
      displacements.at(j) = scale * this->giveDof(i)->
       giveUnknown (type,VM_Total,tStep);
    }
   }
   // perform transformation for desired displacement
   for (i=1 ; i<= 3; i++) 
    coordinate += displacements.at(i) * T-> at(i,ic);
  }
  return coordinate;
 } else { _error("Can't return updatedCoordinate for non-current timestep");}
 return 0.;
}
