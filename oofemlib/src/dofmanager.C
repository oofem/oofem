/* $Header: /home/cvs/bp/oofem/oofemlib/src/dofmanager.C,v 1.18.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   file dofmanager.C

#include "dofmanager.h"
#include "masterdof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "timestep.h"
#include "load.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "debug.h"
#include "verbose.h"
#include "node.h"
#include "elementside.h"
#include "rigidarmnode.h"
#include "hangingnode.h"
#include "usrdefsub.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#endif

#ifndef __MAKEDEPEND
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#endif

#ifdef __PARALLEL_MODE
#include "remotemasterdof.h"
#include "sharedmasterdof.h"
#include "nulldof.h"
#endif

DofManager :: DofManager (int n, Domain* aDomain)
      : FEMComponent (n,aDomain), loadArray()
   // Constructor. Creates a node with number n, belonging to aDomain.
{
   numberOfDofs  = 0 ;
   dofArray      = NULL ;
  isBoundaryFlag= false;
  hasSlaveDofs  = false;
   // locationArray = NULL ;
#ifdef __PARALLEL_MODE
  partitions.resize(0);
#endif
}




DofManager :: ~DofManager()
   // Destructor.
{
   int i = numberOfDofs ;

   if (numberOfDofs) {
      while (i--)
  delete dofArray[i] ;
  //      delete [numberOfDofs] dofArray ;}
      delete dofArray ;}
   // delete locationArray ;
#ifdef __PARALLEL_MODE
#endif
}


void
DofManager :: computeLoadVectorAt (FloatArray& answer, TimeStep* stepN, ValueModeType mode)
   // Computes the vector of the nodal loads of the receiver.
{
 int        i,n,nLoads ;
 Load  *loadN ;
 FloatArray contribution ;
 
 if (this -> giveLoadArray() -> isEmpty()) {
  answer.resize (0);
  return ;
 } else {
  answer.resize(0) ;
  nLoads = loadArray.giveSize() ;          // the node may be subjected
  for (i=1 ; i<=nLoads ; i++) {             // to more than one load
   n            = loadArray.at(i) ;
   loadN        = domain->giveLoad(n) ;
   if (loadN->giveBCGeoType() != NodalLoadBGT) 
    _error ("computeLoadVectorAt: incompatible load type applied");
   loadN -> computeComponentArrayAt(contribution, stepN, mode) ;      // can be NULL
   answer.add(contribution) ;
   //delete contribution ;
  }
 }

 // transform "local dofs" to master dofs
 if (hasSlaveDofs && answer.isNotEmpty ()) {
   FloatMatrix masterTransf;
   
   // assemble transformation contributions from local dofs
   computeLoadTransformation (masterTransf, NULL, _toNodalCS);
   answer.rotatedWith (masterTransf,'n');
 }


 return ;
 
}




Dof*  DofManager :: giveDof (int i) const
   // Returns the i-th degree of freedom of the receiver. Creates the array
   // containing the dofs of the receiver, if such array does not exist yet.
{

   if (! dofArray) {
//      dofArray = new Dof* [this->giveNumberOfDofs()] ;
//      for (j=0 ; j<numberOfDofs ; j++)
//  dofArray[j] = new Dof(j+1,this) ;}
     _error("giveDof: dof is not defined");
   }
   return dofArray[i-1] ;
}

Dof*  DofManager :: giveDofWithID (int dofID) const
  // Returns the degree of freedom of the receiver with 'dofID'.
{
  int indx = this->findDofWithDofId ((DofID) dofID);
  
  // musi zde byt error - spoleham na to
  if (!indx) _error ("giveDofWithID: dof with given DofID doesnot exists");
  
  return dofArray[indx-1] ;
}

IntArray*  DofManager :: giveLoadArray ()
   // Returns the list containing the number of every nodal loads that act on
   // the receiver. If this list does not exist yet, constructs it. This list
   // is not to be confused with the load vector.
{
   return &loadArray ;
}


void 
DofManager :: giveLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const
   // Returns the location array of the receiver. Creates this array if it
   // does not exist yet. The location array contains the equation number of
   // every  requested degree of freedom of the receiver.
  // In dofIDArray are stored DofID's of requsted DOFs in receiver.
  // The DofID's are determining the physical meaning of particular DOFs
{
  if (!hasSlaveDofs) {

    int i,size,indx ;
    // prevents some size problem when connecting different elements with 
    // different number of dofs
    size = dofIDArry.giveSize();
    locationArray.resize(size) ;
    for (i=1;i<=size;i++) {
      if ((indx=this->findDofWithDofId ((DofID) dofIDArry.at(i))) == 0) {
        _error ("giveLocationArray: incompatible dof requested");
      }
      locationArray.at(i)=this->giveDof(indx)->giveEquationNumber();
    }
  } else {

  
    int i, k, indx;
    IntArray dofArray, mstrEqNmbrs;
    
    this-> giveDofArray(dofIDArry, dofArray);
    locationArray.resize (giveNumberOfPrimaryMasterDofs(dofArray));
    
    for (k=1,i=1; i<=dofArray.giveSize(); i++) {
      indx = dofArray.at(i);
      if (!this->giveDof(indx)->isPrimaryDof()) { // slave DOF
        this->giveDof(indx)->giveEquationNumbers (mstrEqNmbrs);
        locationArray.copySubVector (mstrEqNmbrs, k);
        k += mstrEqNmbrs.giveSize();
      }
      else { // primary DOF
        locationArray.at(k++) = this-> giveDof(indx)-> giveEquationNumber();
      }
    }


  }
  return  ;
}



void  DofManager :: giveCompleteLocationArray (IntArray& locationArray) const
   // Returns the complete location array of the receiver.
   // including all available dofs
{
  if (!hasSlaveDofs) {

    int i ;
    // prevents some size problem when connecting different elements with 
    // different number of dofs
    locationArray.resize (numberOfDofs) ;
    for (i=1;i<=numberOfDofs;i++) {
      locationArray.at(i)=this->giveDof(i)->giveEquationNumber();
    }
  } else {
    giveLocationArray (*giveCompleteGlobalDofIDArray(), locationArray);
  }
 return ;
}

void 
DofManager :: givePrescribedLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const
   // Returns the prescribed equation location array of the receiver. Creates this array if it
   // does not exist yet. The location array contains the prescribed equation number of
   // every  requested degree of freedom of the receiver.
  // In dofIDArray are stored DofID's of requsted DOFs in receiver.
  // The DofID's are determining the physical meaning of particular DOFs
{
  if (!hasSlaveDofs) {
    int i,size,indx ;
    // prevents some size problem when connecting different elements with 
    // different number of dofs
    size = dofIDArry.giveSize();
    locationArray.resize(size) ;
    for (i=1;i<=size;i++) {
      if ((indx=this->findDofWithDofId ((DofID) dofIDArry.at(i))) == 0) {
        _error ("givePrescribedLocationArray: incompatible dof requested");
      }
      locationArray.at(i)=this->giveDof(indx)->givePrescribedEquationNumber();
    }
  } else {

  
    int i, k, indx;
    IntArray dofArray, mstrEqNmbrs;
    
    this-> giveDofArray(dofIDArry, dofArray);
    locationArray.resize (giveNumberOfPrimaryMasterDofs(dofArray));
    
    for (k=1,i=1; i<=dofArray.giveSize(); i++) {
      indx = dofArray.at(i);
      if (!this->giveDof(indx)->isPrimaryDof()) { // slave DOF
        this->giveDof(indx)->givePrescribedEquationNumbers (mstrEqNmbrs);
        locationArray.copySubVector (mstrEqNmbrs, k);
        k += mstrEqNmbrs.giveSize();
      }
      else { // primary DOF
        locationArray.at(k++) = this-> giveDof(indx)-> givePrescribedEquationNumber();
      }
    }
    
  }
  
  return  ;
}


void  DofManager :: giveCompletePrescribedLocationArray (IntArray& locationArray) const
   // Returns the complete location array of the receiver.
   // including all available dofs
{
  if (!hasSlaveDofs) {

    int i ;
    // prevents some size problem when connecting different elements with 
    // different number of dofs
    locationArray.resize (numberOfDofs) ;
    for (i=1;i<=numberOfDofs;i++) {
      locationArray.at(i)=this->giveDof(i)->givePrescribedEquationNumber();
    }
  } else {
    givePrescribedLocationArray (*giveCompleteGlobalDofIDArray(), locationArray);
  }
}


void 
DofManager :: giveDofArray (const IntArray& dofIDArry, IntArray& answer) const
   // Returns the dof index array of the receiver.
   // The location array contains the indexes of particular requsted DOFs
  // In dofIDArray are stored DofID's of requsted DOFs in receiver.
  // The DofID's are determining the physical meaning of particular DOFs
{

 int i,size ;
 // IntArray* answer;
 // prevents some size problem when connecting different elements with 
 // different number of dofs
 size = dofIDArry.giveSize();
 //answer = new IntArray(size) ;
 answer.resize (size);
 for (i=1;i<=size;i++) {
   if ((answer.at(i)=this->findDofWithDofId ((DofID) dofIDArry.at(i))) == 0) {
     _error ("giveDofArray : incompatible dof requested");
   }
 }
 return  ;
}

int
DofManager :: findDofWithDofId (DofID dofID) const
{
 // finds index of DOF in receivers node with dofID
 // if such DOF does not exists, returns zero value
 int i;
  for (i=1 ; i<=numberOfDofs ; i++) {
  if (this->giveDof(i)->giveDofID() == dofID) return i;
 }
 // nothing found
 return 0;
}




int  DofManager :: giveNumberOfDofs () const
   // Returns the number of degrees of freedom of the receiver.
{
   if (! numberOfDofs)
      _error ("giveNumberOfDofs: NumberOfDofs is not known yet");

   return numberOfDofs ;
}

int
DofManager :: giveNumberOfPrimaryMasterDofs (IntArray& dofArray) const
{
  if (!hasSlaveDofs) return this->giveNumberOfDofs();
  
  int i, answer=0;
  
  for (i=1; i<=dofArray.giveSize(); i++)
    if (!this->giveDof(dofArray.at(i))->isPrimaryDof())
      answer += ((SlaveDof*)this->giveDof(dofArray.at(i)))->giveNumberOfPrimaryMasterDofs();
    else
      answer += 1;
  
  return answer;
}



IRResultType DofManager::  resolveDofIDArray (InputRecord* ir, IntArray& dofIDArry)
{
 const char *__keyword, *__proc = "resolveDofIDArray";
 IRResultType result; 

 numberOfDofs = 0;
 __keyword = "ndofs"; result = ir->giveOptionalField(numberOfDofs, IFT_DofManager_ndofs, __keyword);
 if (result != IRRT_OK) IR_IOERR (giveClassName(), __proc, IFT_DofManager_ndofs, __keyword, ir, result);

 // returns nonzero if succes
 if(numberOfDofs == 0) {
  numberOfDofs = domain-> giveNumberOfDefaultNodeDofs () ;
  dofIDArry = domain->giveDefaultNodeDofIDArry ();
 } else {
  // if ndofs is prescribed, read the physical meaning of particular dofs
  // for detailed values of DofMask array see cltypes.h file
  // for exaple 1 is for D_u (displacemet in u dir), 2 for D_v, 3 for D_w, ...
  __keyword = "dofidmask"; result = ir->giveField(dofIDArry, IFT_DofManager_dofidmask, __keyword);
  if (result != IRRT_OK) IR_IOERR (giveClassName(), __proc, IFT_DofManager_dofidmask, __keyword, ir, result);
  if (dofIDArry.giveSize() != numberOfDofs) {
    _error ("resolveDofIDArray : DofIDMask size mismatch");
  }
 }
 return IRRT_OK;
}

IRResultType
DofManager :: initializeFrom (InputRecord* ir)
{
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                   // Required by IR_GIVE_FIELD macro

  int j ;
  IntArray dofIDArry;
  IntArray bc, ic, masterMask, dofTypeMask; // termitovo
 
  loadArray.resize(0);
  IR_GIVE_OPTIONAL_FIELD (ir, loadArray, IFT_DofManager_load, "load"); // Macro

  if (this->resolveDofIDArray (ir, dofIDArry) != IRRT_OK) IR_IOERR (giveClassName(), __proc,  IFT_Unknown, "", ir, result);

  // numberOfDofs = domain->giveNumberOfDofs () ;
  bc.resize(0);
  IR_GIVE_OPTIONAL_FIELD (ir, bc, IFT_DofManager_bc, "bc"); // Macro

  ic.resize(0);
  IR_GIVE_OPTIONAL_FIELD (ir, ic, IFT_DofManager_ic, "ic"); // Macro
  // reads master mask - in this array are numbers of master dofManagers
  // to which are connected dofs in receiver.
  // if master mask index is zero then dof is created as master (i.e., having own equation number)
  // othervise slave dof connected to master DofManager is created.
  // by default if masterMask is not specifyed, all dofs are created as masters.
  dofTypeMask.resize(0); // termitovo
  IR_GIVE_OPTIONAL_FIELD (ir, dofTypeMask, IFT_DofManager_doftypemask, "doftype"); // Macro

  // read boundary flag
  if (ir->hasField(IFT_DofManager_boundaryflag, "boundary")) isBoundaryFlag = true;


#ifdef __PARALLEL_MODE
  globalNumber = 0;
  IR_GIVE_OPTIONAL_FIELD (ir, globalNumber, IFT_DofManager_globnum, "globnum"); // Macro

  partitions.resize(0);
  IR_GIVE_OPTIONAL_FIELD (ir, partitions, IFT_DofManager_partitions, "partitions"); // Macro

  if (ir->hasField (IFT_DofManager_sharedflag, "shared")) parallel_mode = DofManager_shared;
  else if (ir->hasField (IFT_DofManager_remoteflag, "remote")) parallel_mode = DofManager_remote;
  else if (ir->hasField (IFT_DofManager_nullflag, "null")) parallel_mode = DofManager_null;
  else parallel_mode = DofManager_local;

  // in parallel mode,  slaves are allowed, because ((Dr. Rypl promissed)
  // masters have to be in same partition as slaves. They can be again Remote copies.
#endif



  int hasIc,hasBc,dofIc=0,dofBc=0, hasTypeinfo=0; 
  dofType dtype;
 
  hasIc = !(ic.giveSize() == 0);
  hasBc = !(bc.giveSize() == 0);
  hasTypeinfo = !(dofTypeMask.giveSize() == 0);

  // check sizes
  if (hasBc) if (bc.giveSize() != this->giveNumberOfDofs()) _error ("initializeFrom: bc size mismatch");
  if (hasIc) if (ic.giveSize() != this->giveNumberOfDofs()) _error ("initializeFrom: ic size mismatch");
  if (hasTypeinfo) if (dofTypeMask.giveSize() != this->giveNumberOfDofs()) 
    _error ("initializeFrom: dofTypeMask size mismatch");

  dofArray = new Dof* [this->giveNumberOfDofs()] ;
  for (j=0 ; j<numberOfDofs ; j++) {
  
    if (hasTypeinfo) dtype = (dofType) dofTypeMask.at(j+1);
    else dtype = DT_master;
    
    if (this->isDofTypeCompatible(dtype)) {
      
      if (dtype == DT_master) {
        if(hasIc) dofIc = ic.at(j+1) ;
        if(hasBc) dofBc = bc.at(j+1) ;
#ifdef __PARALLEL_MODE
        if (parallel_mode == DofManager_remote) 
          dofArray[j] = new RemoteMasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
        else if (parallel_mode == DofManager_shared)
          dofArray[j] = new SharedMasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
        else if (parallel_mode == DofManager_null)
          // ignore applied bc
          dofArray[j] = new NullDof (j+1, this, dofIc, (DofID) dofIDArry.at(j+1)) ;
        else
          dofArray[j] = new MasterDof      (j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
#else
        dofArray[j] = new MasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
#endif
      } else if (dtype == DT_simpleSlave) { // Simple slave dof
        if (masterMask.giveSize() == 0) {
          IR_GIVE_FIELD (ir, masterMask, IFT_DofManager_mastermask, "mastermask"); // Macro
          if (masterMask.giveSize() != numberOfDofs) _error ("initializeFrom: mastermask size mismatch");
        }
        dofArray[j] = new SimpleSlaveDof (j+1,this,masterMask.at(j+1),(DofID) dofIDArry.at(j+1)) ;
      } else if (dtype == DT_slave) { // Slave dof
        dofArray[j] = new SlaveDof (j+1,this,(DofID) dofIDArry.at(j+1)) ;
      } else {
        _error ("initializeFrom: unknown dof type");
      }
    } else {
      _error ("initializeFrom: incompatible dof type");
    }
  }
  return IRRT_OK;
}

/*
IRResultType
DofManager :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 int j ;
 IntArray dofIDArry;
 IntArray bc, ic, masterMask;
 
 loadArray.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, loadArray, IFT_DofManager_load, "load"); // Macro

 if (this->resolveDofIDArray (ir, dofIDArry) != IRRT_OK) IR_IOERR (giveClassName(), __proc,  IFT_Unknown, "", ir, result);

// numberOfDofs = domain->giveNumberOfDofs () ;
 bc.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, bc, IFT_DofManager_bc, "bc"); // Macro

 ic.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, ic, IFT_DofManager_ic, "ic"); // Macro
 // reads master mask - in this array are numbers of master dofManagers
 // to which are connected dofs in receiver.
 // if master mask index is zero then dof is created as master (i.e., having own equation number)
 // othervise slave dof connected to master DofManager is created.
 // by default if masterMask is not specifyed, all dofs are created as masters.
 masterMask.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, masterMask, IFT_DofManager_mastermask, "mastermask"); // Macro

 // read boundary flag
 if (ir->hasField(IFT_DofManager_boundaryflag, "boundary")) isBoundaryFlag = true;


#ifdef __PARALLEL_MODE
 globalNumber = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, globalNumber, IFT_DofManager_globnum, "globnum"); // Macro

 partitions.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, partitions, IFT_DofManager_partitions, "partitions"); // Macro

 if (ir->hasField (IFT_DofManager_sharedflag, "shared")) parallel_mode = DofManager_shared;
 else if (ir->hasField (IFT_DofManager_remoteflag, "remote")) parallel_mode = DofManager_remote;
 else if (ir->hasField (IFT_DofManager_nullflag, "null")) parallel_mode = DofManager_null;
 else parallel_mode = DofManager_local;

 // in parallel mode,  slaves are allowed, because ((Dr. Rypl promissed)
 // masters have to be in same partition as slaves. They can be again Remote copies.
#endif



 int hasIc,hasBc,dofIc=0,dofBc=0, hasSlaveDofs=0, masterManagerIndx=0 ;
 
 hasIc = !(ic.giveSize() == 0);
 hasBc = !(bc.giveSize() == 0);
 hasSlaveDofs = !(masterMask.giveSize() == 0);

 // check sizes
 if (hasBc) if (bc.giveSize() != this->giveNumberOfDofs()) _error ("initializeFrom: bc size mismatch");
 if (hasIc) if (ic.giveSize() != this->giveNumberOfDofs()) _error ("initializeFrom: ic size mismatch");
 if (hasSlaveDofs) if (masterMask.giveSize() != this->giveNumberOfDofs()) 
   _error ("initializeFrom: masterMask size mismatch");

   dofArray = new Dof* [this->giveNumberOfDofs()] ;
   for (j=0 ; j<numberOfDofs ; j++)
     {
       if(hasIc) dofIc = ic.at(j+1) ;
       if(hasBc) dofBc = bc.at(j+1) ;
   
    if (hasSlaveDofs) {
     masterManagerIndx = masterMask.at(j+1);
     if (masterManagerIndx) {
      dofArray[j] = new SimpleSlaveDof (j+1,this,masterManagerIndx, (DofID) dofIDArry.at(j+1)) ;
     } else {
#ifdef __PARALLEL_MODE
      if (parallel_mode == DofManager_remote) 
       dofArray[j] = new RemoteMasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
      else if (parallel_mode == DofManager_shared)
       dofArray[j] = new SharedMasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
      else if (parallel_mode == DofManager_null)
       // ignore applied bc
       dofArray[j] = new NullDof (j+1, this, dofIc, (DofID) dofIDArry.at(j+1)) ;
      else
       dofArray[j] = new MasterDof      (j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
#else
      dofArray[j] = new MasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
#endif
     }
    } else {
#ifdef __PARALLEL_MODE
     if (parallel_mode == DofManager_remote) 
      dofArray[j] = new RemoteMasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
     else if (parallel_mode == DofManager_shared)
      dofArray[j] = new SharedMasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
     else if (parallel_mode == DofManager_null)
      // ignore applied bc
      dofArray[j] = new NullDof (j+1, this, dofIc, (DofID) dofIDArry.at(j+1)) ;
     else 
      dofArray[j] = new MasterDof      (j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
#else
     dofArray[j] = new MasterDof(j+1,this,dofBc,dofIc,(DofID) dofIDArry.at(j+1)) ;
#endif
    }
     }
 // this-> giveLocationArray ();

 return IRRT_OK;
}
*/

void  DofManager :: printOutputAt (FILE* stream, TimeStep* stepN)
{
 EngngModel* emodel = this->giveDomain()->giveEngngModel();
 int  i ;

#ifdef __PARALLEL_MODE
  fprintf (stream,"%-8s%8d [%8d]:\n",this->giveClassName(),this->giveNumber(), this->giveGlobalNumber()) ;
#else
  fprintf (stream,"%-8s%8d:\n",this->giveClassName(), this->giveNumber()) ;
#endif
   for (i=1 ; i<=numberOfDofs ; i++)
   emodel->printDofOutputAt(stream, this -> giveDof(i), stepN) ;
}


void  DofManager :: printYourself ()
   // Prints the receiver on screen.
{
   int    i ;
   // double x,y ;

   printf ("DofManager %d\n",number) ;
   for (i=0 ; i<numberOfDofs ; i++) {
      if (dofArray[i])
  dofArray[i] -> printYourself() ;
      else
  printf ("dof %d is nil \n",i+1) ;}
  loadArray.printYourself() ;
   printf ("\n") ;
}


void  DofManager :: updateYourself (TimeStep* tStep)
   // Updates the receiver at end of step.
{
   int i ;

   for (i=1 ; i<=numberOfDofs ; i++) {
     this -> giveDof(i) -> updateYourself(tStep) ;

  }
 }


contextIOResultType DofManager :: saveContext (DataStream* stream, ContextMode mode, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
  int i;
  contextIOResultType iores;
  
  if ((iores = FEMComponent::saveContext(stream,mode,obj)) != CIO_OK) THROW_CIOERR(iores);
  for (i=1 ; i<=numberOfDofs ; i++) {
    if ((iores = this->giveDof(i)->saveContext(stream,mode,obj)) != CIO_OK) THROW_CIOERR(iores);
  }

  return CIO_OK;
}


contextIOResultType DofManager :: restoreContext (DataStream* stream, ContextMode mode, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
  contextIOResultType iores;
  int i ;

  if ((iores = FEMComponent::restoreContext(stream,mode,obj)) != CIO_OK) THROW_CIOERR(iores);
  for (i=1 ; i<=numberOfDofs ; i++) {
    if ((iores = this->giveDof(i)->restoreContext(stream,mode,obj)) != CIO_OK) THROW_CIOERR(iores);
  }

  return CIO_OK;
}


DofManager*  DofManager :: ofType (char* aClass)
   // Returns a new DofManager, which has the same number than the receiver,
   // but belongs to aClass (Node, ElementSide,..).
{
  DofManager* newDofManager ;
  
  if (! strncasecmp(aClass,"node",4))
    newDofManager = new Node (number,domain) ;
  else if (! strncasecmp(aClass,"elementside",11))
    newDofManager = new ElementSide (number,domain) ; 
  else if (! strncasecmp(aClass,"rigidarmnode",12))
    newDofManager = new RigidArmNode (number,domain) ; 
  else if (! strncasecmp(aClass,"hangingnode",11))
    newDofManager = new HangingNode (number,domain) ; 
  else {   // last resort - call aditional user defined subroutine
    newDofManager = ::CreateUsrDefDofManagerOfType (aClass,number,domain);
    if (newDofManager == NULL) 
      _error2 ("ofType: unknown DofManager type (%s)",aClass) ;
  }
  return newDofManager ;
}


void 
DofManager::giveUnknownVector (FloatArray& answer, const IntArray& dofMask, 
                               EquationID type, ValueModeType mode, TimeStep* stepN)
{
 if (!hasSlaveDofs) {

   int j, size;
   IntArray dofArray;

   answer.resize (size = dofMask.giveSize());
   this-> giveDofArray(dofMask, dofArray);
   
   for (j=1 ; j<=size ; j++)
     //if (this->giveDof(dofArray.at(j))->giveUnknownType() == type)
     answer.at(j) = this->giveDof(dofArray.at(j))->giveUnknown(type, mode, stepN) ;
   //else 
   // _error ("giveUnknownVector :: dofMask compatibility Error");
 } else {
   
   int i, k, indx;
   IntArray dofArray;
   FloatArray mstrUnknwns;
   
   this-> giveDofArray(dofMask, dofArray);
   answer.resize (giveNumberOfPrimaryMasterDofs(dofArray));
   
   for (k=1,i=1; i<=dofArray.giveSize(); i++) {
     indx = dofArray.at(i);
     if (!this->giveDof(indx)->isPrimaryDof()) { // slave DOF
       this->giveDof(indx)->giveUnknowns (mstrUnknwns, type, mode, stepN);
       answer.copySubVector (mstrUnknwns, k);
       k += mstrUnknwns.giveSize();
     }
     else { // primary DOF
       answer.at(k++) = this-> giveDof(indx)-> giveUnknown(type, mode, stepN);
     }
   }
 }  
 
}


void 
DofManager::giveUnknownVector (FloatArray& answer, const IntArray& dofMask, 
                PrimaryField& field, ValueModeType mode, TimeStep* stepN)
{

  if (!hasSlaveDofs) {

    int j, size;
    IntArray dofArray;
    
    answer.resize (size = dofMask.giveSize());
    this-> giveDofArray(dofMask, dofArray);
    
    for (j=1 ; j<=size ; j++)
      //if (this->giveDof(dofArray.at(j))->giveUnknownType() == type)
      answer.at(j) = this->giveDof(dofArray.at(j))->giveUnknown(field, mode, stepN) ;
    //else 
    // _error ("giveUnknownVector :: dofMask compatibility Error");
  } else {

    int i, k, indx;
    IntArray dofArray;
    FloatArray mstrUnknwns;
    
    this-> giveDofArray(dofMask, dofArray);
    answer.resize (giveNumberOfPrimaryMasterDofs(dofArray));
    
    for (k=1,i=1; i<=dofArray.giveSize(); i++) {
      indx = dofArray.at(i);
      if (!this->giveDof(indx)->isPrimaryDof()) { // slave DOF
        this->giveDof(indx)->giveUnknowns (mstrUnknwns, field, mode, stepN);
        answer.copySubVector (mstrUnknwns, k);
        k += mstrUnknwns.giveSize();
      }
      else { // primary DOF
        answer.at(k++) = this-> giveDof(indx)-> giveUnknown(field, mode, stepN);
      }
    }
  }    
  
}

void 
DofManager:: givePrescribedUnknownVector (FloatArray& answer, const IntArray& dofMask, 
                     ValueModeType mode, TimeStep* stepN)
{
  if (!hasSlaveDofs) {

    int j, size;
    IntArray dofArray;
    Dof        *dofJ ;
    
    answer.resize (size = dofMask.giveSize());
    this-> giveDofArray(dofMask, dofArray);
    
    for (j=1 ; j<=size ; j++) {
      dofJ = this->giveDof(dofArray.at(j)) ;
      
      //if (dofJ -> hasBc(stepN) && (dofJ ->giveUnknownType() == type))
      if (dofJ -> hasBc(stepN))
        answer.at(j) = dofJ->giveBcValue(mode,stepN) ;//giveUnknown(u,stepN) ;
      // answer.at(j) = dofJ->giveBcValue(type, mode,stepN) ;//giveUnknown(u,stepN) ;
      else
        answer.at(j) = 0. ;
    }
  } else {

    int i, k, indx;
    IntArray dofArray;
    FloatArray mstrBcVlus;
    Dof *dofJ;
    
    this-> giveDofArray(dofMask, dofArray);
    answer.resize (giveNumberOfPrimaryMasterDofs(dofArray));
    
    for (k=1,i=1; i<=dofArray.giveSize(); i++) {
      indx = dofArray.at(i);
      if (!this->giveDof(indx)->isPrimaryDof()) { // slave DOF
        this->giveDof(indx)->giveBcValues (mstrBcVlus, mode, stepN);
        answer.copySubVector (mstrBcVlus, k);
        k += mstrBcVlus.giveSize();
      }
      else { // primary DOF
        dofJ = this-> giveDof(indx);
        if (dofJ-> hasBc(stepN))
          answer.at(k++) = dofJ-> giveBcValue(mode, stepN);
        else
          answer.at(k++) = 0.0;
      }
    }
  }
}

int
DofManager::hasAnySlaveDofs() 
{
  int i;
  for (i=1; i<=numberOfDofs; i++)
    if (!this->giveDof (i)->isPrimaryDof()) return 1;
  return 0;
}


int
DofManager :: checkConsistency ()
  // Checks internal data consistency in node. 
  // Current implementation checks (when receiver has simple slave dofs) if receiver
  // has the same coordinate system as master dofManager of slave dof.
{
  int i;
  
  hasSlaveDofs = false;
  for (i=1; i<=numberOfDofs; i++) {
    if (this->giveDof(i)->giveClassID() == SlaveDofClass) { 
      hasSlaveDofs = true;
      continue;
    }
  }
  
  return 1;
}


void
DofManager :: computeDofTransformation (FloatMatrix& answer, const IntArray* dofMask, DofManTrasfType mode)
  // computes trasformation matrix of receiver.
  // transformation should include trasformation from global cs to nodal cs,
  // as well as further necessary transformations (for example in case 
  // rigid arms this must include transformation to master dofs).
{
  if (!hasSlaveDofs) {
    int _size = (dofMask == NULL) ? numberOfDofs : dofMask->giveSize();
    answer.resize (_size, _size);
    answer.beUnitMatrix();
  } else {
    this->computeSlaveDofTransformation (answer, dofMask, mode);
  }
}

void
DofManager :: computeLoadTransformation (FloatMatrix& answer, const IntArray* dofMask, DofManTrasfType mode)
  // computes trasformation matrix of receiver.
  // transformation should include trasformation from global cs to nodal cs,
  // as well as further necessary transformations (for example in case 
  // rigid arms this must include transformation to master dofs).
{
  if (mode != _toNodalCS) _error ("computeSlaveLoadTransformation: unsupported mode");
  
  FloatMatrix t;
  
  computeDofTransformation (t, dofMask, _toGlobalCS);
  answer.beTranspositionOf (t);

  /*
  if (!hasSlaveDofs) {
    int _size = (dofMask == NULL) ? numberOfDofs : dofMask->giveSize();
    answer.resize(_size, _size);
    answer.beUnitMatrix ();
  } else {
    this->computeSlaveLoadTransformation (answer, dofMask, mode);
  }
  */
}

void
DofManager :: computeSlaveDofTransformation (FloatMatrix& answer, const IntArray* dofMask, DofManTrasfType mode)
  // computes trasformation matrix of receiver.
  // transformation should include trasformation from global cs to nodal cs,
  // as well as further necessary transformations (for example in case 
  // rigid arms this must include transformation to master dofs).
{
  if (mode != _toGlobalCS) _error ("computeSlaveDofTransformation: unknown mode");
  
  int i, k, indx;
  IntArray dofArray;
  FloatArray mstrContrs;
  
  if (dofMask==NULL) {
    dofArray.resize(numberOfDofs);
    for (i=1; i<=numberOfDofs; i++)  dofArray.at(i) = i;
  }
  else
    this->giveDofArray (*dofMask, dofArray);
  
  answer.resize (dofArray.giveSize(), giveNumberOfPrimaryMasterDofs(dofArray));
  answer.zero ();
  
  for (k=1,i=1; i<=dofArray.giveSize(); i++) {
    indx = dofArray.at(i);
    if (!this->giveDof(indx)->isPrimaryDof()) { // slave DOF
      this->giveDof(indx)->computeDofTransformation (mstrContrs);
      answer.copySubVectorRow (mstrContrs, i, k);
      k += mstrContrs.giveSize();
    }
    else { // primary DOF
      answer.at(i,k++) = 1.0;
    }
  }
}

IntArray*
DofManager :: giveCompleteGlobalDofIDArray (void) const
{
  IntArray* answer = new IntArray(numberOfDofs);
  
  for (int i=1; i<=numberOfDofs; i++)
    answer->at(i) = (int)this->giveDof(i)->giveDofID();
  
  return answer;
}



#ifdef __PARALLEL_MODE
int
DofManager::packDOFsUnknowns (CommunicationBuffer& buff, EquationID type, 
               ValueModeType mode, TimeStep* stepN)
{
 int i, result = 1;
 for (i=1; i<=numberOfDofs; i++)
  result &= this->giveDof (i)->packUnknowns (buff, type, mode, stepN);
 return result;
}

bool 
DofManager :: isLocal () {
  if (parallel_mode == DofManager_local) return true;
  if (parallel_mode == DofManager_shared) {
    // determine if problem is the lowest one sharing the dofman; if yes the receiver is responsible to 
    // deliver number
    int n = partitions.giveSize();
    int myrank = this->giveDomain()->giveEngngModel()->giveRank();
    int minrank = myrank;
    
    for (int j=1; j<=n; j++) minrank = min (minrank, partitions.at(j));
    if (minrank == myrank) return true;
  }
  return false;
}



#endif
