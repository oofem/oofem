/* $Header: /home/cvs/bp/oofem/oofemlib/src/structengngmodel.h,v 1.12.4.1 2004/04/05 15:19:44 bp Exp $ */
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
// Class LinearStatic
//

#ifndef structengngmodel_h
#define structengngmodel_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

#include "engngm.h"

#ifdef __PARALLEL_MODE
#include "problemcomm.h"
#include "processcomm.h"
#endif

class StructuralElement;

/**
   This class implements extension of EngngModel for structural models.
  Its purpose is to declare and implement general methods for computing reaction forces.
*/
class StructuralEngngModel : public EngngModel
{ 
  /*
    This class implements extension of EngngModel for structural models.
    Its purpose is to declare and implement general methods for computing reaction forces.
    
  */

 protected:
  /**
     Contains last time stamp of internal variable update.
     This update is made via various services 
     (like those for computing real internal forces or updating the internal state).
  */
  StateCounterType internalVarUpdateStamp;
  
#ifdef __PARALLEL_MODE
  /// Communicator mode. Determines current strategy used.
  ProblemCommunicator::ProblemCommunicatorMode commMode;
  /// Common Communicator buffer
  CommunicatorBuff* commBuff;
  /// Communicator. 
  ProblemCommunicator* communicator;

  /// Flag indicating if nonlocal extension active
  int nonlocalExt;
  /// NonLocal Communicator. Necessary when nonlocal constitutive models are used.
  ProblemCommunicator* nonlocCommunicator;
#endif

/**
  Computes and prints reaction forces, computed from nodal internal forces. Assumes, that real 
  stresses corresponding to reached state are already computed (uses giveInternalForcesVector 
  structural element service with useUpdatedGpRecord = 1 parameter). Only the dof managers selected for
  output (OutputManager) are handled.
  @see StructuralElement::giveInternalForcesVector
  @see OutputManager
  @param tStep time step
  @param id domain number
 */
void printReactionForces (TimeStep *, int id);
/**
  Buids the reaction force table. For each prescribed equation number it will find
  corresponding node and dof number. The entries in the restrDofMans, restrDofs, and eqn
  arrays are sorted with incresing dofman number and with increasing dof number as
  a second minor criterion.
  @param restrDofMans contains numbers of restrained Dofmanagers, with size equal to total number of prescribed equations.
  @param restrDofs contains numbers of restrained Dofs, with size equal to total number of prescribed equations.
  @param eqn contains the corresponding restrained equation numbers.
  @param tStep time step
   @param di domain number
*/
void buildReactionTable (IntArray& restrDofMans, IntArray& restrDofs, IntArray& eqn, TimeStep* tStep, int di);



/**
 Computes the contribution of internal element forces to reaction forces in given domain.
 @param reactions contains the comuted contributions
 @param tStep solution step
 @param domian number
*/
void computeInternalForceReactionContribution (FloatArray& reactions, TimeStep* tStep, int di);
/**
 Computes the contribution external loading to reaction forces in given domain. Default implementations adds the 
 contibution from computeElementLoadReactionContribution and computeElementLoadReactionContribution methods.
 @param reactions contains the comuted contributions
 @param tStep solution step
 @param domian number
*/
virtual void computeExternalLoadReactionContribution (FloatArray& reactions, TimeStep* tStep, int di);
/**
 Computes the contribution of element load to reaction forces in given domain.
 @param reactions contains the comuted contributions
 @param tStep solution step
 @param domian number
*/
void computeElementLoadReactionContribution (FloatArray& reactions, TimeStep* tStep, int di);
/**
 Computes the contribution of nodal load to reaction forces in given domain.
 @param reactions contains the comuted contributions
 @param tStep solution step
 @param domian number
*/
void computeNodalLoadReactionContribution (FloatArray& reactions, TimeStep* tStep, int di);


/*
 Compute element equivForces in dof managers (nodes and sides).
 It substracts part corresponding to non-nodal loading from internal forces vector.
 @param answer returned equivalent forces
 @param tStep time step
 @param ielem element number
 @param mode ValueModeType (TotalMode should be used)
void computeElementEquivForces (FloatArray& answer, TimeStep *tStep, StructuralElement* ielem, ValueModeType mode) const;
*/
 /**
  Updates nodal values
  (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
  if model supports changes of static system). The element internal state update is also forced using
  updateInternalState service.
  */
 void updateInternalState (TimeStep *);
public:
 /// Constructor - creates new StructuralEngngModel with number i, associated to domain d.
  StructuralEngngModel (int i, EngngModel* _master = NULL) : EngngModel (i, _master)
    {internalVarUpdateStamp = 0;}
  /// Destructor
  ~StructuralEngngModel () ;

 // identification
 /// Returns "StructuralEngngModel" - class name of the receiver.
  const char* giveClassName () const { return "StructuralEngngModel";}
 /// Returns StructuralEngngModelClass - classType id of receiver.
  classType giveClassID ()      const { return StructuralEngngModelClass;}

#ifdef __PARALLEL_MODE
  /**
  Packing function for internal forces of DofManagers. Pascks internal forces of shared DofManagers
  into send communication buffer of given process communicator.
  @param processComm task communicator for which to pack forces
  @src source vector
  @return nonzero if successfull.
   */
  int packInternalForces (FloatArray* src, ProcessCommunicator& processComm);
  /**
  Unpacking function for internal forces of DofManagers . Unpacks internal forces of shared DofManagers
  from  receive communication buffer of given process communicator.
  @param processComm task communicator for which to unpack forces
  @dest destination vector
  @return nonzero if successfull.
  */
  int unpackInternalForces (FloatArray* dest, ProcessCommunicator& processComm);
  /**
  Packing function for  reactions of DofManagers. Pascks reactions of shared DofManagers
  into send communication buffer of given process communicator.
  @param processComm task communicator for which to pack forces
  @src source vector
  @return nonzero if successfull.
   */
  int packReactions (FloatArray* src, ProcessCommunicator& processComm);
  /**
  Unpacking function for rections  of DofManagers . Unpacks rections of shared DofManagers
  from  receive communication buffer of given process communicator.
  @param processComm task communicator for which to unpack forces
  @dest destination vector
  @return nonzero if successfull.
  */
  int unpackReactions (FloatArray* dest, ProcessCommunicator& processComm);
  /**
  Packing function for load vector. Packs load vector values of shared/remote DofManagers
  into send communication buffer of given process communicator.
  @param processComm task communicator for which to pack load
  @src source vector
  @return nonzero if successfull.
   */
  int packLoad (FloatArray* src, ProcessCommunicator& processComm);
  /**
  Unpacking function for load vector values of DofManagers . Unpacks load vector of shared/remote DofManagers
  from  receive communication buffer of given process communicator.
  @param processComm task communicator for which to unpack load
  @dest destination vector
  @return nonzero if successfull.
  */
  int unpackLoad (FloatArray* dest, ProcessCommunicator& processComm);
  /**
  Packs data of local element to be received by their remote counterpart on remote partitions.
  Remote elements are introduced when nonlocal constitutive models are used, in order to 
  allow local averaging procedure (remote elements, which are involved in averaging on local partition are
  mirrored on this local partition) instead of implementing inefficient fine-grain communication.
  Remote element data are exchanged only if necessary and once for all of them.
  Current implementationn calls packUnknowns service for all elements listed in 
  given process communicator send map.
  @param processComm corresponding process communicator.
  @return nonzero if successfull.
  */
  int packRemoteElementData (ProcessCommunicator& processComm);
  /**
  Unpacks data for remote eleemnts (which are mirrors of remote partition's local elements).
  Remote elements are introduced when nonlocal constitutive models are used, in order to 
  allow local averaging procedure (remote elements, which are involved in averaging on local partition are
  mirrored on this local partition) instead of implementing inefficient fine-grain communication.
  Remote element data are exchanged only if necessary and once for all of them.
  Current implementation calls unpackAndUpdateUnknowns service for all elements listed in
  given process communicator receive map.
  @param processComm corresponding process communicator.
  @return nonzero if successfull.
  */
  int unpackRemoteElementData (ProcessCommunicator& processComm);

#endif
#ifdef __PETSC_MODULE
   /**
      Creates Petsc contexts. Must be implemented by derived classes since the governing equation type is reqired 
      for context creation.
    */
  virtual void initPetscContexts ();
#endif

#ifdef __OOFEG   
  /**
  Shows the sparse structure of required matrix, type == 1 stiffness.
  */
  void               showSparseMtrxStructure (int type, oofegGraphicContext& context, TimeStep* atTime);
#endif
};

#endif // structengngmodel_h
