/* $Header: /home/cvs/bp/oofem/oofemlib/src/rigidarmnode.h,v 1.12.4.1 2004/04/05 15:19:43 bp Exp $ */
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


#ifndef rigidarmnode_h

#include "node.h"
#include "domain.h"
#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

class Dof ; class NodalLoad ; class TimeStep ;
class FloatArray ; class IntArray ;

/**
 Class implementing node connected to other node (master) using rigid arm in finite element mesh. 
 Rigid arm node posses no degrees of freedom - all dofs are mapped to master dofs.
 ++ added on August 28, the Rigid arm node supports not only slave dofs mapped to master
 but also some dofs can be primary doofs. Introduced masterDofMask allowing to
 distinguish between primary and mapped (slave) dofs. The primary DOFs can have their own BCs, ICs.

 The introduction of rigid arm connected nodes allows to avoid wery stiff elements used
 for modelling the rigid-arm connection. The rigid arm node maps its dofs to master dofs
 using simple transformations (small rotations are assumed). Therefore, the contribution
 to rigid arm node are localized directly to master related equations.
 The rigid arm node can not have its own boundary or initial conditions,
 they are determined completely from master dof conditions. 
 The local coordinate system in slave is not supported in current implementation, the global lcs applies.
 On the other hand, rigid arm node can be loaded independently of master.
 The transformation for DOFs and load is not ortogonal - the inverse transformation can 
 not be constructed by transposition. Because of time consuming inversion, methods 
 can generally compute both transformations for dofs as well as loads.
*/
class RigidArmNode : public Node
{


private:
 /// number of  master dof (Master DofManager)
 int          masterDofMngr;
  /// mask distinguishing between primary and mapped (slave) dofs (nonzero value)
 IntArray slaveDofMask;

public:

 /**
  Constructor. Creates a rigid-arm node belonging to domain.
  @param n node number in domain aDomain
  @param aDomain domain to which node belongs
  */
 RigidArmNode (int n,Domain* aDomain) ;                        // constructor
  /// Destructor.
 ~RigidArmNode () ;                                            // destructor


   /** Returns location array (array containing for each requested dof related equation number)
  of receiver. Equation number of dof with active boundary condition is zero.
  @param dofIDArry array containing dof mask. This mask containing DofIDItem values 
  (they describe physical meaning of dofs, see cltypes.h) is used to extract only 
  required values. If dof with requested physical meaning dos not exist in receiver,
  an error is generated and execution exits.
  @param locationArray - return parameter containing required equation numbers.
  @see Element::giveDofManDofIDMask function.
  */
 void         giveLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const;
   /** Returns location array of prescribed equations 
   (array containing for each requested dof related equation number)
   of receiver. Equation number of dof with no or inactive boundary condition is zero.
   @param dofIDArry array containing dof mask. This mask containing DofIDItem values 
   (they describe physical meaning of dofs, see cltypes.h) is used to extract only 
   required values. If dof with requested physical meaning dos not exist in receiver,
   an error is generated and execution exits.
   @param locationArray - return parameter containing required equation numbers.
   @see Element::giveDofManDofIDMask function.
  */
 void         givePrescribedLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const;
  /** Returns full location array of receiver containing equation numbers of all dofs 
  of receiver. Their order is specific to every DofManager. Mainly used at EngngModel level 
  to assemble DofManager contribution (typically load vector).
  @param locationArray complete location array of receiver.
  */
 void         giveCompleteLocationArray (IntArray& locationArray) const;
  /** Returns full location array of prescribed equations of receiver containing equation numbers of all dofs 
  of receiver. Their order is specific to every DofManager. Mainly used at EngngModel level 
  to assemble DofManager contribution (typically load vector).
  @param locationArray complete location array of receiver.
  */
 void         giveCompletePrescribedLocationArray (IntArray& locationArray) const;
  /**
  Assembles the vector of unknowns in nodal c.s for given dofs of receiver.
  This vector may have size different from number of dofs requested,
  because some dofs may depend on other dofs. 
  This implementation returns the vector of all unknowns of master node.
  @param answer result (in nodal cs.)
  @param dofMask  dofIDArry array containing dof mask. This mask containing DofIDItem values 
  (they describe physical meaning of dofs, see cltypes.h) is used to extract only 
  required values. If dof with requested physical meaning dos not exist in receiver,
  an error is generated and execution exits.
  @param type physical meaning of  unknown.
  @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
  @stepN time step when unknown requested. See documentation of particular EngngModel 
  class for valid StepN values (most implementaion can return only values for current 
  and possibly for previous time step).
  @see Dof::giveUnknown service.
  */
    void giveUnknownVector (FloatArray& answer, const IntArray& dofMask, 
              EquationID type, ValueModeType mode, TimeStep* stepN);
  /**
  Assembles the vector of unknowns of given filed in nodal c.s for given dofs of receiver.
  This vector may have size different from number of dofs requested,
  because some dofs may depend on other dofs. Default implementation uses 
  Dof::giveUnknown service.
  @param answer result (in nodal cs.)
  @param dofMask  dofIDArry array containing dof mask. This mask containing DofIDItem values 
  (they describe physical meaning of dofs, see cltypes.h) is used to extract only 
  required values. If dof with requested physical meaning dos not exist in receiver,
  an error is generated and execution exits.
  @param field primary filed
  @stepN time step when unknown requested. See documentation of particular EngngModel 
  class for valid StepN values (most implementaion can return only values for current 
  and possibly for previous time step).
  @see Dof::giveUnknown service.
  */
   virtual void giveUnknownVector (FloatArray& answer, const IntArray& dofMask, 
                  PrimaryField& field, ValueModeType mode, TimeStep* stepN);
  /**
  Assembles the vector of prescribed unknowns in nodal c.s for given dofs of receiver.
  This vector may have size different from number of dofs requested,
  because some dofs may depend on other dofs. 
  This implementation returns vector of all prescribed master DOFs.
  @param answer result (in nodal cs.)
  @param dofMask  dofIDArry array containing dof mask. This mask containing DofIDItem values 
  (they describe physical meaning of dofs, see cltypes.h) is used to extract only 
  required values. If dof with requested physical meaning dos not exist in receiver,
  an error is generated and execution exits.
  @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
  @stepN time step when unknown requested. See documentation of particular EngngModel 
  class for valid StepN values (most implementaion can return only values for current 
  and possibly for previous time step).
  @see Dof::giveBcValue, Dof::hasBc service.
  */
    void givePrescribedUnknownVector (FloatArray& answer, const IntArray& dofMask, 
                   ValueModeType mode, TimeStep* stepN);
    //void givePrescribedUnknownVector (FloatArray& answer, IntArray& dofMask, 
  //                 UnknownType type, ValueModeType mode, TimeStep* stepN);
  /** Computes the load vector of receiver in given time.
  @param answer load vector.
  @param stepN time step when answer is computed.
  @param mode determines response mode.
 */
 void         computeLoadVectorAt (FloatArray& answer, TimeStep* stepN, ValueModeType mode)  ;

 // local coordinate system
  /** 
  Returns nonzero if node has prescribed  local coordinate system. 
  Current implementation does not support lcs in slave.
  */
 int          hasLocalCS () { return domain->giveNode (masterDofMngr)->hasLocalCS();}
  /** Returns pointer to local coordinate triplet in node. 
  If not defined, returns NULL. Simply forwards this message to master.
  @return Triplet defining the local coordinate system in node.
  Value at position (i,j) represents angle between e'(i) and e(j),
  where e' is base vector of local coordinate system and e is
  base vector of global c.s.
  */
 FloatMatrix* giveLocalCoordinateTriplet () {return domain->giveNode (masterDofMngr)->giveLocalCoordinateTriplet();}
  /** Computes receiver DOF transformation matrix from global cs. to dofManager specific
  coordinate system - if mode == _toNodalCS, otherwise reverse transformation is computed.
  (in which governing equations are assembled, for example the
  local coordinate system in node). This transformation is not ortogonal.
  @param answer computed transformation matrix. It has generally dofIDArry.size rows and
  if loc is obtained using giveLocationArray(dofIDArry, loc) call, loc.giveSize() columns.
  This is because this transformation should generally include not only transformation to 
  dof manager local coordinate system, but receiver dofs can be expressed using 
  dofs of another dofManager (In this case, squre answer is produced anly if all 
  dof transformation is required).
  @param dofIDArry array containing DofIDItem-type values (this is enumeration 
  identifying physical meaning of particular DOF, see cltypes.h) for which transfromation mtrx is
  assembled. if dofIDArry is NULL, then all receiver dofs are assumed.
  */
  void computeDofTransformation (FloatMatrix& answer, const IntArray* dofIDArry, DofManTrasfType mode);
  /** Computes receiver LOAD transformation matrix from global cs. to dofManager specific
  coordinate system - if mode == _toNodalCS, otherwise reverse transformation is computed.
  (In the dofManager specifis cs the governing equations are assembled, for example the
  local coordinate system in node). This transformation may not be ortogonal.
  @param answer computed transformation matrix. It has generally dofIDArry.size rows and
  if loc is obtained using giveLocationArray(dofIDArry, loc) call, loc.giveSize() columns.
  This is because this transformation should generally include not only transformation to 
  dof manager local coordinate system, but receiver dofs can be expressed using 
  dofs of another dofManager (In this case, squre answer is produced anly if all 
  dof transformation is required).
  @param dofIDArry array containing DofIDItem-type values (this is enumeration 
  identifying physical meaning of particular DOF, see cltypes.h) for which transfromation mtrx is
  assembled. if dofIDArry is NULL, then all receiver dofs are assumed.
  */
  void computeLoadTransformation (FloatMatrix& answer, const IntArray* dofIDArry, DofManTrasfType mode);
  /**
  Indicates, whether dofManager requires the transformation from global c.s. to 
  dof manager specific coordinate system. 
  @return nonzero - the transformation is necessary.
  */
  int requiresTransformation () {return 1;}
 // time step termination
  /** Updates receiver at end of time step (i.e. after equilibrium has been reached).
  If EngngModel formulation ( see giveFormulation() member function) returns actualized
  Lagrange mode, node updates its coordinates according to solution. Implementation is empty 
  function - it is sufficient to update master.
  @see EngngModel::giveFormulation().
  */
 void         updateYourself (TimeStep*) {}
  /// Returns true if receiver contains slave dofs
 virtual int  hasSlaveDofs () {if (masterDofMngr) return 1; else return 0;}

 // miscellaneous
  /// Returns class name of the receiver.
 const char* giveClassName () const      { return "RigidArmNode" ;}
 /** Returns classType id of receiver.
  @see FEMComponent::giveClassID 
  */
 classType    giveClassID () const {return RigidArmNodeClass;}
 ///Initializes receiver acording to object description stored in input record.
 IRResultType initializeFrom (InputRecord* ir);
  //virtual IntArray* ResolveDofIDArray (char* initString);
  /**
  Checks internal data consistency in node. 
  Current implementation checks (when receiver has slave dofs) if receiver has the same 
  coordinate system as master dofManager of slave dof.
  @return nonzero if receiver check is o.k.
  */
  virtual int    checkConsistency () ;

 /// returns reference to master dof. Public because RigidArmSlaveDof need to access.
  Node*    giveMasterDofMngr () const {return domain->giveNode (masterDofMngr);}
 /// returns master node number. 
  int      giveMasterDofMngrNum () const {return (masterDofMngr);}
  /** returns value of slaveDofMask for given dof
  @param i dof number
  @return nonzero if dof is linked (slave), zero othervise
  */
 int      giveSlaveDofMask (int i) const {return slaveDofMask.at(i);}
  /// returns slaveDofMask of receiver
  void     giveSlaveDofMask (IntArray& mask) const { mask = slaveDofMask;}
  // Returns total number of dofs managed by receiver 
  //int          giveNumberOfDofs () const {return this->giveMasterDofMngr()->giveNumberOfDofs();}
  /**
  Returns updated ic-th coordinate of receiver. Return value is computed 
  as coordinate + scale * displacement, where corresponding displacement is obtained
  from cooresponding nodal DOF. Local coodinate system is taken into account.
  Usefull mainly for postprocessing.
  */
 double       giveUpdatedCoordinate(int ic, TimeStep* tStep, 
                    EquationID type, double scale = 1.);

} ;


#define rigidarmnode_h
#endif






