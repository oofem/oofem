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
/*
  Contributed by Ladislav Svoboda
*/

#ifndef hangingnode_h
#define hangingnode_h

#include "node.h"
#include "domain.h"
#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
class Dof ; class NodalLoad ; class TimeStep ;
class FloatArray ; class IntArray ;

/**
   Class implementing hanging node connected to other nodes (masters) using interpolation.
   Hanging node posses no degrees of freedom	- all values are interpolated from corresponding master dofs.
   
   The introduction of hanging nodes allows, for example, to include reinforcing bar elements inside 
   arbitrary fe mesh of concrete specimen or facilitates the local refinment of fe-mesh.

   The contributions of hanging node are localized directly to master related equations.
   The hanging node can not have its own boundary or initial conditions,
   they are determined completely from master dof conditions. 
   The local coordinate system in slave is not supported in current implementation, the global lcs applies.
   On the other hand, hanging node can be loaded independently of master.

   To do: Implement evaluation of natural coordinates using Interpaolation classes, instead of
   using local formulas or supplying then on input.

*/
class HangingNode : public Node
{
  
 private:
  /**
     type of interpolation from hangingnodes = 100*(number of master nodes) + 10*(order of polynomial approximation) + dimension
     211 - linear truss         321 - quadratic truss
     312 - linear triangle      622 - quadratic triangle
     412 - linear rectangle     822 - quadratic rectangle
     413 - linear tetrahedron  1023 - quadratic tetrahedron
     813 - linear hexahedron   2023 - quadratic hexahedron
  */
  int type;
  /// count of Master DofManagers
  int countOfMasterDofMngr;
  /// array of Master DofManager numbers
  IntArray masterDofMngr;
  /// natural(triangular)~local coordinates
  FloatArray locoords;
  //  double ksi,eta,dzeta;
  
  
 public:
  /**
     Constructor. Creates a hanging node belonging to domain.
     @param n node number in domain aDomain
     @param aDomain domain to which node belongs
  */
  HangingNode (int n,Domain* aDomain) ;
  /**
     Destructor.
  */
  ~HangingNode () ;
  /**
     Initializes receiver acording to object description stored in input record.
  */
  IRResultType initializeFrom (InputRecord* ir);
  /**
     Computes the natural coordinates of receiver from global coordinates of master and receiver, 
     taking into account selected interpolation order.
   */
  void compute_naturalcoord ();
  /**
     Checks internal data consistency in node. 
     Current implementation checks (when receiver has slave dofs) if receiver has the same 
     coordinate system as master dofManager of slave dof.
     @return nonzero if receiver check is o.k.
  */
  int checkConsistency () ;
  /**
     Assembles the vector of unknowns in nodal c.s for given dofs of receiver.
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
     @see DofManager::giveUnknownVector
  */
  void giveUnknownVector (FloatArray& answer, const IntArray& dofMask, EquationID type, ValueModeType mode, TimeStep* stepN);
  /**
     Assembles the vector of prescribed unknowns in nodal c.s for given dofs of receiver.
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
  void givePrescribedUnknownVector (FloatArray& answer, const IntArray& dofMask, ValueModeType mode, TimeStep* stepN);
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
  /** Computes the load vector of receiver in given time.
      @param answer load vector.
      @param stepN time step when answer is computed.
      @param mode determines response mode.
  */
  void computeLoadVectorAt (FloatArray& answer, TimeStep* stepN, ValueModeType mode)  ;
  /** Returns location array (array containing for each requested dof related equation number)
      of receiver. Equation number of dof with active boundary condition is zero.
      @param dofIDArry array containing dof mask. This mask containing DofIDItem values 
      (they describe physical meaning of dofs, see cltypes.h)	is used to extract only 
      required values. If dof with requested physical meaning dos not exist in receiver,
      an error is generated and execution exits.
      @param locationArray - return parameter containing required equation numbers.
      @see Element::giveDofManDofIDMask function.
  */
  void giveLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const;
  /** Returns location array of prescribed equations 
      (array containing for each requested dof related equation number)
      of receiver. Equation number of dof with no or inactive boundary condition is zero.
      @param dofIDArry array containing dof mask. This mask containing DofIDItem values 
      (they describe physical meaning of dofs, see cltypes.h)	is used to extract only 
      required values. If dof with requested physical meaning dos not exist in receiver,
      an error is generated and execution exits.
      @param locationArray - return parameter containing required equation numbers.
      @see Element::giveDofManDofIDMask function.
  */
  void givePrescribedLocationArray (const IntArray& dofIDArry, IntArray& locationArray) const;
  /** Returns full location array of receiver containing equation numbers of all dofs 
      of receiver. Their order is specific to every DofManager. Mainly used at EngngModel level 
      to assemble DofManager contribution (typically load vector).
      @param locationArray complete location array of receiver.
  */
  void giveCompleteLocationArray (IntArray& locationArray) const;
  /** Returns full location array of prescribed equations of receiver containing equation numbers of all dofs 
      of receiver. Their order is specific to every DofManager. Mainly used at EngngModel level 
      to assemble DofManager contribution (typically load vector).
      @param locationArray complete location array of receiver.
  */
  void giveCompletePrescribedLocationArray (IntArray& locationArray) const;
  /** 
      Returns nonzero if node has prescribed  local coordinate system. 
      Current implementation does not support lcs in slave.
  */
  int hasLocalCS (void) const { return 0 ; }


  /**
     Indicates, whether dofManager requires the transformation from global c.s. to 
     dof manager specific coordinate system. 
     @return nonzero - the transformation is necessary.
  */
  int requiresTransformation () {return 1;}

  // time step termination
  /**
     Updates receiver at end of time step (i.e. after equilibrium has been reached).
     If EngngModel formulation ( see giveFormulation() member function) returns actualized
     Lagrange mode, node updates its coordinates according to solution. Implementation is empty 
     function - it is sufficient to update master.
     @see EngngModel::giveFormulation().
  */
  //void         updateYourself (TimeStep*) {}                                              ????????????????????????????

  /**
     Returns class name of the receiver.
  */
  const char* giveClassName () const { return "HangingNode" ;}
  /** 
      Returns classType id of receiver.
      @see FEMComponent::giveClassID 
  */
  classType giveClassID () const {return HangingNodeClass;}

  /**
     returns reference to master dof. Public because RigidArmSlaveDof need to access.
  */
  Node* giveMasterDofMngr (int i) const {return domain->giveNode (masterDofMngr.at(i));}
  /**
     returns count of Master DofManagers
   */
  int giveCountOfMasterDofMngr (void) const { return countOfMasterDofMngr ; }

} ;


#endif
