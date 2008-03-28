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

#ifndef hangingnode_h
#define hangingnode_h

#include "node.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif


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
   
   !!!  Pri pouziti slavemask je treba si to dobre rozmyslet, aby nevznikaly nejake nespojitosti.
   !!!  Uzlove zatizeni je treba si dobre rozmyslet vzdy.
   
*/
class HangingNode : public Node
{
 protected:
  /**
     type of interpolation from hangingnodes = 100*(number of master nodes) + 10*(order of polynomial approximation) + dimension
     211 - linear truss         321 - quadratic truss
     312 - linear triangle      622 - quadratic triangle
     412 - linear rectangle     822 - quadratic rectangle
     413 - linear tetrahedron  1023 - quadratic tetrahedron
     813 - linear hexahedron   2023 - quadratic hexahedron
  */
  int typeOfContrib;
  /// natural(triangular)~local coordinates
  FloatArray* locoords;
  
  /// count of Master Nodes
  int countOfMasterNodes;
  /// array of numbers of master DofManagers
  IntArray* masterDofMngr;
  /// array of pointers to master Nodes
  Node** masterNode;
  /// vector of master contribution coefficients - SUMA of contributions == 1.0 ; if node has LCS, contribs are specified in this local coordinate system.
  FloatArray* masterContribution;
  
 private:
  void allocAuxArrays (void);
  void deallocAuxArrays (void);
  
 public:
  /**
     Constructor. Creates a hanging node with number n, belonging to aDomain.
     @param n node number in domain aDomain
     @param aDomain domain to which node belongs
  */
  HangingNode (int n, Domain* aDomain);
  /**
     Destructor.
  */
  ~HangingNode (void) {}
  
  /**
     Initializes receiver acording to object description stored in input record.
  */
  IRResultType initializeFrom (InputRecord* ir);
  /**
     Checks internal data consistency in node. 
     @return nonzero if receiver check is o.k.
  */
  int checkConsistency ();
  /**
     compute vector of master contribution coefficients - SUMA of contributions == 1.0
  */
  int computeMasterContribution ();
  
  /**
     Returns class name of the receiver.
  */
  const char* giveClassName () const { return "HangingNode" ;}
  /**
     Returns classType id of receiver.
     @see FEMComponent::giveClassID 
  */
  classType giveClassID () const {return HangingNodeClass;}
  /// Returns true if dof of given type is allowed to be associated to receiver
  bool isDofTypeCompatible (dofType type) const {return (type == DT_master || type == DT_slave);} // termitovo
  
};


#endif // hangingnode_h
