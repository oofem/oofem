/* $Header: /home/cvs/bp/oofem/oofemlib/src/remotemasterdof.h,v 1.7 2003/04/06 14:08:25 bp Exp $ */
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

 
//   *******************************
//   *** CLASS REMOTE MASTER DOF ***
//   *******************************
 

#ifndef remotemasterdof_h
#ifdef __PARALLEL_MODE

#include "masterdof.h"
#include "compiler.h"
#include "dictionr.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#endif
#include "cltypes.h"
#include "error.h"

#include "combuff.h"

class Domain ; class DofManager ; class TimeStep ; class BoundaryCondition ; 
class InitialCondition ;

/**
 Class representing remote "master" degree of freedom. Master is degree of freedom, which has 
 its related unknown and corresponding equation number.
 Remote master dof is an exact copy of remote dof.
 Remote degrees of freedom are having local equation numbers.
 */
class RemoteMasterDof  : public MasterDof
{
   private:
   public:
 /**
  Constructor. Creates master dof with number i, belonging to DofManager aNode and with
  physical meaning described by id.
  @param i DOF number.
  @param aNode DofManager which possess DOF.
  @param nbc number of associated boundary condition, zero if none.
  @param nic number of associated initial condition, zero if none.
  @param id Physical meaning type.
  @see cltypes.h, DofID type
  */
      RemoteMasterDof (int i, DofManager* aNode,int nbc,int nic,DofID id):
    MasterDof (i, aNode, nbc, nic, id) {}    // constructor
  /// Destructor.
   ~RemoteMasterDof ()   {}    // destructor.
  /// Returns class name of the receiver.
      const char*    giveClassName () const  {return "RemoteMasterDof";}
  /// Returns classType id of receiver.
      classType           giveClassID () const {return RemoteMasterDofClass;}


  /**
  Unpacks DOF unknown from communication buffer and updates unknown if necessary.
  Engng model unknown of receiver is updated 
  to accomodate remote contribution or "prescribed" remote values.
  The unknown dictionary is not updated, it is engng model job to update 
  all unknowns dictionaries.
  @param buff buffer containing packed message
  @param type physical meaning of  unknown.
  @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
  @stepN time step when unknown requested. See documentation of particular EngngModel 
  class for valid StepN values (most implementaion can return only values for current 
  and possibly for previous time step).
  @return nonzero if succesfull
  @see Dof::unpackAndUpdateDOFsUnknown for description
  */
  virtual int unpackAndUpdateUnknown (CommunicationBuffer& buff, EquationID type, 
                   ValueModeType mode, TimeStep* stepN);
//#endif

protected:
} ;

#define masterdof_h
#endif
#endif
