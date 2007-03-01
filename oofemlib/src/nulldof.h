/* $Header: /home/cvs/bp/oofem/oofemlib/src/nulldof.h,v 1.10 2003/04/06 14:08:25 bp Exp $ */
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



//   **********************
//   *** CLASS NULL DOF ***
//   **********************
 

#ifndef nulldof_h

#include "dof.h"
#include "compiler.h"
#include "dictionr.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#endif
#include "cltypes.h"
#include "error.h"

class Domain ; class DofManager ; class TimeStep ; class BoundaryCondition ; 
class InitialCondition ;

/**
 Class representing "null" degree of freedom. NULL dof always have zero valued boundary 
 condition active (==> no equation number) and always return value prescribed by its initial conditions.
 If no initial condition is applied, zero value is returned when asked for its unknown.
 Almost all operation are empty functions.

 */
class NullDof : public Dof
{
/*
   This class implements a "null" degree of freedom.
 A dof is usually attribute of one node.
 DESCRIPTION

  Null dof is special dof - with no equation number ever assigned (always zero valued boundary condition applied). 
 Also unknown dictionary is not necessary - its values are always prescribed by initial condition.

*/

   protected:

    /// Initial condition number associated to dof.
    int      ic;
   public:
   /**
    Constructor. Creates null dof vith number i, belonging to aNode dof manager.
    @param i dof number
    @param aNode receiver will belong to aNode dof manager
    @param id DofID of master dof (and slave too).
   */
      NullDof (int i, DofManager* aNode, int nic, DofID id) : Dof (i, aNode, id) {ic=nic;}    // constructor
    /// Destructor.
      ~NullDof ()   { }    // destructor.

    /// Returns class name of the receiver.
      const char*    giveClassName () const  {return "NullDof";}
    /// Returns classType id of receiver.
      classType           giveClassID () const {return NullDofClass;}
    /**
     Returns equation number corresponding to receiver. 
     @return equation number - returns zero.
     */
      int                 giveEquationNumber ()  {return 0;}
    /**
     Returns prescribed equation number corresponding to receiver. 
     @return prescribed equation number - returns zero.
     */
      int                 givePrescribedEquationNumber ()  {return 0;}
    /**
     Asks new equation number. Empty function.
     */
      int                 askNewEquationNumber (TimeStep* tStep) {return 1;}
    /**
     Returns the value of the unknown associated with the receiver
     at given time step. Always returns value prescribed by initial condition.
     @see MasterDof::giveUnknown function
    */
    double              giveUnknown (EquationID , ValueModeType, TimeStep*)  ;
    /**
     Returns the value of the unknown of the receiver
     at given time step associated to given field.
     Always returns value prescribed by initial condition.
     */
    double              giveUnknown (PrimaryField& field, ValueModeType, TimeStep* stepN);
    /**
     Returns boundary condition of dof if it is precsribed.
     @return always return nonzero indicating active BC.
     */
    int                 hasBc (TimeStep* tStep) {return 1;}
    /**
     Test if Dof has initial condition.
     @return nonzero if IC exists, zero otherwise.
    */
    int                 hasIc () ;
   /**
    Test if Dof has initial condition of required ValueModeType.
    @param u type of required IC
    @return nonzero if IC exists, zero otherwise.
    @see ValueModeType.
    */
   int                 hasIcOn (ValueModeType) ;
    /** Returns the id of associated boundary condition, if there is any.
     Used only for printing purposes. In general, id culd not be used
     to decide whether bc is active. Use appropriate services instead.
     @param id of associated Boubdaray condition, zero otherwise
     */
    int giveBcIdValue () {return 0;}

   /// Stores receiver state to output stream.  
    contextIOResultType    saveContext (FILE* stream, void *obj = NULL)  {return CIO_OK;}
   /// Restores the receiver state previously written in stream.
    contextIOResultType    restoreContext(FILE* stream, void *obj = NULL) {return CIO_OK;}

protected:
    /**
     Returns boundary condition of dof if it is precsribed. 
     Slave simply forwards this mesage to master. 
     @return returns NULL if no BC applied, otherwise pointer to correcpondig BC.
     */
    BoundaryCondition*  giveBc () {return NULL;}
    /**
     Returns initial condition of dof if it is precsribed. 
     Slave forwards this message to master.
     @return returns NULL if no IC applied, otherwise pointer to correcpondig IC.
     */
    InitialCondition*   giveIc ();


} ;

#define nulldof_h
#endif
