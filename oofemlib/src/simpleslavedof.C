/* $Header: /home/cvs/bp/oofem/oofemlib/src/slavedof.C,v 1.8.4.1 2004/05/14 13:45:27 bp Exp $ */
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

//   file SLAVEDOF.CC

#include "simpleslavedof.h"
#include "dofmanager.h"
#include "domain.h"
#include "timestep.h"
#include "boundary.h"
#include "initial.h"

#include "flotarry.h"
#include "dictionr.h"

#include "debug.h"
#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#endif


SimpleSlaveDof :: SimpleSlaveDof (int i, DofManager* aNode, int master, DofID id) : Dof (i, aNode, id)
   // Constructor. Creates a new d.o.f., with number i, belonging
   // to aNode with bc=nbc, ic=nic
{
  masterDofMngr = master;
  masterDofIndx = -1;
/*   unknowns       = new Dictionary() ;           // unknown size ?
   pastUnknowns   = NULL ; */
}


Dof* SimpleSlaveDof::giveMasterDof () 
{
  // returns reference to master dof
  // checks dof compatibility and slave to slave references

  if (this->masterDofIndx == -1) {
  this->masterDofIndx = dofManager->giveDomain()->giveDofManager(masterDofMngr)
  ->findDofWithDofId(this->dofID);

  if (this->masterDofIndx) {

  classType masterDofCT = dofManager->giveDomain()->giveDofManager(masterDofMngr)->
    giveDof(masterDofIndx)->giveClassID();
  if ((masterDofCT != MasterDofClass)&&(masterDofCT != SharedMasterDofClass)&&(masterDofCT != RemoteMasterDofClass)) {
    
    _error ("giveMasterDof: slaveDof to slaveDof reference not allowed");
  }
  } else {
  _error ("giveMasterDof: no dof with dofID in master found");
  }
  }

  return dofManager->giveDomain()->giveDofManager(masterDofMngr)->
  giveDof(masterDofIndx);
}

BoundaryCondition*  SimpleSlaveDof :: giveBc () 
   // Returns the boundary condition the receiver is subjected to.
{
   return  this->giveMasterDof()->giveBc();
}


int  SimpleSlaveDof :: giveEquationNumber ()
   // Returns the number of the equation in the governing system of equations that corres-
   // ponds to the receiver. The equation number is 0 if the receiver is
   // subjected to a boundary condition, else it is n+1, where n is the
   // equation number of the most recently numbered degree of freedom.
{

  return  this->giveMasterDof()->giveEquationNumber();
}

int  SimpleSlaveDof :: givePrescribedEquationNumber ()
   // Returns the number of the equation in the governing system of equations that corres-
   // ponds to the receiver. The equation number is 0 if the receiver is
   // subjected to a boundary condition, else it is n+1, where n is the
   // equation number of the most recently numbered degree of freedom.
{

  return  this->giveMasterDof()->givePrescribedEquationNumber();
}

InitialCondition*  SimpleSlaveDof :: giveIc () 
   // Returns the initial condition on the receiver. Not used.
{
  return this->giveMasterDof()->giveIc ();
}


double  SimpleSlaveDof :: giveUnknown (EquationID type, ValueModeType mode, TimeStep* stepN)
   // The key method of class Dof. Returns the value of the unknown 'u'
   // (e.g., the displacement) of the receiver, at stepN. This value may,
   // or may not, be already available. It may depend on a boundary (if it
   // is not a predicted unknown) or initial condition. stepN is not the
   // current time step n, it is assumed to be the previous one (n-1).
{
  return this->giveMasterDof()->giveUnknown(type, mode, stepN);
}

double  SimpleSlaveDof :: giveUnknown (PrimaryField& field, ValueModeType mode, TimeStep* stepN)
{
  return this->giveMasterDof()->giveUnknown(field, mode, stepN);
}


int  SimpleSlaveDof :: hasBc (TimeStep* tStep) 
   // Returns True if the receiver is subjected to a boundary condition, else
   // returns False. If necessary, reads the answer in the data file.
{
  return this->giveMasterDof()->hasBc (tStep);
}


int  SimpleSlaveDof :: hasIc () 
   // Returns True if the receiver is subjected to an initial condition,
   // else returns False.
{
  return this->giveMasterDof()->hasIc ();
}


int  SimpleSlaveDof :: hasIcOn (ValueModeType u)
   // Returns True if the unknown 'u' (e.g., the displacement 'd') of the
   // receiver is subjected to an initial condition, else returns False.
{
  return this->giveMasterDof()->hasIcOn (u);
}

int SimpleSlaveDof :: giveBcIdValue ()
{
 return this->giveMasterDof()->giveBcIdValue ();
}


contextIOResultType SimpleSlaveDof :: saveContext (DataStream* stream, ContextMode mode, void *obj) 
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
  return Dof::saveContext(stream, mode, obj);
}


contextIOResultType SimpleSlaveDof :: restoreContext (DataStream* stream, ContextMode mode, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
  return Dof::restoreContext (stream, mode, obj);
}
