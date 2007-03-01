/* $Header: /home/cvs/bp/oofem/oofemlib/src/nulldof.C,v 1.6 2003/04/06 14:08:25 bp Exp $ */
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


// file: nulldof.C

#include "nulldof.h"
#include "dofmanager.h"
#include "domain.h"
#include "timestep.h"
#include "boundary.h"
#include "initial.h"
#include "debug.h"
#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <ctype.h>
#endif


InitialCondition*  NullDof :: giveIc () 
   // Returns the initial condition on the receiver. Not used.
{
   return  (InitialCondition*) (dofManager -> giveDomain() -> giveLoad(ic)) ;
}


double  NullDof :: giveUnknown (EquationID type, ValueModeType mode, TimeStep* stepN) 
   // The key method of class Dof. Returns the value of the unknown 'u'
   // (e.g., the displacement) of the receiver, at stepN. This value may,
   // or may not, be already available. It may depend on a boundary (if it
   // is not a predicted unknown) or initial condition. stepN is not the
   // current time step n, it is assumed to be the previous one (n-1).
{
 double value;

#ifdef DEBUG
// if (type != this->giveUnknownType ())
//  _error ("giveUnknown: Noncompatible Request");  
#endif 

 // return always IC value
 if (this->hasIcOn (mode))
  value = this -> giveIc() -> give(mode) ;
 else
  value = 0.;
   
 return value;
}


double  NullDof :: giveUnknown (PrimaryField& field, ValueModeType mode, TimeStep* stepN)
{
 double value;

#ifdef DEBUG
// if (type != this->giveUnknownType ())
//  _error ("giveUnknown: Noncompatible Request");  
#endif 

 // return always IC value
 if (this->hasIcOn (mode))
  value = this -> giveIc() -> give(mode) ;
 else
  value = 0.;
 
 return value;
}


int  NullDof :: hasIc ()
   // Returns True if the receiver is subjected to an initial condition,
   // else returns False.
{
   if (ic == -1) {
   _error ("hasIc:  does not know yet if has InitCond or not \n") ;
   exit(0) ;}

   return ic ;
}

int  NullDof :: hasIcOn (ValueModeType u) 
   // Returns True if the unknown 'u' (e.g., the displacement 'd') of the
   // receiver is subjected to an initial condition, else returns False.
{
   if (this->hasIc())
      return this->giveIc()->hasConditionOn(u) ;
   else
      return FALSE ;
}



