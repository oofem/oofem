/* $Header: /home/cvs/bp/oofem/oofemlib/src/remotemasterdof.C,v 1.6 2003/04/06 14:08:25 bp Exp $ */
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


//   file REMOTEMASTERDOF.CC
 
#ifdef __PARALLEL_MODE


#include "remotemasterdof.h"
#include "dofmanager.h"
#include "domain.h"
#include "timestep.h"
#include "boundary.h"
#include "initial.h"
#include "engngm.h"

#include "flotarry.h"
#include "dictionr.h"

#include "debug.h"
#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#endif

int 
RemoteMasterDof :: unpackAndUpdateUnknown (CommunicationBuffer& buff, EquationID type, 
                      ValueModeType mode, TimeStep* stepN)
{
 int result;
 double value;

 result = buff.unpackDouble (value);

 // if dof belonging to remote DofManager, engng model unknowns are updated 
 // to accomodate remote contribution or "prescribed" remote values.
 // The unknown dictionary is not updated, it is engng model job to update 
 // all unknowns dictionaries.
 dofManager->giveDomain()->giveEngngModel()->
  updateUnknownComponent(type, mode, stepN, this->giveEquationNumber(), 
              value, EngngModel::EngngModel_SET_Mode);

 return result;
}

#endif




