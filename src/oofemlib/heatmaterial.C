/* $Header: /home/cvs/bp/oofem/oofemlib/src/heatmaterial.C,v 1.4 2003/04/06 14:08:24 bp Exp $ */
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


// file heatMaterial.C

#include "heatmaterial.h"
#include "crosssection.h"
#include "domain.h"
#include "verbose.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "cltypes.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <math.h>
#endif


void
HeatMaterial :: giveCharacteristicMatrix (FloatMatrix& answer,MatResponseForm form, 
       MatResponseMode rMode,
       GaussPoint* gp,
       TimeStep* atTime)
  //
// Returns characteristic StructuralMaterial stiffness matrix of the receiver
//
{
 MaterialMode mMode = gp->giveMaterialMode();
 switch(mMode) {
 case _2dHeat:
  this->give2dHeatConductivityCharMtrx (answer, form,rMode,gp,atTime);
  break;
 default:
   _error2 ("giveCharacteristicMatrix : unknown mode (%s)",  __MaterialModeToString(mMode));
  return ;
 }
 return ;
}

int
HeatMaterial :: hasMaterialModeCapability (MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
 if (mode == _2dHeat)
  return  1;
 return 0;
}
