/* $Header: /home/cvs/bp/oofem/tm/src/isoheatmat.C,v 1.2 2003/04/18 10:38:54 bp Exp $ */
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

#include "isoheatmat.h"
#include "domain.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#endif

IRResultType
IsotropicHeatTransferMaterial :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro
  // double value ;
 
 this -> Material::initializeFrom(ir);

 IR_GIVE_FIELD (ir, k, IFT_IsotropicHeatTransferMaterial_k, "k"); // Macro// conductivity
 IR_GIVE_FIELD (ir, c, IFT_IsotropicHeatTransferMaterial_c, "c"); // Macro// specific heat

 return IRRT_OK;
}


double
IsotropicHeatTransferMaterial :: give (int aProperty)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
 if ((aProperty == 'k')) return k;
 else if ((aProperty == HeatCapaCoeff) || (aProperty == 'c')) return (c * this->give('d'));
 return this -> Material :: give(aProperty);
}



void
IsotropicHeatTransferMaterial :: giveCharacteristicMatrix (FloatMatrix& answer,
                                                           MatResponseForm form,
                                                           MatResponseMode mode,
                                                           GaussPoint* gp,
                                                           TimeStep* atTime)
{
/*
  returns constitutive matrix of receiver
*/
 MaterialMode mMode = gp->giveMaterialMode();
 switch  (mMode) {
 case _2dHeat:
  answer.resize (2,2);
  answer.at(1,1) = k;
  answer.at(2,2) = k;
  return;
 case _3dHeat:
  answer.resize (3,3);
  answer.at(1,1) = k;
  answer.at(2,2) = k;
  answer.at(3,3) = k;
  return;
 default:
  _error ("giveCharacteristicMatrix : unknown mode");
 }
}


double
IsotropicHeatTransferMaterial:: giveCharacteristicValue (MatResponseMode mode,
                                                         GaussPoint* gp,
                                                         TimeStep* atTime)
{
  if (mode == Capacity) return (c * this->give('d'));
  else
    _error ("giveCharacteristicValue : unknown mode");
  return 0; // to make compiler happy
}
