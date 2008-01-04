/* $Header: /home/cvs/bp/oofem/oofemlib/src/isolinearheatmat.C,v 1.5 2003/04/06 14:08:24 bp Exp $ */
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

//   file isolinearheatmaterial.C

#include "isolinearheatmat.h"
#include "domain.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#endif

IRResultType
IsotropicLinearHeatMaterial :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro
  // double value ;
 
 this -> Material::initializeFrom(ir);

 IR_GIVE_FIELD (ir, k, IFT_IsotropicLinearHeatMaterial_k, "k"); // Macro// conductivity
 IR_GIVE_FIELD (ir, c, IFT_IsotropicLinearHeatMaterial_c, "c"); // Macro// specific heat

 return IRRT_OK;
}


double
IsotropicLinearHeatMaterial :: give (int aProperty)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
 if ((aProperty == 'k')) return k;
 if ((aProperty == 'c')) return c;
 return this -> Material :: give(aProperty);
}



FloatMatrix* 
IsotropicLinearHeatMaterial :: Give2dHeatConductivityCharMtrx (MatResponseForm form,
                     MatResponseMode mode,
                     GaussPoint* gp,
                     FloatArray* strainIncrement,
                     TimeStep* atTime)
{
/*
  returns constitutive matrix of receiver
*/
  FloatMatrix* constitutiveMatrix ;

  if (form == FullForm) {
  constitutiveMatrix = new FloatMatrix(3,3);
  
  constitutiveMatrix->at(1,1) = k;
  constitutiveMatrix->at(2,2) = k;
  } else {
  constitutiveMatrix = new FloatMatrix(2,2);
  
  constitutiveMatrix->at(1,1) = k;
  constitutiveMatrix->at(2,2) = k;
 }  

  return constitutiveMatrix ;
}
