/* $Header: /home/cvs/bp/oofem/oofemlib/src/heatcrosssection.C,v 1.6 2003/04/06 14:08:24 bp Exp $ */
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


// file: heatcrosssection.C

#include "heatcrosssection.h"
#include "gausspnt.h"
#include "heatmaterial.h"
#include "flotarry.h"
#include "verbose.h"

void
HeatCrossSection::giveCharMaterialConductivityMatrix  (FloatMatrix& answer, 
             MatResponseMode rMode, GaussPoint* gp, 
             TimeStep*tStep) 
{
  HeatMaterial* mat = (HeatMaterial*) gp->giveElement()->giveMaterial();
  mat->giveCharacteristicMatrix (answer, ReducedForm, rMode, gp,tStep) ;
 return;
}


double  
HeatCrossSection :: give (int aProperty)
   // Returns the value of the property aProperty (e.g. the area
   // 'A') of the receiver.
{
   double  value = 0.0;
 
 if (aProperty == 'A') return this->give(THICKNESS)*this->give(WIDTH);
 if (aProperty == 't') return this->give(THICKNESS);
   if (propertyDictionary -> includes(aProperty))
      value = propertyDictionary -> at(aProperty) ;
   else {                                         
//      value = this -> read(aProperty) ;
//      propertyDictionary -> add(aProperty,value) ;}
     _error ("give: property not defined");
   }
   return value ;
}




IRResultType
HeatCrossSection :: initializeFrom (InputRecord* ir)
//
// instanciates receiver from input record
//
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro
 double value;

 this -> CrossSection :: initializeFrom(ir);

 IR_GIVE_FIELD (ir, value, IFT_HeatCrossSection_thick, "thick"); // Macro
 propertyDictionary -> add(THICKNESS,value);

 IR_GIVE_FIELD (ir, value, IFT_HeatCrossSection_width, "width"); // Macro
 propertyDictionary -> add(WIDTH,value);

 return IRRT_OK;
} 
