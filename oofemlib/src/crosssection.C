/* $Header: /home/cvs/bp/oofem/oofemlib/src/crosssection.C,v 1.10.4.1 2004/04/05 15:19:43 bp Exp $ */
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

//   file CROSSSECTION.CC

#include "crosssection.h"
#include "simplecrosssection.h"
//#include "layeredcrosssection.h"
#include "structuralelement.h"
#include "heatcrosssection.h"
#include "emptycs.h"
#include "gausspnt.h"
#include "material.h"
#include "flotarry.h"
#include "verbose.h"
#include "usrdefsub.h"
#include "compiler.h"

#ifndef __MAKEDEPEND
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#endif

IRResultType 
CrossSection::initializeFrom (InputRecord* ir)
//
// instanciates receiver from input record
//
{

 return IRRT_OK;
}


void  
CrossSection :: printYourself ()
   // Prints the receiver on screen.
{
   printf ("Cross Section with properties : \n") ;
   propertyDictionary -> printYourself() ;
}


contextIOResultType
CrossSection :: saveContext (DataStream* stream, ContextMode mode, void *obj)
//
// saves full material context (saves state variables, that completely describe
// current state)
//
{
 contextIOResultType iores;
 GaussPoint *gp = (GaussPoint*) obj;
 Material *mat = gp->giveMaterial();

 if ((iores = FEMComponent::saveContext(stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = mat->saveContext (stream, mode, (void*) gp)) != CIO_OK) THROW_CIOERR(iores);

 return CIO_OK;
}


contextIOResultType
CrossSection :: restoreContext (DataStream* stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
 contextIOResultType iores;
 GaussPoint *gp = (GaussPoint*) obj;
 Material *mat = gp->giveMaterial();

 if ((iores = FEMComponent::restoreContext(stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = mat->restoreContext (stream, mode, (void*) gp)) != CIO_OK) THROW_CIOERR(iores);

 return CIO_OK;
}


CrossSection* 
CrossSection :: ofType (char* aClass)
   // Returns a new cross section, which has the same number than the receiver,
   // but belongs to aClass (simpleCrossSection, LayeredCrossSection, FibredCS, ...).
{
   CrossSection* newCrossSection ;

   if (! strncasecmp(aClass,"simplecs",8))
     newCrossSection = new SimpleCrossSection (this->giveNumber(),domain) ;
   // else if (! strncasecmp(aClass,"layeredcs",14))
   //  newCrossSection = new LayeredCrossSection (number,domain) ;
   else if (! strncasecmp(aClass,"heatcs",14))
     newCrossSection = new HeatCrossSection (number,domain) ;
   /*   else if (! strncasecmp(aClass,"steel1",6))
        newMaterial = new Steel1 (number,domain) ;*/
   else if (! strncasecmp(aClass,"emptycs",7))
     newCrossSection = new EmptyCS (number,domain) ;
   else {
     // last resort - call aditional user defined subroutine
     newCrossSection = ::CreateUsrDefCrossSectionOfType (aClass,number,domain);
     if (newCrossSection == NULL) {
       _error ("ofType:  unknown cross section type \n") ;
       exit(0) ;
     }
   }
   return newCrossSection ;
}


double  
CrossSection :: give (int aProperty)
   // Returns the value of the property aProperty (e.g. the area
   // 'A') of the receiver.
{
 
   if (propertyDictionary -> includes(aProperty))
   return propertyDictionary -> at(aProperty) ;
   else {                                         
     _error ("give: property not defined");
   }
   return 0.0 ;
}


bool 
CrossSection::isCharacteristicMtrxSymmetric (MatResponseMode rMode, int mat) {
  return domain->giveMaterial(mat)->isCharacteristicMtrxSymmetric(rMode);
}
