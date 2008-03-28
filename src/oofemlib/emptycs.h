/* $Header: /home/cvs/bp/oofem/oofemlib/src/crosssection.h,v 1.11 2003/04/06 14:08:23 bp Exp $ */
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


//   ***************************
//   *** CLASS CROSSSECTOION ***
//   ***************************

#ifndef emptycs_h 
#define emptycs_h 

#include "crosssection.h"
#include "material.h"
//#include "perfectlyplasticmaterial.h"
#include "gausspnt.h"
//#include "element.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"

/**
   Empty cross section model, passes all requests to material driver
 */
class EmptyCS : public CrossSection
{
/*
   This abstract class implements a cross section in a finite element problem. A cross
 section  is an attribute of a domain. It is usually also attribute of many 
 elements. This class is base class for SimpleCrossSection, LayeredCrossSection,
 FibredCrossSection, ....

 DESCRIPTION
   The attribute 'propertyDictionary' contains all the properties of a 
 cross section, like its area or thickness.

 TASK
 - Returning a properties of cross section like thickness or area.
*/
protected:
public:
 /**
  Constructor. Creates cross section with number n belonging to domain d.
  @param n cross section number
  @param d domain
  */
  EmptyCS (int n,Domain* d) ;
 /// Destructor.
  ~EmptyCS ()  ;

 // identification and auxiliary functions
 /// Returns class name of the receiver.
 const char*    giveClassName () const      {return "EmptyCS" ;}
 /// Returns classType id of receiver.
 classType giveClassID ()         const      {return EmptyCrossSectionClass;}
} ;


#endif // emptycs_h 

