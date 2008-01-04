/* $Header: /home/cvs/bp/oofem/oofemlib/src/heatcrosssection.h,v 1.8 2003/04/06 14:08:24 bp Exp $ */
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

//   *********************************************
//   *** CLASS CROSSSECTION FOR HEAT TRANSFER ****
//   *********************************************

#ifndef heatcrosssection_h 

#include "femcmpnn.h"
#include "crosssection.h"
#include "material.h"
#include "gausspnt.h"
#include "element.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"

class HeatCrossSection : public CrossSection
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
   HeatCrossSection (int n,Domain* d) : CrossSection(n,d) {}
   ~HeatCrossSection ()                {}
 

   virtual void
     giveCharMaterialConductivityMatrix  (FloatMatrix& answer, MatResponseMode rMode, GaussPoint*, 
       TimeStep*tStep) ;
   double   give (int) ;
   
 // identification and auxiliary functions
 const char*    giveClassName () const      {return "HeatCrossSection" ;}
 classType giveClassID ()         const      {return HeatCrossSectionClass;}
 IRResultType initializeFrom (InputRecord* ir);
  //int hasHeatCapability       () {return 1;}
  int testCrossSectionExtension (CrossSectExtension ext) {return ((ext==CS_HeatCapability)?1:0);}
   
 
 
protected:
} ;


#define heatcrosssection_h
#endif

