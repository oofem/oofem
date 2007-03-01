/* $Header: /home/cvs/bp/oofem/oofemlib/src/isolinearheatmat.h,v 1.9 2003/04/06 14:08:24 bp Exp $ */
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


//   ************************************************
//   *** CLASS ISOTROPIC MATERIAL FOR HEAT 
//   ************************************************

#ifndef isotropiclinearheatmat_h 

#include "heatmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"

class GaussPoint ;


class IsotropicLinearHeatMaterial : public  HeatMaterial
{
/*
   This class implements a isotropic linear heat  material in a finite 
 element problem. A material
   is an attribute of a domain. It is usually also attribute of many elements.

 DESCRIPTION
   ISOTROPIC Linear Heat Material
 
 TASK

*/
protected:
 
  double k, c;

public:
 
  IsotropicLinearHeatMaterial (int n,Domain* d) : HeatMaterial (n,d) {}
  ~IsotropicLinearHeatMaterial () {}
  
  // identification and auxiliary functions
  const char*    giveClassName () const {return "IsotropicLinearHeatMaterial" ;}
  classType giveClassID ()         const {return IsotropicLinearHeatMaterialClass;}
  IRResultType initializeFrom (InputRecord* ir);
  
  // non-standart - returns time independent material constant
  double   give (int) ;
  
  virtual FloatMatrix* Give2dHeatConductivityCharMtrx (MatResponseForm form,
                                                       MatResponseMode mode,
                                                       GaussPoint* gp,
                                                       FloatArray* strainIncrement,
                                                       TimeStep* atTime);

  virtual MaterialStatus* CreateStatus (GaussPoint* gp) const {return NULL;}
protected:
} ;


#define isotropiclinearheatmat_h
#endif
