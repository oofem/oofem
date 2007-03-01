/* $Header: /home/cvs/bp/oofem/oofemlib/src/heatmaterial.h,v 1.7 2003/04/06 14:08:24 bp Exp $ */
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


//   ***********************************************
//   *** CLASS LINEAR  MATERIAL for HEAT TRANSFER***
//   ***********************************************

#ifndef HeatMaterial_h 

  
#include "material.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"

class GaussPoint ;


class HeatMaterial : public Material
{
/*
   This class implements a linear heat material in a finite element problem. A material
   is an attribute of a domain. It is usually also attribute of many elements.

 DESCRIPTION
 
 TASK
 - Returning standard material conductivity  marices for 3d-case.
   according to current state determined by using data stored 
   in Gausspoint.
 - Returning a material property (method 'give'). Only for non-standard elements.
 - Returning real stress state vector(tensor) at gauss point for 3d - case.
*/
public:
 
  HeatMaterial (int n,Domain* d) : Material (n,d) {}
  ~HeatMaterial () {}
 
  // standart matrial stiffness matrices
  virtual void      giveCharacteristicMatrix (FloatMatrix& answer,
                                              MatResponseForm form,
                                              MatResponseMode mode,
                                              GaussPoint* gp,
                                              TimeStep* atTime);

  virtual void  give2dHeatConductivityCharMtrx (FloatMatrix& answer,
                                                MatResponseForm form,
                                                MatResponseMode mode,
                                                GaussPoint* gp,
                                                TimeStep* atTime) 
    {_error ("give2dHeatConductivityCharMtrx: not implemented yet");}
 
  // identification and auxiliary functions
  int hasNonLinearBehaviour ()   { return 0 ;}
  int testMaterialExtension      (MaterialExtension ext) {return ((ext == Material_HeatCapability)?1:0);}
  //int hasHeatCapability       () { return 1;}
  virtual int hasMaterialModeCapability (MaterialMode ) ;
  const char*    giveClassName () const      { return "HeatMaterial" ;}
  classType giveClassID ()         const      {return HeatMaterialClass;}
 
protected:
  
} ;


#define HeatMaterial_h
#endif
