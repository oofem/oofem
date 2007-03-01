/* $Header: /home/cvs/bp/oofem/sm/src/deformationtheorymaterial.h,v 1.3 2003/04/06 14:08:30 bp Exp $ */
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

//   *******************************************
//   *** CLASS MATERIAL - DEFORMATION THEORY ***
//   *******************************************

#ifndef deftheorymaterial_h 

#include "femcmpnn.h"
#include "cltypes.h"
#include "structuralmaterial.h"

class DeformationTheoryMaterial : public StructuralMaterial
{
/*
   This class implements an abstract class material, which behaves
 according to deformation theory.

 DESCRIPTION
 The attribute 'propertyDictionary' possibly contains all the properties of a mate-
   rial, like its Young modulus, its mass density or poisson ratio.
 TASK
 - Returning standard material stiffness and flexibility marices for 3d-case.
   according to current state determined by using data stored 
   in Gausspoint.
 - Returning standard material stiffness for other stress states for point
   in 3d continua - (2dPlanaStress, 2dPlaneStrain, 1dStress);
 - Returning a material property (method 'give'). Only for non-standard elements.
 - Returning real stress state vector(tensor) at gauss point for 3d - case.
 - Imposing constrains according to stressStrain mode in gp to 3dstiffmessMatrix
   (function reduceTo).
 - storing / restoring (possible) Status stored (if defined) in gp->matStatusDict
     with key = (int)obj->giveClassID(), where
   obj is instance of material, yield crit which is associated with this gp.
   if no material status is necessary material doesn'n derive it own instance of
   material Status class and also then saving & restoring context is 
   function with no code.
 */
protected:
 
public:
 
 DeformationTheoryMaterial (int n,Domain* d) : StructuralMaterial(n,d) {}
 ~DeformationTheoryMaterial ()  { }
 
 // identification and auxiliary functions
 virtual int hasNonLinearBehaviour ()   { return 1 ;}
 const char* giveClassName ()  const     { return "DeformationTheoryMaterial" ;}
 classType giveClassID ()          const     {return DeformationTheoryMaterialClass;}

} ;


#define deftheorymaterial_h
#endif
