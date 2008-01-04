/* $Header: /home/cvs/bp/oofem/sm/src/microplanematerial_bazant.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   *************************************************************************
//   *** ABSTRACT CLASS MICROPLANE MATERIAL ACCORDING TO BAZANT'S APPROACH ***
//   *************************************************************************

#ifndef microplanematerial_bazant_h

#include "structuralms.h"
#include "microplanematerial.h"


/**
 Abstract base class for all microplane models according to Bazant's approach. 
 Micro strains on microplane are described using magnitude of normal strain,
 volumetric normal component and by two orthogonal shear components 
 (m and l direction) in microplane.

*/
class MicroplaneMaterial_Bazant : public MicroplaneMaterial
{

protected:

public:
 /**
  Constructor. Creates Abstract Bazant's Microplane Material belonging 
  to domain d, with number n.
  @param n material number
  @param d domain to which newly created material belongs
  */
  MicroplaneMaterial_Bazant (int n,Domain* d) ;
  /// Destructor.
  ~MicroplaneMaterial_Bazant ()                {}
 

/**
 Computes real macro stress in corresponding macro integration point for
 given macroscopic strain  vector.
 */
virtual void giveRealStressVector (FloatArray& answer, MatResponseForm, GaussPoint*, 
                  const FloatArray& , TimeStep*) ;



/**
   Updates the volumetric stress component after computing real stress microplane vectors.
*/
virtual  void  updateVolumetricStressTo (Microplane* mPlane, double sigv) = 0;

/// Returns class name of the receiver.
 const char*    giveClassName () const { return "MicroplaneMaterial_Bazant" ;}
/// Returns classType id of receiver.
 classType giveClassID ()         const {return MicroplaneMaterial_BazantClass;}

 virtual MaterialStatus* CreateStatus (GaussPoint* gp) const {return new StructuralMaterialStatus (1, domain,gp);}

protected:

} ;

#define microplanematerial_bazant_h
#endif
