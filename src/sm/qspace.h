/* $Header: /home/cvs/bp/oofem/sm/src/qspace.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   ************************************
//   *** CLASS QUADRATIC 3D Element   ***
//   ************************************
 

#include "lspace.h"

class QSpace  : public StructuralElement
{
/*
   This class implements an Quadratic 3d  20 - node 
   elasticity finite element. Each node has 3 degrees of freedom.
 DESCRIPTION :
   One single additional attribute is needed for Gauss integration purpose :
   'jacobianMatrix'. This 3x3 matrix contains polynomials.
 TASKS :
   - calculating its Gauss points ;
   - calculating its B,D,N matrices and dV.
*/

protected:
  int           numberOfGaussPoints;

public:
 QSpace (int,Domain*) ;                          // constructor
 ~QSpace ()  {}                                  // destructor

 void    computeJacobianMatrixAt (FloatMatrix& answer, FloatArray*) ;
 virtual int            computeNumberOfDofs (EquationID ut) {return 60;}
 virtual void           giveDofManDofIDMask  (int inode, EquationID, IntArray& ) const;
  double             computeVolumeAround (GaussPoint*) ;
 virtual int computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) ;
 double        giveCharacteristicLenght (GaussPoint* gp, const FloatArray &normalToCrackPlane);

// 
// definition & identification
//
  const char* giveClassName () const { return "QSpace" ;}
  classType            giveClassID () const           { return QSpaceClass; } 
 IRResultType initializeFrom (InputRecord* ir);
  Interface* giveInterface (InterfaceType) {return NULL;}

protected:
 void          computeBmatrixAt (GaussPoint*, FloatMatrix&, int=1,int=ALL_STRAINS) ;
  void          computeNmatrixAt (GaussPoint*, FloatMatrix&) ;
  void       computeGaussPoints () ;
 integrationDomain  giveIntegrationDomain () {return _Cube;}

} ;

   








