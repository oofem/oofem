/* $Header: /home/cvs/bp/oofem/sm/src/ltrelemppde.h,v 1.4 2003/04/06 14:08:30 bp Exp $ */
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

//   **************************
//   *** CLASS Linear Triangular Element for Parabolic Partial Differential eqs.
//   **************************
 
#ifndef ltrelelemppde_h

#include "ppdeelement.h"
#include "gaussintegrationrule.h"

class LTrElementPPDE : public PPdeElement
{
/*
   This class implements an triangular three-node linear
   partial diferential equation (like heat transfer) finite element. 
 Each node has 1 degree of freedom.
 DESCRIPTION :

 TASKS :

   - calculating its B,D,N matrices and dV.
*/

protected:
  double area;
  int numberOfGaussPoints;

public:
  LTrElementPPDE (int,Domain*) ;                          // constructor
  ~LTrElementPPDE ()  {}                                  // destructor
  
  void               printOutputAt (FILE*, TimeStep*) ;

  virtual int            computeNumberOfDofs (EquationID ut) {return 3;}
  virtual void           giveDofManDofIDMask  (int inode, EquationID, IntArray& ) const;
  double             computeVolumeAround (GaussPoint*) ;

#ifdef __OOFEG
      void          drawRawGeometry (oofegGraphicContext&);
#endif
// 
// definition & identification
//
  const char* giveClassName () const { return "LTrElementPPDE" ;}
  classType             giveClassID () const { return LTrElementPPDEClass; } 
  IRResultType initializeFrom (InputRecord* ir);

protected: 
 void               computeBmatrixAt (GaussPoint*, FloatMatrix&) ;
  void               computeNmatrixAt (GaussPoint*, FloatMatrix& );
  void               computeGaussPoints () ;
 integrationDomain  giveIntegrationDomain () {return _Triangle;}

  virtual double     giveArea ();
  virtual FloatArray*        GivebCoeff ();
  virtual FloatArray*        GivecCoeff ();
} ;

#define ltrelelemppde_h
#endif







