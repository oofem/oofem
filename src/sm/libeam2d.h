/* $Header: /home/cvs/bp/oofem/sm/src/libeam2d.h,v 1.6 2003/04/06 14:08:30 bp Exp $ */
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

//   **********************
//   *** CLASS LIBeam2d ***
//   **********************

#ifndef libeam2d_h
#define libeam2d_h

#include "structuralelement.h"
#include "layeredcrosssection.h"

class LIBeam2d : public  StructuralElement, public LayeredCrossSectionInterface
{
/*
   This class implements a 2-dimensional Linear Isoparametric
   Mindlin theory beam element, with reduced integration.
*/
  
public:
  double pitch, length;


    LIBeam2d (int,Domain*) ;                       // constructor
    ~LIBeam2d ()  {}                               // destructor

  // FloatMatrix*  ComputeConstitutiveMatrixAt (GaussPoint*) ;
  // FloatArray*   ComputeResultingBodyForceAt (TimeStep*) ;
  void          computeLumpedMassMatrix (FloatMatrix& answer, TimeStep* tStep) ;
  void          computeMassMatrix (FloatMatrix& answer, TimeStep* tStep) 
  {computeLumpedMassMatrix (answer, tStep);}
  void          computeStiffnessMatrix (FloatMatrix& answer, 
                    MatResponseMode rMode, TimeStep* tStep) ;
  int           computeGtoLRotationMatrix (FloatMatrix&);  // giveRotationMatrix () ;
//
// layered cross section support functions
//
 void          computeStrainVectorInLayer (FloatArray& answer, GaussPoint* masterGp, 
                      GaussPoint* slaveGp, TimeStep* tStep);
 /**
  Computes the global coordinates from given element's local coordinates.
  @returns nonzero if successful
  */
 virtual int computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) ;

  /** Interface requesting service */
  Interface* giveInterface (InterfaceType);

 virtual int            computeNumberOfDofs (EquationID ut) {return 6;}
 virtual void           giveDofManDofIDMask  (int inode, EquationID, IntArray& ) const;
  double        computeVolumeAround (GaussPoint*) ;

// 
// definition & identification
//
  const char* giveClassName () const { return "LIBeam2d" ;}
  classType             giveClassID ()          const { return LIBeam2dClass; } 
  IRResultType initializeFrom (InputRecord* ir);

protected:
  // edge load support
  void  computeEgdeNMatrixAt (FloatMatrix& answer, GaussPoint*) ;
  void  giveEdgeDofMapping (IntArray& answer, int) const;
  double        computeEdgeVolumeAround (GaussPoint*, int);
  void          computeEdgeIpGlobalCoords (FloatArray& answer, GaussPoint* gp, int iEdge) 
  {computeGlobalCoordinates(answer, *(gp->giveCoordinates()));}
  int   computeLoadLEToLRotationMatrix (FloatMatrix&, int, GaussPoint*) ;
 int  computeLoadGToLRotationMtrx (FloatMatrix& answer) ;
  //void          computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint*, TimeStep*, ValueModeType mode);
  void          computeBmatrixAt (GaussPoint*, FloatMatrix&, int=1, int=ALL_STRAINS) ;
  void          computeNmatrixAt (GaussPoint*, FloatMatrix &) ;
  void          computeGaussPoints () ;
 integrationDomain  giveIntegrationDomain () {return _Line;}
  double        giveLength () ;
  double        givePitch () ;

} ;

#endif // libeam2d_h
