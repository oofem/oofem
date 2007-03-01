/* $Header: /home/cvs/bp/oofem/sm/src/truss3d.h,v 1.7 2003/04/06 14:08:32 bp Exp $ */
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

//   ********************* 
//   *** CLASS TRUSS3D ***
//   *********************


#ifndef truss3d_h
#define truss3d_h


#include "nlstructuralelement.h"
#include "gaussintegrationrule.h"
#include "directerrorindicatorrc.h"

class Truss3d : public NLStructuralElement
, public DirectErrorIndicatorRCInterface
{
/*
   This class implements a two-node truss bar element for three-dimensional
   analysis.
*/

protected:
  double        length ;
  //     int           referenceNode;
public:
  Truss3d (int,Domain*) ;                       // constructor
  ~Truss3d ()   {}                              // destructor
  
  // FloatArray*   ComputeBodyLoadVectorAt (TimeStep*) ;
  void          computeLumpedMassMatrix (FloatMatrix& answer, TimeStep* tStep) ;
  void          computeMassMatrix (FloatMatrix& answer, TimeStep* tStep) 
    {computeLumpedMassMatrix (answer, tStep);}
  // FloatArray*   ComputeResultingBodyForceAt (TimeStep*) ;
  // FloatMatrix*  computeStiffnessMatrix () ;
  // FloatArray*   ComputeStrainVector (GaussPoint*,TimeStep*) ;
  // FloatArray*   ComputeInitialStrainVector (TimeStep* );
  int              giveLocalCoordinateSystem (FloatMatrix& answer);
  /**
    Computes the global coordinates from given element's local coordinates.
    Required by nonlocal material models. 
    @returns nonzero if successful
    */
  virtual int computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) ;

  virtual int            computeNumberOfDofs (EquationID ut) {return 6;}
  virtual void           giveDofManDofIDMask  (int inode, EquationID, IntArray& ) const;
  
  // characteristic length in gp (for some material models)
  double        giveCharacteristicLenght (GaussPoint*, const FloatArray&) 
    {return this->giveLength();}
  
  double        computeVolumeAround (GaussPoint*) ;
  
  virtual int testElementExtension (ElementExtension ext) {return ((ext==Element_EdgeLoadSupport)?1:0);}
  //int    hasEdgeLoadSupport () {return 1;}

  /** Interface requesting service */
  Interface* giveInterface (InterfaceType);

#ifdef __OOFEG
      void          drawRawGeometry (oofegGraphicContext&);
      void          drawDeformedGeometry(oofegGraphicContext&, UnknownType);
#endif
// 
// definition & identification
//
  const char* giveClassName () const { return "Truss3d" ;}
  classType            giveClassID () const { return Truss3dClass; } 
  IRResultType initializeFrom (InputRecord* ir);

/**
 @name The element interface required by SPRNodalRecoveryModelInterface
*/
//@{
 /*
   Determines the characteristic size of element. This quantity is defined as follows:
   For 1D it is the element length, for 2D it is the square root of element area.
  */
 virtual double DirectErrorIndicatorRCI_giveCharacteristicSize () {return this->giveLength();}
//@}

protected:
  // edge load support
  void  computeEgdeNMatrixAt (FloatMatrix& answer, GaussPoint*) ;
  void  giveEdgeDofMapping (IntArray& answer, int) const;
  double        computeEdgeVolumeAround (GaussPoint*, int);
  void          computeEdgeIpGlobalCoords (FloatArray& answer, GaussPoint* gp, int iEdge) 
  {computeGlobalCoordinates(answer, *(gp->giveCoordinates()));}
  int   computeLoadLEToLRotationMatrix (FloatMatrix&, int, GaussPoint*) ;
  void          computeBmatrixAt (GaussPoint*, FloatMatrix&, int=1, int=ALL_STRAINS) ;
  void          computeNLBMatrixAt(FloatMatrix& answer, GaussPoint*, int);
  void          computeNmatrixAt (GaussPoint*, FloatMatrix&) ;
//  int           computeGtoNRotationMatrix (FloatMatrix&);
  void          computeGaussPoints () ;
 integrationDomain  giveIntegrationDomain () {return _Line;}

  double        giveLength () ;
  double        givePitch () ;
  int           giveApproxOrder () {return 1;}
      

} ;

#endif








