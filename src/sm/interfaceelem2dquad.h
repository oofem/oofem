/* $Header: /home/cvs/bp/oofem/sm/src/Attic/interfaceelem2dquad.h,v 1.1.2.1 2004/04/05 15:19:47 bp Exp $ */
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

#ifndef interfaceelem2dquad_h
#define interfaceelem2dquad_h


#include "structuralelement.h"
#include "gaussintegrationrule.h"

class InterfaceElem2dQuad : public StructuralElement
{
/*
   This class implements a two dimensional interface element.
   Even if geometry approx is quadratic, the element is assumed straight
   If not straight, the rotation matrix depends on actual integration point
   and stiffness and strain computations should be modified.
*/

protected:
public:
  InterfaceElem2dQuad (int,Domain*) ;                       // constructor
  ~InterfaceElem2dQuad ()   {}                              // destructor
  
  /**
     Computes the global coordinates from given element's local coordinates.
     Required by nonlocal material models. 
     @returns nonzero if successful

  */
  virtual int computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) ;
  /**
     Computes the element local coordinates from given global coordinates.
     @returns nonzero if successful (if point is inside element); zero otherwise
  */
  virtual int computeLocalCoordinates (FloatArray& answer, const FloatArray& gcoords) ;
  
  virtual int            computeNumberOfDofs (EquationID ut) {return 12;}
  virtual void           giveDofManDofIDMask  (int inode, EquationID, IntArray& ) const;
  
  double        computeVolumeAround (GaussPoint*) ;


  virtual int testElementExtension (ElementExtension ext) {return 0;}
  
  /** Interface requesting service */
  Interface* giveInterface (InterfaceType) {return NULL;}
  
#ifdef __OOFEG
  void          drawRawGeometry (oofegGraphicContext&);
  void          drawDeformedGeometry(oofegGraphicContext&, UnknownType);
  void          drawScalar   (oofegGraphicContext& context);
#endif
  // 
  // definition & identification
  //
  const char* giveClassName () const { return "InterfaceElem2dQuad" ;}
  classType            giveClassID () const { return Truss1dClass; } 
  IRResultType initializeFrom (InputRecord* ir);
  Element_Geometry_Type giveGeometryType() const {return EGT_line_2;}
  
protected:
  void          computeBmatrixAt (GaussPoint*, FloatMatrix&, int=1, int=ALL_STRAINS) ;
  void          computeNmatrixAt (GaussPoint*, FloatMatrix&) {}
  void          computeGaussPoints () ;
  integrationDomain  giveIntegrationDomain () {return _Line;}
  
  int           giveApproxOrder () {return 1;}
  /*
    internal
    Note: this is not overloaded computeGtoLRotationMatrix
    the gp parameter is needed, since element geometry can be (in general) curved
  
  //void          computeGtoLRotationMatrix(FloatMatrix& answer, GaussPoint* gp);
  */
  int          computeGtoLRotationMatrix(FloatMatrix& answer);
} ;

#endif // interfaceelem2dquad_h








