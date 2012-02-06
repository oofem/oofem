/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2011   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef lattice2d_h
#define lattice2d_h

#include "latticestructuralelement.h"

namespace oofem {
/**
 * This class implements a 2-dimensional lattice element
 */

class Lattice2d : public  LatticeStructuralElement
{

protected:
  double kappa,pitch, length;

  double width, thickness;
  FloatArray gpCoords;
  
public:
  Lattice2d (int,Domain*) ;                       // constructor
  ~Lattice2d () ;                                 // destructor

  int           giveLocalCoordinateSystem (FloatMatrix& answer);

 /**
  Computes the global coordinates from given element's local coordinates.
  @returns returns nonzero if successful
  */
 virtual int computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) ;

 int giveCrackFlag();
 
 void giveCrossSectionCoordinates(FloatArray & coords);

 double giveCrackWidth();

 // characteristic length in gp (for some material models)
  double        giveCharacteristicLenght (GaussPoint*, const FloatArray&) 
  {return this->giveLength();}
  
  double        giveLength () ;

  void giveGPCoords(FloatArray& coords){coords = this->gpCoords;} 
  double giveArea(){return this->width*this->thickness;}
  double giveDissipation();
  double giveDeltaDissipation();

  Element_Geometry_Type giveGeometryType() const { return EGT_line_2; }


  virtual int testElementExtension (ElementExtension ext) {return ((ext==Element_EdgeLoadSupport)?1:0);}


 virtual int            computeNumberOfDofs (EquationID ut) {return 6;}
 virtual void           giveDofManDofIDMask  (int inode, EquationID, IntArray& ) const;
  double        computeVolumeAround (GaussPoint*) ;

// 
// definition & identification
//
  const char* giveClassName () const { return "Lattice2d" ;}
  classType            giveClassID () const { return Lattice2dClass; } 
  IRResultType initializeFrom (InputRecord* ir);

#ifdef __OOFEG

  void          drawYourself (oofegGraphicContext& context);
  void          drawRawGeometry (oofegGraphicContext&);
  //  void          drawRawVoronoi (oofegGraphicContext&);
  void          drawDeformedGeometry(oofegGraphicContext&, UnknownType);
  //  void          drawDeformedVoronoi(oofegGraphicContext&, UnknownType);

  //  void drawScalar   (oofegGraphicContext& context);

  void drawSpecial(oofegGraphicContext& gc);
#endif

protected:

  void          computeBmatrixAt (GaussPoint*, FloatMatrix&, int=1, int=ALL_STRAINS) ;
  void          computeNmatrixAt (GaussPoint*, FloatMatrix &){;}  ;
  bool           computeGtoLRotationMatrix (FloatMatrix&);

  void  computeStiffnessMatrix (FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep);
    
  int giveNumberOfCrossSectionNodes(){return 2;}

  double        givePitch () ;
  void          computeGaussPoints () ;
  integrationDomain  giveIntegrationDomain () {return _Line;}
  
} ;
} // end namespace oofem
#endif
