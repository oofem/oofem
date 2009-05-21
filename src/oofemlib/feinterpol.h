/* $Header: /home/cvs/bp/oofem/oofemlib/src/feinterpol.h,v 1.1 2003/04/06 14:08:24 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *****************************
//   *** CLASS FEInterpolation ***
//   *****************************


#ifndef feinterpol_h
#define feinterpol_h

#include "flotarry.h"
#include "intarray.h"
#include "domain.h"

class Element;
/**
 * Class representing a general abstraction for cell geometry.
 * The motivation for this class is that the interpolation classes require to pass underlying cell geometry.
 * The aim here is to hide and encapsulate as much as possible from actual cell geometry specification,
 * elements describe its geometry using nodes, which are independent objects, some cells may be
 * directly specified using vertices, etc.
 * 
 */
class FEICellGeometry {
 public:
  virtual int giveNumberOfVertices () const = 0;
  virtual const FloatArray* giveVertexCoordinates(int i) const = 0;
};

/**
 * wrapper aroun element definition to provide FEICellGeometry interface
 */
class FEIElementGeometry : public FEICellGeometry {
 protected:
  Element* elem;
 public:
  FEIElementGeometry (Element* elem) {this->elem = elem;}
  int giveNumberOfVertices () const;
  const FloatArray* giveVertexCoordinates(int i) const;
};


/**
 * Class representing a general abstraction for finite element interpolation class.
 */
class FEInterpolation
{
protected:
    int order;

public:
    FEInterpolation(int o) { order = o; }
    virtual ~FEInterpolation() { }
    /**
     * Returns the interpolation order
     */
    int giveInterpolationOrder() { return order; }
    /**
     * Evaluates the array of interpolation functions (shape functions) at given point.
     * @param answer contains resulting array of evaluated interpolation functions
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time) = 0;
    /**
     * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param matrix contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time) = 0;
    /**
     * Evaluates global coordinates from given local ones
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer contains resulting global coordinates
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time) = 0;
    /**
     * Evaluates local coordinates from given global ones. Returns nonzero if local coordinates are interpolating,
     * zero if extrapolating (nonzero is returned if point is within the element geometry, zero otherwise).
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer contains evaluated local coordinates
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     * @return nonzero is returned if point is within the element geometry, zero otherwise
     */
    virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time) = 0;
    /**
     * Evaluates the jacobian of transformation between local and global coordinates.
     */
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time) = 0;
    /**
       Returns indices (zero based) of nonzero basis functions for given knot span 
       The knot span identifies the sub-region of the finite element
       @returns nonzero if mask is provided, zero otherwise meaning that all 
       basis functions are generally nonzero
    */
    ///Initializes receiver acording to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir) {return IRRT_OK;}

    /**@name Methods to support interpolation defined on patch by patch basis*/
    //@{
    virtual int giveKnotBasisFuncMask (const IntArray& knotSpan, IntArray& mask) {return 0;}
    /** Returns the number of nonzero basis functions at individual knot span,
        @param returns zero in case of all basis functions generally nonzero, answer otherwise */
    virtual int  giveNumberOfKnotBasisFunctions () {return 0;}
    /**
       Returns true, if receiver is formulated on sub-patch basis
     */
    virtual bool hasSubPatchFormulation() {return false;}
    /**
     * 
     * Returns the subdivision of patch parametric space
     */
    virtual double** const giveKnotVector() {return NULL;}
    /**
     * Returns the number of knot spans of the receiver
     */
    virtual int giveNumberOfKnotSpans(int dim) {return 0;}
    /**
     * Returns the knot multiplicity of the receiver
     */
    virtual int** const giveKnotMultiplicity() {return NULL;}

    //@}

};


  


#endif // feinterpol_h






