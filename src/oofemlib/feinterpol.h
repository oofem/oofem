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

namespace oofem {
class Element;
/**
 * Class representing a general abstraction for cell geometry.
 * The motivation for this class is that the interpolation classes require to pass underlying cell geometry.
 * The aim here is to hide and encapsulate as much as possible from actual cell geometry specification,
 * elements describe its geometry using nodes, which are independent objects, some cells may be
 * directly specified using vertices, etc.
 *
 */
class FEICellGeometry
{
public:
    FEICellGeometry() {}
    virtual int giveNumberOfVertices() const = 0;
    virtual const FloatArray *giveVertexCoordinates(int i) const = 0;
};


/**
 * void cell geometry wrapper. Alows to use some interpolation servises not needing the reference to cell geometry.
 */
class FEIVoidCellGeometry : public FEICellGeometry
{
public:
    FEIVoidCellGeometry() : FEICellGeometry() {}
    int giveNumberOfVertices() const { OOFEM_ERROR("FEIVoidCellGeometry: no reference geometry");
                                       return 0; }
    const FloatArray *giveVertexCoordinates(int i) const {
        OOFEM_ERROR("FEIVoidCellGeometry: no reference geometry");
        return NULL;
    }
};

/**
 * wrapper around element definition to provide FEICellGeometry interface
 */
class FEIElementGeometryWrapper : public FEICellGeometry
{
protected:
    Element *elem;
public:
    FEIElementGeometryWrapper(Element *elem) : FEICellGeometry() { this->elem = elem; }
    int giveNumberOfVertices() const;
    const FloatArray *giveVertexCoordinates(int i) const;
};


/**
 * Wrapper around cell with vertex coordinates stored in FloatArray**
 */
class FEIVertexListGeometryWrapper : public FEICellGeometry
{
protected:
    const FloatArray **coords;
    int nvertices;
public:
    FEIVertexListGeometryWrapper(int nvertices, const FloatArray **coords) : FEICellGeometry()
    { this->nvertices = nvertices;
      this->coords = coords; }
    int giveNumberOfVertices() const { return this->nvertices; }
    const FloatArray *giveVertexCoordinates(int i) const { return this->coords [ i - 1 ]; }
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
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time) = 0;
    /**
     * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param matrix contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time) = 0;
    /**
     * Evaluates the matrix of second derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param matrix contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
     * @param nodes array of node numbers defining the interpolation geometry
     * @param lcoords array containing (local) coordinates
     * @param time time
     */
    virtual void evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time) {
        OOFEM_ERROR("FEInterpolation::evald2Ndx2: not implemented");
    }
    /**
     * Evaluates global coordinates from given local ones
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer contains resulting global coordinates
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time) = 0;
    /**
     * Evaluates local coordinates from given global ones. Returns nonzero if local coordinates are interpolating,
     * zero if extrapolating (nonzero is returned if point is within the element geometry, zero otherwise).
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer contains evaluated local coordinates
     * @param gcoords array containing global coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     * @return nonzero is returned if point is within the element geometry, zero otherwise
     */
    virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time) = 0;
    /**
     * Evaluates local coordinates from given global ones. Returns nonzero if local coordinates are interpolating,
     * zero if extrapolating (nonzero is returned if point is within the element geometry, zero otherwise).
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer contains evaluated local coordinates
     * @param coords coordinates of nodes defining the interpolation geometry
     * @param gcoords array containing global coordinates
     * @param time time
     * @return nonzero is returned if point is within the element geometry, zero otherwise
     */
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time) = 0;

    /**
     * Returns indices (zero based) of nonzero basis functions for given knot span
     * The knot span identifies the sub-region of the finite element
     * @returns nonzero if mask is provided, zero otherwise meaning that all
     * basis functions are generally nonzero
     */
    ///Initializes receiver acording to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    /**@name Methods to support interpolation defined on patch by patch basis*/
    //@{
    virtual int giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask) { return 0; }
    /** Returns the number of nonzero basis functions at individual knot span,
     *  @param returns zero in case of all basis functions generally nonzero, answer otherwise */
    virtual int giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan) { return 0; }
    /**
     * Returns true, if receiver is formulated on sub-patch basis
     */
    virtual bool hasSubPatchFormulation() { return false; }
    /**
     *
     * Returns the subdivision of patch parametric space
     */
    virtual double **const giveKnotVector() { return NULL; }
    /**
     * Returns the number of knot spans of the receiver
     */
    virtual int giveNumberOfKnotSpans(int dim) { return 0; }
    /**
     * Returns the knot values of the receiver
     */
    virtual FloatArray *const giveKnotValues(int dim) { return NULL; }
    /**
     * Returns the knot multiplicity of the receiver
     */
    virtual IntArray *const giveKnotMultiplicity(int dim) { return NULL; }
    /**
     * Returns number of spatial dimensions
     */
    virtual int const giveNsd() = 0;

    //@}
};
} // end namespace oofem
#endif // feinterpol_h






