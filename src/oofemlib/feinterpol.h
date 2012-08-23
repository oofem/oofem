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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef feinterpol_h
#define feinterpol_h

#include "error.h"
#include "inputrecord.h"
#include "intarray.h"

namespace oofem {
class Element;
class FloatArray;
class FloatMatrix;
class IntArray;

/**
 * Class representing a general abstraction for cell geometry.
 * The motivation for this class is that the interpolation classes require to pass underlying cell geometry.
 * The aim here is to hide and encapsulate as much as possible from actual cell geometry specification,
 * elements describe its geometry using nodes, which are independent objects, some cells may be
 * directly specified using vertices, etc.
 */
class FEICellGeometry
{
public:
    FEICellGeometry() {}
    virtual ~FEICellGeometry() {}
    virtual int giveNumberOfVertices() const = 0;
    virtual const FloatArray *giveVertexCoordinates(int i) const = 0;
};


/**
 * Void cell geometry wrapper.
 * Allows to use some interpolation services not needing the reference to cell geometry.
 */
class FEIVoidCellGeometry : public FEICellGeometry
{
public:
    FEIVoidCellGeometry() : FEICellGeometry() {}
    virtual ~FEIVoidCellGeometry() {}
    int giveNumberOfVertices() const { OOFEM_ERROR("FEIVoidCellGeometry: no reference geometry");
                                       return 0; }
    const FloatArray *giveVertexCoordinates(int i) const {
        OOFEM_ERROR("FEIVoidCellGeometry: no reference geometry");
        return NULL;
    }
};

/**
 * Wrapper around element definition to provide FEICellGeometry interface.
 */
class FEIElementGeometryWrapper : public FEICellGeometry
{
protected:
    const Element *elem;
public:
    FEIElementGeometryWrapper(const Element *elem) : FEICellGeometry() { this->elem = elem; }
    virtual ~FEIElementGeometryWrapper() {}
    int giveNumberOfVertices() const;
    const FloatArray *giveVertexCoordinates(int i) const;
};


/**
 * Wrapper around cell with vertex coordinates stored in FloatArray**.
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
    virtual ~FEIVertexListGeometryWrapper() {}
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
     * Returns the interpolation order.
     */
    int giveInterpolationOrder() { return order; }
    /**
     * Evaluates the array of interpolation functions (shape functions) at given point.
     * @param answer Contains resulting array of evaluated interpolation functions.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates the matrix of second derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) {
        OOFEM_ERROR("FEInterpolation::evald2Ndx2: not implemented");
    }
    /**
     * Evaluates global coordinates from given local ones.
     * @param answer Contains resulting global coordinates.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates local coordinates from given global ones.
     * @param answer Contains evaluated local coordinates.
     * @param gcoords Array containing global coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Nonzero is returned if point is within the element geometry, zero otherwise.
     */
    virtual int global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates the determinant of the transformation.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Determinant of the transformation.
     */
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;

    /// Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    
    /**@name General boundary functions */
    //@{
    /**
     * Gives the boundary nodes for requested boundary number.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param answer Array to be filled with the boundary nodes.
     * @param boundary Boundary number.
     */
    virtual void boundaryGiveNodes(IntArray &answer, int boundary) = 0;
    
    /**
     * Evaluates the basis functions on the requested boundary.
     * Only basis functions that are nonzero anywhere on the boundary are given. Ordering can be obtained from giveBoundaryNodes.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param answer Basis functions Array to be filled with the boundary nodes.
     * @param boundary Boundary number.
     */
    virtual void boundaryEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;

    /**
     * Evaluates the determinant of the transformation Jacobian on the requested boundary.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param boundary Boundary number.
     * @return The determinant of the boundary transformation Jacobian.
     */
    virtual double boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    //@}

    /**@name Methods to support interpolation defined on patch by patch basis. */
    //@{
    /**
     * Returns indices (zero based) of nonzero basis functions for given knot span.
     * The knot span identifies the sub-region of the finite element.
     * @return Nonzero if mask is provided, zero otherwise meaning that all
     * basis functions are generally nonzero.
     */
    virtual int giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask) { return 0; }
    /**
     * Returns the number of nonzero basis functions at individual knot span,
     * @return Zero in case of all basis functions generally nonzero, answer otherwise.
     */
    virtual int giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan) { return 0; }
    /**
     * Returns true, if receiver is formulated on sub-patch basis.
     */
    virtual bool hasSubPatchFormulation() { return false; }
    /**
     * Returns the subdivision of patch parametric space
     */
    virtual const double * const * giveKnotVector() { return NULL; }
    /**
     * Returns the number of knot spans of the receiver.
     */
    virtual int giveNumberOfKnotSpans(int dim) { return 0; }
    /**
     * Returns the knot values of the receiver.
     */
    virtual const FloatArray * giveKnotValues(int dim) { return NULL; }
    /**
     * Returns the knot multiplicity of the receiver.
     */
    virtual const IntArray * giveKnotMultiplicity(int dim) { return NULL; }
    /**
     * Returns number of spatial dimensions.
     */
    virtual int giveNsd() = 0;
    //@}

    // Needs the jacobian matrix to determine the condition number.
    friend class MeshQualityErrorEstimator;

protected:
    /**
     * Gives the jacobian matrix at the local coordinates.
     * @param jacobianMatrix The requested matrix.
     * @param lcoords Local coordinates.
     * @param cellgeo Element geometry.
     */
    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { OOFEM_ERROR("FEInterpolation::giveJacobianMatrixAt : Not overloaded."); }
};
} // end namespace oofem
#endif // feinterpol_h

