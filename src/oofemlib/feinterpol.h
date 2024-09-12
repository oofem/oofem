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
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef feinterpol_h
#define feinterpol_h

#include "oofemenv.h"
#include "error.h"
#include "inputrecord.h"
#include "intarray.h"
#include "integrationdomain.h"
#include "elementgeometrytype.h"
#include "materialmode.h"
#include "node.h"
#include "element.h"

namespace oofem {
class Element;
class FloatArray;
class FloatMatrix;
class IntArray;
class IntegrationRule;

template <std::size_t N> class FloatArrayF;
template <std::size_t N, std::size_t M> class FloatMatrixF;

/**
 * Class representing a general abstraction for cell geometry.
 * The motivation for this class is that the interpolation classes require to pass underlying cell geometry.
 * The aim here is to hide and encapsulate as much as possible from actual cell geometry specification,
 * elements describe its geometry using nodes, which are independent objects, some cells may be
 * directly specified using vertices, etc.
 */
class OOFEM_EXPORT FEICellGeometry
{
public:
    FEICellGeometry() { }
    virtual ~FEICellGeometry() { }
    virtual int giveNumberOfVertices() const = 0;
    virtual const FloatArray &giveVertexCoordinates(int i) const = 0;
    virtual const Element_Geometry_Type giveGeometryType() const = 0;
};


/**
 * Void cell geometry wrapper.
 * Allows to use some interpolation services not needing the reference to cell geometry.
 */
class OOFEM_EXPORT FEIVoidCellGeometry : public FEICellGeometry
{
    FloatArray tmp;
public:
    FEIVoidCellGeometry() : FEICellGeometry() { }
    virtual ~FEIVoidCellGeometry() { }
    int giveNumberOfVertices() const override
    {
        OOFEM_ERROR("no reference geometry");
    }
    const FloatArray &giveVertexCoordinates(int i) const override
    {
        OOFEM_ERROR("no reference geometry");
    }
    const Element_Geometry_Type giveGeometryType() const override
    {
        return EGT_unknown;
    }
    std :: string errorInfo(const char *func) const { return func; } ///@todo Class name?
};

/**
 * Wrapper around element definition to provide FEICellGeometry interface.
 */
class OOFEM_EXPORT FEIElementGeometryWrapper : public FEICellGeometry
{
protected:
    const Element *elem;
public:
    FEIElementGeometryWrapper(const Element * elem) :
        FEICellGeometry(), elem(elem) { }
    virtual ~FEIElementGeometryWrapper() { }
    int giveNumberOfVertices() const override;
    const FloatArray &giveVertexCoordinates(int i) const override
    {
        return elem->giveNode(i)->giveCoordinates();
    }
    const Element_Geometry_Type giveGeometryType() const override {
        return elem->giveGeometryType();
    }
};


/**
 * Wrapper around cell with vertex coordinates stored in FloatArray**.
 */
class OOFEM_EXPORT FEIVertexListGeometryWrapper : public FEICellGeometry
{
protected:
    const std::vector< FloatArray > &coords;
    Element_Geometry_Type gtype;

public:
    FEIVertexListGeometryWrapper(const std::vector< FloatArray > &coords, const Element_Geometry_Type gt) : 
        FEICellGeometry(), coords(coords), gtype(gt) { }
    virtual ~FEIVertexListGeometryWrapper() { }
    int giveNumberOfVertices() const override { return (int)this->coords.size(); }
    const FloatArray &giveVertexCoordinates(int i) const override { return this->coords [ i - 1 ]; }
    const Element_Geometry_Type giveGeometryType() const override {return gtype;}
};

/**
 * Class representing a general abstraction for finite element interpolation class.
 * The boundary functions denote the (numbered) region that have 1 spatial dimension (i.e. edges) or 2 spatial dimensions.
 */
class OOFEM_EXPORT FEInterpolation
{
protected:
    int order = 0;

public:
    FEInterpolation(int o) : order(o) { }
    virtual ~FEInterpolation() = default;
    /// Initializes receiver according to object description stored in input record.
    virtual void initializeFrom(InputRecord &ir) { }

    /* @name basic interpolation services */
    //@{
    /**
     * Returns the integration domain of the interpolator.
     */
    virtual integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const = 0;
    /**
     * Returns the geometry type fo the interpolator.
     */
    virtual const Element_Geometry_Type giveGeometryType() const = 0;
    /**
     * Returns the interpolation order.
     */
    int giveInterpolationOrder() const { return order; }
    /**
     * Returns local element node numbers defining the approximation. Typically this corresponds to all element nodes.
     * But for elements with mixed interpolation, we need to select subset of nodes 
     * (quadratic triangle with linear interpolation). This method can query element geometry type (from given element)
     * and compile the nodal set.
     * @brief Returns list of element nodes (and list of internal dof managers) (including on edges and surfaces) defining the approximation
     * @note Required by mpm module
     */
    virtual void giveCellDofMans(IntArray& nodes, IntArray& internalDofMans, Element* elem) const {}
    /**
     * Evaluates the array of interpolation functions (shape functions) at given point.
     * @param answer Contains resulting array of evaluated interpolation functions.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Determinant of the Jacobian.
     */
    virtual double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const  = 0;
    /**
     * Evaluates the matrix of second derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
    {
        OOFEM_ERROR("not implemented");
    }
    /**
     * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are wrt local (parent) coordinate system
     * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxij.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
    {
        OOFEM_ERROR("not implemented");
    }
    /**
     * Returns a matrix containing the local coordinates for each node corresponding to the interpolation
     */
    virtual void giveLocalNodeCoords(FloatMatrix &answer, const Element_Geometry_Type) const
    {
        OOFEM_ERROR("FEInterpolation::giveLocalNodeCoords: not implemented");
    }
    /**
     * Evaluates global coordinates from given local ones.
     * @param answer Contains resulting global coordinates.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Evaluates local coordinates from given global ones.
     * If local coordinates cannot be found (generate elements, or point far outside geometry,
     * then the center coordinate will be used as a last resort, and the return value will be zero.
     * @param answer Contains evaluated local coordinates.
     * @param gcoords Array containing global coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Nonzero is returned if point is within the element geometry, zero otherwise.
     */
    virtual int global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Evaluates the determinant of the transformation.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Determinant of the transformation.
     */
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const;
    /**
     * Gives the jacobian matrix at the local coordinates.
     * @param jacobianMatrix The requested matrix.
     * @param lcoords Local coordinates.
     * @param cellgeo Element geometry.
     */
    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
    { OOFEM_ERROR("Not overloaded."); }

    /**
     * Sets up a suitable integration rule for numerical integrating over volume.
     * The required polynomial order for the determinant of the jacobian is added automatically.
     * @param order Polynomial order of integrand (should NOT including determinant of jacobian).
     */
    virtual std::unique_ptr<IntegrationRule> giveIntegrationRule(int order, const Element_Geometry_Type) const;
    //@}

    /** @name Edge boundary functions.
     * Provide interpolation services for boundary edges (entity of dimension 1)
     */
    //@{
    /**
     * Evaluates the basis functions on the requested boundary.
     * Only basis functions that are nonzero anywhere on the boundary are given. Ordering can be obtained from giveBoundaryNodes.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param answer Basis functions Array to be filled with the boundary nodes.
     * @param boundary Boundary number.
     * @param lcoords The local coordinates (on the boundary local coordinate system).
     * @param cellgeo Underlying cell geometry.
     * @todo
     */
    virtual void boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Evaluates the normal out of the edge at given point.
     * @param answer Contains resulting normal vector.
     * @param isurf Determines the surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Surface mapping jacobian.
     */
    virtual double boundaryEdgeEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Evaluates the determinant of the transformation Jacobian on the requested boundary.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param boundary Boundary number.
     * @param lcoords The local coordinates (on the boundary local coordinate system).
     * @param cellgeo Underlying cell geometry.
     * @return The determinant of the boundary transformation Jacobian.
     */
    virtual double boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Maps the local boundary coordinates to global.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param answer Global coordinates.
     * @param boundary Boundary number.
     * @param lcoords The local coordinates (on the boundary local coordinate system).
     * @param cellgeo Underlying cell geometry.
     */
    virtual void boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /// Returns boundary integration domain
    virtual integrationDomain giveBoundaryEdgeIntegrationDomain(int boundary, const Element_Geometry_Type) const = 0;
    /**
     * Sets up a suitable integration rule for integrating over the requested boundary.
     * The required polynomial order for the determinant of the jacobian is added automatically.
     * @param order Polynomial order of the integrand (should NOT including determinant of jacobian).
     * @param boundary Boundary number.
     */
    virtual std::unique_ptr<IntegrationRule> giveBoundaryEdgeIntegrationRule(int order, int boundary, const Element_Geometry_Type) const ;
    /**
     * Gives the boundary nodes for requested boundary number.
     * @param answer Array to be filled with the boundary nodes.
     * @param boundary Boundary number.
     */
    virtual IntArray boundaryEdgeGiveNodes(int boundary, const Element_Geometry_Type) const = 0;
    //@}

    /**@name Surface interpolation services 
     * Provide interpolation services for boundary edges (entities of dimension 2)
     */
    //@{
    /**
     * Evaluates the array of edge interpolation functions (shape functions) at given point.
     * @param answer Contains resulting array of evaluated interpolation functions.
     * @param isurf Surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void boundarySurfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined).
     * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNj/dxi.
     * @param isurf Determines the surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void boundarySurfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Evaluates the normal out of the surface at given point.
     * @param answer Contains resulting normal vector.
     * @param isurf Determines the surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Surface mapping jacobian.
     */
    virtual double boundarySurfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;

    /**
     * Evaluates edge global coordinates from given local ones.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined).
     * @param answer Contains resulting global coordinates.
     * @param isurf Determines the surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void boundarySurfaceLocal2global(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const =0;
    /**
     * Evaluates the edge jacobian of transformation between local and global coordinates.
     * @param isurf Determines the surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Determinant of the transformation.
     */
    virtual double boundarySurfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /// Returns boundary integration domain
    virtual integrationDomain giveBoundarySurfaceIntegrationDomain(int boundary, const Element_Geometry_Type) const = 0;
    /**
     * Sets up a suitable integration rule for integrating over the requested boundary.
     * The required polynomial order for the determinant of the jacobian is added automatically.
     * @param order Polynomial order of the integrand (should NOT including determinant of jacobian).
     * @param boundary Boundary number.
     */
    virtual std::unique_ptr<IntegrationRule> giveBoundarySurfaceIntegrationRule(int order, int boundary, const Element_Geometry_Type) const ;
    /**
     * Gives the boundary nodes for requested boundary number.
     * @param answer Array to be filled with the boundary nodes.
     * @param boundary Boundary number.
     */
    virtual IntArray boundarySurfaceGiveNodes(int boundary, const Element_Geometry_Type) const = 0;
    //@}

    /** @name General boundary interpolation functions.
     * Provide interpolation servises for boundary entities with one dimension lower than the receiver interpolation.
     * Typically these are mapped to boundaryEdge and boundarySurface methods depending on dimension.
     *
     */
    //@{
    /**
     * Gives the boundary nodes for requested boundary number.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param answer Array to be filled with the boundary nodes.
     * @param boundary Boundary number.
     */
    virtual IntArray boundaryGiveNodes(int boundary, const Element_Geometry_Type) const = 0;
    /**
     * Evaluates the basis functions on the requested boundary.
     * Only basis functions that are nonzero anywhere on the boundary are given. Ordering can be obtained from giveBoundaryNodes.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param answer Basis functions Array to be filled with the boundary nodes.
     * @param boundary Boundary number.
     * @param lcoords The local coordinates (on the boundary local coordinate system).
     * @param cellgeo Underlying cell geometry.
     */
    virtual void boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Evaluates the normal on the requested boundary.
     * @param answer The evaluated normal.
     * @param boundary Boundary number.
     * @param lcoords The local coordinates (on the boundary local coordinate system).
     * @param cellgeo Underlying cell geometry.
     * @return The boundary transformation Jacobian.
     */
    virtual double boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Evaluates the determinant of the transformation Jacobian on the requested boundary.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param boundary Boundary number.
     * @param lcoords The local coordinates (on the boundary local coordinate system).
     * @param cellgeo Underlying cell geometry.
     * @return The determinant of the boundary transformation Jacobian.
     */
    virtual double boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Maps the local boundary coordinates to global.
     * Boundaries are defined as the corner nodes for 1D geometries, edges for 2D geometries and surfaces for 3D geometries.
     * @param answer Global coordinates.
     * @param boundary Boundary number.
     * @param lcoords The local coordinates (on the boundary local coordinate system).
     * @param cellgeo Underlying cell geometry.
     */
    virtual void boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const = 0;
    /**
     * Computes the integral @f$ \int_S n \cdot x \mathrm{d}s @f$.
     * @param boundary Boundary number.
     * @param cellgeo Underlying cell geometry.
     * @return Evaluated integral.
     */
    virtual double evalNXIntegral(int boundary, const FEICellGeometry &cellgeo) const 
    {
        OOFEM_ERROR("Not implemented");
    }
    /// Returns boundary integration domain
    virtual integrationDomain giveBoundaryIntegrationDomain(int boundary, const Element_Geometry_Type) const = 0;
    /**
     * Sets up a suitable integration rule for integrating over the requested boundary.
     * The required polynomial order for the determinant of the jacobian is added automatically.
     * @param order Polynomial order of the integrand (should NOT including determinant of jacobian).
     * @param boundary Boundary number.
     */
    virtual std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int order, int boundary, const Element_Geometry_Type) const;
    //@}

    /**@name Methods to support interpolation defined on patch by patch basis. */
    //@{
    /**
     * Returns indices (zero based) of nonzero basis functions for given knot span.
     * The knot span identifies the sub-region of the finite element.
     * @return Nonzero if mask is provided, zero otherwise meaning that all
     * basis functions are generally nonzero.
     */
    virtual int giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask) const { return 0; }
    /**
     * Returns the number of nonzero basis functions at individual knot span,
     * @return Zero in case of all basis functions generaedgeEvaldNdslly nonzero, answer otherwise.
     */
    virtual int giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan) const { return 0; }
    /**
     * Returns true, if receiver is formulated on sub-patch basis.
     */
    virtual bool hasSubPatchFormulation() const  { return false; }
    /**
     * Returns the subdivision of patch parametric space
     */
    virtual const FloatArray *giveKnotVector() const { return nullptr; }
    /**
     * Returns the number of knot spans of the receiver.
     */
    virtual int giveNumberOfKnotSpans(int dim) const { return 0; }
    /**
     * Returns the knot values of the receiver.
     */
    virtual const FloatArray *giveKnotValues(int dim) const { return nullptr; }
    /**
     * Returns the knot multiplicity of the receiver.
     */
    virtual const IntArray *giveKnotMultiplicity(int dim) const { return nullptr; }
    /**
     * Returns number of spatial dimensions.
     */
    virtual int giveNsd(const Element_Geometry_Type) const = 0;
    /**
     * Returns number of edges.
     */
    virtual int giveNumberOfEdges(const Element_Geometry_Type) const 
    { OOFEM_ERROR("FEInterpolation :: giveNumberOfEdges : Not overloaded."); }
    //@}

    /**
     * Returns the number of geometric nodes of the receiver.
     */
    virtual int giveNumberOfNodes(const Element_Geometry_Type) const
    { OOFEM_ERROR("giveNumberOfNodes: Not overloaded."); }
    //@}

    std :: string errorInfo(const char *func) const { return func; } ///@todo Class name?
};
} // end namespace oofem
#endif // feinterpol_h
