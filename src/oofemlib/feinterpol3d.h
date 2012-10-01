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

#ifndef feinterpol3d_h
#define feinterpol3d_h

#include "feinterpol.h"

namespace oofem {
/**
 * Class representing a general abstraction for surface finite element interpolation class.
 */
class FEInterpolation3d : public FEInterpolation
{
public:
    FEInterpolation3d(int o) : FEInterpolation(o) { };
    virtual int giveNsd() { return 3; }

    /**
     * Computes the exact volume.
     * @param cellgeo Cell geometry for the element.
     * @return Volume of geometry.
     */
    virtual double giveVolume(const FEICellGeometry &cellgeo) const
    { OOFEM_ERROR("FEInterpolation3d :: giveVolume - Not implemented in subclass."); return 0; }
    
    virtual void boundaryGiveNodes(IntArray &answer, int boundary)
    { this->computeLocalSurfaceMapping(answer, boundary); }
    virtual void boundaryEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { this->surfaceEvalN(answer, lcoords, cellgeo); }
    virtual double boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { return this->surfaceEvalNormal(answer, boundary, lcoords, cellgeo); }
    virtual double boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { return this->surfaceGiveTransformationJacobian(boundary, lcoords, cellgeo); }

    /**@name Edge interpolation services */
    //@{
    /**
     * Evaluates the array of edge interpolation functions (shape functions) at given point.
     * @param answer Contains resulting array of evaluated interpolation functions.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNj/dxi.
     * @param iedge Determines the edge number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void edgeEvaldNdx(FloatMatrix &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates edge global coordinates from given local ones.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined).
     * @param answer Contains resulting global coordinates.
     * @param iedge Determines edge number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void edgeLocal2global(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates the edge jacobian of transformation between local and global coordinates.
     * @param iedge Determines edge number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Determinant of the transformation.
     */
    virtual double edgeGiveTransformationJacobian(int iedge,
                                                  const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo) = 0;

    virtual void computeLocalEdgeMapping(IntArray &edgeNodes, int iedge) = 0;
    void computeEdgeMapping(IntArray &edgeNodes, IntArray &elemNodes, int iedge) {
        int i, size;
        IntArray ln;
        this->computeLocalEdgeMapping(ln, iedge);
        size = ln.giveSize();
        edgeNodes.resize(size);
        for ( i = 1; i <= size; i++ ) { edgeNodes.at(i) = elemNodes.at( ln.at(i) ); }
    }
    //@}

    /**@name Surface interpolation services */
    //@{
    /**
     * Evaluates the array of edge interpolation functions (shape functions) at given point.
     * @param answer Contains resulting array of evaluated interpolation functions.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void surfaceEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined).
     * @param answer Contains resulting matrix of derivatives, the member at i,j position contains value of dNj/dxi.
     * @param isurf Determines the surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void surfaceEvaldNdx (FloatMatrix&answer, int isurf,
            const FloatArray& lcoords, const FEICellGeometry& cellgeo)
    {
        OOFEM_ERROR("FEInterpolation3D :: surfaceEvaldNdx - Not implemented");
    }
    /**
     * Evaluates the normal out of the surface at given point.
     * @param answer Contains resulting normal vector.
     * @param isurf Determines the surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Surface mapping jacobian.
     */
    virtual double surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords,
            const FEICellGeometry &cellgeo)
    {
        OOFEM_ERROR("FEInterpolation3D :: surfaceEvalNormal - Not implemented");
        return -1.0;
    }

    /**
     * Evaluates edge global coordinates from given local ones.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined).
     * @param answer Contains resulting global coordinates.
     * @param isurf Determines the surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void surfaceLocal2global(FloatArray &answer, int isurf,
                                     const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates the edge jacobian of transformation between local and global coordinates.
     * @param isurf Determines the surface number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Determinant of the transformation.
     */
    virtual double surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                     const FEICellGeometry &cellgeo) = 0;

    virtual void computeLocalSurfaceMapping(IntArray &surfNodes, int isurf) = 0;
    void computeSurfaceMapping(IntArray &surfNodes, IntArray &elemNodes, int isurf) {
        int i, size;
        IntArray ln;
        this->computeLocalSurfaceMapping(ln, isurf);
        size = ln.giveSize();
        surfNodes.resize(size);
        for ( i = 1; i <= size; i++ ) { surfNodes.at(i) = elemNodes.at( ln.at(i) ); }
    }
    //@}
};
} // end namespace oofem
#endif // feinterpol3d_h






