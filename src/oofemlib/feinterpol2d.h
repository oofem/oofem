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

#ifndef feinterpol2d_h
#define feinterpol2d_h

#include "flotarry.h"

#include "feinterpol.h"

namespace oofem {
/**
 * Class representing a general abstraction for surface finite element interpolation class.
 */
class FEInterpolation2d : public FEInterpolation
{
public:
    FEInterpolation2d(int o) : FEInterpolation(o) { }

    virtual int giveNsd() { return 2; }

    /**
     * Computes the exact area.
     * @param cellgeo Cell geometry for the element.
     * @return Area of geometry.
     */
    virtual double giveArea(const FEICellGeometry &cellgeo) const
    { OOFEM_ERROR("FEInterpolation2d :: giveArea - Not implemented in subclass."); return 0; }


    virtual void boundaryGiveNodes(IntArray &answer, int boundary)  { this->computeLocalEdgeMapping(answer, boundary); }
    virtual void boundaryEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) { this->edgeEvalN(answer, lcoords, cellgeo); }
    virtual double boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) { return this->edgeEvalNormal(answer, boundary, lcoords, cellgeo); }
    virtual double boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) { return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo); }
    
    /**@name Edge interpolation services. */
    //@{
    virtual void computeLocalEdgeMapping(IntArray &edgeNodes, int iedge) = 0;
    void computeEdgeMapping(IntArray &edgeNodes, IntArray &elemNodes, int iedge) {
        int i, size;
        IntArray ln;
        this->computeLocalEdgeMapping(ln, iedge);
        size = ln.giveSize();
        edgeNodes.resize(size);
        for ( i = 1; i <= size; i++ ) { edgeNodes.at(i) = elemNodes.at( ln.at(i) ); }
    }
    /**
     * Evaluates the array of edge interpolation functions (shape functions) at given point.
     * @param answer Contains resulting array of evaluated interpolation functions.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates the normal on the given edge.
     * @param answer Contains the evaluated normal.
     * @param iedge Determines the edge number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual double edgeEvalNormal(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) = 0;
    /**
     * Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined).
     * @param answer Contains resulting array of derivatives, the member at i position contains value of @f$ \frac{\mathrm{d}N_i}{\mathrm{d}s} @f$.
     * @param iedge Determines the edge number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     */
    virtual void edgeEvaldNds(FloatArray &answer, int iedge,
                              const FloatArray &lcoords,
                              const FEICellGeometry &cellgeo) = 0;
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
     * Evaluates the edge Jacobian of transformation between local and global coordinates.
     * @param iedge Determines edge number.
     * @param lcoords Array containing (local) coordinates.
     * @param cellgeo Underlying cell geometry.
     * @return Determinant of the mapping on the given edge.
     */
    virtual double edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo) {
        FloatArray normal;
        return this->edgeEvalNormal(normal, iedge, lcoords, cellgeo);
    }
    //@}
};
} // end namespace oofem
#endif // feinterpol2d_h






