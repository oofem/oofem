/* $Header: /home/cvs/bp/oofem/oofemlib/src/element.h,v 1.27 2003/04/06 14:08:24 bp Exp $ */
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

#ifndef feinurbs_h
#define feinurbs_h

/*
 * oofem nodes - control points (coordinates ) + dofs
 * oofem elements - NURBS patches as well as integration elements
 *
 *
 * NURBS PATCH:
 * knot vector - store knot coordinates + multiplicity
 * patch integration rule - keep list of elements
 *
 */

#include "feibspline.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "mathfem.h"

namespace oofem {

class NURBSInterpolation : public BSplineInterpolation
{
public:
    NURBSInterpolation(int nsd) : BSplineInterpolation(nsd) { }
    ~NURBSInterpolation();

    /**
     * Evaluates the array of interpolation functions (shape functions) at given point.
     * @param answer contains resulting array of evaluated interpolation functions
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     *
     * see also giveNonzeroBasisFunctMask method of BSplineInterpolation
     */
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time);
    /**
     * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time);
    /**
     * Evaluates global coordinates from given local ones
     * @param answer contains resulting global coordinates
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time);
    virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time) {
        OOFEM_ERROR("Not yet inplemented, contact lazy dr for implementation");
        return 0;
    }
    /**
     * Evaluates the jacobian of transformation between local and global coordinates.
     */
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time);

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "NURBSInterpolation"; }
}; // end of NURBSInterpolation class definition

} // end namespace oofem
#endif //feinurbs_h
