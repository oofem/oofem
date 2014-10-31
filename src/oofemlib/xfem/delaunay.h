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

#ifndef delaunay_h
#define delaunay_h

#include "oofemcfg.h"

#include <vector>

namespace oofem {
class FloatArray;
class Triangle;

/**
 * O(n4) algorithm, only for testing purposes.
 * Yes, but 4th order in n. For the xfem element subdivision, n does not increase when the mesh is refined.
 * Time will tell if it is too slow ... /ES
 * @author chamrova
 * @author Erik Svenning
 */
class OOFEM_EXPORT Delaunay
{
public:
    Delaunay() : mTol(1.0e-12) { }

    bool colinear(const FloatArray &iP1, const FloatArray &iP2, const FloatArray &iP3) const;
    void printTriangles(std :: vector< Triangle > &triangles);
    bool isInsideCC(const FloatArray &iP, const FloatArray &iP1, const FloatArray &iP2, const FloatArray &iP3) const;
    void triangulate(const std :: vector< FloatArray > &iVertices, std :: vector< Triangle > &oTriangles) const;

private:
    const double mTol; /// Tolerance used when checking if points are colinear.
};
} // end namespace oofem
#endif  // delaunay_h
