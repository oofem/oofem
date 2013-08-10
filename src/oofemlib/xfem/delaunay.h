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

#ifndef delaunay_h
#define delaunay_h

#include "alist.h"
#include <vector>

namespace oofem {
class FloatArray;
class Triangle;

/**
 * O(n4) algorithm, only for testing purposes.
 * @author chamrova
 *
 *  // Yes, but 4th order in n. For the xfem element subdivision, n does not increase when the mesh is refined. Time
 *	// will tell if it is too slow ... /ES
 * @author Erik Svenning
 */
class Delaunay
{
public:
    bool colinear(FloatArray *p1, FloatArray *p2, FloatArray *p3);
    void printTriangles(AList< Triangle > *triangles);
    bool isInsideCC(FloatArray *p, FloatArray *p1, FloatArray *p2, FloatArray *p3);
    void triangulate(const std :: vector< FloatArray > &iVertices, AList< Triangle > *triangles);
};
} // end namespace oofem
#endif  // delaunay_h
