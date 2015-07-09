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

#include "edge2d.h"
#include "delaunaytriangle.h"

namespace oofem {
Edge2D :: Edge2D(int node1, int node2) :
    nodeNumbers{node1, node2}
{}

Edge2D :: ~Edge2D()
{ }

bool
Edge2D :: operator==(const Edge2D &right)
{
    return ( this->nodeNumbers.first == right.nodeNumbers.first && this->nodeNumbers.second == right.nodeNumbers.second ) ||
           ( this->nodeNumbers.first == right.nodeNumbers.second && this->nodeNumbers.second == right.nodeNumbers.first );
}



AlphaEdge2D :: AlphaEdge2D(int node1, int node2, double _length) :
    Edge2D(node1, node2)
    , isOnConvexHull(true)
    , outerAlphaBound(0.0)
    , innerAlphaBound(0.0)
    , length(_length)
{
    sharedByTriangles [ 0 ] = NULL;
    sharedByTriangles [ 1 ] = NULL;
}

AlphaEdge2D :: ~AlphaEdge2D()
{ }

void
AlphaEdge2D :: setSharing(int n, DelaunayTriangle *pTE)
{
    sharedByTriangles [ n - 1 ] = pTE;
}
}
