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

#ifndef edge2d_h
#define edge2d_h

#include "intarray.h"

namespace oofem {
class DelaunayTriangle;
class IntArray;

/**
 * Edge class for Delaunay triangulation
 */
class Edge2D
{
private:
    /// Global node numbers
    IntArray nodes;             //at 1 - nodeA; at 2 -nodeB

public:
    /// Constructor
    Edge2D(int node1, int node2);
    /// Destructor
    ~Edge2D();
    void setNode(int n, int node) { nodes.at(n) = node; }
    int giveNode(int n) { return nodes.at(n); }
    virtual bool operator == ( const Edge2D & right );
};

/**
 * Class for the boundary recognition method - alpha shape.
 * Case 1 - convex hull edge: alpha greater than outerAlphaBound, edge is part of alpha shape
 * Case 2 - shared edge: alpha greater than innerAlphaBound, edge lies within alpha shape;
 * Alpha between bounds, edge is exactly on alpha shape
 * Alpha smaller than outerAlphaBound, edge lies outside of alpha shape
 */
class AlphaEdge2D : public Edge2D
{
public:
    /// Constructor
    AlphaEdge2D(int node1, int node2, double _length);
    /// Destructor
    ~AlphaEdge2D();

    void setOuterAlphaBound(double alphaMin) { outerAlphaBound = alphaMin; }
    void setInnerAlphaBound(double alphaMax) { innerAlphaBound = alphaMax; }
    double giveOuterAlphaBound() { return outerAlphaBound; }
    double giveInnerAlphaBound() { return innerAlphaBound; }
    void setHullFlag(bool flag) { isOnConvexHull = flag; }
    bool giveHullFlag() { return isOnConvexHull; }
    void setSharing(int n, DelaunayTriangle *pTE);
    DelaunayTriangle *giveShared(int n) { return sharedByTriangles [ n - 1 ]; }
    double giveLength() { return length; }

private:
    /// Convex hull flag means edge is not shared by two triangle
    bool isOnConvexHull;
    double outerAlphaBound;
    double innerAlphaBound;
    /// Triangles which share the alphaEdge
    DelaunayTriangle *sharedByTriangles [ 2 ];
    /// Length of edge is stored in order to allow variable alpha value
    double length;
};
}
#endif
