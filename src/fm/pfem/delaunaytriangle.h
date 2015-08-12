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

//   *************************
//   *** DELAUNAY TRIANGLE ***
//   *************************

// used for Delaunay meshing


#ifndef delaunaytriangle_h
#define delaunaytriangle_h

#include <list>
#include "intarray.h"
#include "floatarray.h"

namespace oofem {
class Domain;
class TimeStep;
class IntArray;
class FloatArray;
template< class T >class LocalInsertionData;

/**
 * Delaunay triangle for the triangulation of a set of nodes.
 * According the definition a Delaunay triangle has an empty circumscribed circle.
 *
 * @author David Krybus
 */
class DelaunayTriangle
{
public:
    /// Constructor
    DelaunayTriangle(Domain *d, int node1, int node2, int node3);
    /// Destructor
    ~DelaunayTriangle();

    /// Gives the x coordinate of the center of the circumscribed circle
    double giveXCenterCoordinate() const { return circumCircle.at(1); }
    /// Gives the y coordinate of the center of the circumscribed circle
    double giveYCenterCoordinate() const { return circumCircle.at(2); }
    /// Gives the radius of the circumscribed circle
    double giveCircumRadius() const { return circumCircle.at(3); }
    /// Calculates the distance of a passed point to the center of the circumscribed circle
    double giveDistanceToCenter(const FloatArray &coords);

    /// Gives the i-node of the triangle
    int giveNode(int i) { return nodes.at(i); }
    /// Returns a list of octree cells and with iterator position in their member lists
    std :: list< LocalInsertionData< DelaunayTriangle * > > *giveListOfCellsAndPosition();
    /// Sets the flag whether Delaunay condition is fulfilled
    void setValidFlag(bool newFlag) { validFlag = newFlag; }
    /// Gives true if the delaunay triangle is valid
    bool giveValidFlag() { return validFlag; }
    /// Gives the length of the shortest triangle edge
    double giveShortestEdgeLength();
    /// Gives the length of the edge between two nodes
    double giveEdgeLength(int nodeA, int nodeB);


private:
    /// Calculates the parameters of the circumscribed circle
    void computeCircumcircle();
    /// Sets up the parameters of the calculated circumscribed circle
    void setCircumCircle(double x, double y, double r);
    /// Domain where the nodes are defined
    Domain *domain;
    /// Nodes defining the triangle
    IntArray nodes;
    /// Parameters of the circumscribed circle: coordinates of center (x,y) and its radius
    FloatArray circumCircle;            //x, y, r
    /// Flag for Delaunay property
    bool validFlag;
    /// In order to allow fast search in octree, every triangle stores list of octree cells where its circumscribed circle is contained.
    std :: list< LocalInsertionData< DelaunayTriangle * > >listOfCellsContainedInAndPosition;
};
} // end namespace oofem
#endif // delaunaytriangle_h
