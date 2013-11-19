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
template< class T > class LocalInsertionData;

/**
 * Delaunay triangle for triangulation of set of nodes.
 * According the definition a Delaunay triangle has an empty circumscribed circle.
 */
class DelaunayTriangle
{
public:
    /// Constructor
    DelaunayTriangle(Domain * d, int node1, int node2, int node3);
    /// Destructor
    ~DelaunayTriangle();

    double giveXCenterCoordinate() const { return circumCircle.at(1); }
    double giveYCenterCoordinate() const { return circumCircle.at(2); }
    double giveCircumRadius() const { return circumCircle.at(3); }
    double giveDistanceToCenter(const FloatArray &coords);

    int giveNode(int i) { return nodes.at(i); }
    std :: list< LocalInsertionData< DelaunayTriangle * > > *giveListOfCellsAndPosition();
    void setValidFlag(bool newFlag) { validFlag = newFlag; }
    bool giveValidFlag() { return validFlag; }
    double giveShortestEdgeLength();
    double giveEdgeLength(int nodeA, int nodeB);


private:
    void computeCircumcircle();
    void setCircumCircle(double x, double y, double r);

    Domain *domain;
    IntArray nodes;
    /// Parameters of the circumscribed circle: coordinates of center (x,y) and its radius
    FloatArray circumCircle;            //x, y, r
    /// Flag for Delaunay property
    bool validFlag;
    /// In order to allow fast search in octree, every triangle stores list of octree cells where its circumscribed circle is contained.
    std :: list< LocalInsertionData< DelaunayTriangle * > > *listOfCellsContainedInAndPosition;
};
} // end namespace oofem
#endif // delaunaytriangle_h
