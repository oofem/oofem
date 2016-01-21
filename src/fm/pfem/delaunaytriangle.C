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

#include "delaunaytriangle.h"
#include "dofmanager.h"
#include "math.h"
#include "octreelocalizert.h"
#include "mathfem.h"

namespace oofem {
DelaunayTriangle :: DelaunayTriangle(Domain *d, int node1, int node2, int node3) :
    domain(d),
    nodes(3),

    circumCircle(3),
    validFlag(true)
{
    nodes.at(1) = node1;
    nodes.at(2) = node2;
    nodes.at(3) = node3;

    computeCircumcircle();
}

DelaunayTriangle :: ~DelaunayTriangle()
{
}

void
DelaunayTriangle :: setCircumCircle(double x, double y, double r)
{
    this->circumCircle.at(1) = x;
    this->circumCircle.at(2) = y;
    this->circumCircle.at(3) = r;
}

double
DelaunayTriangle :: giveDistanceToCenter(const FloatArray &coords)
{
    double x = coords.at(1);
    double y = coords.at(2);

    double xC = giveXCenterCoordinate();
    double yC = giveYCenterCoordinate();

    return ( sqrt( ( xC - x ) * ( xC - x ) + ( yC - y ) * ( yC - y ) ) );
}

void
DelaunayTriangle :: computeCircumcircle()
{
    double x1, x2, x3;
    double y1, y2, y3;

    double a, bx, by, c;

    DofManager *dmanA, *dmanB, *dmanC;

    dmanA = domain->giveDofManager( giveNode(1) );
    x1 = dmanA->giveCoordinate(1);
    y1 = dmanA->giveCoordinate(2);

    dmanB = domain->giveDofManager( giveNode(2) );
    x2 = dmanB->giveCoordinate(1);
    y2 = dmanB->giveCoordinate(2);

    dmanC = domain->giveDofManager( giveNode(3) );
    x3 = dmanC->giveCoordinate(1);
    y3 = dmanC->giveCoordinate(2);

    a = x1 * y2 + y1 * x3 + x2 * y3 - 1.0 * ( x1 * y3 + y1 * x2 + y2 * x3 );
    bx = -1.0 * ( ( ( x1 * x1 + y1 * y1 ) * y2 + y1 * ( x3 * x3 + y3 * y3 ) + ( x2 * x2 + y2 * y2 ) * y3 )
                 - 1.0 * ( ( x1 * x1 + y1 * y1 ) * y3 + y1 * ( x2 * x2 + y2 * y2 ) + y2 * ( x3 * x3 + y3 * y3 ) ) );
    by = ( ( ( x1 * x1 + y1 * y1 ) * x2 + x1 * ( x3 * x3 + y3 * y3 ) + ( x2 * x2 + y2 * y2 ) * x3 )
          - 1.0 * ( ( x1 * x1 + y1 * y1 ) * x3 + x1 * ( x2 * x2 + y2 * y2 ) + x2 * ( x3 * x3 + y3 * y3 ) ) );
    c = ( ( ( x1 * x1 + y1 * y1 ) * x2 * y3 + x1 * y2 * ( x3 * x3 + y3 * y3 ) + y1 * ( x2 * x2 + y2 * y2 ) * x3 )
         - 1.0 * ( ( x1 * x1 + y1 * y1 ) * y2 * x3 + x1 * ( x2 * x2 + y2 * y2 ) * y3 + y1 * x2 * ( x3 * x3 + y3 * y3 ) ) );

    double xCenterCoordinate = ( -1.0 * bx / ( 2 * a ) );
    double yCenterCoordinate = ( -1.0 * by / ( 2 * a ) );
    double absA = a < 0 ? -1.0 * a : a;

    double radius = ( ( sqrt(bx * bx + by * by + 4.0 * a * c) ) / ( 2.0 * absA ) );

    setCircumCircle(xCenterCoordinate, yCenterCoordinate, radius);
}

std :: list< LocalInsertionData< DelaunayTriangle * > > *
DelaunayTriangle :: giveListOfCellsAndPosition()
{
    return & listOfCellsContainedInAndPosition;
}

double
DelaunayTriangle :: giveShortestEdgeLength()
{
    double length1 = giveEdgeLength(1, 2);
    double length2 = giveEdgeLength(2, 3);
    double length3 = giveEdgeLength(3, 1);
    return min( length1, min(length2, length3) );
}

double
DelaunayTriangle :: giveEdgeLength(int nodeA, int nodeB)
{
    DofManager *dmanA = domain->giveDofManager( giveNode(nodeA) );
    DofManager *dmanB = domain->giveDofManager( giveNode(nodeB) );

    return dmanA->giveCoordinates()->distance( dmanB->giveCoordinates() );
}
} // end namespace oofem

