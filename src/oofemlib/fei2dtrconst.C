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

#include "fei2dtrconst.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {
void
FEI2dTrConst :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(1);
    answer.at(1) = 1;
}

void
FEI2dTrConst :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(1, 2);
    answer.zero();
}

void
FEI2dTrConst :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l1, l2, l3;
    answer.resize(2);

    l1 = lcoords.at(1);
    l2 = lcoords.at(2);
    l3 = 1.0 - l1 - l2;

    answer.at(1) = ( l1 * cellgeo.giveVertexCoordinates(1)->at(xind) +
                    l2 * cellgeo.giveVertexCoordinates(2)->at(xind) +
                    l3 * cellgeo.giveVertexCoordinates(3)->at(xind) );
    answer.at(2) = ( l1 * cellgeo.giveVertexCoordinates(1)->at(yind) +
                    l2 * cellgeo.giveVertexCoordinates(2)->at(yind) +
                    l3 * cellgeo.giveVertexCoordinates(3)->at(yind) );
}

#define POINT_TOL 1.e-3

int
FEI2dTrConst :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    double area, x1, x2, x3, y1, y2, y3;
    answer.resize(3);

    x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    x3 = cellgeo.giveVertexCoordinates(3)->at(xind);

    y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    y3 = cellgeo.giveVertexCoordinates(3)->at(yind);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );


    answer.at(1) = ( ( x2 * y3 - x3 * y2 ) + ( y2 - y3 ) * coords.at(xind) + ( x3 - x2 ) * coords.at(yind) ) / 2. / area;
    answer.at(2) = ( ( x3 * y1 - x1 * y3 ) + ( y3 - y1 ) * coords.at(xind) + ( x1 - x3 ) * coords.at(yind) ) / 2. / area;
    answer.at(3) = ( ( x1 * y2 - x2 * y1 ) + ( y1 - y2 ) * coords.at(xind) + ( x2 - x1 ) * coords.at(yind) ) / 2. / area;

    // check if point is inside
    int i;
    for ( i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            return 0;
        }

        if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            return 0;
        }
    }

    return 1;
}


double
FEI2dTrConst :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double area, x1, x2, x3, y1, y2, y3;

    x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    x3 = cellgeo.giveVertexCoordinates(3)->at(xind);

    y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    y3 = cellgeo.giveVertexCoordinates(3)->at(yind);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );
    return 2.0 * area;
}


void
FEI2dTrConst :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(1);
    answer.at(1) = 1.0;
}

void
FEI2dTrConst :: edgeEvaldNds(FloatArray &answer, int iedge,
                           const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(2);
    answer.zero();
}

void
FEI2dTrConst :: edgeLocal2global(FloatArray &answer, int iedge,
                               const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray n(2);
    n.at(1) = (1 - lcoords(0))*0.5;
    n.at(2) = (1 + lcoords(0))*0.5;
    this->computeLocalEdgeMapping(edgeNodes, iedge);

    answer.resize(2);
    answer.at(1) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(xind) +
                     n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(xind) );
    answer.at(2) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(yind) +
                     n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(yind) );
}

double
FEI2dTrConst :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    return 0.5 * this->edgeComputeLength(edgeNodes, cellgeo);
}


void
FEI2dTrConst :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    int aNode = 0, bNode = 0;
    edgeNodes.resize(2);

    if ( iedge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        aNode = 3;
        bNode = 1;
    } else {
        OOFEM_ERROR2("FEI2dTrConst :: computeEdgeMapping: wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;

    //OOFEM_ERROR("FEI2dTrConst :: computeLocalEdgeMapping: not implemented");

}

double
FEI2dTrConst :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    double dx, dy;
    int nodeA, nodeB;

    nodeA   = edgeNodes.at(1);
    nodeB   = edgeNodes.at(2);

    dx      = cellgeo.giveVertexCoordinates(nodeB)->at(xind) - cellgeo.giveVertexCoordinates(nodeA)->at(xind);
    dy      = cellgeo.giveVertexCoordinates(nodeB)->at(yind) - cellgeo.giveVertexCoordinates(nodeA)->at(yind);
    return ( sqrt(dx * dx + dy * dy) );
}
} // end namespace oofem
