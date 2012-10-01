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

#include "fei2dtrlin.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {
void
FEI2dTrLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(3);

    answer.at(1) = lcoords.at(1);
    answer.at(2) = lcoords.at(2);
    answer.at(3) = 1. - lcoords.at(1) - lcoords.at(2);
}

void
FEI2dTrLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x1, x2, x3, y1, y2, y3, detJ;
    answer.resize(3, 2);

    x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    x3 = cellgeo.giveVertexCoordinates(3)->at(xind);

    y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    y3 = cellgeo.giveVertexCoordinates(3)->at(yind);

    detJ = x1*(y2 - y3) + x2*(-y1 + y3) + x3*(y1 - y2);

    answer.at(1, 1) = ( y2 - y3 ) / detJ;
    answer.at(1, 2) = ( x3 - x2 ) / detJ;

    answer.at(2, 1) = ( y3 - y1 ) / detJ;
    answer.at(2, 2) = ( x1 - x3 ) / detJ;

    answer.at(3, 1) = ( y1 - y2 ) / detJ;
    answer.at(3, 2) = ( x2 - x1 ) / detJ;
}

void
FEI2dTrLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
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
FEI2dTrLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    double detJ, x1, x2, x3, y1, y2, y3;
    answer.resize(3);

    x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    x3 = cellgeo.giveVertexCoordinates(3)->at(xind);

    y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    y3 = cellgeo.giveVertexCoordinates(3)->at(yind);

    detJ = x1*(y2 - y3) + x2*(-y1 + y3) + x3*(y1 - y2);

    answer.at(1) = ( ( x2 * y3 - x3 * y2 ) + ( y2 - y3 ) * coords.at(xind) + ( x3 - x2 ) * coords.at(yind) ) / detJ;
    answer.at(2) = ( ( x3 * y1 - x1 * y3 ) + ( y3 - y1 ) * coords.at(xind) + ( x1 - x3 ) * coords.at(yind) ) / detJ;
    //answer.at(3) = ( ( x1 * y2 - x2 * y1 ) + ( y1 - y2 ) * coords.at(xind) + ( x2 - x1 ) * coords.at(yind) ) / detJ;

    // check if point is inside
    bool inside = true;
    for ( int i = 1; i <= 2; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            answer.at(i) = 0.;
            inside = false;
        } else if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            answer.at(i) = 1.;
            inside = false;
        }
    }
    
    answer.at(3) = 1. - answer.at(1) - answer.at(2);

    return inside;
}


double
FEI2dTrLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x1, x2, x3, y1, y2, y3;

    x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    x3 = cellgeo.giveVertexCoordinates(3)->at(xind);

    y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    y3 = cellgeo.giveVertexCoordinates(3)->at(yind);

    return ( x1*(y2 - y3) + x2*(-y1 + y3) + x3*(y1 - y2) );
}


void
FEI2dTrLin :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

void
FEI2dTrLin :: edgeEvaldNds(FloatArray &answer, int iedge,
                           const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l;
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    l = this->edgeComputeLength(edgeNodes, cellgeo);

    answer.resize(2);
    answer.at(1) = -1.0 / l;
    answer.at(2) =  1.0 / l;
}

double FEI2dTrLin :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    normal.resize(2);
    normal.at(1) = cellgeo.giveVertexCoordinates(edgeNodes.at(2))->at(xind) - cellgeo.giveVertexCoordinates(edgeNodes.at(1))->at(xind);
    normal.at(2) = -(cellgeo.giveVertexCoordinates(edgeNodes.at(2))->at(yind) - cellgeo.giveVertexCoordinates(edgeNodes.at(1))->at(yind));
    return normal.normalize()*0.5;
}

void
FEI2dTrLin :: edgeLocal2global(FloatArray &answer, int iedge,
                               const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray n;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(n, lcoords, cellgeo);

    answer.resize(2);
    answer.at(1) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(xind) +
                     n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(xind) );
    answer.at(2) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(yind) +
                     n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(yind) );
}


void
FEI2dTrLin :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
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
        OOFEM_ERROR2("FEI2dTrLin :: computeEdgeMapping: wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;
}

double
FEI2dTrLin :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    double dx, dy;
    int nodeA, nodeB;

    nodeA = edgeNodes.at(1);
    nodeB = edgeNodes.at(2);

    dx = cellgeo.giveVertexCoordinates(nodeB)->at(xind) - cellgeo.giveVertexCoordinates(nodeA)->at(xind);
    dy = cellgeo.giveVertexCoordinates(nodeB)->at(yind) - cellgeo.giveVertexCoordinates(nodeA)->at(yind);
    return sqrt(dx * dx + dy * dy);
}

double FEI2dTrLin :: giveArea(const FEICellGeometry &cellgeo) const
{
    const FloatArray *p;
    double x1, x2, x3, y1, y2, y3;

    p = cellgeo.giveVertexCoordinates(1);
    x1 = p->at(1);
    y1 = p->at(2);
    p = cellgeo.giveVertexCoordinates(2);
    x2 = p->at(1);
    y2 = p->at(2);
    p = cellgeo.giveVertexCoordinates(3);
    x3 = p->at(1);
    y3 = p->at(2);

    return 0.5 * ( x1*(y2-y3) + x2*(-y1+y3) + x3*(y1-y2) ); ///@todo Absolute value or not?
}
} // end namespace oofem
