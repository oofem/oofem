/* $Header: /home/cvs/bp/oofem/oofemlib/src/fei2dtrlin.C,v 1.1.4.1 2004/04/05 15:19:43 bp Exp $ */
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

#include "fei2dtrlin.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"

void
FEI2dTrLin :: evalN(FloatArray &answer, const FloatArray &lcoords, double time)
{
    answer.resize(3);

    answer.at(1) = lcoords.at(1);
    answer.at(2) = lcoords.at(2);
    answer.at(3) = 1. - lcoords.at(1) - lcoords.at(2);

    return;
}

void
FEI2dTrLin :: evaldNdx(FloatMatrix &answer, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    double x1, x2, x3, y1, y2, y3, area;
    answer.resize(3, 2);

    x1 = coords [ 0 ]->at(xind);
    x2 = coords [ 1 ]->at(xind);
    x3 = coords [ 2 ]->at(xind);

    y1 = coords [ 0 ]->at(yind);
    y2 = coords [ 1 ]->at(yind);
    y3 = coords [ 2 ]->at(yind);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    answer.at(1, 1) = ( y2 - y3 ) / ( 2. * area );
    answer.at(1, 2) = ( x3 - x2 ) / ( 2. * area );

    answer.at(2, 1) = ( y3 - y1 ) / ( 2. * area );
    answer.at(2, 2) = ( x1 - x3 ) / ( 2. * area );

    answer.at(3, 1) = ( y1 - y2 ) / ( 2. * area );
    answer.at(3, 2) = ( x2 - x1 ) / ( 2. * area );
}

void
FEI2dTrLin :: local2global(FloatArray &answer, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    double l1, l2, l3;
    answer.resize(2);

    l1 = lcoords.at(1);
    l2 = lcoords.at(2);
    l3 = 1.0 - l1 - l2;

    answer.at(1) = ( l1 * coords [ 0 ]->at(xind) +
                    l2 * coords [ 1 ]->at(xind) +
                    l3 * coords [ 2 ]->at(xind) );
    answer.at(2) = ( l1 * coords [ 0 ]->at(yind) +
                    l2 * coords [ 1 ]->at(yind) +
                    l3 * coords [ 2 ]->at(yind) );
}

#define POINT_TOL 1.e-3

int
FEI2dTrLin :: global2local(FloatArray &answer, const FloatArray **nc, const FloatArray &coords, double time)
{
    double area, x1, x2, x3, y1, y2, y3;
    answer.resize(3);

    x1 = nc [ 0 ]->at(xind);
    x2 = nc [ 1 ]->at(xind);
    x3 = nc [ 2 ]->at(xind);

    y1 = nc [ 0 ]->at(yind);
    y2 = nc [ 1 ]->at(yind);
    y3 = nc [ 2 ]->at(yind);

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
FEI2dTrLin :: giveTransformationJacobian(const FloatArray **coords, const FloatArray &lcoords, double time)
{
    double area, x1, x2, x3, y1, y2, y3;

    x1 = coords [ 0 ]->at(xind);
    x2 = coords [ 1 ]->at(xind);
    x3 = coords [ 2 ]->at(xind);

    y1 = coords [ 0 ]->at(yind);
    y2 = coords [ 1 ]->at(yind);
    y3 = coords [ 2 ]->at(yind);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );
    return 2.0 * area;
}


void
FEI2dTrLin :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, double time)
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

void
FEI2dTrLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                           const FloatArray **coords, const FloatArray &lcoords, double time)
{
    double l;
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    l = this->edgeComputeLength(edgeNodes, coords);

    answer.resize(2, 1);
    answer.at(1, 1) = -1.0 / l;
    answer.at(2, 1) =  1.0 / l;
}

void
FEI2dTrLin :: edgeLocal2global(FloatArray &answer, int iedge,
                               const FloatArray **coords, const FloatArray &lcoords, double time)
{
    IntArray edgeNodes;
    FloatArray n;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(n, lcoords, time);

    answer.resize(2);
    answer.at(1) = ( n.at(1) * coords [ edgeNodes.at(1) - 1 ]->at(xind) +
                    n.at(2) * coords [ edgeNodes.at(2) - 1 ]->at(xind) );
    answer.at(2) = ( n.at(1) * coords [ edgeNodes.at(1) - 1 ]->at(yind) +
                    n.at(2) * coords [ edgeNodes.at(2) - 1 ]->at(yind) );
}


double
FEI2dTrLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    return 0.5 * this->edgeComputeLength(edgeNodes, coords);
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
FEI2dTrLin :: edgeComputeLength(IntArray &edgeNodes, const FloatArray **coords)
{
    double dx, dy;
    int nodeA, nodeB;

    nodeA   = edgeNodes.at(1) - 1;
    nodeB   = edgeNodes.at(2) - 1;

    dx      = coords [ nodeB ]->at(xind) - coords [ nodeA ]->at(xind);
    dy      = coords [ nodeB ]->at(yind) - coords [ nodeA ]->at(yind);
    return ( sqrt(dx * dx + dy * dy) );
}
