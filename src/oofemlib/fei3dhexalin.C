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

#include "fei3dhexalin.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {
void
FEI3dHexaLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x, y, z;
    answer.resize(8);

    x = lcoords.at(1);
    y = lcoords.at(2);
    z = lcoords.at(3);

    answer.at(1)  = 0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. + z );
    answer.at(2)  = 0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. + z );
    answer.at(3)  = 0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. + z );
    answer.at(4)  = 0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. + z );
    answer.at(5)  = 0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. - z );
    answer.at(6)  = 0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. - z );
    answer.at(7)  = 0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. - z );
    answer.at(8)  = 0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. - z );
}

void
FEI3dHexaLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    int i;
    FloatMatrix jacobianMatrix(3, 3), inv(3, 3);
    FloatArray dx(8), dy(8), dz(8);
    double u, v, w;

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    inv.beInverseOf(jacobianMatrix);

    this->giveDerivativeKsi(dx, v, w);
    this->giveDerivativeEta(dy, u, w);
    this->giveDerivativeDzeta(dz, u, v);

    answer.resize(8, 3);

    for ( i = 1; i <= 8; i++ ) {
        answer.at(i, 1) = dx.at(i) * inv.at(1, 1) + dy.at(i) * inv.at(1, 2) + dz.at(i) * inv.at(1, 3);
        answer.at(i, 2) = dx.at(i) * inv.at(2, 1) + dy.at(i) * inv.at(2, 2) + dz.at(i) * inv.at(2, 3);
        answer.at(i, 3) = dx.at(i) * inv.at(3, 1) + dy.at(i) * inv.at(3, 2) + dz.at(i) * inv.at(3, 3);
    }
}


void
FEI3dHexaLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    int i;
    double x, y, z;
    FloatArray n(8);

    x = lcoords.at(1);
    y = lcoords.at(2);
    z = lcoords.at(3);

    n.at(1)  = 0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. + z );
    n.at(2)  = 0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. + z );
    n.at(3)  = 0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. + z );
    n.at(4)  = 0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. + z );
    n.at(5)  = 0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. - z );
    n.at(6)  = 0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. - z );
    n.at(7)  = 0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. - z );
    n.at(8)  = 0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. - z );

    answer.resize(3);
    answer.zero();
    for ( i = 1; i <= 8; i++ ) {
        answer.at(1) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(1);
        answer.at(2) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(2);
        answer.at(3) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(3);
    }
}

#define POINT_TOL 1.e-3

int
FEI3dHexaLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    double x1, x2, x3, x4, x5, x6, x7, x8, a1, a2, a3, a4, a5, a6, a7, a8;
    double y1, y2, y3, y4, y5, y6, y7, y8, b1, b2, b3, b4, b5, b6, b7, b8;
    double z1, z2, z3, z4, z5, z6, z7, z8, c1, c2, c3, c4, c5, c6, c7, c8;
    double xp, yp, zp, u, v, w;
    FloatMatrix p(3, 3);
    FloatArray r(3), delta;
    int nite = 0;

    x1 = cellgeo.giveVertexCoordinates(1)->at(1);
    x2 = cellgeo.giveVertexCoordinates(2)->at(1);
    x3 = cellgeo.giveVertexCoordinates(3)->at(1);
    x4 = cellgeo.giveVertexCoordinates(4)->at(1);
    x5 = cellgeo.giveVertexCoordinates(5)->at(1);
    x6 = cellgeo.giveVertexCoordinates(6)->at(1);
    x7 = cellgeo.giveVertexCoordinates(7)->at(1);
    x8 = cellgeo.giveVertexCoordinates(8)->at(1);

    y1 = cellgeo.giveVertexCoordinates(1)->at(2);
    y2 = cellgeo.giveVertexCoordinates(2)->at(2);
    y3 = cellgeo.giveVertexCoordinates(3)->at(2);
    y4 = cellgeo.giveVertexCoordinates(4)->at(2);
    y5 = cellgeo.giveVertexCoordinates(5)->at(2);
    y6 = cellgeo.giveVertexCoordinates(6)->at(2);
    y7 = cellgeo.giveVertexCoordinates(7)->at(2);
    y8 = cellgeo.giveVertexCoordinates(8)->at(2);

    z1 = cellgeo.giveVertexCoordinates(1)->at(3);
    z2 = cellgeo.giveVertexCoordinates(2)->at(3);
    z3 = cellgeo.giveVertexCoordinates(3)->at(3);
    z4 = cellgeo.giveVertexCoordinates(4)->at(3);
    z5 = cellgeo.giveVertexCoordinates(5)->at(3);
    z6 = cellgeo.giveVertexCoordinates(6)->at(3);
    z7 = cellgeo.giveVertexCoordinates(7)->at(3);
    z8 = cellgeo.giveVertexCoordinates(8)->at(3);

    xp = coords.at(1);
    yp = coords.at(2);
    zp = coords.at(3);

    a1 =  x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8;
    a2 = -x1 - x2 + x3 + x4 - x5 - x6 + x7 + x8;
    a3 = -x1 + x2 + x3 - x4 - x5 + x6 + x7 - x8;
    a4 =  x1 + x2 + x3 + x4 - x5 - x6 - x7 - x8;
    a5 =  x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8;
    a6 = -x1 - x2 + x3 + x4 + x5 + x6 - x7 - x8;
    a7 = -x1 + x2 + x3 - x4 + x5 - x6 - x7 + x8;
    a8 =  x1 - x2 + x3 - x4 - x5 + x6 - x7 + x8;

    b1 =  y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8;
    b2 = -y1 - y2 + y3 + y4 - y5 - y6 + y7 + y8;
    b3 = -y1 + y2 + y3 - y4 - y5 + y6 + y7 - y8;
    b4 =  y1 + y2 + y3 + y4 - y5 - y6 - y7 - y8;
    b5 =  y1 - y2 + y3 - y4 + y5 - y6 + y7 - y8;
    b6 = -y1 - y2 + y3 + y4 + y5 + y6 - y7 - y8;
    b7 = -y1 + y2 + y3 - y4 + y5 - y6 - y7 + y8;
    b8 =  y1 - y2 + y3 - y4 - y5 + y6 - y7 + y8;

    c1 =  z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8;
    c2 = -z1 - z2 + z3 + z4 - z5 - z6 + z7 + z8;
    c3 = -z1 + z2 + z3 - z4 - z5 + z6 + z7 - z8;
    c4 =  z1 + z2 + z3 + z4 - z5 - z6 - z7 - z8;
    c5 =  z1 - z2 + z3 - z4 + z5 - z6 + z7 - z8;
    c6 = -z1 - z2 + z3 + z4 + z5 + z6 - z7 - z8;
    c7 = -z1 + z2 + z3 - z4 + z5 - z6 - z7 + z8;
    c8 =  z1 - z2 + z3 - z4 - z5 + z6 - z7 + z8;

    // setup initial guess
    answer.resize(3);
    answer.zero();

    // apply Newton-Raphson to solve the problem
    do {
        if ( ( ++nite ) > 10 ) {
            // fprintf(stderr, "FEI3dHexaLin :: global2local: no convergence after 10 iterations");
            return 0;
        }

        u = answer.at(1);
        v = answer.at(2);
        w = answer.at(3);

        // compute the residual
        r.at(1) = a1 + u * a2 + v * a3 + w * a4 + u * v * a5 + u * w * a6 + v * w * a7 + u * v * w * a8 - 8.0 * xp;
        r.at(2) = b1 + u * b2 + v * b3 + w * b4 + u * v * b5 + u * w * b6 + v * w * b7 + u * v * w * b8 - 8.0 * yp;
        r.at(3) = c1 + u * c2 + v * c3 + w * c4 + u * v * c5 + u * w * c6 + v * w * c7 + u * v * w * c8 - 8.0 * zp;

        // check for convergence
        if ( r.computeSquaredNorm() < 1.e-20 ) {
            break;                                  // sqrt(1.e-20) = 1.e-10
        }

        p.at(1, 1) = a2 + v * a5 + w * a6 + v * w * a8;
        p.at(1, 2) = a3 + u * a5 + w * a7 + u * w * a8;
        p.at(1, 3) = a4 + u * a6 + v * a7 + u * v * a8;

        p.at(2, 1) = b2 + v * b5 + w * b6 + v * w * b8;
        p.at(2, 2) = b3 + u * b5 + w * b7 + u * w * b8;
        p.at(2, 3) = b4 + u * b6 + v * b7 + u * v * b8;

        p.at(3, 1) = c2 + v * c5 + w * c6 + v * w * c8;
        p.at(3, 2) = c3 + u * c5 + w * c7 + u * w * c8;
        p.at(3, 3) = c4 + u * c6 + v * c7 + u * v * c8;

        // solve for corrections
        p.solveForRhs(r, delta);

        // update guess
        answer.subtract(delta);
    } while ( 1 );

    // test if inside
    for ( int i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( -1. - POINT_TOL ) ) {
            return 0;
        }

        if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            return 0;
        }
    }

    return 1;
}


double
FEI3dHexaLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(3, 3);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant();
}


void
FEI3dHexaLin :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

void
FEI3dHexaLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                             const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l;
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    l = this->edgeComputeLength(edgeNodes, cellgeo);

    answer.resize(2, 1);
    answer.at(1, 1) = -1.0 / l;
    answer.at(2, 1) =  1.0 / l;
}

void
FEI3dHexaLin :: edgeLocal2global(FloatArray &answer, int iedge,
                                 const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray n;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(n, lcoords, cellgeo);

    answer.resize(3);
    answer.at(1) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(1) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(1) );
    answer.at(2) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(2) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(2) );
    answer.at(3) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(3) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(3) );
}


double
FEI3dHexaLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    return 0.5 * this->edgeComputeLength(edgeNodes, cellgeo);
}


void
FEI3dHexaLin :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    int aNode = 0, bNode = 0;
    edgeNodes.resize(2);

    if ( iedge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
    } else if ( iedge == 3 ) { // edge between nodes 3 4
        aNode = 3;
        bNode = 4;
    } else if ( iedge == 4 ) { // edge between nodes 4 1
        aNode = 4;
        bNode = 1;
    } else if ( iedge == 5 ) { // edge between nodes 1 5
        aNode = 1;
        bNode = 5;
    } else if ( iedge == 6 ) { // edge between nodes 2 6
        aNode = 2;
        bNode = 6;
    } else if ( iedge == 7 ) { // edge between nodes 3 7
        aNode = 3;
        bNode = 7;
    } else if ( iedge == 8 ) { // edge between nodes 4 8
        aNode = 4;
        bNode = 8;
    } else if ( iedge == 9 ) { // edge between nodes 5 6
        aNode = 5;
        bNode = 6;
    } else if ( iedge == 10 ) { // edge between nodes 6 7
        aNode = 6;
        bNode = 7;
    } else if ( iedge == 11 ) { // edge between nodes 7 8
        aNode = 7;
        bNode = 8;
    } else if ( iedge == 12 ) { // edge between nodes 8 5
        aNode = 8;
        bNode = 5;
    } else {
        OOFEM_ERROR2("FEI3dHexaLin :: computeEdgeMapping: wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;
}

double
FEI3dHexaLin :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    double dx, dy, dz;
    int nodeA, nodeB;

    nodeA   = edgeNodes.at(1);
    nodeB   = edgeNodes.at(2);

    dx      = cellgeo.giveVertexCoordinates(nodeB)->at(1) - cellgeo.giveVertexCoordinates(nodeA)->at(1);
    dy      = cellgeo.giveVertexCoordinates(nodeB)->at(2) - cellgeo.giveVertexCoordinates(nodeA)->at(2);
    dz      = cellgeo.giveVertexCoordinates(nodeB)->at(3) - cellgeo.giveVertexCoordinates(nodeA)->at(3);
    return ( sqrt(dx * dx + dy * dy + dz * dz) );
}

void
FEI3dHexaLin :: surfaceEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi, eta;
    answer.resize(4);

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    answer.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
    answer.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
    answer.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
    answer.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25;
}

void
FEI3dHexaLin :: surfaceLocal2global(FloatArray &answer, int iedge,
                                    const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray nodes(4);
    double ksi, eta, n1, n2, n3, n4;

    computeLocalSurfaceMapping(nodes, iedge);

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    n1 = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
    n2 = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
    n3 = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
    n4 = ( 1. + ksi ) * ( 1. - eta ) * 0.25;

    answer.resize(3);
    answer.at(1) = n1 * cellgeo.giveVertexCoordinates( nodes.at(1) )->at(1) + n2 *cellgeo.giveVertexCoordinates( nodes.at(2) )->at(1) +
                   n3 *cellgeo.giveVertexCoordinates( nodes.at(3) )->at(1) + n4 *cellgeo.giveVertexCoordinates( nodes.at(4) )->at(1);
    answer.at(2) = n1 * cellgeo.giveVertexCoordinates( nodes.at(1) )->at(2) + n2 *cellgeo.giveVertexCoordinates( nodes.at(2) )->at(2) +
                   n3 *cellgeo.giveVertexCoordinates( nodes.at(3) )->at(2) + n4 *cellgeo.giveVertexCoordinates( nodes.at(4) )->at(2);
    answer.at(3) = n1 * cellgeo.giveVertexCoordinates( nodes.at(1) )->at(3) + n2 *cellgeo.giveVertexCoordinates( nodes.at(2) )->at(3) +
                   n3 *cellgeo.giveVertexCoordinates( nodes.at(3) )->at(3) + n4 *cellgeo.giveVertexCoordinates( nodes.at(4) )->at(3);
}

double
FEI3dHexaLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo)
{
    // only plane surface is supported !!!

    // FEI2dQuadLin cannot be used without modifying coordinates of surface nodes
    // therefore local calculation is done without using FEI2dQuadLin :: giveJacobianMatrixAt

    int n1, n2, n3, n4, n;
    IntArray snodes(4);
    FloatMatrix jacobianMatrix(2, 2);
    FloatArray sn [ 4 ];
    double length;
    int i;

    this->computeLocalSurfaceMapping(snodes, isurf);

    // check whether all surface nodes are in plane
    n1 = snodes.at(1);
    n2 = snodes.at(2);
    n3 = snodes.at(3);
    n4 = snodes.at(4);

    // get the normal using nodes 1 2 3
    FloatArray a(3), b(3), c(3);
    for ( i = 1; i <= 3; i++ ) {
        b.at(i) = cellgeo.giveVertexCoordinates(n2)->at(i) - cellgeo.giveVertexCoordinates(n1)->at(i);
        a.at(i) = cellgeo.giveVertexCoordinates(n3)->at(i) - cellgeo.giveVertexCoordinates(n1)->at(i);
    }

    n = n4;
    c.beVectorProductOf(b, a);
    length = c.computeNorm();

    if ( length < 1.0e-10 ) {
        // try nodes 1 3 4
        for ( i = 1; i <= 3; i++ ) {
            b.at(i) = cellgeo.giveVertexCoordinates(n4)->at(i) - cellgeo.giveVertexCoordinates(n1)->at(i);
        }

        n = n2;
        c.beVectorProductOf(a, b);
        length = c.computeNorm();

        if ( length < 1.0e-10 ) {
            OOFEM_ERROR("FEI3dHexaLin :: surfaceGiveTransformationJacobian: degenerated surface");
        }
    }

    c.times(1.0 / length);
    for ( i = 1; i <= 3; i++ ) {
        b.at(i) = cellgeo.giveVertexCoordinates(n)->at(i) - cellgeo.giveVertexCoordinates(n1)->at(i);
    }

    // check distance of the 4th node n
    if ( fabs( b.dotProduct(c) ) > 1.0e-6 ) {
        OOFEM_ERROR("FEI3dHexaLin :: surfaceGiveTransformationJacobian: not planar surface");
    }

    // map nodes to the surface (x,y) plane
    // get x and y unit vectors a and b (a is already computed as n3-n1)
    length = a.computeNorm();
    a.times(1.0 / length);
    b.beVectorProductOf(c, a);

    for ( i = 0; i < 4; i++ ) {
        sn [ i ].resize(2);
    }

    sn [ 0 ].at(1) = 0.0;
    sn [ 0 ].at(2) = 0.0;

    sn [ 2 ].at(1) = length;
    sn [ 2 ].at(2) = 0.0;

    for ( i = 1; i <= 3; i++ ) {
        c.at(i) = cellgeo.giveVertexCoordinates(n2)->at(i) - cellgeo.giveVertexCoordinates(n1)->at(i);
    }

    sn [ 1 ].at(1) = c.dotProduct(a);
    sn [ 1 ].at(2) = c.dotProduct(b);

    for ( i = 1; i <= 3; i++ ) {
        c.at(i) = cellgeo.giveVertexCoordinates(n4)->at(i) - cellgeo.giveVertexCoordinates(n1)->at(i);
    }

    sn [ 3 ].at(1) = c.dotProduct(a);
    sn [ 3 ].at(2) = c.dotProduct(b);

    double ksi, eta, x, y;
    FloatArray nx(4), ny(4);

    jacobianMatrix.resize(2, 2);
    jacobianMatrix.zero();

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    /*
     * FEI2dQuadLin interpolation(1, 2);
     *
     * interpolation.giveDerivativeKsi(nx, eta);
     * interpolation.giveDerivativeEta(ny, ksi);
     */

    nx.at(1) =  0.25 * ( 1. + eta );
    nx.at(2) = -0.25 * ( 1. + eta );
    nx.at(3) = -0.25 * ( 1. - eta );
    nx.at(4) =  0.25 * ( 1. - eta );

    ny.at(1) =  0.25 * ( 1. + ksi );
    ny.at(2) =  0.25 * ( 1. - ksi );
    ny.at(3) = -0.25 * ( 1. - ksi );
    ny.at(4) = -0.25 * ( 1. + ksi );


    for ( i = 1; i <= 4; i++ ) {
        x = sn [ i - 1 ].at(1);
        y = sn [ i - 1 ].at(2);

        jacobianMatrix.at(1, 1) += nx.at(i) * x;
        jacobianMatrix.at(1, 2) += nx.at(i) * y;
        jacobianMatrix.at(2, 1) += ny.at(i) * x;
        jacobianMatrix.at(2, 2) += ny.at(i) * y;
    }

    return jacobianMatrix.giveDeterminant();
}

void
FEI3dHexaLin :: computeLocalSurfaceMapping(IntArray &surfNodes, int isurf)
{
    int aNode = 0, bNode = 0, cNode = 0, dNode = 0;
    surfNodes.resize(4);

    if ( isurf == 1 ) { // surface 1 - nodes 1 4 3 2
        aNode = 1;
        bNode = 4;
        cNode = 3;
        dNode = 2;
    } else if ( isurf == 2 ) { // surface 2 - nodes 5 6 7 8
        aNode = 5;
        bNode = 6;
        cNode = 7;
        dNode = 8;
    } else if ( isurf == 3 ) { // surface 3  - nodes 1 2 6 5
        aNode = 1;
        bNode = 2;
        cNode = 6;
        dNode = 5;
    } else if ( isurf == 4 ) { // surface 4 - nodes 2 3 7 6
        aNode = 2;
        bNode = 3;
        cNode = 7;
        dNode = 6;
    } else if ( isurf == 5 ) { // surface 5 - nodes 3 4 8 7
        aNode = 3;
        bNode = 4;
        cNode = 8;
        dNode = 7;
    } else if ( isurf == 6 ) { // surface 6 - nodes 4 1 5 8
        aNode = 4;
        bNode = 1;
        cNode = 5;
        dNode = 8;
    } else {
        OOFEM_ERROR2("FEI3dHexaLin :: computeSurfaceMapping: wrong surface number (%d)", isurf);
    }

    surfNodes.at(1) = ( aNode );
    surfNodes.at(2) = ( bNode );
    surfNodes.at(3) = ( cNode );
    surfNodes.at(4) = ( dNode );
}


void
FEI3dHexaLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
// Computes it if it does not exist yet.
{
    int i;
    double x, y, z, u, v, w;
    FloatArray dx(8), dy(8), dz(8);

    jacobianMatrix.resize(3, 3);
    jacobianMatrix.zero();

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    this->giveDerivativeKsi(dx, v, w);
    this->giveDerivativeEta(dy, u, w);
    this->giveDerivativeDzeta(dz, u, v);

    for ( i = 1; i <= 8; i++ ) {
        x = cellgeo.giveVertexCoordinates(i)->at(1);
        y = cellgeo.giveVertexCoordinates(i)->at(2);
        z = cellgeo.giveVertexCoordinates(i)->at(3);

        jacobianMatrix.at(1, 1) += dx.at(i) * x;
        jacobianMatrix.at(1, 2) += dx.at(i) * y;
        jacobianMatrix.at(1, 3) += dx.at(i) * z;
        jacobianMatrix.at(2, 1) += dy.at(i) * x;
        jacobianMatrix.at(2, 2) += dy.at(i) * y;
        jacobianMatrix.at(2, 3) += dy.at(i) * z;
        jacobianMatrix.at(3, 1) += dz.at(i) * x;
        jacobianMatrix.at(3, 2) += dz.at(i) * y;
        jacobianMatrix.at(3, 3) += dz.at(i) * z;
    }
}


void
FEI3dHexaLin :: giveDerivativeKsi(FloatArray &dx, double v, double w)
{
    dx.at(1)  = -0.125 * ( 1. - v ) * ( 1. + w );
    dx.at(2)  = -0.125 * ( 1. + v ) * ( 1. + w );
    dx.at(3)  =  0.125 * ( 1. + v ) * ( 1. + w );
    dx.at(4)  =  0.125 * ( 1. - v ) * ( 1. + w );
    dx.at(5)  = -0.125 * ( 1. - v ) * ( 1. - w );
    dx.at(6)  = -0.125 * ( 1. + v ) * ( 1. - w );
    dx.at(7)  =  0.125 * ( 1. + v ) * ( 1. - w );
    dx.at(8)  =  0.125 * ( 1. - v ) * ( 1. - w );
}

void
FEI3dHexaLin :: giveDerivativeEta(FloatArray &dy, double u, double w)
{
    dy.at(1)  = -0.125 * ( 1. - u ) * ( 1. + w );
    dy.at(2)  =  0.125 * ( 1. - u ) * ( 1. + w );
    dy.at(3)  =  0.125 * ( 1. + u ) * ( 1. + w );
    dy.at(4)  = -0.125 * ( 1. + u ) * ( 1. + w );
    dy.at(5)  = -0.125 * ( 1. - u ) * ( 1. - w );
    dy.at(6)  =  0.125 * ( 1. - u ) * ( 1. - w );
    dy.at(7)  =  0.125 * ( 1. + u ) * ( 1. - w );
    dy.at(8)  = -0.125 * ( 1. + u ) * ( 1. - w );
}

void
FEI3dHexaLin :: giveDerivativeDzeta(FloatArray &dz, double u, double v)
{
    dz.at(1)  =  0.125 * ( 1. - u ) * ( 1. - v );
    dz.at(2)  =  0.125 * ( 1. - u ) * ( 1. + v );
    dz.at(3)  =  0.125 * ( 1. + u ) * ( 1. + v );
    dz.at(4)  =  0.125 * ( 1. + u ) * ( 1. - v );
    dz.at(5)  = -0.125 * ( 1. - u ) * ( 1. - v );
    dz.at(6)  = -0.125 * ( 1. - u ) * ( 1. + v );
    dz.at(7)  = -0.125 * ( 1. + u ) * ( 1. + v );
    dz.at(8)  = -0.125 * ( 1. + u ) * ( 1. - v );
}
} // end namespace oofem
