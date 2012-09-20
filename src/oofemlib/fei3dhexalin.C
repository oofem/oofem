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
    FloatMatrix jacobianMatrix, inv, dNduvw, coords;
    
    this->giveLocalDerivative(dNduvw, lcoords);
    coords.resize(3, 8);
    for ( int i = 1; i <= 8; i++ ) {
        coords.setColumn(*cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
    inv.beInverseOf(jacobianMatrix);
    
    answer.beProductOf(dNduvw, inv);
    //return detJ;
}

void
FEI3dHexaLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);

    answer.resize(0);
    for (int i = 1; i <= 8; i++ ) {
        answer.add(n.at(i), *cellgeo.giveVertexCoordinates(i));
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
            answer.zero();
            return false;
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
    bool inside = true;
    for ( int i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( -1. - POINT_TOL ) ) {
            answer.at(i) = -1.;
            inside = false;
        } else if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            answer.at(i) = 1.;
            inside = false;
        }
    }

    return inside;
}


double
FEI3dHexaLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix;
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

    edgeNodes.setValues(2, aNode, bNode);
}

double
FEI3dHexaLin :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    return cellgeo.giveVertexCoordinates(edgeNodes.at(2))->distance(cellgeo.giveVertexCoordinates(edgeNodes.at(1)));
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
    IntArray nodes;
    FloatArray n;

    this->computeLocalSurfaceMapping(nodes, iedge);
    this->surfaceEvalN(n, lcoords, cellgeo);

    answer.resize(3);
    answer.at(1) = n.at(1) * cellgeo.giveVertexCoordinates( nodes.at(1) )->at(1) + n.at(2) *cellgeo.giveVertexCoordinates( nodes.at(2) )->at(1) +
                   n.at(3) * cellgeo.giveVertexCoordinates( nodes.at(3) )->at(1) + n.at(4) *cellgeo.giveVertexCoordinates( nodes.at(4) )->at(1);
    answer.at(2) = n.at(1) * cellgeo.giveVertexCoordinates( nodes.at(1) )->at(2) + n.at(2) *cellgeo.giveVertexCoordinates( nodes.at(2) )->at(2) +
                   n.at(3) * cellgeo.giveVertexCoordinates( nodes.at(3) )->at(2) + n.at(4) *cellgeo.giveVertexCoordinates( nodes.at(4) )->at(2);
    answer.at(3) = n.at(1) * cellgeo.giveVertexCoordinates( nodes.at(1) )->at(3) + n.at(2) *cellgeo.giveVertexCoordinates( nodes.at(2) )->at(3) +
                   n.at(3) * cellgeo.giveVertexCoordinates( nodes.at(3) )->at(3) + n.at(4) *cellgeo.giveVertexCoordinates( nodes.at(4) )->at(3);
}

double
FEI3dHexaLin :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray a, b, dNdksi(4), dNdeta(4);
    double ksi, eta;
    IntArray snodes;
    
    this->computeLocalSurfaceMapping(snodes, isurf);

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    // No need to divide by 1/4, we'll normalize anyway;
    dNdksi.at(1) =  ( 1. + eta );
    dNdksi.at(2) = -( 1. + eta );
    dNdksi.at(3) = -( 1. - eta );
    dNdksi.at(4) =  ( 1. - eta );

    dNdeta.at(1) =  ( 1. + ksi );
    dNdeta.at(2) =  ( 1. - ksi );
    dNdeta.at(3) = -( 1. - ksi );
    dNdeta.at(4) = -( 1. + ksi );

    for (int i = 1; i <= 4; ++i) {
        a.add(dNdksi.at(i), *cellgeo.giveVertexCoordinates(snodes.at(i)));
        b.add(dNdeta.at(i), *cellgeo.giveVertexCoordinates(snodes.at(i)));
    }
    
    answer.beVectorProductOf(a, b);
    return answer.normalize()*0.0625;
}

double
FEI3dHexaLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo)
{
    FloatArray normal;
    return this->surfaceEvalNormal(normal, isurf, lcoords, cellgeo);
}

void
FEI3dHexaLin :: computeLocalSurfaceMapping(IntArray &surfNodes, int isurf)
{
    int aNode = 0, bNode = 0, cNode = 0, dNode = 0;

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

    surfNodes.setValues(4, aNode, bNode, cNode, dNode);
}

void
FEI3dHexaLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
// Computes it if it does not exist yet.
{
    FloatMatrix dNduvw, coords;
    this->giveLocalDerivative(dNduvw, lcoords);
    coords.resize(3, 8);
    for ( int i = 1; i <= 8; i++ ) {
        coords.setColumn(*cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
}

void
FEI3dHexaLin :: giveLocalDerivative(FloatMatrix &dN, const FloatArray &lcoords)
{
    double u, v, w;
    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);
    
    dN.resize(8, 3);
    
    dN.at(1, 1) = -0.125 * ( 1. - v ) * ( 1. + w );
    dN.at(2, 1) = -0.125 * ( 1. + v ) * ( 1. + w );
    dN.at(3, 1) =  0.125 * ( 1. + v ) * ( 1. + w );
    dN.at(4, 1) =  0.125 * ( 1. - v ) * ( 1. + w );
    dN.at(5, 1) = -0.125 * ( 1. - v ) * ( 1. - w );
    dN.at(6, 1) = -0.125 * ( 1. + v ) * ( 1. - w );
    dN.at(7, 1) =  0.125 * ( 1. + v ) * ( 1. - w );
    dN.at(8, 1) =  0.125 * ( 1. - v ) * ( 1. - w );

    dN.at(1, 2) = -0.125 * ( 1. - u ) * ( 1. + w );
    dN.at(2, 2) =  0.125 * ( 1. - u ) * ( 1. + w );
    dN.at(3, 2) =  0.125 * ( 1. + u ) * ( 1. + w );
    dN.at(4, 2) = -0.125 * ( 1. + u ) * ( 1. + w );
    dN.at(5, 2) = -0.125 * ( 1. - u ) * ( 1. - w );
    dN.at(6, 2) =  0.125 * ( 1. - u ) * ( 1. - w );
    dN.at(7, 2) =  0.125 * ( 1. + u ) * ( 1. - w );
    dN.at(8, 2) = -0.125 * ( 1. + u ) * ( 1. - w );

    dN.at(1, 3) =  0.125 * ( 1. - u ) * ( 1. - v );
    dN.at(2, 3) =  0.125 * ( 1. - u ) * ( 1. + v );
    dN.at(3, 3) =  0.125 * ( 1. + u ) * ( 1. + v );
    dN.at(4, 3) =  0.125 * ( 1. + u ) * ( 1. - v );
    dN.at(5, 3) = -0.125 * ( 1. - u ) * ( 1. - v );
    dN.at(6, 3) = -0.125 * ( 1. - u ) * ( 1. + v );
    dN.at(7, 3) = -0.125 * ( 1. + u ) * ( 1. + v );
    dN.at(8, 3) = -0.125 * ( 1. + u ) * ( 1. - v );
}

} // end namespace oofem
