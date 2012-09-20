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

#include "fei2dquadlin.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {

double
FEI2dQuadLin :: giveArea(const FEICellGeometry &cellgeo) const
{
    const FloatArray *node1 = cellgeo.giveVertexCoordinates(1);
    const FloatArray *node2 = cellgeo.giveVertexCoordinates(2);
    const FloatArray *node3 = cellgeo.giveVertexCoordinates(3);
    const FloatArray *node4 = cellgeo.giveVertexCoordinates(4);

    double x13 = node1->at(xind) - node3->at(xind);
    double y13 = node1->at(yind) - node3->at(yind);
    double x24 = node2->at(xind) - node4->at(xind);
    double y24 = node2->at(yind) - node4->at(yind);

    return fabs(0.5*(x13*y24 - x24*y13));
}

void
FEI2dQuadLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
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
FEI2dQuadLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(2, 2), inv(2, 2);
    FloatArray nx(4), ny(4);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    inv.beInverseOf(jacobianMatrix);

    this->giveDerivativeKsi( nx, lcoords.at(2) );
    this->giveDerivativeEta( ny, lcoords.at(1) );

    answer.resize(4, 2);

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(i, 1) = nx.at(i) * inv.at(1, 1) + ny.at(i) * inv.at(1, 2);
        answer.at(i, 2) = nx.at(i) * inv.at(2, 1) + ny.at(i) * inv.at(2, 2);
    }
}

void
FEI2dQuadLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi, eta, n1, n2, n3, n4;

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    n1 = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
    n2 = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
    n3 = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
    n4 = ( 1. + ksi ) * ( 1. - eta ) * 0.25;

    answer.resize(2);
    answer.at(1) = n1 * cellgeo.giveVertexCoordinates(1)->at(xind) + n2 *cellgeo.giveVertexCoordinates(2)->at(xind) +
                   n3 * cellgeo.giveVertexCoordinates(3)->at(xind) + n4 *cellgeo.giveVertexCoordinates(4)->at(xind);
    answer.at(2) = n1 * cellgeo.giveVertexCoordinates(1)->at(yind) + n2 *cellgeo.giveVertexCoordinates(2)->at(yind) +
                   n3 * cellgeo.giveVertexCoordinates(3)->at(yind) + n4 *cellgeo.giveVertexCoordinates(4)->at(yind);
}

#define POINT_TOL 1.e-6

int
FEI2dQuadLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    double x1, x2, x3, x4, y1, y2, y3, y4, a1, a2, a3, a4, b1, b2, b3, b4;
    double a, b, c, ksi1, ksi2, ksi3, eta1 = 0.0, eta2 = 0.0, denom;
    int nroot;

    answer.resize(2);

    x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    x3 = cellgeo.giveVertexCoordinates(3)->at(xind);
    x4 = cellgeo.giveVertexCoordinates(4)->at(xind);

    y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    y3 = cellgeo.giveVertexCoordinates(3)->at(yind);
    y4 = cellgeo.giveVertexCoordinates(4)->at(yind);

    a1 = x1 + x2 + x3 + x4;
    a2 = x1 - x2 - x3 + x4;
    a3 = x1 + x2 - x3 - x4;
    a4 = x1 - x2 + x3 - x4;

    b1 = y1 + y2 + y3 + y4;
    b2 = y1 - y2 - y3 + y4;
    b3 = y1 + y2 - y3 - y4;
    b4 = y1 - y2 + y3 - y4;

    a = a2 * b4 - b2 * a4;
    b = a1 * b4 + a2 * b3 - a3 * b2 - b1 * a4 - b4 * 4.0 * coords.at(xind) + a4 * 4.0 * coords.at(yind);
    c = a1 * b3 - a3 * b1 - 4.0 * coords.at(xind) * b3 + 4.0 * coords.at(yind) * a3;

    // solve quadratic equation for ksi
    cubic(0.0, a, b, c, & ksi1, & ksi2, & ksi3, & nroot);

    if ( nroot == 0 ) {
        answer.zero();
        return false;
    }

    if ( nroot ) {
        denom = ( b3 + ksi1 * b4 );
        if ( fabs(denom) <= 1.0e-10 ) {
            eta1 = ( 4.0 * coords.at(xind) - a1 - ksi1 * a2 ) / ( a3 + ksi1 * a4 );
        } else {
            eta1 = ( 4.0 * coords.at(yind) - b1 - ksi1 * b2 ) / denom;
        }
    }

    if ( nroot > 1 ) {
        double diff_ksi1, diff_eta1, diff_ksi2, diff_eta2, diff1, diff2;

        denom = b3 + ksi2 * b4;
        if ( fabs(denom) <= 1.0e-10 ) {
            eta2 = ( 4.0 * coords.at(xind) - a1 - ksi2 * a2 ) / ( a3 + ksi2 * a4 );
        } else {
            eta2 = ( 4.0 * coords.at(yind) - b1 - ksi2 * b2 ) / denom;
        }

        // choose the one which seems to be closer to the parametric space (square <-1;1>x<-1;1>)
        diff_ksi1 = 0.0;
        if ( ksi1 > 1.0 ) {
            diff_ksi1 = ksi1 - 1.0;
        }

        if ( ksi1 < -1.0 ) {
            diff_ksi1 = ksi1 + 1.0;
        }

        diff_eta1 = 0.0;
        if ( eta1 > 1.0 ) {
            diff_eta1 = eta1 - 1.0;
        }

        if ( eta1 < -1.0 ) {
            diff_eta1 = eta1 + 1.0;
        }

        diff_ksi2 = 0.0;
        if ( ksi2 > 1.0 ) {
            diff_ksi2 = ksi2 - 1.0;
        }

        if ( ksi2 < -1.0 ) {
            diff_ksi2 = ksi2 + 1.0;
        }

        diff_eta2 = 0.0;
        if ( eta2 > 1.0 ) {
            diff_eta2 = eta2 - 1.0;
        }

        if ( eta2 < -1.0 ) {
            diff_eta2 = eta2 + 1.0;
        }

        diff1 = diff_ksi1 * diff_ksi1 + diff_eta1 * diff_eta1;
        diff2 = diff_ksi2 * diff_ksi2 + diff_eta2 * diff_eta2;

        // ksi2, eta2 seems to be closer
        if ( diff1 > diff2 ) {
            ksi1 = ksi2;
            eta1 = eta2;
        }
    }

    answer.at(1) = ksi1;
    answer.at(2) = eta1;

    // test if inside
    bool inside = true;
    for ( int i = 1; i <= 2; i++ ) {
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
FEI2dQuadLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(2, 2);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant();
}

void
FEI2dQuadLin :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords,  const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

void
FEI2dQuadLin :: edgeEvaldNds(FloatArray &answer, int iedge,
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

void
FEI2dQuadLin :: edgeLocal2global(FloatArray &answer, int iedge,
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


double
FEI2dQuadLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    return 0.5 * this->edgeComputeLength(edgeNodes, cellgeo);
}


void
FEI2dQuadLin :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
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
    } else {
        OOFEM_ERROR2("FEI2dQuadLin :: computeEdgeMapping: wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;
}

double
FEI2dQuadLin :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    double dx, dy;
    int nodeA, nodeB;

    nodeA = edgeNodes.at(1);
    nodeB = edgeNodes.at(2);

    dx = cellgeo.giveVertexCoordinates(nodeB)->at(xind) - cellgeo.giveVertexCoordinates(nodeA)->at(xind);
    dy = cellgeo.giveVertexCoordinates(nodeB)->at(yind) - cellgeo.giveVertexCoordinates(nodeA)->at(yind);
    return sqrt(dx * dx + dy * dy);
}

void
FEI2dQuadLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
// Computes it if it does not exist yet.
{
    double ksi, eta, x, y;
    FloatArray nx(4), ny(4);

    jacobianMatrix.resize(2, 2);
    jacobianMatrix.zero();

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    this->giveDerivativeKsi(nx, eta);
    this->giveDerivativeEta(ny, ksi);

    for ( int i = 1; i <= 4; i++ ) {
        x = cellgeo.giveVertexCoordinates(i)->at(xind);
        y = cellgeo.giveVertexCoordinates(i)->at(yind);

        jacobianMatrix.at(1, 1) += nx.at(i) * x;
        jacobianMatrix.at(1, 2) += nx.at(i) * y;
        jacobianMatrix.at(2, 1) += ny.at(i) * x;
        jacobianMatrix.at(2, 2) += ny.at(i) * y;
    }
}


void
FEI2dQuadLin :: giveDerivativeKsi(FloatArray &n, double eta)
{
    n.resize(4);

    n.at(1) =  0.25 * ( 1. + eta );
    n.at(2) = -0.25 * ( 1. + eta );
    n.at(3) = -0.25 * ( 1. - eta );
    n.at(4) =  0.25 * ( 1. - eta );
}

void
FEI2dQuadLin :: giveDerivativeEta(FloatArray &n, double ksi)
{
    n.resize(4);

    n.at(1) =  0.25 * ( 1. + ksi );
    n.at(2) =  0.25 * ( 1. - ksi );
    n.at(3) = -0.25 * ( 1. - ksi );
    n.at(4) = -0.25 * ( 1. + ksi );
}
} // end namespace oofem
