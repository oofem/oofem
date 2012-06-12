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

#include "fei2dquadquad.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {

double
FEI2dQuadQuad :: giveArea(const FEICellGeometry &cellgeo) const
{
    double x1, x2, x3, x4, y1, y2, y3, y4;
    double x85, x56, x67, x78, y85, y56, y67, y78;

    const FloatArray *node1 = cellgeo.giveVertexCoordinates(1);
    const FloatArray *node2 = cellgeo.giveVertexCoordinates(2);
    const FloatArray *node3 = cellgeo.giveVertexCoordinates(3);
    const FloatArray *node4 = cellgeo.giveVertexCoordinates(4);
    const FloatArray *node5 = cellgeo.giveVertexCoordinates(5);
    const FloatArray *node6 = cellgeo.giveVertexCoordinates(6);
    const FloatArray *node7 = cellgeo.giveVertexCoordinates(7);
    const FloatArray *node8 = cellgeo.giveVertexCoordinates(8);

    x1 = node1->at(xind);
    x2 = node2->at(xind);
    x3 = node3->at(xind);
    x4 = node4->at(xind);

    y1 = node1->at(yind);
    y2 = node2->at(yind);
    y3 = node3->at(yind);
    y4 = node4->at(yind);

    x85 = node8->at(xind) - node5->at(xind);
    x56 = node5->at(xind) - node6->at(xind);
    x67 = node6->at(xind) - node7->at(xind);
    x78 = node7->at(xind) - node8->at(xind);

    y85 = node8->at(yind) - node5->at(yind);
    y56 = node5->at(yind) - node6->at(yind);
    y67 = node6->at(yind) - node7->at(yind);
    y78 = node7->at(yind) - node8->at(yind);

    double p1 = (x2-x4)*(y1-y3) - (x1-x3)*(y2-y4);
    double p2 = y1*x85 + y2*x56 + y3*x67 + y4*x78 - x1*y85 - x2*y56 - x3*y67 - x4*x78;

    return fabs(p1 + p2*4.0)/6.; // Expression derived with mathematica, but not verified in any computations
}

void
FEI2dQuadQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    /* Local Node Numbering
     *
     *  4----7--- 3
     *  |         |
     *  |         |
     *  8         6
     *  |         |
     *  |         |
     *  1----5----2
     */

    double ksi, eta;

    answer.resize(8);

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    answer.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25 * ( ksi + eta - 1. );
    answer.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25 * ( -ksi + eta - 1. );
    answer.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25 * ( -ksi - eta - 1. );
    answer.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25 * ( ksi - eta - 1. );
    answer.at(5) = 0.5 * ( 1. - ksi * ksi ) * ( 1. + eta );
    answer.at(6) = 0.5 * ( 1. - ksi ) * ( 1. - eta * eta );
    answer.at(7) = 0.5 * ( 1. - ksi * ksi ) * ( 1. - eta );
    answer.at(8) = 0.5 * ( 1. + ksi ) * ( 1. - eta * eta );
}

void
FEI2dQuadQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(8, 2);
    int i;
    FloatMatrix jacobianMatrix(2, 2), inv(2, 2);
    FloatArray dndxi(8), dndeta(8);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    inv.beInverseOf(jacobianMatrix);

    this->giveDerivativeXi(dndxi, lcoords);
    this->giveDerivativeEta(dndeta, lcoords);

    for ( i = 1; i <= 8; i++ ) {
        answer.at(i, 1) = dndxi.at(i) * inv.at(1, 1) + dndeta.at(i) * inv.at(1, 2);
        answer.at(i, 2) = dndxi.at(i) * inv.at(2, 1) + dndeta.at(i) * inv.at(2, 2);
    }
}

void
FEI2dQuadQuad :: local2global(FloatArray &answer, const FloatArray &lcoords,  const FEICellGeometry &cellgeo)
{
    int i;
    FloatArray n(8);
    answer.resize(2);
    answer.zero();

    this->evalN(n, lcoords, cellgeo);

    for ( i = 1; i <= 8; i++ ) {
        answer.at(1) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(xind);
        answer.at(2) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(yind);
    }
}

double FEI2dQuadQuad :: giveCharacteristicLength(const FEICellGeometry &cellgeo) const
{
    const FloatArray *n1 = cellgeo.giveVertexCoordinates(1);
    const FloatArray *n2 = cellgeo.giveVertexCoordinates(3);
    return n1->distance(n2);
}

#define POINT_TOL 1.e-3

int
FEI2dQuadQuad :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    FloatArray res, delta, guess;
    FloatMatrix jac;
    double convergence_limit, error;

    // find a suitable convergence limit
    convergence_limit = 1e-6 * this->giveCharacteristicLength(cellgeo);

    // setup initial guess
    answer.resize(gcoords.giveSize());
    answer.zero();

    // apply Newton-Raphson to solve the problem
    for (int nite = 0; nite < 10; nite++) {
        // compute the residual
        this->local2global(guess, answer, cellgeo);
        res.beDifferenceOf(gcoords, guess);

        // check for convergence
        error = res.computeNorm();
        if ( error < convergence_limit ) {
            break;
        }

        // compute the corrections
        this->giveJacobianMatrixAt(jac, answer, cellgeo);
        jac.solveForRhs(res, delta, true);

        // update guess
        answer.add(delta);
    }
    if ( error > convergence_limit) { // Imperfect, could give false negatives.
        //OOFEM_ERROR ("global2local: no convergence after 10 iterations");
        return false;
    }

    // check limits for each local coordinate [-1,1] for quadrilaterals. (different for other elements, typically [0,1]).
    for (int i = 1; i <= answer.giveSize(); i++ ) {
        if ( fabs( answer.at(i) ) > ( 1. + POINT_TOL ) ) {
            return false;
        }
    }

    return true;
}


double
FEI2dQuadQuad :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(2, 2);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant();
}


void
FEI2dQuadQuad :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    /*
     *   1-------3-------2
     */

    double n3, ksi = lcoords.at(1);
    answer.resize(3);

    n3 = 1. - ksi * ksi;
    answer.at(1) = ( 1. - ksi ) * 0.5 - 0.5 * n3;
    answer.at(2) = ( 1. + ksi ) * 0.5 - 0.5 * n3;
    answer.at(3) = n3;
}

void
FEI2dQuadQuad :: edgeEvaldNds(FloatArray &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI2dQuadQuad :: edgeEvaldNds: not implemented");
}

void
FEI2dQuadQuad :: edgeLocal2global(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray n;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(n, lcoords, cellgeo);

    answer.resize(2);
    answer.at(1) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(xind) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(xind) +
                    n.at(3) * cellgeo.giveVertexCoordinates( edgeNodes.at(3) )->at(xind) );
    answer.at(2) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(yind) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(yind) +
                    n.at(3) * cellgeo.giveVertexCoordinates( edgeNodes.at(3) )->at(yind) );
}


void
FEI2dQuadQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    int aNode = 0, bNode = 0, cNode = 0;
    edgeNodes.resize(3);

    if ( iedge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
        cNode = 5;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
        cNode = 6;
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        aNode = 3;
        bNode = 4;
        cNode = 7;
    } else if ( iedge == 4 ) { // edge between nodes 3 4
        aNode = 4;
        bNode = 1;
        cNode = 8;
    } else {
        OOFEM_ERROR2("FEI2dQuadQuad :: computeEdgeMapping: wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;
    edgeNodes.at(3) = cNode;
}

void FEI2dQuadQuad :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    double xi = lcoords(0);
    double dN1dxi = -0.5 + xi;
    double dN2dxi =  0.5 + xi;
    double dN3dxi = -2.0 * xi;

    normal.resize(2);

    normal.at(1) =-dN1dxi*cellgeo.giveVertexCoordinates(edgeNodes.at(1))->at(yind) +
                  -dN2dxi*cellgeo.giveVertexCoordinates(edgeNodes.at(2))->at(yind) +
                  -dN3dxi*cellgeo.giveVertexCoordinates(edgeNodes.at(3))->at(yind);

    normal.at(2) = dN1dxi*cellgeo.giveVertexCoordinates(edgeNodes.at(1))->at(xind) +
                   dN2dxi*cellgeo.giveVertexCoordinates(edgeNodes.at(2))->at(xind) +
                   dN3dxi*cellgeo.giveVertexCoordinates(edgeNodes.at(3))->at(xind);

    normal.normalize();
}

double
FEI2dQuadQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    double xi = lcoords.at(1);
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    FloatArray dNdxi(3);
    dNdxi.at(1) = xi-0.5;
    dNdxi.at(2) = xi+0.5;
    dNdxi.at(3) = -2*xi;

    FloatArray dxdxi(2);
    dxdxi.at(1) = dNdxi.at(1)*cellgeo.giveVertexCoordinates(edgeNodes.at(1))->at(xind) +
                  dNdxi.at(2)*cellgeo.giveVertexCoordinates(edgeNodes.at(2))->at(xind) +
                  dNdxi.at(3)*cellgeo.giveVertexCoordinates(edgeNodes.at(3))->at(xind);
    dxdxi.at(2) = dNdxi.at(1)*cellgeo.giveVertexCoordinates(edgeNodes.at(1))->at(yind) +
                  dNdxi.at(2)*cellgeo.giveVertexCoordinates(edgeNodes.at(2))->at(yind) +
                  dNdxi.at(3)*cellgeo.giveVertexCoordinates(edgeNodes.at(3))->at(yind);

    return dxdxi.computeNorm();
}


void
FEI2dQuadQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
// Computes it if it does not exist yet.
{
    int i;
    double x, y;
    FloatArray dxi, deta;

    jacobianMatrix.resize(2, 2);
    jacobianMatrix.zero();

    this->giveDerivativeXi(dxi, lcoords);
    this->giveDerivativeEta(deta, lcoords);

    for ( i = 1; i <= 8; i++ ) {
        x = cellgeo.giveVertexCoordinates(i)->at(xind);
        y = cellgeo.giveVertexCoordinates(i)->at(yind);

        jacobianMatrix.at(1, 1) += dxi.at(i) * x;
        jacobianMatrix.at(1, 2) += dxi.at(i) * y;
        jacobianMatrix.at(2, 1) += deta.at(i) * x;
        jacobianMatrix.at(2, 2) += deta.at(i) * y;
    }
}


void
FEI2dQuadQuad :: giveDerivativeXi(FloatArray &n, const FloatArray &lc)
{
    double ksi, eta;
    ksi = lc.at(1);
    eta = lc.at(2);
    n.resize(8);


    n.at(1) =  0.25 * ( 1. + eta ) * ( 2.0 * ksi + eta );
    n.at(2) = -0.25 * ( 1. + eta ) * ( -2.0 * ksi + eta );
    n.at(3) = -0.25 * ( 1. - eta ) * ( -2.0 * ksi - eta );
    n.at(4) =  0.25 * ( 1. - eta ) * ( 2.0 * ksi - eta );
    n.at(5) = -ksi * ( 1. + eta );
    n.at(6) = -0.5 * ( 1. - eta * eta );
    n.at(7) = -ksi * ( 1. - eta );
    n.at(8) =  0.5 * ( 1. - eta * eta );
}

void
FEI2dQuadQuad :: giveDerivativeEta(FloatArray &n, const FloatArray &lc)
{
    double ksi, eta;
    ksi = lc.at(1);
    eta = lc.at(2);
    n.resize(8);

    n.at(1) =  0.25 * ( 1. + ksi ) * ( 2.0 * eta + ksi );
    n.at(2) =  0.25 * ( 1. - ksi ) * ( 2.0 * eta - ksi );
    n.at(3) = -0.25 * ( 1. - ksi ) * ( -2.0 * eta - ksi );
    n.at(4) = -0.25 * ( 1. + ksi ) * ( -2.0 * eta + ksi );
    n.at(5) =  0.5 * ( 1. - ksi * ksi );
    n.at(6) = -eta * ( 1. - ksi );
    n.at(7) = -0.5 * ( 1. - ksi * ksi );
    n.at(8) = -eta * ( 1. + ksi );
}
} // end namespace oofem
