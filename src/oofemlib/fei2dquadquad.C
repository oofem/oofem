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

#include "fei2dquadquad.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

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

    double p1 = ( x2 - x4 ) * ( y1 - y3 ) - ( x1 - x3 ) * ( y2 - y4 );
    double p2 = y1 * x85 + y2 * x56 + y3 * x67 + y4 * x78 - x1 * y85 - x2 * y56 - x3 * y67 - x4 * y78;

    return fabs(p1 + p2 * 4.0) / 6.; // Expression derived with mathematica, but not verified in any computations
}

void
FEI2dQuadQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi, eta;

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    answer = {
        ( 1. + ksi ) * ( 1. + eta ) * 0.25 * ( ksi + eta - 1. ),
        ( 1. - ksi ) * ( 1. + eta ) * 0.25 * ( -ksi + eta - 1. ),
        ( 1. - ksi ) * ( 1. - eta ) * 0.25 * ( -ksi - eta - 1. ),
        ( 1. + ksi ) * ( 1. - eta ) * 0.25 * ( ksi - eta - 1. ),
        0.5 * ( 1. - ksi * ksi ) * ( 1. + eta ),
        0.5 * ( 1. - ksi ) * ( 1. - eta * eta ),
        0.5 * ( 1. - ksi * ksi ) * ( 1. - eta ),
        0.5 * ( 1. + ksi ) * ( 1. - eta * eta )
    };
}

double
FEI2dQuadQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(2, 2), inv, dn;

    this->giveDerivatives(dn, lcoords);
    for ( int i = 1; i <= dn.giveNumberOfRows(); i++ ) {
        double x = cellgeo.giveVertexCoordinates(i)->at(xind);
        double y = cellgeo.giveVertexCoordinates(i)->at(yind);

        jacobianMatrix.at(1, 1) += dn.at(i, 1) * x;
        jacobianMatrix.at(1, 2) += dn.at(i, 1) * y;
        jacobianMatrix.at(2, 1) += dn.at(i, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(i, 2) * y;
    }
    inv.beInverseOf(jacobianMatrix);

    answer.beProductTOf(dn, inv);
    return jacobianMatrix.giveDeterminant();
}

void
FEI2dQuadQuad :: local2global(FloatArray &answer, const FloatArray &lcoords,  const FEICellGeometry &cellgeo)
{
    FloatArray n;

    this->evalN(n, lcoords, cellgeo);

    answer.resize(2);
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
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
    double convergence_limit, error = 0.0;

    // find a suitable convergence limit
    convergence_limit = 1e-6 * this->giveCharacteristicLength(cellgeo);

    // setup initial guess
    answer.resize( gcoords.giveSize() );
    answer.zero();

    // apply Newton-Raphson to solve the problem
    for ( int nite = 0; nite < 10; nite++ ) {
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
    if ( error > convergence_limit ) { // Imperfect, could give false negatives.
        //OOFEM_ERROR("no convergence after 10 iterations");
        answer.zero();
        return false;
    }

    // check limits for each local coordinate [-1,1] for quadrilaterals
    bool inside = true;
    for ( int i = 1; i <= answer.giveSize(); i++ ) {
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

void
FEI2dQuadQuad :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    // 1-------3-------2

    double n3, ksi = lcoords.at(1);
    n3 = 1. - ksi * ksi;

    answer = { ( 1. - ksi - n3 ) * 0.5, ( 1. + ksi - n3 ) * 0.5, n3 };
}

void
FEI2dQuadQuad :: edgeEvaldNds(FloatArray &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer = { ksi - 0.5, ksi + 0.5, ksi * 2.0 };
}

void
FEI2dQuadQuad :: edgeLocal2global(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray n;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

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
        OOFEM_ERROR("wrong edge number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;
    edgeNodes.at(3) = cNode;
}

double FEI2dQuadQuad :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    double xi = lcoords(0);
    double dN1dxi = -0.5 + xi;
    double dN2dxi =  0.5 + xi;
    double dN3dxi = -2.0 * xi;

    normal.resize(2);

    normal.at(1) = dN1dxi * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(yind) +
    dN2dxi *cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(yind) +
    dN3dxi *cellgeo.giveVertexCoordinates( edgeNodes.at(3) )->at(yind);

    normal.at(2) = -dN1dxi *cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(xind) +
    - dN2dxi *cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(xind) +
    - dN3dxi *cellgeo.giveVertexCoordinates( edgeNodes.at(3) )->at(xind);

    return normal.normalize();
}

void
FEI2dQuadQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
{
    double x, y;
    FloatMatrix dn;

    jacobianMatrix.resize(2, 2);
    jacobianMatrix.zero();

    this->giveDerivatives(dn, lcoords);

    for ( int i = 1; i <= dn.giveNumberOfRows(); i++ ) {
        x = cellgeo.giveVertexCoordinates(i)->at(xind);
        y = cellgeo.giveVertexCoordinates(i)->at(yind);

        jacobianMatrix.at(1, 1) += dn.at(i, 1) * x;
        jacobianMatrix.at(1, 2) += dn.at(i, 1) * y;
        jacobianMatrix.at(2, 1) += dn.at(i, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(i, 2) * y;
    }
}


void
FEI2dQuadQuad :: giveDerivatives(FloatMatrix &dn, const FloatArray &lc)
{
    double ksi, eta;
    ksi = lc.at(1);
    eta = lc.at(2);
    dn.resize(8, 2);

    // dn/dxi
    dn.at(1, 1) =  0.25 * ( 1. + eta ) * ( 2.0 * ksi + eta );
    dn.at(2, 1) = -0.25 * ( 1. + eta ) * ( -2.0 * ksi + eta );
    dn.at(3, 1) = -0.25 * ( 1. - eta ) * ( -2.0 * ksi - eta );
    dn.at(4, 1) =  0.25 * ( 1. - eta ) * ( 2.0 * ksi - eta );
    dn.at(5, 1) = -ksi * ( 1. + eta );
    dn.at(6, 1) = -0.5 * ( 1. - eta * eta );
    dn.at(7, 1) = -ksi * ( 1. - eta );
    dn.at(8, 1) =  0.5 * ( 1. - eta * eta );

    dn.at(1, 2) =  0.25 * ( 1. + ksi ) * ( 2.0 * eta + ksi );
    dn.at(2, 2) =  0.25 * ( 1. - ksi ) * ( 2.0 * eta - ksi );
    dn.at(3, 2) = -0.25 * ( 1. - ksi ) * ( -2.0 * eta - ksi );
    dn.at(4, 2) = -0.25 * ( 1. + ksi ) * ( -2.0 * eta + ksi );
    dn.at(5, 2) =  0.5 * ( 1. - ksi * ksi );
    dn.at(6, 2) = -eta * ( 1. - ksi );
    dn.at(7, 2) = -0.5 * ( 1. - ksi * ksi );
    dn.at(8, 2) = -eta * ( 1. + ksi );
}


double
FEI2dQuadQuad :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo)
{
    IntArray eNodes;
    const FloatArray *node;
    double x1, x2, x3, y1, y2, y3;

    this->computeLocalEdgeMapping(eNodes, iEdge);

    node = cellgeo.giveVertexCoordinates( eNodes.at(1) );
    x1 = node->at(xind);
    y1 = node->at(yind);

    node = cellgeo.giveVertexCoordinates( eNodes.at(2) );
    x2 = node->at(xind);
    y2 = node->at(yind);

    node = cellgeo.giveVertexCoordinates( eNodes.at(3) );
    x3 = node->at(xind);
    y3 = node->at(yind);

    return -( x1 * y2 - x2 * y1 + 4 * ( x3 * ( y1 - y2 ) + y3 * ( x2 - x1 ) ) ) / 3.0;
}


IntegrationRule *
FEI2dQuadQuad :: giveIntegrationRule(int order)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Square, order + 4);
    iRule->SetUpPointsOnSquare(points, _Unknown);
    return iRule;
}
} // end namespace oofem
