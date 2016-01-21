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

#include "fei2dtrquad.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
void
FEI2dTrQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1. - l1 - l2;

    answer = {
        ( 2. * l1 - 1. ) * l1,
        ( 2. * l2 - 1. ) * l2,
        ( 2. * l3 - 1. ) * l3,
        4. * l1 * l2,
        4. * l2 * l3,
        4. * l3 * l1
    };
}

double
FEI2dTrQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
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
FEI2dTrQuad :: evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x1, x2, x3, y1, y2, y3, y23, x32, y31, x13, area;

    answer.resize(6, 3);

    x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    x3 = cellgeo.giveVertexCoordinates(3)->at(xind);

    y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    y3 = cellgeo.giveVertexCoordinates(3)->at(yind);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    y23 = ( y2 - y3 ) / ( 2. * area );
    x32 = ( x3 - x2 ) / ( 2. * area );

    y31 = ( y3 - y1 ) / ( 2. * area );
    x13 = ( x1 - x3 ) / ( 2. * area );

    answer.at(1, 1) = 4 * y23 * y23;
    answer.at(1, 2) = 4 * x32 * x32;
    answer.at(1, 3) = 4 * y23 * x32;

    answer.at(2, 1) = 4 * y31 * y31;
    answer.at(2, 2) = 4 * x13 * x13;
    answer.at(2, 3) = 4 * y31 * x13;

    answer.at(3, 1) = 4 * y23 * y23 + 8 * y31 * y23 + 4 * y31 * y31;
    answer.at(3, 2) = 4 * x32 * x32 + 8 * x13 * x32 + 4 * x13 * x13;
    answer.at(3, 3) = 4 * y23 * x32 + 4 * y31 * x32 + 4 * y23 * x13 + 4 * y31 * x13;

    answer.at(4, 1) = 8 * y31 * y23;
    answer.at(4, 2) = 8 * x13 * x32;
    answer.at(4, 3) = 4 * y31 * x32 + 4 * y23 * x13;

    answer.at(5, 1) = ( -8 ) * y31 * y23 + ( -8 ) * y31 * y31;
    answer.at(5, 2) = ( -8 ) * x13 * x32 + ( -8 ) * x13 * x13;
    answer.at(5, 3) = ( -4 ) * y31 * x32 + ( -4 ) * y23 * x13 + ( -8 ) * y31 * x13;

    answer.at(6, 1) = ( -8 ) * y23 * y23 + ( -8 ) * y31 * y23;
    answer.at(6, 2) = ( -8 ) * x32 * x32 + ( -8 ) * x13 * x32;
    answer.at(6, 3) = ( -8 ) * y23 * x32 + ( -4 ) * y31 * x32 + ( -4 ) * y23 * x13;
}



void
FEI2dTrQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);

    answer.resize(2);
    answer.zero();
    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(xind);
        answer.at(2) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(yind);
    }
}


double FEI2dTrQuad :: giveCharacteristicLength(const FEICellGeometry &cellgeo) const
{
    return sqrt( this->giveArea(cellgeo) );
}


#define POINT_TOL 1.e-3

int
FEI2dTrQuad :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    FloatArray res, delta, guess, lcoords_guess;
    FloatMatrix jac;
    double convergence_limit, error = 0.0;

    // find a suitable convergence limit
    convergence_limit = 1e-6 * this->giveCharacteristicLength(cellgeo);

    // setup initial guess
    lcoords_guess.resize( 2 );
    lcoords_guess.zero();

    // apply Newton-Raphson to solve the problem
    for ( int nite = 0; nite < 10; nite++ ) {
        // compute the residual
        this->local2global(guess, lcoords_guess, cellgeo);
        res = {gcoords(0) - guess(0), gcoords(1) - guess(1)};

        // check for convergence
        error = res.computeNorm();
        if ( error < convergence_limit ) {
            break;
        }

        // compute the corrections
        this->giveJacobianMatrixAt(jac, lcoords_guess, cellgeo);
        jac.solveForRhs(res, delta, true);

        // update guess
        lcoords_guess.add(delta);
    }
    if ( error > convergence_limit ) { // Imperfect, could give false negatives.
        //OOFEM_WARNING("Failed convergence");
        answer = {1. / 3., 1. / 3., 1. / 3.};
        return false;
    }

    answer = { lcoords_guess(0), lcoords_guess(1), 1.0 - lcoords_guess(0) - lcoords_guess(1) };

    bool inside = true;
    for ( int i = 0; i < 3; i++ ) {
        if ( answer(i) < ( 0. - POINT_TOL ) ) {
            answer(i) = 0.;
            inside = false;
        } else if ( answer(i) > ( 1. + POINT_TOL ) ) {
            answer(i) = 1.;
            inside = false;
        }
    }

    return inside;
}

void
FEI2dTrQuad :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double n3, ksi = lcoords.at(1);
    n3 = 1. - ksi * ksi;

    answer = { ( 1. - ksi - n3 ) * 0.5, ( 1. + ksi - n3 ) * 0.5, n3 };
}

void
FEI2dTrQuad :: edgeEvaldNds(FloatArray &answer, int iedge,
                            const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    // I think it at least should return both dNds and J. Both are almost always needed.
    // In fact, dxdxi is also needed sometimes (surface tension)
#if 0
    IntArray edgeNodes;
    FloatArray dNdxi(3);
    FloatArray dxdxi(2);
    double xi = lcoords.at(1);
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    dNdxi.at(1) = xi - 0.5;
    dNdxi.at(2) = xi + 0.5;
    dNdxi.at(3) = -2 * xi;

    dxdxi.at(1) = dNdxi.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(xind) +
    dNdxi.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(xind) +
    dNdxi.at(3) * cellgeo.giveVertexCoordinates( edgeNodes.at(3) )->at(xind);
    dxdxi.at(2) = dNdxi.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(yind) +
    dNdxi.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(yind) +
    dNdxi.at(3) * cellgeo.giveVertexCoordinates( edgeNodes.at(3) )->at(yind);

    double J = dxdxi.computeNorm();
    answer = dNdxi;
    answer.times(1 / J);
    return J;

#endif
    double xi = lcoords.at(1);
    double J = edgeGiveTransformationJacobian(iedge, lcoords, cellgeo);
    answer = {
        ( xi - 0.5 ) / J,
        ( xi + 0.5 ) / J,
        -2 * xi / J
    };
}

double FEI2dTrQuad :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
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
FEI2dTrQuad :: edgeLocal2global(FloatArray &answer, int iedge,
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
FEI2dTrQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    int aNode = 0, bNode = 0, cNode = 0;
    edgeNodes.resize(3);

    if ( iedge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
        cNode = 4;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
        cNode = 5;
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        aNode = 3;
        bNode = 1;
        cNode = 6;
    } else {
        OOFEM_ERROR("wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;
    edgeNodes.at(3) = cNode;
}


void
FEI2dTrQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
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
FEI2dTrQuad :: giveDerivatives(FloatMatrix &dn, const FloatArray &lc)
{
    double l1, l2, l3;
    l1 = lc.at(1);
    l2 = lc.at(2);
    l3 = 1.0 - l1 - l2;

    dn.resize(6, 2);

    dn.at(1, 1) =  4.0 * l1 - 1.0;
    dn.at(2, 1) =  0.0;
    dn.at(3, 1) = -1.0 * ( 4.0 * l3 - 1.0 );
    dn.at(4, 1) =  4.0 * l2;
    dn.at(5, 1) = -4.0 * l2;
    dn.at(6, 1) =  4.0 * l3 - 4.0 * l1;

    dn.at(1, 2) =  0.0;
    dn.at(2, 2) =  4.0 * l2 - 1.0;
    dn.at(3, 2) = -1.0 * ( 4.0 * l3 - 1.0 );
    dn.at(4, 2) =  4.0 * l1;
    dn.at(5, 2) =  4.0 * l3 - 4.0 * l2;
    dn.at(6, 2) = -4.0 * l1;
}

double
FEI2dTrQuad :: giveArea(const FEICellGeometry &cellgeo) const
{
    const FloatArray *p;
    double x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6;

    p = cellgeo.giveVertexCoordinates(1);
    x1 = p->at(1);
    y1 = p->at(2);
    p = cellgeo.giveVertexCoordinates(2);
    x2 = p->at(1);
    y2 = p->at(2);
    p = cellgeo.giveVertexCoordinates(3);
    x3 = p->at(1);
    y3 = p->at(2);
    p = cellgeo.giveVertexCoordinates(4);
    x4 = p->at(1);
    y4 = p->at(2);
    p = cellgeo.giveVertexCoordinates(5);
    x5 = p->at(1);
    y5 = p->at(2);
    p = cellgeo.giveVertexCoordinates(6);
    x6 = p->at(1);
    y6 = p->at(2);

    return fabs( ( 4 * ( -( x4 * y1 ) + x6 * y1 + x4 * y2 - x5 * y2 + x5 * y3 - x6 * y3 ) + x2 * ( y1 - y3 - 4 * y4 + 4 * y5 ) +
            x1 * ( -y2 + y3 + 4 * y4 - 4 * y6 ) + x3 * ( -y1 + y2 - 4 * y5 + 4 * y6 ) ) / 6 );
}

double
FEI2dTrQuad :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo)
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
FEI2dTrQuad :: giveIntegrationRule(int order)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order + 2);
    iRule->SetUpPointsOnTriangle(points, _Unknown);
    return iRule;
}
} // end namespace oofem
