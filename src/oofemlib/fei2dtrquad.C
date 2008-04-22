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

#include "fei2dtrquad.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"

void
FEI2dTrQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, double time)
{
    double l1, l2, l3;
    answer.resize(6);

    answer.at(1) = ( 2. * l1 - 1. ) * l1;
    answer.at(2) = ( 2. * l2 - 1. ) * l2;
    answer.at(3) = ( 2. * l3 - 1. ) * l3;
    answer.at(4) = 4. * l1 * l2;
    answer.at(5) = 4. * l2 * l3;
    answer.at(6) = 4. * l3 * l1;

    return;
}

void
FEI2dTrQuad :: evaldNdx(FloatMatrix &answer, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    answer.resize(6, 2);
    int i;
    FloatMatrix jacobianMatrix(2, 2), inv(2, 2);
    FloatArray nx(6), ny(6);

    this->giveJacobianMatrixAt(jacobianMatrix, coords, lcoords);
    inv.beInverseOf(jacobianMatrix);

    this->giveDerivativeXi(nx, lcoords);
    this->giveDerivativeEta(ny, lcoords);

    for ( i = 1; i <= 6; i++ ) {
        answer.at(i, 1) = nx.at(i) * inv.at(1, 1) + ny.at(i) * inv.at(1, 2);
        answer.at(i, 2) = nx.at(i) * inv.at(2, 1) + ny.at(i) * inv.at(2, 2);
    }
}

void
FEI2dTrQuad :: local2global(FloatArray &answer, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    int i;
    FloatArray n(6);
    answer.resize(2);
    answer.zero();

    this->evalN(n, lcoords, time);

    for ( i = 1; i <= 6; i++ ) {
        answer.at(1) += n.at(i) * coords [ i - 1 ]->at(xind);
        answer.at(2) += n.at(i) * coords [ i - 1 ]->at(yind);
    }
}

#define POINT_TOL 1.e-3

int
FEI2dTrQuad :: global2local(FloatArray &answer, const FloatArray **nc, const FloatArray &coords, double time)
{
    FloatArray lc(2);
    FloatArray r(2), n(6), dksi, deta, delta;
    FloatMatrix p(2, 2);
    double l1, l2, l3, x, y;
    int i, nite = 0;

    // setup initial guess
    lc.resize(2);
    lc.at(1) = lc.at(2) = 1. / 3.;

    // apply Newton-Raphson to solve the problem
    do {
        if ( ( ++nite ) > 10 ) {
            //_error ("computeLocalCoords: no convergence after 10 iterations");
            return 0;
        }

        // compute the residual
        l1 = lc.at(1);
        l2 = lc.at(2);
        l3 = 1.0 - l1 - l2;

        n.at(1) = ( 2. * l1 - 1. ) * l1;
        n.at(2) = ( 2. * l2 - 1. ) * l2;
        n.at(3) = ( 2. * l3 - 1. ) * l3;
        n.at(4) = 4. * l1 * l2;
        n.at(5) = 4. * l2 * l3;
        n.at(6) = 4. * l3 * l1;

        r.at(1) = coords.at(1);
        r.at(2) = coords.at(2);
        for ( i = 1; i <= 6; i++ ) {
            r.at(1) -= n.at(i) * nc [ i - 1 ]->at(xind);
            r.at(2) -= n.at(i) * nc [ i - 1 ]->at(yind);
        }

        // check for convergence
        if ( sqrt( dotProduct(r, r, 2) ) < 1.e-10 ) {
            break;
        }

        // compute the corrections
        this->giveDerivativeXi(dksi, lc);
        this->giveDerivativeEta(deta, lc);

        p.zero();
        for ( i = 1; i <= 6; i++ ) {
            x = nc [ i - 1 ]->at(xind);
            y = nc [ i - 1 ]->at(yind);

            p.at(1, 1) += dksi.at(i) * x;
            p.at(1, 2) += deta.at(i) * x;
            p.at(2, 1) += dksi.at(i) * y;
            p.at(2, 2) += deta.at(i) * y;
        }

        // solve for corrections
        p.solveForRhs(r, delta);
        // update guess
        lc.add(delta);
    } while ( 1 );

    answer.resize(3);
    answer.at(1) = lc(1);
    answer.at(2) = lc(2);
    answer.at(3) = 1.0 - lc.at(1) - lc.at(2);

    for ( i = 1; i <= 2; i++ ) {
        if ( lc.at(i) < ( 0. - POINT_TOL ) ) {
            return 0;
        }

        if ( lc.at(i) > ( 1. + POINT_TOL ) ) {
            return 0;
        }
    }

    return 1;
}


double
FEI2dTrQuad :: giveTransformationJacobian(const FloatArray **coords, const FloatArray &lcoords, double time)
{
    FloatMatrix jacobianMatrix(2, 2);

    this->giveJacobianMatrixAt(jacobianMatrix, coords, lcoords);
    return jacobianMatrix.giveDeterminant();
}


void
FEI2dTrQuad :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, double time)
{
    double n3, ksi = lcoords.at(1);
    answer.resize(3);

    n3 = 1. - ksi * ksi;
    answer.at(1) = ( 1. - ksi ) * 0.5 - 0.5 * n3;
    answer.at(2) = ( 1. + ksi ) * 0.5 - 0.5 * n3;
    answer.at(3) = n3;
}

void
FEI2dTrQuad :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                            const FloatArray **coords, const FloatArray &lcoords, double time)
{
    OOFEM_ERROR("FEI2dTrQuad :: edgeEvaldNdx: not implemented");
}

void
FEI2dTrQuad :: edgeLocal2global(FloatArray &answer, int iedge,
                                const FloatArray **coords, const FloatArray &lcoords, double time)
{
    IntArray edgeNodes;
    FloatArray n;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(n, lcoords, time);

    answer.resize(2);
    answer.at(1) = ( n.at(1) * coords [ edgeNodes.at(1) - 1 ]->at(xind) +
                    n.at(2) * coords [ edgeNodes.at(2) - 1 ]->at(xind) +
                    n.at(3) * coords [ edgeNodes.at(3) - 1 ]->at(xind) );
    answer.at(2) = ( n.at(1) * coords [ edgeNodes.at(1) - 1 ]->at(yind) +
                    n.at(2) * coords [ edgeNodes.at(2) - 1 ]->at(yind) +
                    n.at(3) * coords [ edgeNodes.at(3) - 1 ]->at(yind) );
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
        OOFEM_ERROR2("FEI2dTrQuad :: computeEdgeMapping: wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;
    edgeNodes.at(3) = cNode;
}

double
FEI2dTrQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    OOFEM_ERROR("FEI2dTrQuad :: edgeGiveTransformationJacobian: not implemented");
    return 0.0;
}



void
FEI2dTrQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray **coords, const FloatArray &lcoords)
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

    for ( i = 1; i <= 6; i++ ) {
        x = coords [ i - 1 ]->at(xind);
        y = coords [ i - 1 ]->at(yind);

        jacobianMatrix.at(1, 1) += dxi.at(i) * x;
        jacobianMatrix.at(1, 2) += dxi.at(i) * y;
        jacobianMatrix.at(2, 1) += deta.at(i) * x;
        jacobianMatrix.at(2, 2) += deta.at(i) * y;
    }
}


void
FEI2dTrQuad :: giveDerivativeXi(FloatArray &n, const FloatArray &lc)
{
    double l1, l2, l3;

    l1 = lc.at(1);
    l2 = lc.at(2);
    l3 = 1.0 - l1 - l2;

    n.resize(6);

    n.at(1) = 4.0 * l1 - 1.0;
    n.at(2) = 0.0;
    n.at(3) = -1.0 * ( 4.0 * l3 - 1.0 );
    n.at(4) = 4.0 * l2;
    n.at(5) = -4.0 * l2;
    n.at(6) = 4.0 * l3 - 4.0 * l1;
}

void
FEI2dTrQuad :: giveDerivativeEta(FloatArray &n, const FloatArray &lc)
{
    double l1, l2, l3;

    l1 = lc.at(1);
    l2 = lc.at(2);
    l3 = 1.0 - l1 - l2;

    n.resize(6);

    n.at(1) = 0.0;
    n.at(2) = 4.0 * l2 - 1.0;
    n.at(3) = -1.0 * ( 4.0 * l3 - 1.0 );
    n.at(4) = 4.0 * l1;
    n.at(5) = 4.0 * l3 - 4.0 * l2;
    n.at(6) = -4.0 * l1;
}
