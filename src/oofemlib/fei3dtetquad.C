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

#include "fei3dtetquad.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"

namespace oofem {
void
FEI3dTetQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x1 = lcoords(0);
    double x2 = lcoords(1);
    double x3 = lcoords(2);
    double x4 = 1.0 - x1 - x2 - x3;

    answer.resize(10);
    answer(0) = x1*(2*x1-1);
    answer(1) = x2*(2*x2-1);
    answer(2) = x3*(2*x3-1);
    answer(3) = x4*(2*x4-1);

    answer(4) = 4*x1*x2;
    answer(5) = 4*x2*x3;
    answer(6) = 4*x3*x1;
    answer(7) = 4*x1*x4;
    answer(8) = 4*x2*x4;
    answer(9) = 4*x3*x4;
}

void
FEI3dTetQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix;
    FloatMatrix dNdxi;
    FloatMatrix coords(10,3);
    this->evaldNdxi(dNdxi, lcoords);
    for (int i = 1; i <= 10; ++i) {
        const FloatArray *c = cellgeo.giveVertexCoordinates(i);
        coords.at(i,1) = c->at(1);
        coords.at(i,2) = c->at(2);
        coords.at(i,3) = c->at(3);
    }
    jacobianMatrix.beProductOf(dNdxi, coords);
    jacobianMatrix.solveForRhs(answer, dNdxi);
}

void
FEI3dTetQuad :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords)
{
    double x1 = lcoords(0);
    double x2 = lcoords(1);
    double x3 = lcoords(2);
    double x4 = 1.0 - x1 - x2 - x3;

    answer.resize(3,10);

    // dNj/dx1
    answer(0,0) = 4*x1 - 1;
    answer(0,1) = 0;
    answer(0,2) = 0;
    answer(0,3) = -4*x4 + 1;
    answer(0,4) = 4*x2;
    answer(0,5) = 0;
    answer(0,6) = 4*x3;
    answer(0,7) = 4*(x4-x1);
    answer(0,8) = -4*x2;
    answer(0,9) = -4*x3;

    // dNj/dx2
    answer(1,0) = 0;
    answer(1,1) = 4*x2 - 1;
    answer(1,2) = 0;
    answer(1,3) = -4*x4 + 1;
    answer(1,4) = 4*x1;
    answer(1,5) = 4*x3;
    answer(1,6) = 0;
    answer(1,7) = -4*x1;
    answer(1,8) = 4*(x4-x2);
    answer(1,9) = -4*x3;

    // dNj/dx3
    answer(2,0) = 0;
    answer(2,1) = 0;
    answer(2,2) = 4*x3 - 1;
    answer(2,3) = -4*x4 + 1;
    answer(2,4) = 0;
    answer(2,5) = 4*x2;
    answer(2,6) = 4*x1;
    answer(2,7) = -4*x1;
    answer(2,8) = -4*x2;
    answer(2,9) = 4*(x4-x3);
}


void
FEI3dTetQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray N;
    this->evalN(N,lcoords, cellgeo, time);
    answer.resize(0);
    for ( int i = 1; i <= N.giveSize(); i++ ) {
        answer.add( N(i), *cellgeo.giveVertexCoordinates(i));
    }
}

#define POINT_TOL 1e-6

int
FEI3dTetQuad :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    FloatArray res, delta, guess, lcoords_guess;
    FloatMatrix jac;
    double convergence_limit, error;

    // find a suitable convergence limit
    convergence_limit = 1e-6 * this->giveCharacteristicLength(cellgeo);

    // setup initial guess
    lcoords_guess.resize(gcoords.giveSize());
    lcoords_guess.zero();

    // apply Newton-Raphson to solve the problem
    for (int nite = 0; nite < 10; nite++) {
        // compute the residual
        this->local2global(guess, lcoords_guess, cellgeo, time);
        res.beDifferenceOf(gcoords, guess);

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
    if ( error > convergence_limit) { // Imperfect, could give false negatives.
        //OOFEM_WARNING("FEI2dTrQuad :: global2local - Failed convergence");
        answer.resize(0);
        return false;
    }

    answer.resize(4);
    answer(0) = lcoords_guess(0);
    answer(1) = lcoords_guess(1);
    answer(2) = lcoords_guess(2);
    answer(3) = 1.0 - lcoords_guess(0) - lcoords_guess(1) - lcoords_guess(2);

    for (int  i = 0; i < 4; i++ ) {
        if ( answer(i) < ( 0. - POINT_TOL ) ) {
            return false;
        } else if ( answer(i) > ( 1. + POINT_TOL ) ) {
            return false;
        }
    }
    return true;
}


double
FEI3dTetQuad :: giveCharacteristicLength(const FEICellGeometry &cellgeo) const
{
    return cellgeo.giveVertexCoordinates(1)->distance(cellgeo.giveVertexCoordinates(2));
}


double
FEI3dTetQuad :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix;
    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant();
}


void
FEI3dTetQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix dNdxi;
    FloatMatrix coords;
    this->evaldNdxi(dNdxi, lcoords);
    jacobianMatrix.resize(3,3);
    coords.resize(10,3);
    for (int i = 1; i <= 10; ++i) {
        const FloatArray *c = cellgeo.giveVertexCoordinates(i);
        coords.at(i,1) = c->at(1);
        coords.at(i,1) = c->at(2);
        coords.at(i,1) = c->at(3);
    }
    jacobianMatrix.beProductOf(dNdxi, coords);
}


void
FEI3dTetQuad :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double xi = lcoords.at(1);
    answer.resize(3);
    answer(0) = 0.5*(xi-1.0)*xi;
    answer(1) = 0.5*(xi+1.0)*xi;
    answer(2) = 1.0-xi*xi;
}

void
FEI3dTetQuad :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                           const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    ///@todo Implement this
    OOFEM_ERROR("FEI3dTetQuad :: edgeEvaldNdx - Not supported");
}

void
FEI3dTetQuad :: edgeLocal2global(FloatArray &answer, int iedge,
                               const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray N;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(N, lcoords, cellgeo, time);

    answer.resize(0);
    for (int i = 0; i < N.giveSize(); ++i) {
        answer.add( N(i), *cellgeo.giveVertexCoordinates( edgeNodes(i)) );
    }
}


double
FEI3dTetQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    ///@todo Implement this
    OOFEM_ERROR("FEI3dTetQuad :: edgeGiveTransformationJacobian - Not supported");
    return -1;
}


void
FEI3dTetQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    edgeNodes.resize(3);

    if ( iedge == 1 ) { // edge between nodes 1 2
        edgeNodes(0) = 1;
        edgeNodes(1) = 2;
        edgeNodes(2) = 5;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        edgeNodes(0) = 2;
        edgeNodes(1) = 3;
        edgeNodes(2) = 6;
    } else if ( iedge == 3 ) { // edge between nodes 3 1
        edgeNodes(0) = 3;
        edgeNodes(1) = 1;
        edgeNodes(2) = 7;
    } else if ( iedge == 4 ) { // edge between nodes 1 4
        edgeNodes(0) = 1;
        edgeNodes(1) = 4;
        edgeNodes(2) = 8;
    } else if ( iedge == 5 ) { // edge between nodes 2 4
        edgeNodes(0) = 2;
        edgeNodes(1) = 4;
        edgeNodes(2) = 9;
    } else if ( iedge == 6 ) { // edge between nodes 3 4
        edgeNodes(0) = 3;
        edgeNodes(1) = 4;
        edgeNodes(2) = 10;
    } else {
        OOFEM_ERROR2("FEI3dTetQuad :: computeEdgeMapping: wrong egde number (%d)", iedge);
    }
}

double
FEI3dTetQuad :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    ///@todo Implement this
    OOFEM_ERROR("FEI3dTetQuad :: edgeComputeLength - Not supported");
    return -1;
}

void
FEI3dTetQuad :: surfaceEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1. - l1 - l2;

    answer.resize(6);

    answer.at(1) = ( 2. * l1 - 1. ) * l1;
    answer.at(2) = ( 2. * l2 - 1. ) * l2;
    answer.at(3) = ( 2. * l3 - 1. ) * l3;
    answer.at(4) = 4. * l1 * l2;
    answer.at(5) = 4. * l2 * l3;
    answer.at(6) = 4. * l3 * l1;
}

void
FEI3dTetQuad :: surfaceLocal2global(FloatArray &answer, int isurf,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray nodes;
    FloatArray N;
    this->computeLocalSurfaceMapping(nodes, isurf);
    this->surfaceEvalN(N, lcoords, cellgeo, 0.0);

    answer.resize(0);
    for (int i = 0; i < N.giveSize(); ++i) {
        answer.add( N(i), *cellgeo.giveVertexCoordinates(nodes(i)) );
    }
}

void
FEI3dTetQuad :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    // Translate the local surface coordinate to the volume coordinates and compute the gradient there.
    double a, b, c;
    a = lcoords.at(1);
    b = lcoords.at(2);
    c = 1.0 - a - b;
    FloatArray lcoords_tet(4);
    lcoords_tet.at(isurf) = 0.0;
    if (isurf == 1) {
        lcoords_tet.at(4) = a;
        lcoords_tet.at(3) = b;
        lcoords_tet.at(2) = c;
    } else if (isurf == 2) {
        lcoords_tet.at(1) = a;
        lcoords_tet.at(3) = b;
        lcoords_tet.at(4) = c;
    } else if (isurf == 3) {
        lcoords_tet.at(1) = a;
        lcoords_tet.at(4) = b;
        lcoords_tet.at(2) = c;
    } else if (isurf == 4) {
        lcoords_tet.at(2) = a;
        lcoords_tet.at(3) = b;
        lcoords_tet.at(1) = c;
    }
    this->evaldNdx(answer, lcoords_tet, cellgeo, 0.0);
}

double
FEI3dTetQuad :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray snodes(3);
    FloatArray a,b;
    this->computeLocalSurfaceMapping(snodes, isurf);

    double l1, l2, l3;
    l1 = lcoords.at(1);
    l2 = lcoords.at(2);
    l3 = 1.0 - l1 - l2;

    FloatArray dNdxi(6), dNdeta(6);

    dNdxi(0) = 4.0 * l1 - 1.0;
    dNdxi(1) = 0.0;
    dNdxi(2) = -1.0 * ( 4.0 * l3 - 1.0 );
    dNdxi(3) = 4.0 * l2;
    dNdxi(4) = -4.0 * l2;
    dNdxi(5) = 4.0 * l3 - 4.0 * l1;

    dNdeta(0) = 0.0;
    dNdeta(1) = 4.0 * l2 - 1.0;
    dNdeta(2) = -1.0 * ( 4.0 * l3 - 1.0 );
    dNdeta(3) = 4.0 * l1;
    dNdeta(4) = 4.0 * l3 - 4.0 * l2;
    dNdeta(5) = -4.0 * l1;

    for (int i = 0; i < 6; ++i) {
        a.add(dNdxi(i),  *cellgeo.giveVertexCoordinates(snodes(i)));
        b.add(dNdeta(i), *cellgeo.giveVertexCoordinates(snodes(i)));
    }
    answer.beVectorProductOf(a, b);
    double J = answer.computeNorm();
    answer.times(1/J);
    return J;
}

double
FEI3dTetQuad :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                const FEICellGeometry &cellgeo)
{
    FloatArray normal;
    return this->surfaceEvalNormal(normal, isurf, lcoords, cellgeo);
}

void
FEI3dTetQuad :: computeLocalSurfaceMapping(IntArray &surfNodes, int isurf)
{
    int aNode = 0, bNode = 0, cNode = 0, dNode = 0, eNode = 0, fNode = 0;
    surfNodes.resize(3);

    if ( isurf == 1 ) {
        aNode = 1;
        bNode = 3;
        cNode = 2;
        dNode = 7;
        eNode = 6;
        fNode = 5;
    } else if ( isurf == 2 ) {
        aNode = 1;
        bNode = 2;
        cNode = 4;
        dNode = 5;
        eNode = 9;
        fNode = 8;
    } else if ( isurf == 3 ) {
        aNode = 2;
        bNode = 3;
        cNode = 4;
        dNode = 6;
        eNode = 10;
        fNode = 9;
    } else if ( isurf == 4 ) {
        aNode = 1;
        bNode = 4;
        cNode = 3;
        dNode = 8;
        eNode = 10;
        fNode = 7;
    } else {
        OOFEM_ERROR2("FEI3dTetQuad :: computeLocalSurfaceMapping: wrong surface number (%d)", isurf);
    }

    surfNodes.at(1) = aNode;
    surfNodes.at(2) = bNode;
    surfNodes.at(3) = cNode;
    surfNodes.at(4) = dNode;
    surfNodes.at(5) = eNode;
    surfNodes.at(6) = fNode;

}
} // end namespace oofem
