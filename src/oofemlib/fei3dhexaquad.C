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

#include "fei3dhexaquad.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "gaussintegrationrule.h"

namespace oofem {
void
FEI3dHexaQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double u, v, w;
    answer.resize(20);

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    answer.at(1)  = 0.125 * ( 1.0 - u ) * ( 1.0 - v ) * ( 1.0 + w ) * ( -u - v + w - 2.0 );
    answer.at(2)  = 0.125 * ( 1.0 - u ) * ( 1.0 + v ) * ( 1.0 + w ) * ( -u + v + w - 2.0 );
    answer.at(3)  = 0.125 * ( 1.0 + u ) * ( 1.0 + v ) * ( 1.0 + w ) * ( u + v + w - 2.0 );
    answer.at(4)  = 0.125 * ( 1.0 + u ) * ( 1.0 - v ) * ( 1.0 + w ) * ( u - v + w - 2.0 );
    answer.at(5)  = 0.125 * ( 1.0 - u ) * ( 1.0 - v ) * ( 1.0 - w ) * ( -u - v - w - 2.0 );
    answer.at(6)  = 0.125 * ( 1.0 - u ) * ( 1.0 + v ) * ( 1.0 - w ) * ( -u + v - w - 2.0 );
    answer.at(7)  = 0.125 * ( 1.0 + u ) * ( 1.0 + v ) * ( 1.0 - w ) * ( u + v - w - 2.0 );
    answer.at(8)  = 0.125 * ( 1.0 + u ) * ( 1.0 - v ) * ( 1.0 - w ) * ( u - v - w - 2.0 );

    answer.at(9)  = 0.25 * ( 1.0 - v * v ) * ( 1.0 - u ) * ( 1.0 + w );
    answer.at(10)  = 0.25 * ( 1.0 - u * u ) * ( 1.0 + v ) * ( 1.0 + w );
    answer.at(11)  = 0.25 * ( 1.0 - v * v ) * ( 1.0 + u ) * ( 1.0 + w );
    answer.at(12)  = 0.25 * ( 1.0 - u * u ) * ( 1.0 - v ) * ( 1.0 + w );

    answer.at(13)  = 0.25 * ( 1.0 - v * v ) * ( 1.0 - u ) * ( 1.0 - w );
    answer.at(14)  = 0.25 * ( 1.0 - u * u ) * ( 1.0 + v ) * ( 1.0 - w );
    answer.at(15)  = 0.25 * ( 1.0 - v * v ) * ( 1.0 + u ) * ( 1.0 - w );
    answer.at(16)  = 0.25 * ( 1.0 - u * u ) * ( 1.0 - v ) * ( 1.0 - w );

    answer.at(17)  = 0.25 * ( 1.0 - u ) * ( 1.0 - v ) * ( 1.0 - w * w );
    answer.at(18)  = 0.25 * ( 1.0 - u ) * ( 1.0 + v ) * ( 1.0 - w * w );
    answer.at(19)  = 0.25 * ( 1.0 + u ) * ( 1.0 + v ) * ( 1.0 - w * w );
    answer.at(20)  = 0.25 * ( 1.0 + u ) * ( 1.0 - v ) * ( 1.0 - w * w );
}

double
FEI3dHexaQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix, inv, dNduvw, coords;
    this->giveLocalDerivative(dNduvw, lcoords);
    coords.resize( 3, dNduvw.giveNumberOfRows() );
    for ( int i = 1; i <= dNduvw.giveNumberOfRows(); i++ ) {
        coords.setColumn(* cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
    inv.beInverseOf(jacobianMatrix);

    answer.beProductOf(dNduvw, inv);
    return jacobianMatrix.giveDeterminant();
}

void
FEI3dHexaQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;

    this->evalN(n, lcoords, cellgeo);
    answer.clear();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates(i) );
    }
}

double FEI3dHexaQuad :: giveCharacteristicLength(const FEICellGeometry &cellgeo) const
{
    const FloatArray *n1 = cellgeo.giveVertexCoordinates(1);
    const FloatArray *n2 = cellgeo.giveVertexCoordinates(7);
    ///@todo Change this so that it is not dependent on node order.
    return n1->distance(n2);
}


#define POINT_TOL 1.e-3

int
FEI3dHexaQuad :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
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
    for ( int nite = 0; nite < 40; nite++ ) {
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
        jac.solveForRhs(res, delta);

        // update guess
        answer.add(delta);
    }
    if ( error > convergence_limit ) { // Imperfect, could give false negatives.
        //OOFEM_ERROR("no convergence after 10 iterations");
        answer.zero();
        return false;
    }

    // check limits for each local coordinate [-1,1] for quadrilaterals. (different for other elements, typically [0,1]).
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

void FEI3dHexaQuad :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double u = lcoords.at(1);
    answer.resize(3);
    answer.at(1) = 0.5 * ( u - 1. ) * u;
    answer.at(2) = 0.5 * ( u + 1. ) * u;
    answer.at(3) = 1. - u * u;
}

void FEI3dHexaQuad :: edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray eNodes;
    double u = lcoords.at(1);
    this->computeLocalEdgeMapping(eNodes, iedge);
    answer.clear();
    answer.add( 0.5 * ( u - 1. ) * u, * cellgeo.giveVertexCoordinates( eNodes.at(1) ) );
    answer.add( 0.5 * ( u - 1. ) * u, * cellgeo.giveVertexCoordinates( eNodes.at(2) ) );
    answer.add( 1. - u * u,       * cellgeo.giveVertexCoordinates( eNodes.at(3) ) );
}

void FEI3dHexaQuad :: edgeEvaldNdx(FloatMatrix &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray eNodes;
    FloatArray dNdu;
    double u = lcoords.at(1);
    this->computeLocalEdgeMapping(eNodes, iedge);
    dNdu.add( u - 0.5, * cellgeo.giveVertexCoordinates( eNodes.at(1) ) );
    dNdu.add( u + 0.5, * cellgeo.giveVertexCoordinates( eNodes.at(2) ) );
    dNdu.add( -2. * u, * cellgeo.giveVertexCoordinates( eNodes.at(3) ) );
    // Why matrix output?
    answer.resize(3, 1);
    answer.setColumn(dNdu, 1);
}

double FEI3dHexaQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray eNodes;
    FloatArray dNdu;
    double u = lcoords.at(1);
    this->computeLocalEdgeMapping(eNodes, iedge);
    dNdu.add( u - 0.5, * cellgeo.giveVertexCoordinates( eNodes.at(1) ) );
    dNdu.add( u + 0.5, * cellgeo.giveVertexCoordinates( eNodes.at(2) ) );
    dNdu.add( -2. * u, * cellgeo.giveVertexCoordinates( eNodes.at(3) ) );
    return dNdu.computeNorm();
}

void
FEI3dHexaQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    if ( iedge == 1 ) {
        edgeNodes = { 1, 2,  9};
    } else if ( iedge == 2 ) {
        edgeNodes = { 2, 3, 10};
    } else if ( iedge == 3 ) {
        edgeNodes = { 3, 4, 11};
    } else if ( iedge == 4 ) {
        edgeNodes = { 4, 1, 12};
    } else if ( iedge == 5 ) {
        edgeNodes = { 5, 6, 13};
    } else if ( iedge == 6 ) {
        edgeNodes = { 6, 7, 14};
    } else if ( iedge == 7 ) {
        edgeNodes = { 7, 8, 15};
    } else if ( iedge == 8 ) {
        edgeNodes = { 8, 5, 16};
    } else if ( iedge == 9 ) {
        edgeNodes = { 1, 5, 17};
    } else if ( iedge == 10 ) {
        edgeNodes = { 2, 6, 18};
    } else if ( iedge == 11 ) {
        edgeNodes = { 3, 7, 19};
    } else if ( iedge == 12 ) {
        edgeNodes = { 4, 8, 20};
    } else {
        OOFEM_ERROR("wrong edge number (%d)", iedge);
    }
}

void
FEI3dHexaQuad :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
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


double
FEI3dHexaQuad :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray a, b, dNdksi(8), dNdeta(8);
    double ksi, eta;
    IntArray snodes;

    this->computeLocalSurfaceMapping(snodes, isurf);

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    dNdksi.at(1) =  0.25 * ( 1. + eta ) * ( 2.0 * ksi + eta );
    dNdksi.at(2) = -0.25 * ( 1. + eta ) * ( -2.0 * ksi + eta );
    dNdksi.at(3) = -0.25 * ( 1. - eta ) * ( -2.0 * ksi - eta );
    dNdksi.at(4) =  0.25 * ( 1. - eta ) * ( 2.0 * ksi - eta );
    dNdksi.at(5) = -ksi * ( 1. + eta );
    dNdksi.at(6) = -0.5 * ( 1. - eta * eta );
    dNdksi.at(7) = -ksi * ( 1. - eta );
    dNdksi.at(8) =  0.5 * ( 1. - eta * eta );

    dNdeta.at(1) =  0.25 * ( 1. + ksi ) * ( 2.0 * eta + ksi );
    dNdeta.at(2) =  0.25 * ( 1. - ksi ) * ( 2.0 * eta - ksi );
    dNdeta.at(3) = -0.25 * ( 1. - ksi ) * ( -2.0 * eta - ksi );
    dNdeta.at(4) = -0.25 * ( 1. + ksi ) * ( -2.0 * eta + ksi );
    dNdeta.at(5) =  0.5 * ( 1. - ksi * ksi );
    dNdeta.at(6) = -eta * ( 1. - ksi );
    dNdeta.at(7) = -0.5 * ( 1. - ksi * ksi );
    dNdeta.at(8) = -eta * ( 1. + ksi );

    for ( int i = 1; i <= 8; ++i ) {
        a.add( dNdksi.at(i), * cellgeo.giveVertexCoordinates( snodes.at(i) ) );
        b.add( dNdeta.at(i), * cellgeo.giveVertexCoordinates( snodes.at(i) ) );
    }

    answer.beVectorProductOf(a, b);
    return answer.normalize();
}

void
FEI3dHexaQuad :: surfaceLocal2global(FloatArray &answer, int isurf,
                                     const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray nodes;
    FloatArray n;

    this->computeLocalSurfaceMapping(nodes, isurf);

    this->surfaceEvalN(n, isurf, lcoords, cellgeo);

    answer.clear();
    for ( int i = 1; i <= n.giveSize(); ++i ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates( nodes.at(i) ) );
    }
}


double
FEI3dHexaQuad :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                   const FEICellGeometry &cellgeo)
{
    FloatArray normal;
    return this->surfaceEvalNormal(normal, isurf, lcoords, cellgeo);
}


void
FEI3dHexaQuad :: computeLocalSurfaceMapping(IntArray &nodes, int isurf)
{
    nodes.resize(8);

    // the actual numbering  has a positive normal pointing outwards from the element  - (LSpace compatible)

    if ( isurf == 1 ) {
        nodes = { 3, 2, 1, 4, 10,  9, 12, 11};
    } else if ( isurf == 2 ) {
        nodes = { 7, 8, 5, 6, 15, 16, 13, 14};
    } else if ( isurf == 3 ) {
        nodes = { 2, 6, 5, 1, 18, 13, 17,  9};
    } else if ( isurf == 4 ) {
        nodes = { 3, 7, 6, 2, 19, 14, 18, 10};
    } else if ( isurf == 5 ) {
        nodes = { 3, 4, 8, 7, 11, 20, 15, 19};
    } else if ( isurf == 6 ) {
        nodes = { 4, 1, 5, 8, 12, 17, 16, 20};
    } else {
        OOFEM_ERROR("wrong surface number (%d)", isurf);
    }

#if 0
    // this commented numbering is symmetrical with respect to local coordinate axes

    if      ( isurf == 1 ) { // surface 1 - nodes   3 2 1 4  10  9 12 11
        nodes.at(1) =  3;
        nodes.at(2) =  2;
        nodes.at(3) =  1;
        nodes.at(4) =  4;
        nodes.at(5) = 10;
        nodes.at(6) =  9;
        nodes.at(7) = 12;
        nodes.at(8) = 11;
    } else if ( iSurf == 2 ) { // surface 2 - nodes   7 6 5 8  14 13 16 15
        nodes.at(1) =  7;
        nodes.at(2) =  6;
        nodes.at(3) =  5;
        nodes.at(4) =  8;
        nodes.at(5) = 14;
        nodes.at(6) = 13;
        nodes.at(7) = 16;
        nodes.at(8) = 15;
    } else if ( iSurf == 3 ) { // surface 3 - nodes   2 1 5 6   9 17 13 18
        nodes.at(1) =  2;
        nodes.at(2) =  1;
        nodes.at(3) =  5;
        nodes.at(4) =  6;
        nodes.at(5) =  9;
        nodes.at(6) = 17;
        nodes.at(7) = 13;
        nodes.at(8) = 18;
    } else if ( isurf == 4 ) { // surface 4 - nodes   3 7 6 2  19 14 18 10
        nodes.at(1) =  3;
        nodes.at(2) =  7;
        nodes.at(3) =  6;
        nodes.at(4) =  2;
        nodes.at(5) = 19;
        nodes.at(6) = 14;
        nodes.at(7) = 18;
        nodes.at(8) = 10;
    } else if ( isurf == 5 ) { // surface 5 - nodes   3 4 8 7  11 20 15 19
        nodes.at(1) =  3;
        nodes.at(2) =  4;
        nodes.at(3) =  8;
        nodes.at(4) =  7;
        nodes.at(5) = 11;
        nodes.at(6) = 20;
        nodes.at(7) = 15;
        nodes.at(8) = 19;
    } else if ( iSurf == 6 ) { // surface 6 - nodes   4 8 5 1  20 16 17 12
        nodes.at(1) =  4;
        nodes.at(2) =  8;
        nodes.at(3) =  5;
        nodes.at(4) =  1;
        nodes.at(5) = 20;
        nodes.at(6) = 16;
        nodes.at(7) = 17;
        nodes.at(8) = 12;
    } else {
        OOFEM_ERROR("wrong surface number (%d)", isurf);
    }
#endif
}


void
FEI3dHexaQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
{
    FloatMatrix dNduvw, coords;
    this->giveLocalDerivative(dNduvw, lcoords);
    coords.resize( 3, dNduvw.giveNumberOfRows() );
    for ( int i = 1; i <= dNduvw.giveNumberOfRows(); i++ ) {
        coords.setColumn(* cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
}


void
FEI3dHexaQuad :: giveLocalDerivative(FloatMatrix &dN, const FloatArray &lcoords)
{
    double u, v, w;
    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    dN.resize(20, 3);

    dN.at(1, 1) = 0.125 * ( 1.0 - v ) * ( 1.0 + w ) * ( 2.0 * u + v - w + 1.0 );
    dN.at(2, 1) = 0.125 * ( 1.0 + v ) * ( 1.0 + w ) * ( 2.0 * u - v - w + 1.0 );
    dN.at(3, 1) = 0.125 * ( 1.0 + v ) * ( 1.0 + w ) * ( 2.0 * u + v + w - 1.0 );
    dN.at(4, 1) = 0.125 * ( 1.0 - v ) * ( 1.0 + w ) * ( 2.0 * u - v + w - 1.0 );
    dN.at(5, 1) = 0.125 * ( 1.0 - v ) * ( 1.0 - w ) * ( 2.0 * u + v + w + 1.0 );
    dN.at(6, 1) = 0.125 * ( 1.0 + v ) * ( 1.0 - w ) * ( 2.0 * u - v + w + 1.0 );
    dN.at(7, 1) = 0.125 * ( 1.0 + v ) * ( 1.0 - w ) * ( 2.0 * u + v - w - 1.0 );
    dN.at(8, 1) = 0.125 * ( 1.0 - v ) * ( 1.0 - w ) * ( 2.0 * u - v - w - 1.0 );
    dN.at(9, 1) = -0.25 * ( 1.0 - v * v ) * ( 1.0 + w );
    dN.at(10, 1) =  -0.5 * u * ( 1.0 + v ) * ( 1.0 + w );
    dN.at(11, 1) =  0.25 * ( 1.0 - v * v ) * ( 1.0 + w );
    dN.at(12, 1) =  -0.5 * u * ( 1.0 - v ) * ( 1.0 + w );
    dN.at(13, 1) = -0.25 * ( 1.0 - v * v ) * ( 1.0 - w );
    dN.at(14, 1) =  -0.5 * u * ( 1.0 + v ) * ( 1.0 - w );
    dN.at(15, 1) =  0.25 * ( 1.0 - v * v ) * ( 1.0 - w );
    dN.at(16, 1) =  -0.5 * u * ( 1.0 - v ) * ( 1.0 - w );
    dN.at(17, 1) = -0.25 * ( 1.0 - v ) * ( 1.0 - w * w );
    dN.at(18, 1) = -0.25 * ( 1.0 + v ) * ( 1.0 - w * w );
    dN.at(19, 1) =  0.25 * ( 1.0 + v ) * ( 1.0 - w * w );
    dN.at(20, 1) =  0.25 * ( 1.0 - v ) * ( 1.0 - w * w );

    dN.at(1, 2) = 0.125 * ( 1.0 - u ) * ( 1.0 + w ) * ( 2.0 * v + u - w + 1.0 );
    dN.at(2, 2) = 0.125 * ( 1.0 - u ) * ( 1.0 + w ) * ( 2.0 * v - u + w - 1.0 );
    dN.at(3, 2) = 0.125 * ( 1.0 + u ) * ( 1.0 + w ) * ( 2.0 * v + u + w - 1.0 );
    dN.at(4, 2) = 0.125 * ( 1.0 + u ) * ( 1.0 + w ) * ( 2.0 * v - u - w + 1.0 );
    dN.at(5, 2) = 0.125 * ( 1.0 - u ) * ( 1.0 - w ) * ( 2.0 * v + u + w + 1.0 );
    dN.at(6, 2) = 0.125 * ( 1.0 - u ) * ( 1.0 - w ) * ( 2.0 * v - u - w - 1.0 );
    dN.at(7, 2) = 0.125 * ( 1.0 + u ) * ( 1.0 - w ) * ( 2.0 * v + u - w - 1.0 );
    dN.at(8, 2) = 0.125 * ( 1.0 + u ) * ( 1.0 - w ) * ( 2.0 * v - u + w + 1.0 );
    dN.at(9, 2) =  -0.5 * v * ( 1.0 - u ) * ( 1.0 + w );
    dN.at(10, 2) =  0.25 * ( 1.0 - u * u ) * ( 1.0 + w );
    dN.at(11, 2) =  -0.5 * v * ( 1.0 + u ) * ( 1.0 + w );
    dN.at(12, 2) = -0.25 * ( 1.0 - u * u ) * ( 1.0 + w );
    dN.at(13, 2) =  -0.5 * v * ( 1.0 - u ) * ( 1.0 - w );
    dN.at(14, 2) =  0.25 * ( 1.0 - u * u ) * ( 1.0 - w );
    dN.at(15, 2) =  -0.5 * v * ( 1.0 + u ) * ( 1.0 - w );
    dN.at(16, 2) = -0.25 * ( 1.0 - u * u ) * ( 1.0 - w );
    dN.at(17, 2) = -0.25 * ( 1.0 - u ) * ( 1.0 - w * w );
    dN.at(18, 2) =  0.25 * ( 1.0 - u ) * ( 1.0 - w * w );
    dN.at(19, 2) =  0.25 * ( 1.0 + u ) * ( 1.0 - w * w );
    dN.at(20, 2) = -0.25 * ( 1.0 + u ) * ( 1.0 - w * w );

    dN.at(1, 3) = 0.125 * ( 1.0 - u ) * ( 1.0 - v ) * ( 2.0 * w - u - v - 1.0 );
    dN.at(2, 3) = 0.125 * ( 1.0 - u ) * ( 1.0 + v ) * ( 2.0 * w - u + v - 1.0 );
    dN.at(3, 3) = 0.125 * ( 1.0 + u ) * ( 1.0 + v ) * ( 2.0 * w + u + v - 1.0 );
    dN.at(4, 3) = 0.125 * ( 1.0 + u ) * ( 1.0 - v ) * ( 2.0 * w + u - v - 1.0 );
    dN.at(5, 3) = 0.125 * ( 1.0 - u ) * ( 1.0 - v ) * ( 2.0 * w + u + v + 1.0 );
    dN.at(6, 3) = 0.125 * ( 1.0 - u ) * ( 1.0 + v ) * ( 2.0 * w + u - v + 1.0 );
    dN.at(7, 3) = 0.125 * ( 1.0 + u ) * ( 1.0 + v ) * ( 2.0 * w - u - v + 1.0 );
    dN.at(8, 3) = 0.125 * ( 1.0 + u ) * ( 1.0 - v ) * ( 2.0 * w - u + v + 1.0 );
    dN.at(9, 3) =  0.25 * ( 1.0 - v * v ) * ( 1.0 - u );
    dN.at(10, 3) =  0.25 * ( 1.0 - u * u ) * ( 1.0 + v );
    dN.at(11, 3) =  0.25 * ( 1.0 - v * v ) * ( 1.0 + u );
    dN.at(12, 3) =  0.25 * ( 1.0 - u * u ) * ( 1.0 - v );
    dN.at(13, 3) = -0.25 * ( 1.0 - v * v ) * ( 1.0 - u );
    dN.at(14, 3) = -0.25 * ( 1.0 - u * u ) * ( 1.0 + v );
    dN.at(15, 3) = -0.25 * ( 1.0 - v * v ) * ( 1.0 + u );
    dN.at(16, 3) = -0.25 * ( 1.0 - u * u ) * ( 1.0 - v );
    dN.at(17, 3) =  -0.5 * ( 1.0 - u ) * ( 1.0 - v ) * w;
    dN.at(18, 3) =  -0.5 * ( 1.0 - u ) * ( 1.0 + v ) * w;
    dN.at(19, 3) =  -0.5 * ( 1.0 + u ) * ( 1.0 + v ) * w;
    dN.at(20, 3) =  -0.5 * ( 1.0 + u ) * ( 1.0 - v ) * w;
}


double
FEI3dHexaQuad :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo)
{
    IntArray fNodes;
    this->computeLocalSurfaceMapping(fNodes, iEdge);

    const FloatArray &c1 = * cellgeo.giveVertexCoordinates( fNodes.at(1) );
    const FloatArray &c2 = * cellgeo.giveVertexCoordinates( fNodes.at(2) );
    const FloatArray &c3 = * cellgeo.giveVertexCoordinates( fNodes.at(3) );
    const FloatArray &c4 = * cellgeo.giveVertexCoordinates( fNodes.at(4) );
    const FloatArray &c5 = * cellgeo.giveVertexCoordinates( fNodes.at(5) );
    const FloatArray &c6 = * cellgeo.giveVertexCoordinates( fNodes.at(6) );
    const FloatArray &c7 = * cellgeo.giveVertexCoordinates( fNodes.at(7) );
    const FloatArray &c8 = * cellgeo.giveVertexCoordinates( fNodes.at(8) );

    // Generated with Mathematica (rather unwieldy expression, tried to simplify it as good as possible, but it could probably be better)
    return (
               c1(2) * ( c2(1) * ( -3 * c3(0) - 3 * c4(0) - 12 * c5(0) + 14 * c6(0) + 14 * c8(0) ) +
                        c3(1) * ( 3 * c2(0) - 3 * c4(0) - 6 * c5(0) - 6 * c6(0) + 6 * c7(0) + 6 * c8(0) ) +
                        c4(1) * ( 3 * c2(0) + 3 * c3(0) - 14 * c5(0) - 14 * c7(0) + 12 * c8(0) ) +
                        c5(1) * ( 12 * c2(0) + 6 * c3(0) + 14 * c4(0) - 4 * c6(0) - 8 * c7(0) - 60 * c8(0) ) +
                        c6(1) * ( -14 * c2(0) + 6 * c3(0) + 4 * c5(0) + 12 * c7(0) - 8 * c8(0) ) +
                        c7(1) * ( -6 * c3(0) + 14 * c4(0) + 8 * c5(0) - 12 * c6(0) - 4 * c8(0) ) +
                        c8(1) * ( -14 * c2(0) - 6 * c3(0) - 12 * c4(0) + 60 * c5(0) + 8 * c6(0) + 4 * c7(0) ) ) +
               c2(2) * ( c1(1) * ( 3 * c3(0) + 3 * c4(0) + 12 * c5(0) - 14 * c6(0) - 14 * c8(0) ) +
                        c3(1) * ( -3 * c1(0) - 3 * c4(0) + 14 * c5(0) - 12 * c6(0) + 14 * c7(0) ) +
                        c4(1) * ( -3 * c1(0) + 3 * c3(0) + 6 * c5(0) - 6 * c6(0) - 6 * c7(0) + 6 * c8(0) ) +
                        c5(1) * ( -12 * c1(0) - 14 * c3(0) - 6 * c4(0) + 60 * c6(0) + 8 * c7(0) + 4 * c8(0) ) +
                        c6(1) * ( 14 * c1(0) + 12 * c3(0) + 6 * c4(0) - 60 * c5(0) - 4 * c7(0) - 8 * c8(0) ) +
                        c7(1) * ( -14 * c3(0) + 6 * c4(0) - 8 * c5(0) + 4 * c6(0) + 12 * c8(0) ) +
                        c8(1) * ( 14 * c1(0) - 6 * c4(0) - 4 * c5(0) + 8 * c6(0) - 12 * c7(0) ) ) +
               c3(2) * ( c1(1) * ( -3 * c2(0) + 3 * c4(0) + 6 * c5(0) + 6 * c6(0) - 6 * c7(0) - 6 * c8(0) ) +
                        c2(1) * ( 3 * c1(0) + 3 * c4(0) - 14 * c5(0) + 12 * c6(0) - 14 * c7(0) ) +
                        c4(1) * ( -3 * c1(0) - 3 * c2(0) + 14 * c6(0) - 12 * c7(0) + 14 * c8(0) ) +
                        c5(1) * ( -6 * c1(0) + 14 * c2(0) - 4 * c6(0) + 8 * c7(0) - 12 * c8(0) ) +
                        c6(1) * ( -6 * c1(0) - 12 * c2(0) - 14 * c4(0) + 4 * c5(0) + 60 * c7(0) + 8 * c8(0) ) +
                        c7(1) * ( 6 * c1(0) + 14 * c2(0) + 12 * c4(0) - 8 * c5(0) - 60 * c6(0) - 4 * c8(0) ) +
                        c8(1) * ( 6 * c1(0) - 14 * c4(0) + 12 * c5(0) - 8 * c6(0) + 4 * c7(0) ) ) +
               c4(2) * ( c1(1) * ( -3 * c2(0) - 3 * c3(0) + 14 * c5(0) + 14 * c7(0) - 12 * c8(0) ) +
                        c2(1) * ( 3 * c1(0) - 3 * c3(0) - 6 * c5(0) + 6 * c6(0) + 6 * c7(0) - 6 * c8(0) ) +
                        c3(1) * ( 3 * c1(0) + 3 * c2(0) - 14 * c6(0) + 12 * c7(0) - 14 * c8(0) ) +
                        c5(1) * ( -14 * c1(0) + 6 * c2(0) + 12 * c6(0) - 8 * c7(0) + 4 * c8(0) ) +
                        c6(1) * ( -6 * c2(0) + 14 * c3(0) - 12 * c5(0) - 4 * c7(0) + 8 * c8(0) ) +
                        c7(1) * ( -14 * c1(0) - 6 * c2(0) - 12 * c3(0) + 8 * c5(0) + 4 * c6(0) + 60 * c8(0) ) +
                        c8(1) * ( 12 * c1(0) + 6 * c2(0) + 14 * c3(0) - 4 * c5(0) - 8 * c6(0) - 60 * c7(0) ) ) +
               c5(2) * ( c1(1) * ( -12 * c2(0) - 6 * c3(0) - 14 * c4(0) + 4 * c6(0) + 8 * c7(0) + 60 * c8(0) ) +
                        c2(1) * ( 12 * c1(0) + 14 * c3(0) + 6 * c4(0) - 60 * c6(0) - 8 * c7(0) - 4 * c8(0) ) +
                        c3(1) * ( 6 * c1(0) - 14 * c2(0) + 4 * c6(0) - 8 * c7(0) + 12 * c8(0) ) +
                        c4(1) * ( 14 * c1(0) - 6 * c2(0) - 12 * c6(0) + 8 * c7(0) - 4 * c8(0) ) +
                        c6(1) * ( -4 * c1(0) + 60 * c2(0) - 4 * c3(0) + 12 * c4(0) - 32 * c7(0) - 32 * c8(0) ) +
                        c7(1) * ( -8 * c1(0) + 8 * c2(0) + 8 * c3(0) - 8 * c4(0) + 32 * c6(0) - 32 * c8(0) ) +
                        c8(1) * ( -60 * c1(0) + 4 * c2(0) - 12 * c3(0) + 4 * c4(0) + 32 * c6(0) + 32 * c7(0) ) ) +
               c6(2) * ( c1(1) * ( 14 * c2(0) - 6 * c3(0) - 4 * c5(0) - 12 * c7(0) + 8 * c8(0) ) +
                        c2(1) * ( -14 * c1(0) - 12 * c3(0) - 6 * c4(0) + 60 * c5(0) + 4 * c7(0) + 8 * c8(0) ) +
                        c3(1) * ( 6 * c1(0) + 12 * c2(0) + 14 * c4(0) - 4 * c5(0) - 60 * c7(0) - 8 * c8(0) ) +
                        c4(1) * ( 6 * c2(0) - 14 * c3(0) + 12 * c5(0) + 4 * c7(0) - 8 * c8(0) ) +
                        c5(1) * ( 4 * c1(0) - 60 * c2(0) + 4 * c3(0) - 12 * c4(0) + 32 * c7(0) + 32 * c8(0) ) +
                        c7(1) * ( 12 * c1(0) - 4 * c2(0) + 60 * c3(0) - 4 * c4(0) - 32 * c5(0) - 32 * c8(0) ) +
                        c8(1) * ( -8 * c1(0) - 8 * c2(0) + 8 * c3(0) + 8 * c4(0) - 32 * c5(0) + 32 * c7(0) ) ) +
               c7(2) * ( c1(1) * ( 6 * c3(0) - 14 * c4(0) - 8 * c5(0) + 12 * c6(0) + 4 * c8(0) ) +
                        c2(1) * ( 14 * c3(0) - 6 * c4(0) + 8 * c5(0) - 4 * c6(0) - 12 * c8(0) ) +
                        c3(1) * ( -6 * c1(0) - 14 * c2(0) - 12 * c4(0) + 8 * c5(0) + 60 * c6(0) + 4 * c8(0) ) +
                        c4(1) * ( 14 * c1(0) + 6 * c2(0) + 12 * c3(0) - 8 * c5(0) - 4 * c6(0) - 60 * c8(0) ) +
                        c5(1) * ( 8 * c1(0) - 8 * c2(0) - 8 * c3(0) + 8 * c4(0) - 32 * c6(0) + 32 * c8(0) ) +
                        c6(1) * ( -12 * c1(0) + 4 * c2(0) - 60 * c3(0) + 4 * c4(0) + 32 * c5(0) + 32 * c8(0) ) +
                        c8(1) * ( -4 * c1(0) + 12 * c2(0) - 4 * c3(0) + 60 * c4(0) - 32 * c5(0) - 32 * c6(0) ) ) +
               c8(2) * ( c1(1) * ( 14 * c2(0) + 6 * c3(0) + 12 * c4(0) - 60 * c5(0) - 8 * c6(0) - 4 * c7(0) ) +
                        c2(1) * ( -14 * c1(0) + 6 * c4(0) + 4 * c5(0) - 8 * c6(0) + 12 * c7(0) ) +
                        c3(1) * ( -6 * c1(0) + 14 * c4(0) - 12 * c5(0) + 8 * c6(0) - 4 * c7(0) ) +
                        c4(1) * ( -12 * c1(0) - 6 * c2(0) - 14 * c3(0) + 4 * c5(0) + 8 * c6(0) + 60 * c7(0) ) +
                        c5(1) * ( 60 * c1(0) - 4 * c2(0) + 12 * c3(0) - 4 * c4(0) - 32 * c6(0) - 32 * c7(0) ) +
                        c6(1) * ( 8 * c1(0) + 8 * c2(0) - 8 * c3(0) - 8 * c4(0) + 32 * c5(0) - 32 * c7(0) ) +
                        c7(1) * ( 4 * c1(0) - 12 * c2(0) + 4 * c3(0) - 60 * c4(0) + 32 * c5(0) + 32 * c6(0) ) )
               ) / 60.0;
}


IntegrationRule *
FEI3dHexaQuad :: giveIntegrationRule(int order)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Cube, order + 9);
    iRule->SetUpPointsOnCube(points, _Unknown);
    return iRule;
}

IntegrationRule *
FEI3dHexaQuad :: giveBoundaryIntegrationRule(int order, int boundary)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Square, order + 4);
    iRule->SetUpPointsOnSquare(points, _Unknown);
    return iRule;
}
} // end namespace oofem
