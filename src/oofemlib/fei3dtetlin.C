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

#include "fei3dtetlin.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
void
FEI3dTetLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(4);

    answer.at(1) = lcoords.at(1);
    answer.at(2) = lcoords.at(2);
    answer.at(3) = lcoords.at(3);
    answer.at(4) = 1. - lcoords.at(1) - lcoords.at(2) - lcoords.at(3);
}

double
FEI3dTetLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, detJ;
    answer.resize(4, 3);

    x1 = cellgeo.giveVertexCoordinates(1)->at(1);
    x2 = cellgeo.giveVertexCoordinates(2)->at(1);
    x3 = cellgeo.giveVertexCoordinates(3)->at(1);
    x4 = cellgeo.giveVertexCoordinates(4)->at(1);

    y1 = cellgeo.giveVertexCoordinates(1)->at(2);
    y2 = cellgeo.giveVertexCoordinates(2)->at(2);
    y3 = cellgeo.giveVertexCoordinates(3)->at(2);
    y4 = cellgeo.giveVertexCoordinates(4)->at(2);

    z1 = cellgeo.giveVertexCoordinates(1)->at(3);
    z2 = cellgeo.giveVertexCoordinates(2)->at(3);
    z3 = cellgeo.giveVertexCoordinates(3)->at(3);
    z4 = cellgeo.giveVertexCoordinates(4)->at(3);

    detJ = ( ( x4 - x1 ) * ( y2 - y1 ) * ( z3 - z1 ) - ( x4 - x1 ) * ( y3 - y1 ) * ( z2 - z1 ) +
            ( x3 - x1 ) * ( y4 - y1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( y4 - y1 ) * ( z3 - z1 ) +
            ( x2 - x1 ) * ( y3 - y1 ) * ( z4 - z1 ) - ( x3 - x1 ) * ( y2 - y1 ) * ( z4 - z1 ) );

    if ( detJ <= 0.0 ) {
        OOFEM_ERROR("negative volume");
    }

    answer.at(1, 1) = -( ( y3 - y2 ) * ( z4 - z2 ) - ( y4 - y2 ) * ( z3 - z2 ) );
    answer.at(2, 1) = ( y4 - y3 ) * ( z1 - z3 ) - ( y1 - y3 ) * ( z4 - z3 );
    answer.at(3, 1) = -( ( y1 - y4 ) * ( z2 - z4 ) - ( y2 - y4 ) * ( z1 - z4 ) );
    answer.at(4, 1) = ( y2 - y1 ) * ( z3 - z1 ) - ( y3 - y1 ) * ( z2 - z1 );

    answer.at(1, 2) = -( ( x4 - x2 ) * ( z3 - z2 ) - ( x3 - x2 ) * ( z4 - z2 ) );
    answer.at(2, 2) = ( x1 - x3 ) * ( z4 - z3 ) - ( x4 - x3 ) * ( z1 - z3 );
    answer.at(3, 2) = -( ( x2 - x4 ) * ( z1 - z4 ) - ( x1 - x4 ) * ( z2 - z4 ) );
    answer.at(4, 2) = ( x3 - x1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( z3 - z1 );

    answer.at(1, 3) = -( ( x3 - x2 ) * ( y4 - y2 ) - ( x4 - x2 ) * ( y3 - y2 ) );
    answer.at(2, 3) = ( x4 - x3 ) * ( y1 - y3 ) - ( x1 - x3 ) * ( y4 - y3 );
    answer.at(3, 3) = -( ( x1 - x4 ) * ( y2 - y4 ) - ( x2 - x4 ) * ( y1 - y4 ) );
    answer.at(4, 3) = ( x2 - x1 ) * ( y3 - y1 ) - ( x3 - x1 ) * ( y2 - y1 );

    answer.times(1. / detJ);

    return detJ;
}

void
FEI3dTetLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);

    answer.clear();
    for ( int i = 1; i <= 4; i++ ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates(i) );
    }
}

#define POINT_TOL 1.e-3

int
FEI3dTetLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, xp, yp, zp, volume;
    answer.resize(4);

    x1 = cellgeo.giveVertexCoordinates(1)->at(1);
    x2 = cellgeo.giveVertexCoordinates(2)->at(1);
    x3 = cellgeo.giveVertexCoordinates(3)->at(1);
    x4 = cellgeo.giveVertexCoordinates(4)->at(1);

    y1 = cellgeo.giveVertexCoordinates(1)->at(2);
    y2 = cellgeo.giveVertexCoordinates(2)->at(2);
    y3 = cellgeo.giveVertexCoordinates(3)->at(2);
    y4 = cellgeo.giveVertexCoordinates(4)->at(2);

    z1 = cellgeo.giveVertexCoordinates(1)->at(3);
    z2 = cellgeo.giveVertexCoordinates(2)->at(3);
    z3 = cellgeo.giveVertexCoordinates(3)->at(3);
    z4 = cellgeo.giveVertexCoordinates(4)->at(3);

    xp = coords.at(1);
    yp = coords.at(2);
    zp = coords.at(3);

    volume = ( ( x4 - x1 ) * ( y2 - y1 ) * ( z3 - z1 ) - ( x4 - x1 ) * ( y3 - y1 ) * ( z2 - z1 ) +
              ( x3 - x1 ) * ( y4 - y1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( y4 - y1 ) * ( z3 - z1 ) +
              ( x2 - x1 ) * ( y3 - y1 ) * ( z4 - z1 ) - ( x3 - x1 ) * ( y2 - y1 ) * ( z4 - z1 ) ) / 6.;

    answer.resize(4);

    answer.at(1) = ( ( x3 - x2 ) * ( yp - y2 ) * ( z4 - z2 ) - ( xp - x2 ) * ( y3 - y2 ) * ( z4 - z2 ) +
                    ( x4 - x2 ) * ( y3 - y2 ) * ( zp - z2 ) - ( x4 - x2 ) * ( yp - y2 ) * ( z3 - z2 ) +
                    ( xp - x2 ) * ( y4 - y2 ) * ( z3 - z2 ) - ( x3 - x2 ) * ( y4 - y2 ) * ( zp - z2 ) ) / 6. / volume;

    answer.at(2) = ( ( x4 - x1 ) * ( yp - y1 ) * ( z3 - z1 ) - ( xp - x1 ) * ( y4 - y1 ) * ( z3 - z1 ) +
                    ( x3 - x1 ) * ( y4 - y1 ) * ( zp - z1 ) - ( x3 - x1 ) * ( yp - y1 ) * ( z4 - z1 ) +
                    ( xp - x1 ) * ( y3 - y1 ) * ( z4 - z1 ) - ( x4 - x1 ) * ( y3 - y1 ) * ( zp - z1 ) ) / 6. / volume;

    answer.at(3) = ( ( x2 - x1 ) * ( yp - y1 ) * ( z4 - z1 ) - ( xp - x1 ) * ( y2 - y1 ) * ( z4 - z1 ) +
                    ( x4 - x1 ) * ( y2 - y1 ) * ( zp - z1 ) - ( x4 - x1 ) * ( yp - y1 ) * ( z2 - z1 ) +
                    ( xp - x1 ) * ( y4 - y1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( y4 - y1 ) * ( zp - z1 ) ) / 6. / volume;

    // test if inside + clamping
    bool inside = true;
    for ( int i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            answer.at(i) = 0.;
            inside = false;
        } else if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            answer.at(i) = 1.;
            inside = false;
        }
    }

    answer.at(4) = 1.0 - answer.at(1) - answer.at(2) - answer.at(3);

    return inside;
}


double
FEI3dTetLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double detJ, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;

    x1 = cellgeo.giveVertexCoordinates(1)->at(1);
    x2 = cellgeo.giveVertexCoordinates(2)->at(1);
    x3 = cellgeo.giveVertexCoordinates(3)->at(1);
    x4 = cellgeo.giveVertexCoordinates(4)->at(1);

    y1 = cellgeo.giveVertexCoordinates(1)->at(2);
    y2 = cellgeo.giveVertexCoordinates(2)->at(2);
    y3 = cellgeo.giveVertexCoordinates(3)->at(2);
    y4 = cellgeo.giveVertexCoordinates(4)->at(2);

    z1 = cellgeo.giveVertexCoordinates(1)->at(3);
    z2 = cellgeo.giveVertexCoordinates(2)->at(3);
    z3 = cellgeo.giveVertexCoordinates(3)->at(3);
    z4 = cellgeo.giveVertexCoordinates(4)->at(3);

    detJ = ( ( x4 - x1 ) * ( y2 - y1 ) * ( z3 - z1 ) - ( x4 - x1 ) * ( y3 - y1 ) * ( z2 - z1 ) +
            ( x3 - x1 ) * ( y4 - y1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( y4 - y1 ) * ( z3 - z1 ) +
            ( x2 - x1 ) * ( y3 - y1 ) * ( z4 - z1 ) - ( x3 - x1 ) * ( y2 - y1 ) * ( z4 - z1 ) );

    return detJ;
}


void
FEI3dTetLin :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

void
FEI3dTetLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                            const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double coeff, l, x1, x2, y1, y2, z1, z2;
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    l = this->edgeComputeLength(edgeNodes, cellgeo);
    coeff = 1.0 / l / l;

    x1 = cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(1);
    y1 = cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(2);
    z1 = cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(3);
    x2 = cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(1);
    y2 = cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(2);
    z2 = cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(3);

    answer.resize(2, 3);
    answer.at(1, 1) = ( x1 - x2 ) * coeff;
    answer.at(1, 2) = ( y1 - y2 ) * coeff;
    answer.at(1, 3) = ( z1 - z2 ) * coeff;

    answer.at(2, 1) = ( x2 - x1 ) * coeff;
    answer.at(2, 2) = ( y2 - y1 ) * coeff;
    answer.at(2, 3) = ( z2 - z1 ) * coeff;
}

void
FEI3dTetLin :: edgeLocal2global(FloatArray &answer, int iedge,
                                const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray n;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

    answer.resize(3);
    answer.at(1) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(1) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(1) );
    answer.at(2) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(2) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(2) );
    answer.at(3) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(3) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(3) );
}


double
FEI3dTetLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    return 0.5 * this->edgeComputeLength(edgeNodes, cellgeo);
}


void
FEI3dTetLin :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    int aNode = 0, bNode = 0;
    edgeNodes.resize(2);

    if ( iedge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
    } else if ( iedge == 3 ) { // edge between nodes 3 1
        aNode = 3;
        bNode = 1;
    } else if ( iedge == 4 ) { // edge between nodes 1 4
        aNode = 1;
        bNode = 4;
    } else if ( iedge == 5 ) { // edge between nodes 2 4
        aNode = 2;
        bNode = 4;
    } else if ( iedge == 6 ) { // edge between nodes 3 4
        aNode = 3;
        bNode = 4;
    } else {
        OOFEM_ERROR("wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = ( aNode );
    edgeNodes.at(2) = ( bNode );
}

double
FEI3dTetLin :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    return cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->distance( cellgeo.giveVertexCoordinates( edgeNodes.at(1) ) );
}

void
FEI3dTetLin :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(3);

    answer.at(1) = lcoords.at(1);
    answer.at(2) = lcoords.at(2);
    answer.at(3) = 1. - lcoords.at(1) - lcoords.at(2);
}

void
FEI3dTetLin :: surfaceLocal2global(FloatArray &answer, int iedge,
                                   const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l1, l2, l3;
    answer.resize(3);
    IntArray nodes(3);

    computeLocalSurfaceMapping(nodes, iedge);

    l1 = lcoords.at(1);
    l2 = lcoords.at(2);
    l3 = 1.0 - l1 - l2;

    answer.at(1) = ( l1 * cellgeo.giveVertexCoordinates( nodes.at(1) )->at(1) +
                    l2 * cellgeo.giveVertexCoordinates( nodes.at(2) )->at(1) +
                    l3 * cellgeo.giveVertexCoordinates( nodes.at(3) )->at(1) );
    answer.at(2) = ( l1 * cellgeo.giveVertexCoordinates( nodes.at(1) )->at(2) +
                    l2 * cellgeo.giveVertexCoordinates( nodes.at(2) )->at(2) +
                    l3 * cellgeo.giveVertexCoordinates( nodes.at(3) )->at(2) );
    answer.at(3) = ( l1 * cellgeo.giveVertexCoordinates( nodes.at(1) )->at(3) +
                    l2 * cellgeo.giveVertexCoordinates( nodes.at(2) )->at(3) +
                    l3 * cellgeo.giveVertexCoordinates( nodes.at(3) )->at(3) );
}

void
FEI3dTetLin :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    // Note, this must be in correct order, not just the correct nodes, therefore we must use snodes;
    IntArray snodes;
    this->computeLocalSurfaceMapping(snodes, isurf);

    FloatArray lcoords_tet(4);
    lcoords_tet.at(snodes.at(1)) = lcoords.at(1);
    lcoords_tet.at(snodes.at(2)) = lcoords.at(2);
    lcoords_tet.at(snodes.at(3)) = 1. - lcoords.at(1) - lcoords.at(2);

    FloatMatrix fullB;
    this->evaldNdx(fullB, lcoords_tet, cellgeo);
    answer.resize(snodes.giveSize(), 3);
    for ( int i = 1; i <= snodes.giveSize(); ++i ) {
        for ( int j = 1; j <= 3; ++j ) {
            answer.at(i, j) = fullB.at(snodes.at(i), j);
        }
    }
}

double
FEI3dTetLin :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray a, b;
    IntArray snodes(3);
    this->computeLocalSurfaceMapping(snodes, isurf);

    a.beDifferenceOf( * cellgeo.giveVertexCoordinates( snodes.at(2) ), * cellgeo.giveVertexCoordinates( snodes.at(1) ) );
    b.beDifferenceOf( * cellgeo.giveVertexCoordinates( snodes.at(3) ), * cellgeo.giveVertexCoordinates( snodes.at(1) ) );
    answer.beVectorProductOf(a, b);

    return answer.normalize();
}

double
FEI3dTetLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                 const FEICellGeometry &cellgeo)
{
    FloatArray c;
    return this->surfaceEvalNormal(c, isurf, lcoords, cellgeo);
}

void
FEI3dTetLin :: computeLocalSurfaceMapping(IntArray &surfNodes, int isurf)
{
    int aNode = 0, bNode = 0, cNode = 0;
    surfNodes.resize(3);

    if ( isurf == 1 ) { // surface 1 - nodes 1 3 2
        aNode = 1;
        bNode = 3;
        cNode = 2;
    } else if ( isurf == 2 ) { // surface 2 - nodes 1 2 4
        aNode = 1;
        bNode = 2;
        cNode = 4;
    } else if ( isurf == 3 ) { // surface 3  - nodes 2 3 4
        aNode = 2;
        bNode = 3;
        cNode = 4;
    } else if ( isurf == 4 ) { // surface 4 - nodes 1 4 3
        aNode = 1;
        bNode = 4;
        cNode = 3;
    } else {
        OOFEM_ERROR("wrong surface number (%d)", isurf);
    }

    surfNodes.at(1) = ( aNode );
    surfNodes.at(2) = ( bNode );
    surfNodes.at(3) = ( cNode );
}

double
FEI3dTetLin :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo)
{
    IntArray fNodes;
    this->computeLocalSurfaceMapping(fNodes, iEdge);

    const FloatArray &c1 = * cellgeo.giveVertexCoordinates( fNodes.at(1) );
    const FloatArray &c2 = * cellgeo.giveVertexCoordinates( fNodes.at(2) );
    const FloatArray &c3 = * cellgeo.giveVertexCoordinates( fNodes.at(3) );

    return ( ( c2.at(1) * c3.at(2) - c3.at(1) * c2.at(2) ) * c1.at(3) +
            ( c3.at(1) * c1.at(2) - c1.at(1) * c3.at(2) ) * c2.at(3) +
            ( c1.at(1) * c2.at(2) - c2.at(1) * c1.at(2) ) * c3.at(3) ) * 0.5;
}

IntegrationRule *
FEI3dTetLin :: giveIntegrationRule(int order)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Tetrahedra, order + 0);
    iRule->SetUpPointsOnTetrahedra(points, _Unknown);
    return iRule;
}

IntegrationRule *
FEI3dTetLin :: giveBoundaryIntegrationRule(int order, int boundary)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order + 0);
    iRule->SetUpPointsOnTriangle(points, _Unknown);
    return iRule;
}
} // end namespace oofem
