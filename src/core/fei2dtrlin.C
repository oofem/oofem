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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "fei2dtrlin.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "gaussintegrationrule.h"

namespace oofem {

FloatArrayF<3>
FEI2dTrLin :: evalN(const FloatArrayF<2> &lcoords)
{
    return {lcoords[0], lcoords[1], 1. - lcoords[0] - lcoords[1]};
}

void
FEI2dTrLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    answer = Vec3(
        lcoords.at(1),
        lcoords.at(2),
        1. - lcoords.at(1) - lcoords.at(2)
    );
}

std::pair<double, FloatMatrixF<2,3>>
FEI2dTrLin :: evaldNdx(const FEICellGeometry &cellgeo) const
{
    double x1 = cellgeo.giveVertexCoordinates(1).at(xind);
    double x2 = cellgeo.giveVertexCoordinates(2).at(xind);
    double x3 = cellgeo.giveVertexCoordinates(3).at(xind);
    double y1 = cellgeo.giveVertexCoordinates(1).at(yind);
    double y2 = cellgeo.giveVertexCoordinates(2).at(yind);
    double y3 = cellgeo.giveVertexCoordinates(3).at(yind);
    double detJ = x1 * ( y2 - y3 ) + x2 * ( y3 - y1 ) + x3 * ( y1 - y2 );

    FloatMatrixF<2,3> ans = {
        y2 - y3, x3 - x2,
        y3 - y1, x1 - x3,
        y1 - y2, x2 - x1,
    };

    return {detJ, ans * (1./detJ)};
}


double
FEI2dTrLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double x1, x2, x3, y1, y2, y3, detJ;

    x1 = cellgeo.giveVertexCoordinates(1).at(xind);
    x2 = cellgeo.giveVertexCoordinates(2).at(xind);
    x3 = cellgeo.giveVertexCoordinates(3).at(xind);

    y1 = cellgeo.giveVertexCoordinates(1).at(yind);
    y2 = cellgeo.giveVertexCoordinates(2).at(yind);
    y3 = cellgeo.giveVertexCoordinates(3).at(yind);

    detJ = x1 * ( y2 - y3 ) + x2 * ( -y1 + y3 ) + x3 * ( y1 - y2 );

    answer.resize(3, 2);
    answer.at(1, 1) = ( y2 - y3 ) / detJ;
    answer.at(1, 2) = ( x3 - x2 ) / detJ;

    answer.at(2, 1) = ( y3 - y1 ) / detJ;
    answer.at(2, 2) = ( x1 - x3 ) / detJ;

    answer.at(3, 1) = ( y1 - y2 ) / detJ;
    answer.at(3, 2) = ( x2 - x1 ) / detJ;

    return detJ;
}

void
FEI2dTrLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1.0 - l1 - l2;

    answer.resize(3);
    answer.zero();
    answer.at(1) = l1 * cellgeo.giveVertexCoordinates(1).at(xind) +
                   l2 * cellgeo.giveVertexCoordinates(2).at(xind) +
                   l3 * cellgeo.giveVertexCoordinates(3).at(xind);
    answer.at(2) = l1 * cellgeo.giveVertexCoordinates(1).at(yind) +
                   l2 * cellgeo.giveVertexCoordinates(2).at(yind) +
                   l3 * cellgeo.giveVertexCoordinates(3).at(yind);
}

#define POINT_TOL 1.e-3

int
FEI2dTrLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo) const
{
    double x1 = cellgeo.giveVertexCoordinates(1).at(xind);
    double x2 = cellgeo.giveVertexCoordinates(2).at(xind);
    double x3 = cellgeo.giveVertexCoordinates(3).at(xind);

    double y1 = cellgeo.giveVertexCoordinates(1).at(yind);
    double y2 = cellgeo.giveVertexCoordinates(2).at(yind);
    double y3 = cellgeo.giveVertexCoordinates(3).at(yind);

    double detJ = x1 * ( y2 - y3 ) + x2 * ( -y1 + y3 ) + x3 * ( y1 - y2 );

    answer.resize(3);
    answer.at(1) = ( ( x2 * y3 - x3 * y2 ) + ( y2 - y3 ) * coords.at(xind) + ( x3 - x2 ) * coords.at(yind) ) / detJ;
    answer.at(2) = ( ( x3 * y1 - x1 * y3 ) + ( y3 - y1 ) * coords.at(xind) + ( x1 - x3 ) * coords.at(yind) ) / detJ;

    // check if point is inside
    bool inside = true;
    for ( int i = 1; i <= 2; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            answer.at(i) = 0.;
            inside = false;
        } else if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            answer.at(i) = 1.;
            inside = false;
        }
    }

    if( ( answer.at(1) + answer.at(2)) > 1.0 ) {
        const double temp = 0.5*( answer.at(1) + answer.at(2) - 1.);
        answer.at(1) -= temp;
        answer.at(2) -= temp;
        inside = false;
    }

    answer.at(3) = 1. - answer.at(1) - answer.at(2);

    return inside;
}


double
FEI2dTrLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double x1 = cellgeo.giveVertexCoordinates(1).at(xind);
    double x2 = cellgeo.giveVertexCoordinates(2).at(xind);
    double x3 = cellgeo.giveVertexCoordinates(3).at(xind);

    double y1 = cellgeo.giveVertexCoordinates(1).at(yind);
    double y2 = cellgeo.giveVertexCoordinates(2).at(yind);
    double y3 = cellgeo.giveVertexCoordinates(3).at(yind);

    return ( x1 * ( y2 - y3 ) + x2 * ( -y1 + y3 ) + x3 * ( y1 - y2 ) );
}

bool FEI2dTrLin :: inside(const FloatArray &lcoords) const
{
	const double point_tol = 1.0e-3;
    bool inside = true;
    for ( int i = 1; i <= 2; i++ ) {
        if ( lcoords.at(i) < - point_tol ) {
            inside = false;
        } else if ( lcoords.at(i) > ( 1. + point_tol ) ) {
            inside = false;
        }
    }

    if ( 1. - lcoords.at(1) - lcoords.at(2) < - point_tol ) {
        inside = false;
    } else if ( 1. - lcoords.at(1) - lcoords.at(2) > ( 1. + point_tol ) ) {
        inside = false;
    }

    return inside;
}

void
FEI2dTrLin :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double ksi = lcoords.at(1);
    answer = Vec2( ( 1. - ksi ) * 0.5, ( 1. + ksi ) * 0.5 );
}

void
FEI2dTrLin :: edgeEvaldNds(FloatArray &answer, int iedge,
                           const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    double l = this->edgeComputeLength(edgeNodes, cellgeo);

    answer = Vec2( -1.0 / l, 1.0 / l );
}

double FEI2dTrLin :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    normal = Vec2(
        cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(yind) - cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(yind),
        cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(xind) - cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(xind)
    );
    return normal.normalize_giveNorm() * 0.5;
}

void
FEI2dTrLin :: edgeLocal2global(FloatArray &answer, int iedge,
                               const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

    answer.resize(2);
    answer.at(1) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(xind) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(xind) );
    answer.at(2) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(yind) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(yind) );
}


IntArray
FEI2dTrLin :: computeLocalEdgeMapping(int iedge) const
{
    if ( iedge == 1 ) { // edge between nodes 1 2
        return {1, 2};
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        return {2, 3};
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        return {3, 1};
    } else {
        throw std::range_error("invalid edge number");
        //return {};
    }
}

double
FEI2dTrLin :: edgeComputeLength(const IntArray &edgeNodes, const FEICellGeometry &cellgeo) const
{
    int nodeA = edgeNodes.at(1);
    int nodeB = edgeNodes.at(2);

    double dx = cellgeo.giveVertexCoordinates(nodeB).at(xind) - cellgeo.giveVertexCoordinates(nodeA).at(xind);
    double dy = cellgeo.giveVertexCoordinates(nodeB).at(yind) - cellgeo.giveVertexCoordinates(nodeA).at(yind);
    return sqrt(dx * dx + dy * dy);
}

double
FEI2dTrLin :: giveArea(const FEICellGeometry &cellgeo) const
{
    const auto &p1 = cellgeo.giveVertexCoordinates(1);
    double x1 = p1.at(xind);
    double y1 = p1.at(yind);
    const auto &p2 = cellgeo.giveVertexCoordinates(2);
    double x2 = p2.at(xind);
    double y2 = p2.at(yind);
    const auto &p3 = cellgeo.giveVertexCoordinates(3);
    double x3 = p3.at(xind);
    double y3 = p3.at(yind);

    return fabs( 0.5 * ( x1 * ( y2 - y3 ) + x2 * ( -y1 + y3 ) + x3 * ( y1 - y2 ) ) );
}

double
FEI2dTrLin :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo) const
{
    const auto &eNodes = this->computeLocalEdgeMapping(iEdge);

    const auto &node1 = cellgeo.giveVertexCoordinates( eNodes.at(1) );
    double x1 = node1.at(xind);
    double y1 = node1.at(yind);

    const auto &node2 = cellgeo.giveVertexCoordinates( eNodes.at(2) );
    double x2 = node2.at(xind);
    double y2 = node2.at(yind);

    return -( x2 * y1 - x1 * y2 );
}

std::unique_ptr<IntegrationRule>
FEI2dTrLin :: giveIntegrationRule(int order, Element_Geometry_Type egt) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order + 0);
    iRule->SetUpPointsOnTriangle(points, _Unknown);
    return std::move(iRule);
}

// FEI2dTrLinAxi element

double
FEI2dTrLinAxi :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray N;
    this->evalN( N, lcoords, cellgeo);

    double r = 0.0;
    for ( int i = 1; i <= 3; i++ ) {
        double x = cellgeo.giveVertexCoordinates(i).at(1);
        r += x * N.at(i);
    }

    return r * FEI2dTrLin::giveTransformationJacobian(lcoords, cellgeo);
}

double
FEI2dTrLinAxi::edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords,
                                                const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

    double r = n.at(1)*cellgeo.giveVertexCoordinates(edgeNodes.at(1)).at(1) + n.at(2)*cellgeo.giveVertexCoordinates(edgeNodes.at(2)).at(1);
    return r * FEI2dTrLin::edgeGiveTransformationJacobian(iedge, lcoords, cellgeo);

}

double
FEI2dTrLinAxi::boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

double
FEI2dTrLinAxi::boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

} // end namespace oofem
