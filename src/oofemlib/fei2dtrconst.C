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

#include "fei2dtrconst.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"

namespace oofem {
void
FEI2dTrConst :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer = FloatArray{1.};
}

double
FEI2dTrConst :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(1, 2);
    answer.zero();

    return 0.0;
}

void
FEI2dTrConst :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1.0 - l1 - l2;

    answer.resize(2);
    answer.at(1) = ( l1 * cellgeo.giveVertexCoordinates(1)->at(xind) +
                    l2 * cellgeo.giveVertexCoordinates(2)->at(xind) +
                    l3 * cellgeo.giveVertexCoordinates(3)->at(xind) );
    answer.at(2) = ( l1 * cellgeo.giveVertexCoordinates(1)->at(yind) +
                    l2 * cellgeo.giveVertexCoordinates(2)->at(yind) +
                    l3 * cellgeo.giveVertexCoordinates(3)->at(yind) );
}

#define POINT_TOL 1.e-3

int
FEI2dTrConst :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    double x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    double x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    double x3 = cellgeo.giveVertexCoordinates(3)->at(xind);

    double y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    double y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    double y3 = cellgeo.giveVertexCoordinates(3)->at(yind);

    double detJ = ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    answer.resize(3);
    answer.at(1) = ( ( x2 * y3 - x3 * y2 ) + ( y2 - y3 ) * coords.at(xind) + ( x3 - x2 ) * coords.at(yind) ) / detJ;
    answer.at(2) = ( ( x3 * y1 - x1 * y3 ) + ( y3 - y1 ) * coords.at(xind) + ( x1 - x3 ) * coords.at(yind) ) / detJ;
    //answer.at(3) = ( ( x1 * y2 - x2 * y1 ) + ( y1 - y2 ) * coords.at(xind) + ( x2 - x1 ) * coords.at(yind) ) / detJ;

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
    answer.at(3) = 1. - answer.at(1) - answer.at(2);

    return inside;
}


double
FEI2dTrConst :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    double x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    double x3 = cellgeo.giveVertexCoordinates(3)->at(xind);

    double y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    double y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    double y3 = cellgeo.giveVertexCoordinates(3)->at(yind);

    return x1 * ( y2 - y3 ) + x2 * ( -y1 + y3 ) + x3 * ( y1 - y2 );
}


void
FEI2dTrConst :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(1);
    answer.at(1) = 1.;
}

double
FEI2dTrConst :: edgeEvalNormal(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not applicable to constant interpolation");
    return 0.;
}

void
FEI2dTrConst :: edgeEvaldNds(FloatArray &answer, int iedge,
                             const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer = { 0., 0. };
}

void
FEI2dTrConst :: edgeLocal2global(FloatArray &answer, int iedge,
                                 const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n = { ( 1 - lcoords(0) ) * 0.5, ( 1 + lcoords(0) ) * 0.5 };
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);

    answer.resize(2);
    answer.at(1) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(xind) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(xind) );
    answer.at(2) = ( n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) )->at(yind) +
                    n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) )->at(yind) );
}

IntArray
FEI2dTrConst :: computeLocalEdgeMapping(int iedge) const
{
    if ( iedge == 1 ) { // edge between nodes 1 2
        return {1, 2};
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        return {2, 3};
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        return {3, 1};
    } else {
        throw std::range_error("invalid edge number");
        return {};
    }
}

double
FEI2dTrConst :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    auto nodeA = edgeNodes.at(1);
    auto nodeB = edgeNodes.at(2);

    double dx = cellgeo.giveVertexCoordinates(nodeB)->at(xind) - cellgeo.giveVertexCoordinates(nodeA)->at(xind);
    double dy = cellgeo.giveVertexCoordinates(nodeB)->at(yind) - cellgeo.giveVertexCoordinates(nodeA)->at(yind);
    return sqrt(dx * dx + dy * dy);
}

std::unique_ptr<IntegrationRule>
FEI2dTrConst :: giveIntegrationRule(int order)
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order + 0);
    iRule->SetUpPointsOnTriangle(points, _Unknown);
    return std::move(iRule);
}
} // end namespace oofem
