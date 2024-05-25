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

#include "fei2dlinelin.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "gaussintegrationrule.h"

namespace oofem {

FloatArrayF<2> FEI2dLineLin :: evalN(double xi)
{
    return {( 1. - xi ) * 0.5, ( 1. + xi ) * 0.5};
}

void FEI2dLineLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double xi = lcoords(0);
    answer.resize(2);
    answer.at(1) = ( 1. - xi ) * 0.5;
    answer.at(2) = ( 1. + xi ) * 0.5;
}

double FEI2dLineLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    // Not meaningful to return anything.
    answer.clear();
    return 0.;
}

void FEI2dLineLin :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    answer.resize(2, 1);
    answer(0, 0) = -0.5;
    answer(1, 0) =  0.5;
}

void FEI2dLineLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);
    answer.resize( max(xind, yind) );
    answer.zero();
    answer.at(xind) = n(0) * cellgeo.giveVertexCoordinates(1).at(xind) +
                      n(1) * cellgeo.giveVertexCoordinates(2).at(xind);
    answer.at(yind) = n(0) * cellgeo.giveVertexCoordinates(1).at(yind) +
                      n(1) * cellgeo.giveVertexCoordinates(2).at(yind);
}

int FEI2dLineLin :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo) const
{
    double x2_x1 = cellgeo.giveVertexCoordinates(2).at(xind) - cellgeo.giveVertexCoordinates(1).at(xind);
    double y2_y1 = cellgeo.giveVertexCoordinates(2).at(yind) - cellgeo.giveVertexCoordinates(1).at(yind);

    // Projection of the global coordinate gives the value interpolated in [0,1].
    double dx = gcoords(0) - cellgeo.giveVertexCoordinates(1).at(xind);
    double dy = gcoords(1) - cellgeo.giveVertexCoordinates(1).at(yind);
    double xi = (x2_x1) ? ( sqrt( dx*dx*(1 + (y2_y1 / x2_x1)*(y2_y1 / x2_x1) )) ) / ( sqrt(x2_x1 * x2_x1 + y2_y1 * y2_y1) ) : sqrt(dy*dy) / ( sqrt(x2_x1 * x2_x1 + y2_y1 * y2_y1) );
    // Map to [-1,1] domain.
    xi = xi * 2 - 1;

    answer.resize(1);
    answer(0) = clamp(xi, -1., 1.);
    return false;
}

void FEI2dLineLin :: edgeEvaldNds(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double xi = lcoords(0);
    answer.resize(2);
    answer(0) = -0.5 * xi;
    answer(1) =  0.5 * xi;

    double es1 = answer(0) * cellgeo.giveVertexCoordinates(1).at(xind) +
                 answer(1) * cellgeo.giveVertexCoordinates(2).at(xind);
    double es2 = answer(0) * cellgeo.giveVertexCoordinates(1).at(yind) +
                 answer(1) * cellgeo.giveVertexCoordinates(2).at(yind);

    double J = sqrt(es1 * es1 + es2 * es2);
    answer.times(1 / J);
    //return J;
}

double FEI2dLineLin :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    normal.resize(2);
    normal.at(1) = cellgeo.giveVertexCoordinates(2).at(yind) - cellgeo.giveVertexCoordinates(1).at(yind);
    normal.at(2) = -( cellgeo.giveVertexCoordinates(2).at(xind) - cellgeo.giveVertexCoordinates(1).at(xind) );
    return normal.normalize() * 0.5;
}

double FEI2dLineLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double x2_x1 = cellgeo.giveVertexCoordinates(2).at(xind) - cellgeo.giveVertexCoordinates(1).at(xind);
    double y2_y1 = cellgeo.giveVertexCoordinates(2).at(yind) - cellgeo.giveVertexCoordinates(1).at(yind);
    return sqrt(x2_x1 * x2_x1 + y2_y1 * y2_y1) / 2.0;
}

IntArray FEI2dLineLin :: boundaryEdgeGiveNodes(int boundary, Element_Geometry_Type egt) const
{
    return {1, 2};
}

IntArray FEI2dLineLin :: computeLocalEdgeMapping(int iedge) const
{
    if ( iedge != 1 ) {
        OOFEM_ERROR("wrong egde number (%d)", iedge);
    }

    return {1, 2};
}

void FEI2dLineLin :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->evalN(answer, lcoords, cellgeo);
}

double FEI2dLineLin :: edgeComputeLength(const IntArray &edgeNodes, const FEICellGeometry &cellgeo) const
{
    double x2_x1 = cellgeo.giveVertexCoordinates(2).at(xind) - cellgeo.giveVertexCoordinates(1).at(xind);
    double y2_y1 = cellgeo.giveVertexCoordinates(2).at(yind) - cellgeo.giveVertexCoordinates(1).at(yind);
    return sqrt(x2_x1 * x2_x1 + y2_y1 * y2_y1);
}

double FEI2dLineLin :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo) const
{
    const auto &node1 = cellgeo.giveVertexCoordinates(1);
    double x1 = node1.at(xind);
    double y1 = node1.at(yind);

    const auto &node2 = cellgeo.giveVertexCoordinates(2);
    double x2 = node2.at(xind);
    double y2 = node2.at(yind);

    return x2 * y1 - x1 * y2;
}

std::unique_ptr<IntegrationRule> FEI2dLineLin :: giveIntegrationRule(int order, Element_Geometry_Type egt) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Line, order + 0);
    iRule->SetUpPointsOnLine(points, _Unknown);
    return std::move(iRule);
}
} // end namespace oofem
