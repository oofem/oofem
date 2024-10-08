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

#include "fei2dlinequad.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
void FEI2dLineQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double xi = lcoords(0);
    answer.resize(3);
    answer(0) = 0.5 * ( xi - 1.0 ) * xi;
    answer(1) = 0.5 * ( xi + 1.0 ) * xi;
    answer(2) = 1.0 - xi * xi;
}

double FEI2dLineQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    // Not meaningful to return anything.
    answer.clear();
    return 0.;
}

void
FEI2dLineQuad :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double xi = lcoords(0);
    answer.resize(3, 1);
    answer.at(1, 1) = xi - 0.5;
    answer.at(2, 1) = xi + 0.5;
    answer.at(3, 1) = -2.0 * xi;
}

void FEI2dLineQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);
    answer.resize(2);
    answer.at(1) = n(0) * cellgeo.giveVertexCoordinates(1).at(xind) +
                   n(1) * cellgeo.giveVertexCoordinates(2).at(xind) +
                   n(2) * cellgeo.giveVertexCoordinates(3).at(xind);
    answer.at(2) = n(0) * cellgeo.giveVertexCoordinates(1).at(yind) +
                   n(1) * cellgeo.giveVertexCoordinates(2).at(yind) +
                   n(2) * cellgeo.giveVertexCoordinates(3).at(yind);
}

int FEI2dLineQuad :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo) const
{
    double x1_x2, y1_y2, px_x3, py_y3, x3_x2_x1, y3_y2_y1;
    double b0, b1, b2, b3;

    x1_x2 = cellgeo.giveVertexCoordinates(1).at(xind) - cellgeo.giveVertexCoordinates(2).at(xind);
    y1_y2 = cellgeo.giveVertexCoordinates(1).at(yind) - cellgeo.giveVertexCoordinates(2).at(yind);
    px_x3 = gcoords.at(1) - cellgeo.giveVertexCoordinates(3).at(xind);
    py_y3 = gcoords.at(2) - cellgeo.giveVertexCoordinates(3).at(yind);
    x3_x2_x1 = 2 * cellgeo.giveVertexCoordinates(3).at(xind) - cellgeo.giveVertexCoordinates(2).at(xind) - cellgeo.giveVertexCoordinates(1).at(xind);
    y3_y2_y1 = 2 * cellgeo.giveVertexCoordinates(3).at(yind) - cellgeo.giveVertexCoordinates(2).at(yind) - cellgeo.giveVertexCoordinates(1).at(yind);

    b0 = 0.50 * ( x1_x2 * px_x3 + y1_y2 * py_y3 );
    b1 = 0.25 * ( x1_x2 * x1_x2 + y1_y2 * y1_y2 ) + x3_x2_x1 * px_x3 + y3_y2_y1 * py_y3;
    b2 = 0.75 * ( x1_x2 * x3_x2_x1 + y1_y2 * y3_y2_y1 );
    b3 = 0.50 * ( x3_x2_x1 * x3_x2_x1 + y3_y2_y1 * y3_y2_y1 );

    double r [ 3 ];
    int roots;
    cubic(b3, b2, b1, b0, & r [ 0 ], & r [ 1 ], & r [ 2 ], & roots);

    // Copy them over, along with boundary cases.
    int points = 2;
    double p [ 5 ] = {
        -1.0, 1.0
    };
    for ( int i = 0; i < roots; i++ ) {
        if ( r [ i ] > -1.0 && r [ i ] < 1.0 ) {
            // The cubic solver has pretty bad accuracy, performing a single newton iteration
            // which typically improves the solution by many many orders of magnitude.
            ///@todo Move this to the cubic solver.
            r [ i ] -= ( b0 + b1 * r [ i ] + b2 * r [ i ] * r [ i ] + b3 * r [ i ] * r [ i ] * r [ i ] ) / ( b1 + 2 * b2 * r [ i ] + 3 * b3 * r [ i ] * r [ i ] );
            p [ points ] = r [ i ];
            points++;
        }
    }

    double min_distance2 = 0.0, min_xi = 0;
    FloatArray f(2);
    answer.resize(1);

    for ( int i = 0; i < points; i++ ) {
        answer(0) = p [ i ];
        this->local2global(f, answer, cellgeo);
        double distance2 = distance_square(f, gcoords);
        if ( i == 0 || distance2 < min_distance2 ) {
            min_distance2 = distance2;
            min_xi = answer(0);
        }
    }

    answer(0) = clamp(min_xi, -1., 1.); // Make sure to only give coordinates within the geometry.

    return false; // It will never land exactly on the geometry.
}

IntArray FEI2dLineQuad :: computeLocalEdgeMapping(int iedge) const
{
    return {1, 2, 3};
}

void FEI2dLineQuad :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->evalN(answer, lcoords, cellgeo);
}

void FEI2dLineQuad :: edgeEvaldNds(FloatArray &answer, int iedge,
                                   const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double xi = lcoords(0);
    answer.resize(3);
    answer(0) = -0.5 + xi;
    answer(1) =  0.5 + xi;
    answer(2) = -2.0 * xi;

    double es1 = answer(0) * cellgeo.giveVertexCoordinates(1).at(xind) +
                 answer(1) * cellgeo.giveVertexCoordinates(2).at(xind) +
                 answer(2) * cellgeo.giveVertexCoordinates(3).at(xind);

    double es2 = answer(0) * cellgeo.giveVertexCoordinates(1).at(yind) +
                 answer(1) * cellgeo.giveVertexCoordinates(2).at(yind) +
                 answer(2) * cellgeo.giveVertexCoordinates(3).at(yind);

    double J = sqrt(es1 * es1 + es2 * es2);
    answer.times(1 / J);
    //return J;
}

double FEI2dLineQuad :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    double xi = lcoords(0);
    double dN1dxi = -0.5 + xi;
    double dN2dxi =  0.5 + xi;
    double dN3dxi = -2.0 * xi;

    normal.resize(2);

    normal.at(1) = dN1dxi * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(yind) +
                   dN2dxi * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(yind) +
                   dN3dxi * cellgeo.giveVertexCoordinates( edgeNodes.at(3) ).at(yind);

    normal.at(2) = - dN1dxi * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(xind) +
                   - dN2dxi * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(xind) +
                   - dN3dxi * cellgeo.giveVertexCoordinates( edgeNodes.at(3) ).at(xind);

    return normal.normalize();
}

void FEI2dLineQuad :: edgeLocal2global(FloatArray &answer, int iedge,
                                       const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->local2global(answer, lcoords, cellgeo);
}

double FEI2dLineQuad :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double xi = lcoords(0);
    double a1 = -0.5 + xi;
    double a2 =  0.5 + xi;
    double a3 = -2.0 * xi;

    double es1 = a1 * cellgeo.giveVertexCoordinates(1).at(xind) +
                 a2 * cellgeo.giveVertexCoordinates(2).at(xind) +
                 a3 * cellgeo.giveVertexCoordinates(3).at(xind);
    double es2 = a1 * cellgeo.giveVertexCoordinates(1).at(yind) +
                 a2 * cellgeo.giveVertexCoordinates(2).at(yind) +
                 a3 * cellgeo.giveVertexCoordinates(3).at(yind);

    return sqrt(es1 * es1 + es2 * es2);
}

void FEI2dLineQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double xi = lcoords(0);
    double dN1dxi = -0.5 + xi;
    double dN2dxi =  0.5 + xi;
    double dN3dxi = -2.0 * xi;

    double es1 = dN1dxi * cellgeo.giveVertexCoordinates(1).at(xind) +
                 dN2dxi * cellgeo.giveVertexCoordinates(2).at(xind) +
                 dN3dxi * cellgeo.giveVertexCoordinates(3).at(xind);
    double es2 = dN1dxi * cellgeo.giveVertexCoordinates(1).at(yind) +
                 dN2dxi * cellgeo.giveVertexCoordinates(2).at(yind) +
                 dN3dxi * cellgeo.giveVertexCoordinates(3).at(yind);

    double J = sqrt(es1 * es1 + es2 * es2);

    // Only used for determined deformation of element, not sure about the transpose here;
    jacobianMatrix.resize(2, 2);
    jacobianMatrix(0, 0) = es1 / J;
    jacobianMatrix(0, 1) = es2 / J;
    jacobianMatrix(1, 0) = -es2 / J;
    jacobianMatrix(1, 1) = es1 / J;
}

double FEI2dLineQuad :: edgeComputeLength(const IntArray &edgeNodes, const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented");
    //return 0.0;
}

double FEI2dLineQuad :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo) const
{
    const auto &node1 = cellgeo.giveVertexCoordinates(1);
    double x1 = node1.at(xind);
    double y1 = node1.at(yind);

    const auto &node2 = cellgeo.giveVertexCoordinates(2);
    double x2 = node2.at(xind);
    double y2 = node2.at(yind);

    const auto &node3 = cellgeo.giveVertexCoordinates(3);
    double x3 = node3.at(xind);
    double y3 = node3.at(yind);

    return ( x1 * y2 - x2 * y1 + 4 * ( x3 * ( y1 - y2 ) + y3 * ( x2 - x1 ) ) ) / 3.0;
}

std::unique_ptr<IntegrationRule> FEI2dLineQuad :: giveIntegrationRule(int order, Element_Geometry_Type egt) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Line, order + 1);
    iRule->SetUpPointsOnLine(points, _Unknown);
    return std::move(iRule);
}
} // end namespace oofem
