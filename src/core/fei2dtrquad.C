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

#include "fei2dtrquad.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "gaussintegrationrule.h"

namespace oofem {

FloatArrayF<6>
FEI2dTrQuad :: evalN(const FloatArrayF<2> &lcoords)
{
    double l1 = lcoords[0];
    double l2 = lcoords[1];
    double l3 = 1. - l1 - l2;

    return {
        ( 2. * l1 - 1. ) * l1,
        ( 2. * l2 - 1. ) * l2,
        ( 2. * l3 - 1. ) * l3,
        4. * l1 * l2,
        4. * l2 * l3,
        4. * l3 * l1
    };
}

void
FEI2dTrQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
#if 0
    answer = evalN(lcoords);
#else
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1. - l1 - l2;

    answer = Vec6(
        ( 2. * l1 - 1. ) * l1,
        ( 2. * l2 - 1. ) * l2,
        ( 2. * l3 - 1. ) * l3,
        4. * l1 * l2,
        4. * l2 * l3,
        4. * l3 * l1
    );
#endif
}


FloatMatrixF<2,6>
FEI2dTrQuad :: evaldNdxi(const FloatArrayF<2> &lcoords) 
{
    double l1 = lcoords[0];
    double l2 = lcoords[1];
    double l3 = 1.0 - l1 - l2;

    return {
        4.0 * l1 - 1.0, // col0
        0.0,
        0.0, // col1
        4.0 * l2 - 1.0,
        -1.0 * ( 4.0 * l3 - 1.0 ), // col2
        -1.0 * ( 4.0 * l3 - 1.0 ),
        4.0 * l2, // col3
        4.0 * l1,
        -4.0 * l2, // col4
        4.0 * l3 - 4.0 * l2,
        4.0 * l3 - 4.0 * l1, // col5
        -4.0 * l1,
    };
}


std::pair<double,FloatMatrixF<2,6>>
FEI2dTrQuad :: evaldNdx(const FloatArrayF<2> &lcoords, const FEICellGeometry &cellgeo) const
{
    auto dn = evaldNdxi(lcoords);
    FloatMatrixF<2,2> jacT;
    for ( std::size_t i = 1; i <= dn.cols(); i++ ) {
        const auto &c = cellgeo.giveVertexCoordinates(i);
        double x = c.at(xind);
        double y = c.at(yind);

        ///@todo check transpose
        jacT(0, 0) += dn.at(1, i) * x;
        jacT(0, 1) += dn.at(1, i) * y;
        jacT(1, 0) += dn.at(2, i) * x;
        jacT(1, 1) += dn.at(2, i) * y;
    }

    return {det(jacT), dot(inv(jacT), dn)};
}


double
FEI2dTrQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
#if 0
    auto tmp = evaldNdx(lcoords, cellgeo);
    answer = tmp.second;
    return tmp.first;
#else
    FloatMatrix jacobianMatrix(2, 2), inv, dn;

    this->evaldNdxi(dn, lcoords, cellgeo);
    for ( int i = 1; i <= dn.giveNumberOfRows(); i++ ) {
        double x = cellgeo.giveVertexCoordinates(i).at(xind);
        double y = cellgeo.giveVertexCoordinates(i).at(yind);

        jacobianMatrix.at(1, 1) += dn.at(i, 1) * x;
        jacobianMatrix.at(1, 2) += dn.at(i, 1) * y;
        jacobianMatrix.at(2, 1) += dn.at(i, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(i, 2) * y;
    }
    inv.beInverseOf(jacobianMatrix);

    answer.beProductTOf(dn, inv);
    return jacobianMatrix.giveDeterminant();
#endif
}

void
FEI2dTrQuad :: evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double x1 = cellgeo.giveVertexCoordinates(1).at(xind);
    double x2 = cellgeo.giveVertexCoordinates(2).at(xind);
    double x3 = cellgeo.giveVertexCoordinates(3).at(xind);

    double y1 = cellgeo.giveVertexCoordinates(1).at(yind);
    double y2 = cellgeo.giveVertexCoordinates(2).at(yind);
    double y3 = cellgeo.giveVertexCoordinates(3).at(yind);

    double area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    double y23 = ( y2 - y3 ) / ( 2. * area );
    double x32 = ( x3 - x2 ) / ( 2. * area );

    double y31 = ( y3 - y1 ) / ( 2. * area );
    double x13 = ( x1 - x3 ) / ( 2. * area );

    answer.resize(6, 3);
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
FEI2dTrQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);

    answer.resize(2);
    answer.zero();
    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1) += n.at(i) * cellgeo.giveVertexCoordinates(i).at(xind);
        answer.at(2) += n.at(i) * cellgeo.giveVertexCoordinates(i).at(yind);
    }
}


void
FEI2dTrQuad :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double ksi = lcoords.at(1);
    double n3 = 1. - ksi * ksi;

    answer = Vec3( ( 1. - ksi - n3 ) * 0.5, ( 1. + ksi - n3 ) * 0.5, n3 );
}

void
FEI2dTrQuad :: edgeEvaldNds(FloatArray &answer, int iedge,
                            const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
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

    dxdxi.at(1) = dNdxi.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(xind) +
    dNdxi.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(xind) +
    dNdxi.at(3) * cellgeo.giveVertexCoordinates( edgeNodes.at(3) ).at(xind);
    dxdxi.at(2) = dNdxi.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(yind) +
    dNdxi.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(yind) +
    dNdxi.at(3) * cellgeo.giveVertexCoordinates( edgeNodes.at(3) ).at(yind);

    double J = dxdxi.computeNorm();
    answer = dNdxi;
    answer.times(1 / J);
    return J;

#endif
    double xi = lcoords.at(1);
    double J = edgeGiveTransformationJacobian(iedge, lcoords, cellgeo);
    answer = Vec3(
        ( xi - 0.5 ) / J,
        ( xi + 0.5 ) / J,
        -2 * xi / J
    );
}

double FEI2dTrQuad :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
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

    normal.at(2) = - dN1dxi *cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(xind) +
                   - dN2dxi *cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(xind) +
                   - dN3dxi *cellgeo.giveVertexCoordinates( edgeNodes.at(3) ).at(xind);

    return normal.normalize_giveNorm();
}

void
FEI2dTrQuad :: edgeLocal2global(FloatArray &answer, int iedge,
                                const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

    answer.resize(2);
    answer.at(1) = n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(xind) +
                   n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(xind) +
                   n.at(3) * cellgeo.giveVertexCoordinates( edgeNodes.at(3) ).at(xind);
    answer.at(2) = n.at(1) * cellgeo.giveVertexCoordinates( edgeNodes.at(1) ).at(yind) +
                   n.at(2) * cellgeo.giveVertexCoordinates( edgeNodes.at(2) ).at(yind) +
                   n.at(3) * cellgeo.giveVertexCoordinates( edgeNodes.at(3) ).at(yind);
}


IntArray
FEI2dTrQuad :: computeLocalEdgeMapping(int iedge) const
{
    if ( iedge == 1 ) { // edge between nodes 1 2
        return {1, 2, 4};
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        return {2, 3, 5};
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        return {3, 1, 6};
    } else {
        throw std::range_error("invalid edge number");
        //return {};
    }
}



void FEI2dTrQuad :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1.0 - l1 - l2;

    answer.resize(6, 2);

    answer.at(1, 1) =  4.0 * l1 - 1.0;
    answer.at(2, 1) =  0.0;
    answer.at(3, 1) = -1.0 * ( 4.0 * l3 - 1.0 );
    answer.at(4, 1) =  4.0 * l2;
    answer.at(5, 1) = -4.0 * l2;
    answer.at(6, 1) =  4.0 * l3 - 4.0 * l1;

    answer.at(1, 2) =  0.0;
    answer.at(2, 2) =  4.0 * l2 - 1.0;
    answer.at(3, 2) = -1.0 * ( 4.0 * l3 - 1.0 );
    answer.at(4, 2) =  4.0 * l1;
    answer.at(5, 2) =  4.0 * l3 - 4.0 * l2;
    answer.at(6, 2) = -4.0 * l1;
}


double
FEI2dTrQuad :: giveArea(const FEICellGeometry &cellgeo) const
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
    const auto &p4 = cellgeo.giveVertexCoordinates(4);
    double x4 = p4.at(xind);
    double y4 = p4.at(yind);
    const auto &p5 = cellgeo.giveVertexCoordinates(5);
    double x5 = p5.at(xind);
    double y5 = p5.at(yind);
    const auto &p6 = cellgeo.giveVertexCoordinates(6);
    double x6 = p6.at(xind);
    double y6 = p6.at(yind);

    return fabs( ( 4 * ( -( x4 * y1 ) + x6 * y1 + x4 * y2 - x5 * y2 + x5 * y3 - x6 * y3 ) + x2 * ( y1 - y3 - 4 * y4 + 4 * y5 ) +
            x1 * ( -y2 + y3 + 4 * y4 - 4 * y6 ) + x3 * ( -y1 + y2 - 4 * y5 + 4 * y6 ) ) / 6 );
}

bool FEI2dTrQuad :: inside(const FloatArray &lcoords) const
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

double
FEI2dTrQuad :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo) const
{
    const auto &eNodes = this->computeLocalEdgeMapping(iEdge);

    const auto &node1 = cellgeo.giveVertexCoordinates( eNodes.at(1) );
    double x1 = node1.at(xind);
    double y1 = node1.at(yind);

    const auto &node2 = cellgeo.giveVertexCoordinates( eNodes.at(2) );
    double x2 = node2.at(xind);
    double y2 = node2.at(yind);

    const auto &node3 = cellgeo.giveVertexCoordinates( eNodes.at(3) );
    double x3 = node3.at(xind);
    double y3 = node3.at(yind);

    return -( x1 * y2 - x2 * y1 + 4 * ( x3 * ( y1 - y2 ) + y3 * ( x2 - x1 ) ) ) / 3.0;
}

std::unique_ptr<IntegrationRule>
FEI2dTrQuad :: giveIntegrationRule(int order, Element_Geometry_Type egt) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order + 2);
    iRule->SetUpPointsOnTriangle(points, _Unknown);
    return std::move(iRule);
}
} // end namespace oofem
