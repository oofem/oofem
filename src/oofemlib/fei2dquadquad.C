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

#include "fei2dquadquad.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "gaussintegrationrule.h"

namespace oofem {
double
FEI2dQuadQuad :: giveArea(const FEICellGeometry &cellgeo) const
{
    double x1, x2, x3, x4, y1, y2, y3, y4;
    double x85, x56, x67, x78, y85, y56, y67, y78;

    const auto &node1 = cellgeo.giveVertexCoordinates(1);
    const auto &node2 = cellgeo.giveVertexCoordinates(2);
    const auto &node3 = cellgeo.giveVertexCoordinates(3);
    const auto &node4 = cellgeo.giveVertexCoordinates(4);
    const auto &node5 = cellgeo.giveVertexCoordinates(5);
    const auto &node6 = cellgeo.giveVertexCoordinates(6);
    const auto &node7 = cellgeo.giveVertexCoordinates(7);
    const auto &node8 = cellgeo.giveVertexCoordinates(8);

    x1 = node1.at(xind);
    x2 = node2.at(xind);
    x3 = node3.at(xind);
    x4 = node4.at(xind);

    y1 = node1.at(yind);
    y2 = node2.at(yind);
    y3 = node3.at(yind);
    y4 = node4.at(yind);

    x85 = node8.at(xind) - node5.at(xind);
    x56 = node5.at(xind) - node6.at(xind);
    x67 = node6.at(xind) - node7.at(xind);
    x78 = node7.at(xind) - node8.at(xind);

    y85 = node8.at(yind) - node5.at(yind);
    y56 = node5.at(yind) - node6.at(yind);
    y67 = node6.at(yind) - node7.at(yind);
    y78 = node7.at(yind) - node8.at(yind);

    double p1 = ( x2 - x4 ) * ( y1 - y3 ) - ( x1 - x3 ) * ( y2 - y4 );
    double p2 = y1 * x85 + y2 * x56 + y3 * x67 + y4 * x78 - x1 * y85 - x2 * y56 - x3 * y67 - x4 * y78;

    return fabs(p1 + p2 * 4.0) / 6.; // Expression derived with mathematica, but not verified in any computations
}


FloatArrayF<8>
FEI2dQuadQuad :: evalN(const FloatArrayF<2> &lcoords)
{
    double ksi = lcoords[0];
    double eta = lcoords[1];
    return {
        ( 1. + ksi ) * ( 1. + eta ) * 0.25 * ( ksi + eta - 1. ),
        ( 1. - ksi ) * ( 1. + eta ) * 0.25 * ( -ksi + eta - 1. ),
        ( 1. - ksi ) * ( 1. - eta ) * 0.25 * ( -ksi - eta - 1. ),
        ( 1. + ksi ) * ( 1. - eta ) * 0.25 * ( ksi - eta - 1. ),
        0.5 * ( 1. - ksi * ksi ) * ( 1. + eta ),
        0.5 * ( 1. - ksi ) * ( 1. - eta * eta ),
        0.5 * ( 1. - ksi * ksi ) * ( 1. - eta ),
        0.5 * ( 1. + ksi ) * ( 1. - eta * eta )
    };
}


void
FEI2dQuadQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double ksi = lcoords.at(1);
    double eta = lcoords.at(2);

    answer = {
        ( 1. + ksi ) * ( 1. + eta ) * 0.25 * ( ksi + eta - 1. ),
        ( 1. - ksi ) * ( 1. + eta ) * 0.25 * ( -ksi + eta - 1. ),
        ( 1. - ksi ) * ( 1. - eta ) * 0.25 * ( -ksi - eta - 1. ),
        ( 1. + ksi ) * ( 1. - eta ) * 0.25 * ( ksi - eta - 1. ),
        0.5 * ( 1. - ksi * ksi ) * ( 1. + eta ),
        0.5 * ( 1. - ksi ) * ( 1. - eta * eta ),
        0.5 * ( 1. - ksi * ksi ) * ( 1. - eta ),
        0.5 * ( 1. + ksi ) * ( 1. - eta * eta )
    };
}


FloatMatrixF<2,8>
FEI2dQuadQuad :: evaldNdxi(const FloatArrayF<2> &lcoords)
{
    double ksi = lcoords.at(1);
    double eta = lcoords.at(2);

    return {
        0.25 * ( 1. + eta ) * ( 2.0 * ksi + eta ), // 1
        0.25 * ( 1. + ksi ) * ( 2.0 * eta + ksi ),
        -0.25 * ( 1. + eta ) * ( -2.0 * ksi + eta ), // 2
        0.25 * ( 1. - ksi ) * ( 2.0 * eta - ksi ),
        -0.25 * ( 1. - eta ) * ( -2.0 * ksi - eta ), // 3
        -0.25 * ( 1. - ksi ) * ( -2.0 * eta - ksi ),
        0.25 * ( 1. - eta ) * ( 2.0 * ksi - eta ), // 4
        -0.25 * ( 1. + ksi ) * ( -2.0 * eta + ksi ),
        -ksi * ( 1. + eta ), // 5
        0.5 * ( 1. - ksi * ksi ),
        -0.5 * ( 1. - eta * eta ), // 6
        -eta * ( 1. - ksi ),
        -ksi * ( 1. - eta ), // 7
        -0.5 * ( 1. - ksi * ksi ),
        0.5 * ( 1. - eta * eta ), // 83
        -eta * ( 1. + ksi )
    };
}


void FEI2dQuadQuad :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
#if 0
    answer = dNdxi(lcoords);
#else
    double ksi = lcoords.at(1);
    double eta = lcoords.at(2);
    answer.resize(8, 2);

    // dn/dxi
    answer.at(1, 1) =  0.25 * ( 1. + eta ) * ( 2.0 * ksi + eta );
    answer.at(2, 1) = -0.25 * ( 1. + eta ) * ( -2.0 * ksi + eta );
    answer.at(3, 1) = -0.25 * ( 1. - eta ) * ( -2.0 * ksi - eta );
    answer.at(4, 1) =  0.25 * ( 1. - eta ) * ( 2.0 * ksi - eta );
    answer.at(5, 1) = -ksi * ( 1. + eta );
    answer.at(6, 1) = -0.5 * ( 1. - eta * eta );
    answer.at(7, 1) = -ksi * ( 1. - eta );
    answer.at(8, 1) =  0.5 * ( 1. - eta * eta );

    answer.at(1, 2) =  0.25 * ( 1. + ksi ) * ( 2.0 * eta + ksi );
    answer.at(2, 2) =  0.25 * ( 1. - ksi ) * ( 2.0 * eta - ksi );
    answer.at(3, 2) = -0.25 * ( 1. - ksi ) * ( -2.0 * eta - ksi );
    answer.at(4, 2) = -0.25 * ( 1. + ksi ) * ( -2.0 * eta + ksi );
    answer.at(5, 2) =  0.5 * ( 1. - ksi * ksi );
    answer.at(6, 2) = -eta * ( 1. - ksi );
    answer.at(7, 2) = -0.5 * ( 1. - ksi * ksi );
    answer.at(8, 2) = -eta * ( 1. + ksi );
#endif
}


std::pair<double, FloatMatrixF<2,8>> 
FEI2dQuadQuad :: evaldNdx(const FloatArrayF<2> &lcoords, const FEICellGeometry &cellgeo) const
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
FEI2dQuadQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
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
FEI2dQuadQuad :: local2global(FloatArray &answer, const FloatArray &lcoords,  const FEICellGeometry &cellgeo) const
{
    FloatArray n;

    this->evalN(n, lcoords, cellgeo);

    answer.resize(2);
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer.at(1) += n.at(i) * cellgeo.giveVertexCoordinates(i).at(xind);
        answer.at(2) += n.at(i) * cellgeo.giveVertexCoordinates(i).at(yind);
    }
}

double FEI2dQuadQuad :: giveCharacteristicLength(const FEICellGeometry &cellgeo) const
{
    const auto &n1 = cellgeo.giveVertexCoordinates(1);
    const auto &n2 = cellgeo.giveVertexCoordinates(3);
    return distance(n1, n2);
}

bool FEI2dQuadQuad :: inside(const FloatArray &lcoords) const
{
	const double point_tol = 1.0e-3;
    bool inside = true;
    for ( int i = 1; i <= 2; i++ ) {
        if ( lcoords.at(i) < ( -1. - point_tol ) ) {
            inside = false;
        } else if ( lcoords.at(i) > ( 1. + point_tol ) ) {
            inside = false;
        }
    }

    return inside;
}


void
FEI2dQuadQuad :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    // 1-------3-------2

    double ksi = lcoords.at(1);
    double n3 = 1. - ksi * ksi;

    answer = { ( 1. - ksi - n3 ) * 0.5, ( 1. + ksi - n3 ) * 0.5, n3 };
}

void
FEI2dQuadQuad :: edgeEvaldNds(FloatArray &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double ksi = lcoords.at(1);
    answer = { ksi - 0.5, ksi + 0.5, ksi * 2.0 };
}

void
FEI2dQuadQuad :: edgeLocal2global(FloatArray &answer, int iedge,
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
FEI2dQuadQuad :: computeLocalEdgeMapping(int iedge) const
{
    if ( iedge == 1 ) { // edge between nodes 1 2
        return {1, 2, 5};
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        return {2, 3, 6};
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        return {3, 4, 7};
    } else if ( iedge == 4 ) { // edge between nodes 3 4
        return {4, 1, 8};
    } else {
        throw std::range_error("invalid edge number");
        return {};
    }
}

double FEI2dQuadQuad :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
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


double
FEI2dQuadQuad :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo) const
{
    auto eNodes = this->computeLocalEdgeMapping(iEdge);

    const FloatArray &node1 = cellgeo.giveVertexCoordinates( eNodes.at(1) );
    double x1 = node1.at(xind);
    double y1 = node1.at(yind);

    const FloatArray &node2 = cellgeo.giveVertexCoordinates( eNodes.at(2) );
    double x2 = node2.at(xind);
    double y2 = node2.at(yind);

    const FloatArray &node3 = cellgeo.giveVertexCoordinates( eNodes.at(3) );
    double x3 = node3.at(xind);
    double y3 = node3.at(yind);

    return -( x1 * y2 - x2 * y1 + 4 * ( x3 * ( y1 - y2 ) + y3 * ( x2 - x1 ) ) ) / 3.0;
}


std::unique_ptr<IntegrationRule> 
FEI2dQuadQuad :: giveIntegrationRule(int order) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Square, order + 4);
    iRule->SetUpPointsOnSquare(points, _Unknown);
    return std::move(iRule);
}


/*
 * FEI2dQuadQuadAxi element
 */
double
FEI2dQuadQuadAxi :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatArray N;
    this->evalN( N, lcoords, cellgeo);

    double r = 0.0;
    for ( int i = 1; i <= 8; i++ ) {
        double x = cellgeo.giveVertexCoordinates(i).at(1);
        r += x * N.at(i);
    }

    return r * FEI2dQuadQuad::giveTransformationJacobian(lcoords, cellgeo);
}

double
FEI2dQuadQuadAxi::edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords,
                                                 const FEICellGeometry &cellgeo) const
{
    FloatArray n;
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

    double r = n.at(1)*cellgeo.giveVertexCoordinates(edgeNodes.at(1)).at(1) +
               n.at(2)*cellgeo.giveVertexCoordinates(edgeNodes.at(2)).at(1) +
               n.at(3)*cellgeo.giveVertexCoordinates(edgeNodes.at(3)).at(1);
    return r * FEI2dQuadQuad::edgeGiveTransformationJacobian(iedge, lcoords, cellgeo);

}

double
FEI2dQuadQuadAxi::boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

double
FEI2dQuadQuadAxi::boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

} // end namespace oofem
