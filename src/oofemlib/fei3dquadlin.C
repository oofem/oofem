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

#include "fei3dquadlin.h"

#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"
#include <stdexcept>

namespace oofem {
void
FEI3dQuadLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->surfaceEvalN(answer, 1, lcoords, cellgeo);
}

double
FEI3dQuadLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dQuadLin :: evaldNdx - Not supported");
    return 0.;
}


void
FEI3dQuadLin :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->surfaceEvaldNdxi(answer, lcoords);
}


void
FEI3dQuadLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);
    answer.resize(0);
    for ( int i = 1; i <= 4; ++i ) {
        answer.add( n.at(i), cellgeo.giveVertexCoordinates(i) );
    }
}

#define POINT_TOL 1.e-3
int
FEI3dQuadLin :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dQuadLin :: global2local - Not supported");
    return -1;

}


void
FEI3dQuadLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dQuadLin :: giveJacobianMatrixAt - Not supported");
}


void
FEI3dQuadLin :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double xi = lcoords.at(1);
    answer.resize(2);
    answer.at(1) = ( 1. - xi ) * 0.5;
    answer.at(2) = ( 1. + xi ) * 0.5;
}



void
FEI3dQuadLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                            const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dQuadLin :: edgeEvaldNdx - Not supported");
}

void
FEI3dQuadLin :: edgeEvaldNdxi(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(2);
    answer(0) = -0.5;
    answer(1) =  0.5;
}

void
FEI3dQuadLin :: edgeLocal2global(FloatArray &answer, int iedge,
                                const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray N;
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    this->edgeEvalN(N, iedge, lcoords, cellgeo);

    answer.resize(0);
    for ( int i = 0; i < N.giveSize(); ++i ) {
        answer.add( N[i], cellgeo.giveVertexCoordinates( edgeNodes[i] ) );
    }
}


double
FEI3dQuadLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    const auto &edgeNodes = this->computeLocalEdgeMapping(iedge);
    ///@todo Implement this
    OOFEM_ERROR("FEI3dQuadLin :: edgeGiveTransformationJacobian - Not supported");
    return -1;
}


IntArray
FEI3dQuadLin :: computeLocalEdgeMapping(int iedge) const
{
    if ( iedge == 1 ) { // edge between nodes 1 2
        return { 1, 2 };
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        return { 2, 3 };
    } else if ( iedge == 3 ) { // edge between nodes 3 4
        return { 3, 4 };
    } else if ( iedge == 4 ) { // edge between nodes 4 1
        return { 4, 1 };
    } else {
        throw std::range_error("invalid edge number");
        return {};
    }
}

double
FEI3dQuadLin :: edgeComputeLength(const IntArray &edgeNodes, const FEICellGeometry &cellgeo) const
{
    ///@todo Implement this
    OOFEM_ERROR("FEI3dQuadLin :: edgeComputeLength - Not supported");
    return -1;
}

void
FEI3dQuadLin :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords[0];
    double eta = lcoords[1];

    answer.resize(4);
    answer.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
    answer.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
    answer.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
    answer.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25;
}

void
FEI3dQuadLin :: surfaceEvaldNdxi(FloatMatrix &answer, const FloatArray &lcoords)
{
    const double ksi = lcoords[0];
    const double eta = lcoords[1];

    answer.resize(4,2);
    // dn/dxi
    answer.at(1, 1) =  0.25 * ( 1. + eta );
    answer.at(2, 1) = -0.25 * ( 1. + eta );
    answer.at(3, 1) = -0.25 * ( 1. - eta );
    answer.at(4, 1) =  0.25 * ( 1. - eta );

    // dn/deta
    answer.at(1, 2) =  0.25 * ( 1. + ksi );
    answer.at(2, 2) =  0.25 * ( 1. - ksi );
    answer.at(3, 2) = -0.25 * ( 1. - ksi );
    answer.at(4, 2) = -0.25 * ( 1. + ksi );
}



void
FEI3dQuadLin :: surfaceLocal2global(FloatArray &answer, int isurf,
                                   const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    //Note: This gives the coordinate in the reference system
    FloatArray N;
    this->surfaceEvalN(N, isurf, lcoords, cellgeo);

    answer.resize(0);
    for ( int i = 1; i <= N.giveSize(); ++i ) {
        answer.add( N.at(i), cellgeo.giveVertexCoordinates(i) );
    }
}

void
FEI3dQuadLin :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    ///@todo Implement this
    OOFEM_ERROR("FEI3dQuadLin :: surfaceEvaldNdx - Not supported");
}

void
FEI3dQuadLin :: surfaceEvalBaseVectorsAt(FloatArray &G1, FloatArray &G2, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    // Note: These are not normalized. Returns the two tangent vectors to the surface.
    FloatMatrix dNdxi;
    this->surfaceEvaldNdxi(dNdxi, lcoords);

    G1.resize(0);
    G2.resize(0);
    for ( int i = 0; i < 4; ++i ) {
        G1.add( dNdxi(i, 1), cellgeo.giveVertexCoordinates(i) );
        G2.add( dNdxi(i, 2), cellgeo.giveVertexCoordinates(i) );
    }
}

double
FEI3dQuadLin :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray G1, G2; // local curvilinear base vectors
    this->surfaceEvalBaseVectorsAt(G1, G2, lcoords, cellgeo);
    answer.beVectorProductOf(G1, G2);
    double J = answer.computeNorm();
    answer.times(1 / J);
    return J;
}

void
FEI3dQuadLin :: surfaceGiveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    // Jacobian matrix consists of the three curvilinear base vectors. The third is taken as the normal to the surface.
    // Note! The base vectors are not normalized except the third (normal)
    FloatArray G1, G2, G3;
    this->surfaceEvalBaseVectorsAt(G1, G2, lcoords, cellgeo);
    G3.beVectorProductOf(G1, G2);

    jacobianMatrix.resize(3, 3);
    jacobianMatrix.at(1, 1) = G1.at(1);
    jacobianMatrix.at(1, 2) = G2.at(1);
    jacobianMatrix.at(1, 3) = G3.at(1);
    jacobianMatrix.at(2, 1) = G1.at(2);
    jacobianMatrix.at(2, 2) = G2.at(2);
    jacobianMatrix.at(2, 3) = G3.at(2);
    jacobianMatrix.at(3, 1) = G1.at(3);
    jacobianMatrix.at(3, 2) = G2.at(3);
    jacobianMatrix.at(3, 3) = G3.at(3);
}

double
FEI3dQuadLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dQuadLin :: surfaceGiveTransformationJacobian - Not supported yet");
    return 0;
}

IntArray
FEI3dQuadLin :: computeLocalSurfaceMapping(int isurf) const
{
    //surfNodes.setValues(3, 1, 2, 3);
    return computeLocalEdgeMapping(isurf);

}

std::unique_ptr<IntegrationRule>
FEI3dQuadLin :: giveIntegrationRule(int order)
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Square, order);
    iRule->SetUpPointsOnSquare(points, _Unknown);
    return std::move(iRule);
}

std::unique_ptr<IntegrationRule>
FEI3dQuadLin :: giveBoundaryIntegrationRule(int order, int boundary)
{
    ///@todo Not sure about what defines boundaries on these elements.
    OOFEM_ERROR("FEI3dQuadLin :: giveBoundaryIntegrationRule - Not supported");
    return nullptr;
}


double
FEI3dQuadLin :: giveArea(const FEICellGeometry &cellgeo) const
{
    ///@todo Not sure about what defines boundaries on these elements.
    OOFEM_ERROR("FEI3dQuadLin :: giveArea - Not supported");
    return 0.;

}
} // end namespace oofem
