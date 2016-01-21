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

#include "fei3dtrlin.h"

#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
void
FEI3dTrLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->surfaceEvalN(answer, 1, lcoords, cellgeo);
}

double
FEI3dTrLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dTrLin :: evaldNdx - Not supported");
    return 0.;
}


void
FEI3dTrLin :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->surfaceEvaldNdxi(answer, lcoords);
}


void
FEI3dTrLin :: giveDerivativeXi(FloatArray &n, const FloatArray &lc)
{
    n.resize(3);
    n.at(1) =  1.0;;
    n.at(2) =  0.0;
    n.at(3) = -1.0;
}

void
FEI3dTrLin :: giveDerivativeEta(FloatArray &n, const FloatArray &lc)
{
    n.resize(3);
    n.at(1) =  0.0;;
    n.at(2) =  1.0;
    n.at(3) = -1.0;
}


void
FEI3dTrLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);
    answer.resize(0);
    for ( int i = 1; i <= 3; ++i ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates(i) );
    }
}

#define POINT_TOL 1.e-3
int
FEI3dTrLin :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dTrLin :: global2local - Not supported");
    return -1;

}


void
FEI3dTrLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dTrLin :: giveJacobianMatrixAt - Not supported");
}


void
FEI3dTrLin :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double xi = lcoords.at(1);
    answer.resize(2);
    answer.at(1) = ( 1. - xi ) * 0.5;
    answer.at(2) = ( 1. + xi ) * 0.5;
}



void
FEI3dTrLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                            const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dTrLin :: edgeEvaldNdx - Not supported");
}

void
FEI3dTrLin :: edgeEvaldNdxi(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(2);
    answer(0) = -0.5;
    answer(1) =  0.5;
}

void
FEI3dTrLin :: edgeLocal2global(FloatArray &answer, int iedge,
                                const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray N;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(N, iedge, lcoords, cellgeo);

    answer.resize(0);
    for ( int i = 0; i < N.giveSize(); ++i ) {
        answer.add( N(i), * cellgeo.giveVertexCoordinates( edgeNodes(i) ) );
    }
}


double
FEI3dTrLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    ///@todo Implement this
    OOFEM_ERROR("FEI3dTrLin :: edgeGiveTransformationJacobian - Not supported");
    return -1;
}


void
FEI3dTrLin :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{

    if ( iedge == 1 ) { // edge between nodes 1 2
        edgeNodes = { 1, 2 };

    } else if ( iedge == 2 ) { // edge between nodes 2 3
        edgeNodes = { 2, 3 };

    } else if ( iedge == 3 ) { // edge between nodes 2 3
        edgeNodes = { 3, 1 };

    } else {
        OOFEM_ERROR("Wrong edge number (%d)", iedge);
    }

}

double
FEI3dTrLin :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    ///@todo Implement this
    OOFEM_ERROR("FEI3dTrLin :: edgeComputeLength - Not supported");
    return -1;
}

void
FEI3dTrLin :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(3);
    answer.at(1) = lcoords.at(1);
    answer.at(2) = lcoords.at(2);
    answer.at(3) = 1. - lcoords.at(1) - lcoords.at(2);
}

void
FEI3dTrLin :: surfaceEvaldNdxi(FloatMatrix &answer, const FloatArray &lcoords)
{
    // Returns matrix with derivatives wrt local coordinates
    answer.resize(3, 2);
    FloatArray dndxi(3), dndeta(3);

    this->giveDerivativeXi(dndxi, lcoords);
    this->giveDerivativeEta(dndeta, lcoords);
    for ( int i = 1; i <= 3; ++i ) {
        answer.at(i, 1) = dndxi.at(i);
        answer.at(i, 2) = dndeta.at(i);
    }
}



void
FEI3dTrLin :: surfaceLocal2global(FloatArray &answer, int isurf,
                                   const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    //Note: This gives the coordinate in the reference system
    FloatArray N;
    this->surfaceEvalN(N, isurf, lcoords, cellgeo);

    answer.resize(0);
    for ( int i = 0; i < N.giveSize(); ++i ) {
        answer.add( N(i), * cellgeo.giveVertexCoordinates(i) );
    }
}

void
FEI3dTrLin :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    ///@todo Implement this
    OOFEM_ERROR("FEI3dTrLin :: surfaceEvaldNdx - Not supported");
}

void
FEI3dTrLin :: surfaceEvalBaseVectorsAt(FloatArray &G1, FloatArray &G2, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    // Note: These are not normalized. Returns the two tangent vectors to the surface.
    FloatMatrix dNdxi;
    this->surfaceEvaldNdxi(dNdxi, lcoords);

    G1.resize(0);
    G2.resize(0);
    for ( int i = 0; i < 3; ++i ) {
        G1.add( dNdxi(i, 1), * cellgeo.giveVertexCoordinates(i) );
        G2.add( dNdxi(i, 2), * cellgeo.giveVertexCoordinates(i) );
    }
}

double
FEI3dTrLin :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray G1, G2; // local curvilinear base vectors
    this->surfaceEvalBaseVectorsAt(G1, G2, lcoords, cellgeo);
    answer.beVectorProductOf(G1, G2);
    double J = answer.computeNorm();
    answer.times(1 / J);
    return J;
}

void
FEI3dTrLin :: surfaceGiveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
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
FEI3dTrLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dTrLin :: surfaceGiveTransformationJacobian - Not supported yet");
    return 0;
}

void
FEI3dTrLin :: computeLocalSurfaceMapping(IntArray &surfNodes, int isurf)
{
    //surfNodes.setValues(3, 1, 2, 3);
    computeLocalEdgeMapping(surfNodes, isurf);

}

IntegrationRule *
FEI3dTrLin :: giveIntegrationRule(int order)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order);
    iRule->SetUpPointsOnTriangle(points, _Unknown);
    return iRule;
}

IntegrationRule *
FEI3dTrLin :: giveBoundaryIntegrationRule(int order, int boundary)
{
    ///@todo Not sure about what defines boundaries on these elements. 2 surfaces + 3 edges? Ask Jim about this.
    OOFEM_ERROR("FEI3dTrLin :: giveBoundaryIntegrationRule - Not supported");
    return NULL;
}


double
FEI3dTrLin :: giveArea(const FEICellGeometry &cellgeo) const
{
    // A = 0.5 * |AB x AC|
    FloatArray AB, AC;
    AB = *cellgeo.giveVertexCoordinates(2) - *cellgeo.giveVertexCoordinates(1);
    AC = *cellgeo.giveVertexCoordinates(3) - *cellgeo.giveVertexCoordinates(1);
    FloatArray temp;
    temp.beVectorProductOf(AB, AC);
    return 0.5 * temp.computeNorm();

}
} // end namespace oofem
