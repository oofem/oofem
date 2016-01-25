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

#include "fei3dtrquad.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
void
FEI3dTrQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->surfaceEvalN(answer, 1, lcoords, cellgeo);
}

double
FEI3dTrQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dTrQuad :: evaldNdx - Not supported");
    return 0.;
}


void
FEI3dTrQuad :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->surfaceEvaldNdxi(answer, lcoords);
}


void
FEI3dTrQuad :: giveDerivativeXi(FloatArray &n, const FloatArray &lc)
{
    double l1, l2, l3;

    l1 = lc.at(1);
    l2 = lc.at(2);
    l3 = 1.0 - l1 - l2;

    n.resize(6);

    n.at(1) = 4.0 * l1 - 1.0;
    n.at(2) = 0.0;
    n.at(3) = -1.0 * ( 4.0 * l3 - 1.0 );
    n.at(4) = 4.0 * l2;
    n.at(5) = -4.0 * l2;
    n.at(6) = 4.0 * l3 - 4.0 * l1;
}

void
FEI3dTrQuad :: giveDerivativeEta(FloatArray &n, const FloatArray &lc)
{
    double l1, l2, l3;

    l1 = lc.at(1);
    l2 = lc.at(2);
    l3 = 1.0 - l1 - l2;

    n.resize(6);

    n.at(1) = 0.0;
    n.at(2) = 4.0 * l2 - 1.0;
    n.at(3) = -1.0 * ( 4.0 * l3 - 1.0 );
    n.at(4) = 4.0 * l1;
    n.at(5) = 4.0 * l3 - 4.0 * l2;
    n.at(6) = -4.0 * l1;
}


void
FEI3dTrQuad :: giveLocalNodeCoords(FloatMatrix &answer)
{

    answer.resize(3,6);
    answer.zero();
    answer.at(1,1) = 1.0;
    answer.at(1,2) = 0.0;
    answer.at(1,3) = 0.0;
    answer.at(1,4) = 0.5;
    answer.at(1,5) = 0.0;
    answer.at(1,6) = 0.5;

    answer.at(2,1) = 0.0;
    answer.at(2,2) = 1.0;
    answer.at(2,3) = 0.0;
    answer.at(2,4) = 0.5;
    answer.at(2,5) = 0.5;
    answer.at(2,6) = 0.0;

}


void
FEI3dTrQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);
    answer.clear();
    for ( int i = 1; i <= 6; ++i ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates(i) );
    }
}

#define POINT_TOL 1.e-3
int
FEI3dTrQuad :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    //OOFEM_ERROR("FEI3dTrQuad :: global2local - Not supported");
    //return -1;

    ///@todo this is for linear triangle
    int xind = 1;
    int yind = 2;
    
    
    double detJ, x1, x2, x3, y1, y2, y3;
    answer.resize(3);

    x1 = cellgeo.giveVertexCoordinates(1)->at(xind);
    x2 = cellgeo.giveVertexCoordinates(2)->at(xind);
    x3 = cellgeo.giveVertexCoordinates(3)->at(xind);

    y1 = cellgeo.giveVertexCoordinates(1)->at(yind);
    y2 = cellgeo.giveVertexCoordinates(2)->at(yind);
    y3 = cellgeo.giveVertexCoordinates(3)->at(yind);

    detJ = x1*(y2 - y3) + x2*(-y1 + y3) + x3*(y1 - y2);

    answer.at(1) = ( ( x2 * y3 - x3 * y2 ) + ( y2 - y3 ) * gcoords.at(xind) + ( x3 - x2 ) * gcoords.at(yind) ) / detJ;
    answer.at(2) = ( ( x3 * y1 - x1 * y3 ) + ( y3 - y1 ) * gcoords.at(xind) + ( x1 - x3 ) * gcoords.at(yind) ) / detJ;
    answer.at(3) = 1. - answer.at(1) - answer.at(2);

    // check if point is inside
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
    
    return inside;

}


void
FEI3dTrQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not supported");
}


void
FEI3dTrQuad :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double xi = lcoords.at(1);
    answer.resize(3);
    answer(0) = 0.5 * ( xi - 1.0 ) * xi;
    answer(1) = 0.5 * ( xi + 1.0 ) * xi;
    answer(2) = 1.0 - xi * xi;
}



void
FEI3dTrQuad :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                            const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not supported");
}

void
FEI3dTrQuad :: edgeEvaldNdxi(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double xi = lcoords.at(1);
    answer.resize(3);
    answer(0) = xi - 0.5;
    answer(1) = xi + 0.5;
    answer(2) = -2 * xi;
}

void
FEI3dTrQuad :: edgeLocal2global(FloatArray &answer, int iedge,
                                const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray N;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(N, iedge, lcoords, cellgeo);

    answer.clear();
    for ( int i = 0; i < N.giveSize(); ++i ) {
        answer.add( N(i), * cellgeo.giveVertexCoordinates( edgeNodes(i) ) );
    }
}


double
FEI3dTrQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    ///@todo Implement this
    OOFEM_ERROR("Not supported");
    return -1;
}


void
FEI3dTrQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    int aNode = 0, bNode = 0, cNode = 0;
    edgeNodes.resize(3);

    if ( iedge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
        cNode = 4;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
        cNode = 5;
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        aNode = 3;
        bNode = 1;
        cNode = 6;
    } else {
        OOFEM_ERROR("Wrong edge number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;
    edgeNodes.at(3) = cNode;
}

double
FEI3dTrQuad :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    ///@todo Implement this
    OOFEM_ERROR("Not supported");
    return -1;
}

void
FEI3dTrQuad :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1. - l1 - l2;

    answer.resize(6);

    answer.at(1) = ( 2. * l1 - 1. ) * l1;
    answer.at(2) = ( 2. * l2 - 1. ) * l2;
    answer.at(3) = ( 2. * l3 - 1. ) * l3;
    answer.at(4) = 4. * l1 * l2;
    answer.at(5) = 4. * l2 * l3;
    answer.at(6) = 4. * l3 * l1;
}

void
FEI3dTrQuad :: surfaceEvaldNdxi(FloatMatrix &answer, const FloatArray &lcoords)
{
    // Returns matrix with derivatives wrt local coordinates
    answer.resize(6, 2);
    FloatArray dndxi(6), dndeta(6);

    this->giveDerivativeXi(dndxi, lcoords);
    this->giveDerivativeEta(dndeta, lcoords);
    for ( int i = 1; i <= 6; ++i ) {
        answer.at(i, 1) = dndxi.at(i);
        answer.at(i, 2) = dndeta.at(i);
    }
}



void
FEI3dTrQuad :: surfaceLocal2global(FloatArray &answer, int isurf,
                                   const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    //Note: This gives the coordinate in the reference system
    FloatArray N;
    this->surfaceEvalN(N, isurf, lcoords, cellgeo);

    answer.clear();
    for ( int i = 0; i < N.giveSize(); ++i ) {
        answer.add( N(i), * cellgeo.giveVertexCoordinates(i) );
    }
}

void
FEI3dTrQuad :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    ///@todo Implement this
    OOFEM_ERROR("Not supported");
}

void
FEI3dTrQuad :: surfaceEvalBaseVectorsAt(FloatArray &G1, FloatArray &G2, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    // Note: These are not normalized. Returns the two tangent vectors to the surface.
    FloatMatrix dNdxi;
    this->surfaceEvaldNdxi(dNdxi, lcoords);

    G1.clear();
    G2.clear();
    for ( int i = 0; i < 6; ++i ) {
        G1.add( dNdxi(i, 1), * cellgeo.giveVertexCoordinates(i) );
        G2.add( dNdxi(i, 2), * cellgeo.giveVertexCoordinates(i) );
    }
}

double
FEI3dTrQuad :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray G1, G2; // local curvilinear base vectors
    this->surfaceEvalBaseVectorsAt(G1, G2, lcoords, cellgeo);
    answer.beVectorProductOf(G1, G2);
    double J = answer.computeNorm();
    answer.times(1 / J);
    return J;
}

void
FEI3dTrQuad :: surfaceGiveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
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
FEI3dTrQuad :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not supported yet");
    return 0;
}

void
FEI3dTrQuad :: computeLocalSurfaceMapping(IntArray &surfNodes, int isurf)
{
    //surfNodes.setValues(6, 1, 2, 3, 4, 5, 6);
    //surfNodes = {1, 2, 3, 4, 5, 6};
    ///@todo - fix wrt xfem
    computeLocalEdgeMapping(surfNodes, isurf);

}

IntegrationRule *
FEI3dTrQuad :: giveIntegrationRule(int order)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order);
    iRule->SetUpPointsOnTriangle(points, _Unknown);
    return iRule;
}

IntegrationRule *
FEI3dTrQuad :: giveBoundaryIntegrationRule(int order, int boundary)
{
    ///@todo Not sure about what defines boundaries on these elements. 2 surfaces + 3 edges? Ask Jim about this.
    OOFEM_ERROR("FEI3dTrQuad :: giveBoundaryIntegrationRule - Not supported");
    return NULL;
}


double
FEI3dTrQuad :: giveArea(const FEICellGeometry &cellgeo) const
{
    ///@todo this only correct for a planar triangle in the xy-plane
    const FloatArray *p;
    double x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6;

    p = cellgeo.giveVertexCoordinates(1);
    x1 = p->at(1);
    y1 = p->at(2);
    p = cellgeo.giveVertexCoordinates(2);
    x2 = p->at(1);
    y2 = p->at(2);
    p = cellgeo.giveVertexCoordinates(3);
    x3 = p->at(1);
    y3 = p->at(2);
    p = cellgeo.giveVertexCoordinates(4);
    x4 = p->at(1);
    y4 = p->at(2);
    p = cellgeo.giveVertexCoordinates(5);
    x5 = p->at(1);
    y5 = p->at(2);
    p = cellgeo.giveVertexCoordinates(6);
    x6 = p->at(1);
    y6 = p->at(2);

    return (4*(-(x4*y1) + x6*y1 + x4*y2 - x5*y2 + x5*y3 - x6*y3) + x2*(y1 - y3 - 4*y4 + 4*y5) +
            x1*(-y2 + y3 + 4*y4 - 4*y6) + x3*(-y1 + y2 - 4*y5 + 4*y6))/6;
}
} // end namespace oofem
