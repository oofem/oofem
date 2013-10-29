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

#include "fei3dwedgelin.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
void
FEI3dWedgeLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double u, v, w;
    answer.resize(6);

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    answer.at(1) = 0.5 * ( 1. - w ) * ( 1. - u - v );
    answer.at(2) = 0.5 * ( 1. - w ) * u;
    answer.at(3) = 0.5 * ( 1. - w ) * v;
    answer.at(4) = 0.5 * ( 1. + w ) * ( 1. - u - v );
    answer.at(5) = 0.5 * ( 1. + w ) * u;
    answer.at(6) = 0.5 * ( 1. + w ) * v;
}


double
FEI3dWedgeLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(2, 2), inv, dNduvw, coords;

    this->giveLocalDerivative(dNduvw, lcoords);
    coords.resize(3, 6);
    for ( int i = 1; i <= 6; i++ ) {
        coords.setColumn(*cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
    inv.beInverseOf(jacobianMatrix);

    answer.beProductOf(dNduvw, inv);
    return jacobianMatrix.giveDeterminant();
}


void
FEI3dWedgeLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;

    this->evalN(n, lcoords, cellgeo);

    answer.resize(0);
    for ( int i = 1; i <= 6; i++ ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates(i) );
    }
}


int
FEI3dWedgeLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
     OOFEM_ERROR("FEI3dHexaQuad :: global2local not implemented");
     return 1;
}


double
FEI3dWedgeLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(3, 3);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant()/2.; ///@todo Should this really be a factor 1/2 here?
}


void
FEI3dWedgeLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
{
    FloatMatrix dNduvw, coords;
    this->giveLocalDerivative(dNduvw, lcoords);
    coords.resize(3, 6);
    for ( int i = 1; i <= 6; i++ ) {
        coords.setColumn(*cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
}


void
FEI3dWedgeLin :: giveLocalDerivative(FloatMatrix &dN, const FloatArray &lcoords)
{
    double u, v, w;
    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    dN.resize(6, 3);

    dN.at(1, 1) = -0.5 * ( 1. - w );
    dN.at(2, 1) =  0.5 * ( 1. - w );
    dN.at(3, 1) =  0.;
    dN.at(4, 1) = -0.5 * ( 1. + w );
    dN.at(5, 1) =  0.5 * ( 1. + w );
    dN.at(6, 1) =  0.;

    dN.at(1, 2) = -0.5 * ( 1. - w );
    dN.at(2, 2) =  0.;
    dN.at(3, 2) =  0.5 * ( 1. - w );
    dN.at(4, 2) =  0.5 * ( 1. + w );
    dN.at(5, 2) =  0.;
    dN.at(6, 2) =  0.5 * ( 1. + w );

    dN.at(1, 3) = -0.5 * ( 1. - u - v );
    dN.at(2, 3) = -0.5 * u;
    dN.at(3, 3) = -0.5 * v;
    dN.at(4, 3) =  0.5 * ( 1. - u - v );
    dN.at(5, 3) =  0.5 * u;
    dN.at(6, 3) =  0.5 * v;
}


void
FEI3dWedgeLin :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer.resize(2);
    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. - ksi ) * 0.5;
}


void
FEI3dWedgeLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dWedgeLin :: edgeEvaldNdx not implemented");
}


void
FEI3dWedgeLin :: edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray nodes;
    FloatArray n;

    this->computeLocalEdgeMapping(nodes, iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

    answer.resize(0);
    for ( int i = 1; i <= n.giveSize(); ++i ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates(nodes.at(i)));
    }
}


void
FEI3dWedgeLin :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    if ( iedge == 1 ) {
        edgeNodes.setValues(2, 1, 2);
    } else if ( iedge == 2 ) {
        edgeNodes.setValues(2, 2, 3);
    } else if ( iedge == 3 ) {
        edgeNodes.setValues(2, 3, 1);
    } else if ( iedge == 4 ) {
        edgeNodes.setValues(2, 4, 5);
    } else if ( iedge == 5 ) {
        edgeNodes.setValues(2, 5, 6);
    } else if ( iedge == 6 ) {
        edgeNodes.setValues(2, 6, 4);
    } else if ( iedge == 7 ) {
        edgeNodes.setValues(2, 1, 4);
    } else if ( iedge == 8 ) {
        edgeNodes.setValues(2, 2, 5);
    } else if ( iedge == 9 ) {
        edgeNodes.setValues(2, 3, 6);
    } else {
        OOFEM_ERROR2("FEI3dWedgeQuad :: computeLocalEdgeMapping - Edge %d doesn't exist.\n", iedge);
    }
}


double
FEI3dWedgeLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dWedgeLin :: edgeGiveTransformationJacobian not implemented");
    return 0.0;
}


void
FEI3dWedgeLin :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    double eta = lcoords.at(2);

    if ( isurf <= 2) {
        answer.resize(3);
        answer.at(1) = ksi;
        answer.at(2) = eta;
        answer.at(3) = 1.0 - ksi - eta;
    } else {
        answer.resize(4);
        answer.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
        answer.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
        answer.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
        answer.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25;
    }
}


void
FEI3dWedgeLin :: surfaceLocal2global(FloatArray &answer, int isurf,
                                     const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray nodes;
    FloatArray n;

    this->computeLocalSurfaceMapping(nodes, isurf);
    this->surfaceEvalN(n, isurf, lcoords, cellgeo);

    answer.resize(0);
    for ( int i = 1; i <= n.giveSize(); ++i ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates(nodes.at(i)));
    }
}


void
FEI3dWedgeLin :: computeLocalSurfaceMapping(IntArray &nodes, int isurf)
{
    if ( isurf == 1 ) {
        nodes.setValues(3, 1, 2, 3);
    } else if ( isurf == 2 ) {
        nodes.setValues(3, 4, 5, 6);
    } else if ( isurf == 3 ) {
        nodes.setValues(4, 1, 2, 5, 4);
    } else if ( isurf == 4 ) {
        nodes.setValues(4, 2, 3, 6, 5);
    } else if ( isurf == 5 ) {
        nodes.setValues(4, 3, 1, 4, 6);
    } else {
        OOFEM_ERROR2("FEI3dWedgeQuad :: computeLocalSurfaceMapping - Surface %d doesn't exist.\n", isurf);
    }
}


double
FEI3dWedgeLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                   const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dWedgeLin :: surfaceGiveTransformationJacobian not implemented");
    return 0;
}


IntegrationRule *
FEI3dWedgeLin :: giveIntegrationRule(int order)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    ///@todo This function below isn't supported for wedges. We must decide how we should do this.
    //int points = iRule->getRequiredNumberOfIntegrationPoints(_Wedge, order);
    OOFEM_WARNING("Warning.. ignoring 'order' argument: FIXME");
    int pointsZeta = 1;
    int pointsTriangle = 1;
    iRule->SetUpPointsOnWedge(pointsTriangle, pointsZeta, _Unknown);
    return iRule;
}


IntegrationRule *
FEI3dWedgeLin :: giveBoundaryIntegrationRule(int order, int boundary)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    if ( boundary <= 2 ) {
        int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order + 0);
        iRule->SetUpPointsOnTriangle(points, _Unknown);
    } else {
        int points = iRule->getRequiredNumberOfIntegrationPoints(_Square, order + 2);
        iRule->SetUpPointsOnSquare(points, _Unknown);
    }
    return iRule;
}

} // end namespace oofem
