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
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "intarray.h"
#include "gaussintegrationrule.h"

namespace oofem {

FloatArrayF<6>
FEI3dWedgeLin :: evalN(const FloatArrayF<3> &lcoords)
{
    double u = lcoords[0];
    double v = lcoords[1];
    double w = lcoords[2];
    double x = 1. - u - v;
    return {
        0.5 * ( 1. - w ) * x,
        0.5 * ( 1. - w ) * u,
        0.5 * ( 1. - w ) * v,
        0.5 * ( 1. + w ) * x,
        0.5 * ( 1. + w ) * u,
        0.5 * ( 1. + w ) * v,
    };
}


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


FloatMatrixF<3,6>
FEI3dWedgeLin :: evaldNdxi(const FloatArrayF<3> &lcoords)
{
    double u = lcoords.at(1);
    double v = lcoords.at(2);
    double w = lcoords.at(3);
    double x = 1. - u - v;

    return {
        -0.5 * ( 1. - w ),
        -0.5 * ( 1. - w ),
        -0.5 * x,
        0.5 * ( 1. - w ),
        0.,
        -0.5 * u,
        0.,
        0.5 * ( 1. - w ),
        -0.5 * v,
        -0.5 * ( 1. + w ),
        -0.5 * ( 1. + w ),
        0.5 * x,
        0.5 * ( 1. + w ),
        0.,
        0.5 * u,
        0.,
        0.5 * ( 1. + w ),
        0.5 * v,
    };
}


void
FEI3dWedgeLin :: evaldNdxi(FloatMatrix &dN, const FloatArray &lcoords, const FEICellGeometry &)
{
    double u = lcoords.at(1);
    double v = lcoords.at(2);
    double w = lcoords.at(3);

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
    dN.at(4, 2) = -0.5 * ( 1. + w );
    dN.at(5, 2) =  0.;
    dN.at(6, 2) =  0.5 * ( 1. + w );

    dN.at(1, 3) = -0.5 * ( 1. - u - v );
    dN.at(2, 3) = -0.5 * u;
    dN.at(3, 3) = -0.5 * v;
    dN.at(4, 3) =  0.5 * ( 1. - u - v );
    dN.at(5, 3) =  0.5 * u;
    dN.at(6, 3) =  0.5 * v;
}


std::pair<double, FloatMatrixF<3,6>>
FEI3dWedgeLin :: evaldNdx(const FloatArrayF<3> &lcoords, const FEICellGeometry &cellgeo)
{
    auto dNduvw = evaldNdxi(lcoords);
    FloatMatrixF<3,6> coords;
    for ( int i = 0; i < 6; i++ ) {
        coords.setColumn(cellgeo.giveVertexCoordinates(i+1), i);
    }
    auto jacT = dotT(dNduvw, coords);
    return {det(jacT), dot(inv(jacT), dNduvw)};
}


double
FEI3dWedgeLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix, inv, dNduvw, coords;

    this->evaldNdxi(dNduvw, lcoords, cellgeo);
    coords.resize(3, 6);
    for ( int i = 1; i <= 6; i++ ) {
        coords.setColumn(cellgeo.giveVertexCoordinates(i), i);
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

    answer.clear();
    for ( int i = 1; i <= 6; i++ ) {
        answer.add( n.at(i), cellgeo.giveVertexCoordinates(i) );
    }
}


double FEI3dWedgeLin :: giveCharacteristicLength(const FEICellGeometry &cellgeo) const
{
    const auto &n1 = cellgeo.giveVertexCoordinates(1);
    const auto &n2 = cellgeo.giveVertexCoordinates(6);
    ///@todo Change this so that it is not dependent on node order.
    return distance(n1, n2);
}

#define POINT_TOL 1.e-3

int
FEI3dWedgeLin :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    FloatArray res, delta, guess;
    FloatMatrix jac;
    double convergence_limit, error = 0.0;

    // find a suitable convergence limit
    convergence_limit = 1e-6 * this->giveCharacteristicLength(cellgeo);

    // setup initial guess
    answer.resize( gcoords.giveSize() );
    answer.zero();

    // apply Newton-Raphson to solve the problem
    for ( int nite = 0; nite < 10; nite++ ) {
        // compute the residual
        this->local2global(guess, answer, cellgeo);
        res.beDifferenceOf(gcoords, guess);

        // check for convergence
        error = res.computeNorm();
        if ( error < convergence_limit ) {
            break;
        }

        // compute the corrections
        this->giveJacobianMatrixAt(jac, answer, cellgeo);
        jac.solveForRhs(res, delta);

        // update guess
        answer.add(delta);
    }
    if ( error > convergence_limit ) { // Imperfect, could give false negatives.
        //OOFEM_ERROR("no convergence after 10 iterations");
        answer = {1. / 3., 1. / 3., 1. / 3.};
        return false;
    }

    // check limits for each local coordinate [-1,1] for quadrilaterals. (different for other elements, typically [0,1]).
    bool inside = true;
    for ( int i = 1; i <= answer.giveSize(); i++ ) {
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


double
FEI3dWedgeLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix;

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant();
}


void
FEI3dWedgeLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
{
    FloatMatrix dNduvw, coords;
    this->evaldNdxi(dNduvw, lcoords, cellgeo);
    coords.resize(3, 6);
    for ( int i = 1; i <= 6; i++ ) {
        coords.setColumn(cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
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
    OOFEM_ERROR("not implemented");
}


void
FEI3dWedgeLin :: edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;

    const auto &nodes = this->computeLocalEdgeMapping(iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

    answer.clear();
    for ( int i = 1; i <= n.giveSize(); ++i ) {
        answer.add( n.at(i), cellgeo.giveVertexCoordinates( nodes.at(i) ) );
    }
}


IntArray
FEI3dWedgeLin :: computeLocalEdgeMapping(int iedge) const
{
    if ( iedge == 1 ) {
        return {1, 2};
    } else if ( iedge == 2 ) {
        return {2, 3};
    } else if ( iedge == 3 ) {
        return {3, 1};
    } else if ( iedge == 4 ) {
        return {4, 5};
    } else if ( iedge == 5 ) {
        return {5, 6};
    } else if ( iedge == 6 ) {
        return {6, 4};
    } else if ( iedge == 7 ) {
        return {1, 4};
    } else if ( iedge == 8 ) {
        return {2, 5};
    } else if ( iedge == 9 ) {
        return {3, 6};
    } else {
        throw std::range_error("invalid edge number");
        return {};
    }
}


double
FEI3dWedgeLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("not implemented");
    return 0.0;
}


void
FEI3dWedgeLin :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    double eta = lcoords.at(2);

    if ( isurf <= 2 ) {
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
    FloatArray n;

    const auto &nodes = this->computeLocalSurfaceMapping(isurf);
    this->surfaceEvalN(n, isurf, lcoords, cellgeo);

    answer.clear();
    for ( int i = 1; i <= n.giveSize(); ++i ) {
        answer.add( n.at(i), cellgeo.giveVertexCoordinates( nodes.at(i) ) );
    }
}


IntArray
FEI3dWedgeLin :: computeLocalSurfaceMapping(int isurf) const
{
    if ( isurf == 1 ) {
        return {1, 2, 3};
    } else if ( isurf == 2 ) {
        return {4, 5, 6};
    } else if ( isurf == 3 ) {
        return {1, 2, 5, 4};
    } else if ( isurf == 4 ) {
        return {2, 3, 6, 5};
    } else if ( isurf == 5 ) {
        return {3, 1, 4, 6};
    } else {
        throw std::range_error("invalid surface number");
        return {};
    }
}


double
FEI3dWedgeLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                   const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("not implemented");
    return 0;
}


std::unique_ptr<IntegrationRule>
FEI3dWedgeLin :: giveIntegrationRule(int order)
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    ///@todo This function below isn't supported for wedges. We must decide how we should do this.
    //int points = iRule->getRequiredNumberOfIntegrationPoints(_Wedge, order);
    OOFEM_WARNING("Warning.. ignoring 'order' argument: FIXME");
    int pointsZeta = 1;
    int pointsTriangle = 1;
    iRule->SetUpPointsOnWedge(pointsTriangle, pointsZeta, _Unknown);
    return std::move(iRule);
}


std::unique_ptr<IntegrationRule>
FEI3dWedgeLin :: giveBoundaryIntegrationRule(int order, int boundary)
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    if ( boundary <= 2 ) {
        int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order + 0);
        iRule->SetUpPointsOnTriangle(points, _Unknown);
    } else {
        int points = iRule->getRequiredNumberOfIntegrationPoints(_Square, order + 2);
        iRule->SetUpPointsOnSquare(points, _Unknown);
    }
    return std::move(iRule);
}
} // end namespace oofem
