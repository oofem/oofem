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

#include "feinterpol2d.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
void FEInterpolation2d :: boundaryEdgeGiveNodes(IntArray &answer, int boundary)
{
    ///@todo What is this?! ...Doesn't seem to be used anyway / JB
    answer.resize(1);
    answer.at(1) = boundary;
}

void FEInterpolation2d :: boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(1);
    answer.at(1) = 1.;
}

double FEInterpolation2d :: boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return 1.;
}

void FEInterpolation2d :: boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(2);
    answer.at(1) = cellgeo.giveVertexCoordinates(boundary)->at(xind);
    answer.at(2) = cellgeo.giveVertexCoordinates(boundary)->at(yind);
}

double FEInterpolation2d :: giveArea(const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented in subclass.");
    return 0;
}

#define POINT_TOL 1.e-3

int FEInterpolation2d :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    FloatArray res, delta, guess, lcoords_guess;
    FloatMatrix jac;
    double convergence_limit, error = 0.0;

    // find a suitable convergence limit
    convergence_limit = 1e-6 * this->giveCharacteristicLength(cellgeo);

    // setup initial guess
    lcoords_guess.resize( 2 );
    lcoords_guess.zero();

    // apply Newton-Raphson to solve the problem
    for ( int nite = 0; nite < 10; nite++ ) {
        // compute the residual
        this->local2global(guess, lcoords_guess, cellgeo);
        res = {gcoords(0) - guess(0), gcoords(1) - guess(1)};

        // check for convergence
        error = res.computeNorm();
        if ( error < convergence_limit ) {
            break;
        }

        // compute the corrections
        this->giveJacobianMatrixAt(jac, lcoords_guess, cellgeo);
        jac.solveForRhs(res, delta);

        // update guess
        lcoords_guess.add(delta);
    }
    if ( error > convergence_limit ) { // Imperfect, could give false negatives.
        OOFEM_WARNING("Failed convergence");
        answer = {1. / 3., 1. / 3.};
        return false;
    }

    answer = { lcoords_guess(0), lcoords_guess(1) };

    return inside(answer);
}

void
FEInterpolation2d :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
{
    double x, y;
    FloatMatrix dn;

    jacobianMatrix.resize(2, 2);
    jacobianMatrix.zero();

    this->evaldNdxi(dn, lcoords, cellgeo );

    for ( int i = 1; i <= dn.giveNumberOfRows(); i++ ) {
        x = cellgeo.giveVertexCoordinates(i)->at(xind);
        y = cellgeo.giveVertexCoordinates(i)->at(yind);

        jacobianMatrix.at(1, 1) += dn.at(i, 1) * x;
        jacobianMatrix.at(2, 1) += dn.at(i, 1) * y;
        jacobianMatrix.at(1, 2) += dn.at(i, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(i, 2) * y;
    }
}

bool FEInterpolation2d ::inside(const FloatArray &lcoords) const
{
	OOFEM_ERROR("Not implemented.")
}

void FEInterpolation2d :: boundaryGiveNodes(IntArray &answer, int boundary)
{
    this->computeLocalEdgeMapping(answer, boundary);
}

void FEInterpolation2d :: boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->edgeEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->edgeEvalNormal(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation2d :: boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->edgeLocal2global(answer, boundary, lcoords, cellgeo);
}

void FEInterpolation2d :: computeEdgeMapping(IntArray &edgeNodes, IntArray &elemNodes, int iedge)
{
    IntArray ln;
    this->computeLocalEdgeMapping(ln, iedge);
    int size = ln.giveSize();
    edgeNodes.resize(size);
    for ( int i = 1; i <= size; i++ ) {
        edgeNodes.at(i) = elemNodes.at( ln.at(i) );
    }
}

double FEInterpolation2d :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray normal;
    return this->edgeEvalNormal(normal, iedge, lcoords, cellgeo);
}

IntegrationRule *FEInterpolation2d :: giveBoundaryIntegrationRule(int order, int boundary)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Line, order + this->order);
    iRule->SetUpPointsOnLine(points, _Unknown);
    return iRule;
}

IntegrationRule *FEInterpolation2d :: giveBoundaryEdgeIntegrationRule(int order, int boundary)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    iRule->SetUpPoint(_Unknown);
    return iRule;
}
} // end namespace oofem
