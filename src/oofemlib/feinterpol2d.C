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
IntArray FEInterpolation2d :: boundaryEdgeGiveNodes(int boundary) const
{
    return this->computeLocalEdgeMapping(boundary);
}

void FEInterpolation2d :: boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->edgeEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation2d :: boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->edgeLocal2global(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: giveArea(const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented in subclass.");
    return 0;
}

#define POINT_TOL 1.e-3

int FEInterpolation2d :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo) const
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
FEInterpolation2d :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
{
    FloatMatrix dn;

    jacobianMatrix.resize(2, 2);
    jacobianMatrix.zero();

    this->evaldNdxi(dn, lcoords, cellgeo );

    for ( int i = 1; i <= dn.giveNumberOfRows(); i++ ) {
        double x = cellgeo.giveVertexCoordinates(i).at(xind);
        double y = cellgeo.giveVertexCoordinates(i).at(yind);

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

IntArray FEInterpolation2d :: boundaryGiveNodes(int boundary) const
{
    return this->computeLocalEdgeMapping(boundary);
}

void FEInterpolation2d :: boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->edgeEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return this->edgeEvalNormal(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation2d :: boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    return this->edgeLocal2global(answer, boundary, lcoords, cellgeo);
}

IntArray FEInterpolation2d :: computeEdgeMapping(const IntArray &elemNodes, int iedge) const
{
    const auto& ln = this->computeLocalEdgeMapping(iedge);
    int size = ln.giveSize();
    IntArray edgeNodes(size);
    for ( int i = 1; i <= size; i++ ) {
        edgeNodes.at(i) = elemNodes.at( ln.at(i) );
    }
    return edgeNodes;
}

double FEInterpolation2d :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    FloatArray normal;
    return this->edgeEvalNormal(normal, iedge, lcoords, cellgeo);
}

void FEInterpolation2d::boundarySurfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
  this->evalN(answer, lcoords, cellgeo);
}

void FEInterpolation2d::boundarySurfaceEvaldNdx(FloatMatrix &answer, int isurf,
                            const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
  this->evaldNdx(answer, lcoords, cellgeo);
}

double FEInterpolation2d::boundarySurfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords,
                            const FEICellGeometry &cellgeo) const
{
    answer = {0, 0, 1};
    return this->giveTransformationJacobian(lcoords, cellgeo);
}

void FEInterpolation2d::boundarySurfaceLocal2global(FloatArray &answer, int isurf,
                            const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->local2global(answer, lcoords, cellgeo);
}

double FEInterpolation2d::boundarySurfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                            const FEICellGeometry &cellgeo) const 
{
  return this->giveTransformationJacobian(lcoords, cellgeo);
}

IntArray FEInterpolation2d::boundarySurfaceGiveNodes(int boundary) const
{
    int nnode = this->giveNumberOfNodes();
    IntArray answer(nnode);
    answer.enumerate(nnode);
    return answer;
}

  
} // end namespace oofem
