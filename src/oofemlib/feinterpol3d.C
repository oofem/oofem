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

#include "feinterpol3d.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
double FEInterpolation3d :: giveVolume(const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented in subclass.");
    return 0;
}

void FEInterpolation3d :: boundaryEdgeGiveNodes(IntArray &answer, int boundary)
{
    this->computeLocalEdgeMapping(answer, boundary);
}

void FEInterpolation3d :: boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->edgeEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation3d :: boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation3d :: boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->edgeLocal2global(answer, boundary, lcoords, cellgeo);
}

void FEInterpolation3d :: boundaryGiveNodes(IntArray &answer, int boundary)
{
    this->computeLocalSurfaceMapping(answer, boundary);
}

void FEInterpolation3d :: boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->surfaceEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation3d :: boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->surfaceEvalNormal(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation3d :: boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->surfaceGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation3d :: boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->surfaceLocal2global(answer, boundary, lcoords, cellgeo);
}

void FEInterpolation3d :: computeEdgeMapping(IntArray &edgeNodes, IntArray &elemNodes, int iedge)
{
    int size;
    IntArray ln;
    this->computeLocalEdgeMapping(ln, iedge);
    size = ln.giveSize();
    edgeNodes.resize(size);
    for ( int i = 1; i <= size; i++ ) {
        edgeNodes.at(i) = elemNodes.at( ln.at(i) );
    }
}

void FEInterpolation3d :: computeSurfaceMapping(IntArray &surfNodes, IntArray &elemNodes, int isurf)
{
    int size;
    IntArray ln;
    this->computeLocalSurfaceMapping(ln, isurf);
    size = ln.giveSize();
    surfNodes.resize(size);
    for ( int i = 1; i <= size; i++ ) {
        surfNodes.at(i) = elemNodes.at( ln.at(i) );
    }
}

IntegrationRule *FEInterpolation3d :: giveBoundaryEdgeIntegrationRule(int order, int boundary)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Line, order + this->order);
    iRule->SetUpPointsOnLine(points, _Unknown);
    return iRule;
}

void FEInterpolation3d :: edgeEvaldNdxi(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not implemented");
}

void FEInterpolation3d :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not implemented");
}

double FEInterpolation3d :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("Not implemented");
    return -1.0;
}
} // end namespace oofem
