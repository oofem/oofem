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

#include "feinterpol1d.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"
#include <stdexcept>

namespace oofem {


void FEInterpolation1d :: boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    answer.resize(1);
    answer.at(1) = 1.;
}

IntArray FEInterpolation1d :: boundaryGiveNodes(int boundary) const
{
    throw std::runtime_error("Not implemented");
}

double FEInterpolation1d :: boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented");
    return 1.;
}

double FEInterpolation1d :: boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    return 1.;
}

void FEInterpolation1d :: boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    answer = cellgeo.giveVertexCoordinates(boundary);
}

std::unique_ptr<IntegrationRule> FEInterpolation1d :: giveIntegrationRule(int order) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Line, order + this->order);
    iRule->SetUpPointsOnLine(points, _Unknown);
    return std::move(iRule);
}

std::unique_ptr<IntegrationRule> FEInterpolation1d :: giveBoundaryIntegrationRule(int order, int boundary)  const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    iRule->SetUpPoint(_Unknown);
    return std::move(iRule);
}

std::unique_ptr<IntegrationRule> FEInterpolation1d :: giveBoundaryEdgeIntegrationRule(int order, int boundary) const
{
    return this->giveIntegrationRule(order);
}
} // end namespace oofem
