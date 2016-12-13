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

#include "feinterpol.h"
#include "element.h"
#include "gaussintegrationrule.h"

namespace oofem {
int FEIElementGeometryWrapper :: giveNumberOfVertices() const { return elem->giveNumberOfNodes(); }

double
FEInterpolation :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix;
    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant();
}


IntegrationRule*
FEInterpolation:: giveIntegrationRule(int order)
{
  integrationDomain id = this->giveIntegrationDomain();
  IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);

  int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
  iRule->SetUpPointsOnLine(points, _Unknown);
  return iRule;
}

IntegrationRule*
FEInterpolation::giveBoundaryIntegrationRule(int order, int boundary)
{
  integrationDomain id = this->giveBoundaryIntegrationDomain(boundary);
  IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);

  int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
  iRule->setUpIntegrationPoints(id, points, _Unknown);
  return iRule;
}

IntegrationRule*
FEInterpolation::giveBoundaryEdgeIntegrationRule(int order, int boundary)
{
  integrationDomain id = this->giveBoundaryEdgeIntegrationDomain(boundary);
  IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);

  int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
  iRule->setUpIntegrationPoints(id, points, _Unknown);
  return iRule;
}

IntegrationRule*
FEInterpolation::giveBoundarySurfaceIntegrationRule(int order, int boundary)
{
  integrationDomain id = this->giveBoundarySurfaceIntegrationDomain(boundary);
  IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);

  int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
  iRule->setUpIntegrationPoints(id, points, _Unknown);
  return iRule;
}

  
  
} // end namespace oofem
