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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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


FEIElementDeformedGeometryWrapper::FEIElementDeformedGeometryWrapper(const Element *elem) : FEICellGeometry() {
    this->elem = elem;
    this->alpha = 1;
    this->tStep = NULL;
}

FEIElementDeformedGeometryWrapper::FEIElementDeformedGeometryWrapper(const Element *elem, TimeStep *tStep) : FEICellGeometry() {
    this->elem = elem;
    this->alpha = 1;
    this->tStep = tStep;
}



int
FEIElementDeformedGeometryWrapper::giveNumberOfVertices() const
{
    return elem->giveNumberOfNodes();
}


const FloatArray 
FEIElementDeformedGeometryWrapper::giveVertexCoordinates(int i) const
{
    FloatArray actualCoords = elem->giveNode(i)->giveCoordinates();
    if ( tStep != NULL ) {
        FloatArray u;
        elem->giveNode(i)->giveUnknownVector(u, { D_u, D_v, D_w }, VM_Total, tStep);
        u.times(alpha);
        actualCoords.add(u);
    }
    return actualCoords;
}



  
double
FEInterpolation :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    FloatMatrix jacobianMatrix;
    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant();
}


std::unique_ptr<IntegrationRule>
FEInterpolation:: giveIntegrationRule(int order, Element_Geometry_Type egt) const
{
    integrationDomain id = this->giveIntegrationDomain(egt);
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);

    int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
    iRule->SetUpPointsOnLine(points, _Unknown);
    return std::move(iRule);
}

std::unique_ptr<IntegrationRule>
FEInterpolation::giveBoundaryIntegrationRule(int order, int boundary, Element_Geometry_Type egt) const 
{
    integrationDomain id = this->giveBoundaryIntegrationDomain(boundary, egt);
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);

    int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
    iRule->setUpIntegrationPoints(id, points, _Unknown);
    return std::move(iRule);
}

std::unique_ptr<IntegrationRule>
FEInterpolation::giveBoundaryEdgeIntegrationRule(int order, int boundary, Element_Geometry_Type egt) const 
{
    integrationDomain id = this->giveBoundaryEdgeIntegrationDomain(boundary, egt);
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);

    int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
    iRule->setUpIntegrationPoints(id, points, _Unknown);
    return std::move(iRule);
}

std::unique_ptr<IntegrationRule>
FEInterpolation::giveBoundarySurfaceIntegrationRule(int order, int boundary, Element_Geometry_Type egt) const
{
    integrationDomain id = this->giveBoundarySurfaceIntegrationDomain(boundary, egt);
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);

    int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
    iRule->setUpIntegrationPoints(id, points, _Unknown);
    return std::move(iRule);
}  
  
} // end namespace oofem
