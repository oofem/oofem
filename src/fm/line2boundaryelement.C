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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "line2boundaryelement.h"
#include "node.h"
#include "fei2dlinequad.h"
#include "gaussintegrationrule.h"

namespace oofem {
FEI2dLineQuad Line2BoundaryElement :: fei(1, 2);

Line2BoundaryElement :: Line2BoundaryElement(int n, Domain *aDomain) : FMElement(n, aDomain)
{
    this->numberOfDofMans = 3;
    this->numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ 1 ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 2, _Unknown);
}

Line2BoundaryElement :: ~Line2BoundaryElement()
{
}

IRResultType Line2BoundaryElement :: initializeFrom(InputRecord *ir)
{
    return FMElement :: initializeFrom(ir);
}

void Line2BoundaryElement :: computeN(FloatArray &answer, const FloatArray &lcoords) const
{
    this->fei.evalN(answer, lcoords, FEIElementGeometryWrapper(this));
}

int Line2BoundaryElement :: computeLocalCoordinates(FloatArray &lcoords, const FloatArray &gcoords)
{
    this->fei.global2local(lcoords, gcoords, FEIElementGeometryWrapper(this));
    return true;
}

int Line2BoundaryElement :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    this->fei.local2global(answer, lcoords, FEIElementGeometryWrapper(this));
    return true;
}

FEInterpolation * Line2BoundaryElement :: giveInterpolation()
{
    return &this->fei;
}

double Line2BoundaryElement :: computeNXIntegral() const
{
    ///@todo Use the FEI classes for this
    //return this->fei.evalNXIntegral(FEIElementGeometryWrapper(this));
    Node *node;
    double x1, x2, x3, y1, y2, y3;

    node = this->giveNode(1);
    x1 = node->giveCoordinate(1);
    y1 = node->giveCoordinate(2);

    node = this->giveNode(2);
    x2 = node->giveCoordinate(1);
    y2 = node->giveCoordinate(2);

    node = this->giveNode(3);
    x3 = node->giveCoordinate(1);
    y3 = node->giveCoordinate(2);

    return (x1*y2 - x2*y1 + 4*(x3*(y1 - y2) + y3*(x2 - x1)))/3.0;
}

void Line2BoundaryElement :: giveDofManDofIDMask(int i, EquationID eid, IntArray &nodeDofIDMask) const
{
    if (eid == EID_MomentumBalance) {
        nodeDofIDMask.setValues(2, V_u, V_v);
    }
    else {
        nodeDofIDMask.resize(0);
    }
}

Interface *Line2BoundaryElement :: giveInterface(InterfaceType it)
{
    switch ( it ) {
        case SpatialLocalizerInterfaceType:
            return static_cast< SpatialLocalizerInterface * >(this);
        case EIPrimaryUnknownMapperInterfaceType:
            return static_cast< EIPrimaryUnknownMapperInterface * >(this);
        default:
            return FMElement :: giveInterface(it);
    }
}

double Line2BoundaryElement :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray c;
    c = *this->giveNode(1)->giveCoordinates();
    c.add(*this->giveNode(2)->giveCoordinates());
    c.times(0.5);
    return c.distance(coords);
}

double Line2BoundaryElement :: SpatialLocalizerI_giveClosestPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &gcoords)
{
    if (!this->computeLocalCoordinates(lcoords, gcoords)) {
        lcoords.resize(0);
        closest.resize(0);
        return -1.0;
    }
    // compute local coordinates already gives closest point.
    this->computeGlobalCoordinates(closest, lcoords);
    return closest.distance(gcoords);
}

int Line2BoundaryElement :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
        TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
{
    FloatArray lcoords, closest;
    double distance = this->SpatialLocalizerI_giveClosestPoint(lcoords, closest, gcoords);
    if (distance < 0) {
        answer.resize(0);
        return false;
    }
    this->EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(mode, tStep, lcoords, answer);
    return true;
}

void Line2BoundaryElement :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    ///@todo Generalize this
    answer.setValues(2, V_u, V_v);
}

void Line2BoundaryElement :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
        TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n;
    this->computeN(n, lcoords);

    IntArray dofIDs;
    this->EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(dofIDs);

    answer.resize(dofIDs.giveSize());
    answer.zero();
    for (int i = 1; i <= n.giveSize(); ++i) {
        for (int j = 1; j <= dofIDs.giveSize(); ++j) {
            answer.at(j) += n.at(i)*this->giveNode(i)->giveDofWithID(dofIDs.at(j))->giveUnknown(EID_MomentumBalance, mode, tStep);
        }
    }
}

} // end namespace oofem
