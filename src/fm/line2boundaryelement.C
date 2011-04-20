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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
#include "timestep.h"
#include "mathfem.h"

#include <stdlib.h>
#include <math.h>

namespace oofem {
Line2BoundaryElement :: Line2BoundaryElement(int n, Domain *aDomain) : FMElement(n, aDomain)
{
    this->numberOfDofMans = 3;
    this->numberOfIntegrationRules = 0;
}

Line2BoundaryElement :: ~Line2BoundaryElement() {}

IRResultType Line2BoundaryElement :: initializeFrom(InputRecord *ir)
{
    return FMElement :: initializeFrom(ir);
}

void Line2BoundaryElement :: computeN(FloatArray &answer, const FloatArray &lcoords) const
{
    double xi = lcoords(0);
    answer.resize(3);
    answer(0) = 0.5*(xi-1.0)*xi;
    answer(1) = 0.5*(xi+1.0)*xi;
    answer(2) = 1.0-xi*xi;
}

int Line2BoundaryElement :: computeLocalCoordinates(FloatArray &lcoords, const FloatArray &gcoords)
{
    FloatArray x1_x2, p_x3, x3_x2_x1;
    x1_x2.beDifferenceOf(*this->giveNode(1)->giveCoordinates(), *this->giveNode(2)->giveCoordinates());
    p_x3.beDifferenceOf(gcoords, *this->giveNode(3)->giveCoordinates());
    x3_x2_x1.beDifferenceOf(*this->giveNode(3)->giveCoordinates(), *this->giveNode(1)->giveCoordinates());
    x3_x2_x1.add(*this->giveNode(3)->giveCoordinates());
    x3_x2_x1.subtract(*this->giveNode(2)->giveCoordinates());
    //x1_x2 = *this->giveNode(1)->giveCoordinates();
    //x1_x2.subtract(*this->giveNode(2)->giveCoordinates());
    //p_x3 = gcoords;
    //p_x3.subtract(*this->giveNode(3)->giveCoordinates());
    //x3_x2_x1 = *this->giveNode(3)->giveCoordinates();
    //x3_x2_x1.add(*this->giveNode(3)->giveCoordinates());
    //x3_x2_x1.subtract(*this->giveNode(2)->giveCoordinates());
    //x3_x2_x1.subtract(*this->giveNode(1)->giveCoordinates());

    double b0 = 0.50*x1_x2.dotProduct(p_x3); // 1/2*(x1 - x2).(p - x3)
    double b1 = 0.25*x1_x2.computeSquaredNorm() + x3_x2_x1.dotProduct(p_x3); // 1/4*(x1 - x2)^2 + (2*x3 - x2 - x1).(p - x3)
    double b2 = 0.75*x1_x2.dotProduct(x3_x2_x1); // 3/4*(x1-x2).(2x3 - x2 - x1)
    double b3 = 0.50*x3_x2_x1.computeSquaredNorm(); // 1/2*(2*x3 - x2 - x3)^2

    double r[3];
    int roots;
    cubic(b3, b2, b1, b0, &r[0], &r[1], &r[2], &roots);

    // Copy them over, along with boundary cases.
    int points = 2;
    double p[5] = {-1.0, 1.0};
    for (int i = 0; i < roots; i++) {
        if (r[i] > -1.0 && r[i] < 1.0) {
            // The cubic solver has pretty bad accuracy, performing a single newton iteration
            // which typically improves the solution by many many orders of magnitude.
            // TODO: Move this to the cubic solver.
            r[i] -= (b0 + b1*r[i] + b2*r[i]*r[i] + b3*r[i]*r[i]*r[i])/(b1 + 2*b2*r[i] + 3*b3*r[i]*r[i]);
            //printf("r[%d] -> %e\n",i, (b0 + b1*r[i] + b2*r[i]*r[i] + b3*r[i]*r[i]*r[i]));
            p[points] = r[i];
            points++;
        }
    }

    double min_distance2 = 0.0, min_xi, distance2;
    FloatArray f(2), xi(1);

    for (int i = 0; i < points; i++) {
        xi(0) = p[i];
        this->computeGlobalCoordinates(f,xi);
        distance2 = f.distance_square(gcoords);
        if ( i == 0 || distance2 < min_distance2 ) {
            min_distance2 = distance2;
            min_xi = xi(0);
        }
    }

    lcoords.resize(1);
    lcoords(0) = min_xi;

    return true;
}

int Line2BoundaryElement :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FloatArray n;
    this->computeN(n, lcoords);
    answer.resize(2);
    answer.zero();
    for (int i = 1; i <= n.giveSize(); ++i) {
        answer.add(n.at(i), *this->giveNode(i)->giveCoordinates());
    }
    return true;
}

double Line2BoundaryElement :: computeNXIntegral() const
{
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
    int ok = this->SpatialLocalizerI_giveClosestPoint(lcoords, closest, gcoords);
    if (!ok) {
        answer.resize(0);
        return false;
    }
    this->EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(mode, tStep, lcoords, answer);
    return true;
}

void Line2BoundaryElement :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    // TODO: Generalize this
    answer.resize(2);
    answer(0) = V_u;
    answer(1) = V_v;
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
