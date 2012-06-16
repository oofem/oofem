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

#include "enrichmentfunction.h"
#include "enrichmentitem.h"
#include "domain.h"
#include "gausspnt.h"
#include "xfemmanager.h"
#include "mathfem.h"
#include "geometry.h"
#include "feinterpol.h"

namespace oofem {
IRResultType EnrichmentFunction :: initializeFrom(InputRecord *ir)
{
    return IRRT_OK;
}

double EnrichmentFunction :: evaluateFunctionAt(GaussPoint *gp, EnrichmentItem *ei)
{
    FloatArray gcoords;
    gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
    return this->evaluateFunctionAt(& gcoords, ei);
}

void EnrichmentFunction :: evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp, EnrichmentItem *ei)
{
    FloatArray gc;
    gp->giveElement()->computeGlobalCoordinates( gc, * gp->giveCoordinates() );
    this->evaluateDerivativeAt(answer, & gc, ei);
}


double DiscontinuousFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentItem *ei)
{
    double dist = ei->giveGeometry()->computeDistanceTo(point);
    return sgn(dist);
}

void DiscontinuousFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem *ei)
{
    answer.resize(2);
    answer.zero();
}

double RampFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentItem *ei)
{
    return fabs( ei->giveGeometry()->computeDistanceTo(point) );
}

void RampFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem *ei)
{
    double dist = ei->giveGeometry()->computeDistanceTo(point);

    answer.resize(2);
    answer.zero();
    answer.at(1) = answer.at(2) = sgn(dist);
}

double RampFunction :: evaluateFunctionAt(GaussPoint *gp, EnrichmentItem *ei)
{
    FloatArray N;
    Element *el = gp->giveElement();
    el->giveInterpolation()->evalN(N, * gp->giveCoordinates(), FEIElementGeometryWrapper(el));
    double dist = 0;
    double absMember = 0;
    double member = 0;
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        dist = ei->giveGeometry()->computeDistanceTo( el->giveDofManager(i)->giveCoordinates() );
        member += N.at(i) * dist;
    }

    return fabs(member);
}

void RampFunction :: evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp, EnrichmentItem *ei)
{
    FloatArray N;
    Element *el = gp->giveElement();
    el->giveInterpolation()->evalN(N, * gp->giveCoordinates(), FEIElementGeometryWrapper(el));
    IntArray dofManArray( el->giveNumberOfDofManagers() );
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        dofManArray.at(i) = el->giveDofManagerNumber(i);
    }

    FloatMatrix dNdx;
    el->giveInterpolation()->evaldNdx(dNdx, * gp->giveCoordinates(), FEIElementGeometryWrapper(el));
    double dist = 0;
    double dfdx = 0;
    double dfdy = 0;
    double phi = 0;
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        dist = ei->giveGeometry()->computeDistanceTo( el->giveDofManager(i)->giveCoordinates() );
        phi += N.at(i) * dist;
        dfdx += dNdx.at(i, 1) * dist;
        dfdy += dNdx.at(i, 2) * dist;
    }

    answer.resize(2);
    answer.zero();
    answer.at(1) = dfdx * sgn(phi);
    answer.at(2) = dfdy * sgn(phi);
}

double BranchFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentItem *ei)
{
    ///@todo change
#if 0
    double ret = 0;
    CrackTip *cr = (CrackTip*) ei;
    FloatArray transfCoor;
    Line *l = (Line*) cr->giveGeometry();
    l->transformIntoPolar(point, transfCoor);
    double rS = sqrt(transfCoor.at(1));
    double theta = transfCoor.at(2);
    switch (number) {
        case 1: ret = rS * sin(0.5 * theta);
        case 2: ret = rS * cos(0.5 * theta);
        case 3: ret = rS * sin(0.5 * theta) * cos(theta);
        case 4: ret = rS * cos(0.5 * theta) * cos(theta);
    }
    return ret;
#endif
    return 0.0;
}

void BranchFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem *ei)
{
    ///@todo change
#if 0
    answer.resize(2);
    CrackTip *cr = (CrackTip*) ei;
    Line *l = (Line*) cr->giveGeometry();
    FloatArray transfCoor;
    l->transformIntoPolar(point, transfCoor);
    double alpha = l->computeInclinationAngle();
    double fac = 0.5 / sqrt(transfCoor.at(1));
    double theta = transfCoor.at(2);
    double dPhidr = 0;
    double dPhidt = 0;
    double drdx = cos(alpha);
    double drdy = sin(alpha);
    double dtdx = (-1) * sin(alpha);
    double dtdy = cos(alpha);

    switch (number) {
        case 1:
            dPhidr = (-1) * fac * sin(0.5 * theta);
            dPhidt = fac * cos(0.5 * theta);
            break;
        case 2:
            dPhidr = fac * cos(0.5 * theta);
            dPhidt = fac * sin(0.5 * theta);
            break;
        case 3:
            dPhidr = fac * sin(0.5 * theta)*(2 * sin(theta) * sin(theta) - cos(theta));
            dPhidt = fac * cos(theta) * cos(1.5 * theta);
            break;
        case 4:
            dPhidr = fac * cos(0.5 * theta)*(cos(theta) + 2 * sin(theta) * sin(theta));
            dPhidt = (-1) * fac * cos(theta) * sin(1.5 * theta);
            break;
    }
    answer.at(1) = dPhidr * drdx + dPhidt*dtdx;
    answer.at(2) = dPhidr * drdy + dPhidt*dtdy;
#endif
}
} // end namespace oofem
