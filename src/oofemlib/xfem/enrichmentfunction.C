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
#include "gausspoint.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "enrichmentdomain.h"
#include "classfactory.h"

namespace oofem {

//REGISTER_EnrichmentFunction( DiscontinuousFunction )
//REGISTER_EnrichmentFunction( BranchFunction )
//REGISTER_EnrichmentFunction( RampFunction )

IRResultType EnrichmentFunction :: initializeFrom(InputRecord *ir)
{
    return IRRT_OK;
}

double EnrichmentFunction :: evaluateFunctionAt(GaussPoint *gp, EnrichmentDomain *ed)
{
    FloatArray gcoords;
    gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
    return this->evaluateFunctionAt(& gcoords, ed);
}

void EnrichmentFunction :: evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp, EnrichmentDomain *ed)
{
    FloatArray gc;
    gp->giveElement()->computeGlobalCoordinates( gc, * gp->giveCoordinates() );
    this->evaluateDerivativeAt(answer, & gc, ed);
}


double DiscontinuousFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentDomain *ed)
{
    if ( EnrichmentDomain_BG *edbg = dynamic_cast< EnrichmentDomain_BG * > (ed) ) {  
        return sgn( edbg->bg->computeDistanceTo( point ) );
    } else {
        OOFEM_ERROR("DiscontinuousFunction :: evaluateFunctionAt - only supports enrichment domains of type Basic Geometry");
        return 0.;
    }
}

void DiscontinuousFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentDomain *ed)
{
    answer.resize(2);
    answer.zero();
}

double RampFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentDomain *ed)
{
    if ( EnrichmentDomain_BG *edbg = dynamic_cast< EnrichmentDomain_BG * > (ed) ) {  
        return fabs( edbg->bg->computeDistanceTo( point ) );
    } else {
        OOFEM_ERROR("RampFunction :: evaluateFunctionAt - only supports enrichment domains of type Basic Geometry");
        return 0.;
    }

}

void RampFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentDomain *ed)
{
    if ( EnrichmentDomain_BG *edbg = dynamic_cast< EnrichmentDomain_BG * > (ed) ) {  
        double dist = edbg->bg->computeDistanceTo( point );
        answer.resize(2);
        answer.zero();
        answer.at(1) = answer.at(2) = sgn(dist);
    } else {
        OOFEM_ERROR("RampFunction :: evaluateDerivativeAt - only supports enrichment domains of type Basic Geometry");
    }
}

double RampFunction :: evaluateFunctionAt(GaussPoint *gp, EnrichmentDomain *ed)
{
    FloatArray N;
    Element *el = gp->giveElement();
    el->giveInterpolation()->evalN(N, * gp->giveCoordinates(), FEIElementGeometryWrapper(el));
    double dist = 0;
    double member = 0;
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        if ( EnrichmentDomain_BG *edbg = dynamic_cast< EnrichmentDomain_BG * > (ed) ) {;  
            dist = edbg->bg->computeDistanceTo( el->giveDofManager(i)->giveCoordinates() );
            member += N.at(i) * dist;
        } else {
            OOFEM_ERROR("RampFunction :: evaluateFunctionAt - only supports enrichment domains of type Basic Geometry");
        }
    }

    return fabs(member);
}

void RampFunction :: evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp, EnrichmentDomain *ed)
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
        if ( EnrichmentDomain_BG *edbg = dynamic_cast< EnrichmentDomain_BG * > (ed) ) {;  
            dist = edbg->bg->computeDistanceTo( el->giveDofManager(i)->giveCoordinates() );
            phi += N.at(i) * dist;
            dfdx += dNdx.at(i, 1) * dist;
            dfdy += dNdx.at(i, 2) * dist;
        } else {
            OOFEM_ERROR("RampFunction :: evaluateFunctionAt - only supports enrichment domains of type Basic Geometry");
        }
    }

    answer.resize(2);
    answer.zero();
    answer.at(1) = dfdx * sgn(phi);
    answer.at(2) = dfdy * sgn(phi);
}

double BranchFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentDomain *ed)
{
    ///@todo change
#if 0
    double ret = 0;
    CrackTip *cr = (CrackTip*) ed;
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

void BranchFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentDomain *ed)
{
    ///@todo change
#if 0
    answer.resize(2);
    CrackTip *cr = (CrackTip*) ed;
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
