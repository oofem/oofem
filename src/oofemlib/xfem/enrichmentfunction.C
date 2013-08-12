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
REGISTER_EnrichmentFunction(DiscontinuousFunction)
REGISTER_EnrichmentFunction(BranchFunction)
REGISTER_EnrichmentFunction(RampFunction)

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
    if ( EnrichmentDomain_BG * edbg = dynamic_cast< EnrichmentDomain_BG * >(ed) ) {
        return sgn( edbg->bg->computeDistanceTo(point) );
    } else {
        OOFEM_ERROR("DiscontinuousFunction :: evaluateFunctionAt - only supports enrichment domains of type Basic Geometry");
        return 0.;
    }
}

void DiscontinuousFunction :: evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, const EnrichmentDomain *ipEnrDom) const
{
    oEnrFunc = sgn(iLevelSet);
}

void DiscontinuousFunction :: evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, const EnrichmentDomain *ipEnrDom) const
{
    oEnrFuncDeriv.resize(2);
    oEnrFuncDeriv.zero();
}

void DiscontinuousFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentDomain *ed)
{
    answer.resize(2);
    answer.zero();
}

double RampFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentDomain *ed)
{
    if ( EnrichmentDomain_BG * edbg = dynamic_cast< EnrichmentDomain_BG * >(ed) ) {
        return fabs( edbg->bg->computeDistanceTo(point) );
    } else {
        OOFEM_ERROR("RampFunction :: evaluateFunctionAt - only supports enrichment domains of type Basic Geometry");
        return 0.;
    }
}

void RampFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentDomain *ed)
{
    if ( EnrichmentDomain_BG * edbg = dynamic_cast< EnrichmentDomain_BG * >(ed) ) {
        double dist = edbg->bg->computeDistanceTo(point);
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
    el->giveInterpolation()->evalN( N, * gp->giveCoordinates(), FEIElementGeometryWrapper(el) );
    double dist = 0;
    double member = 0;
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        if ( EnrichmentDomain_BG * edbg = dynamic_cast< EnrichmentDomain_BG * >(ed) ) {
            ;
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
    el->giveInterpolation()->evalN( N, * gp->giveCoordinates(), FEIElementGeometryWrapper(el) );
    IntArray dofManArray( el->giveNumberOfDofManagers() );
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        dofManArray.at(i) = el->giveDofManagerNumber(i);
    }

    FloatMatrix dNdx;
    el->giveInterpolation()->evaldNdx( dNdx, * gp->giveCoordinates(), FEIElementGeometryWrapper(el) );
    double dist = 0;
    double dfdx = 0;
    double dfdy = 0;
    double phi = 0;
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        if ( EnrichmentDomain_BG * edbg = dynamic_cast< EnrichmentDomain_BG * >(ed) ) {
            ;
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

void RampFunction :: evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, const EnrichmentDomain *ipEnrDom) const
{
    oEnrFunc = fabs(iLevelSet);
}

void RampFunction :: evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, const EnrichmentDomain *ipEnrDom) const
{
    oEnrFuncDeriv.resize(2);
    oEnrFuncDeriv.zero();
    oEnrFuncDeriv.at(1) = iGradLevelSet.at(1) * sgn(iLevelSet);
    oEnrFuncDeriv.at(2) = iGradLevelSet.at(2) * sgn(iLevelSet);
}

double BranchFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentDomain *ed)
{
    ///@todo change
#if 0
    double ret = 0;
    CrackTip *cr = ( CrackTip * ) ed;
    FloatArray transfCoor;
    Line *l = ( Line * ) cr->giveGeometry();
    l->transformIntoPolar(point, transfCoor);
    double rS = sqrt( transfCoor.at(1) );
    double theta = transfCoor.at(2);
    switch ( number ) {
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
    CrackTip *cr = ( CrackTip * ) ed;
    Line *l = ( Line * ) cr->giveGeometry();
    FloatArray transfCoor;
    l->transformIntoPolar(point, transfCoor);
    double alpha = l->computeInclinationAngle();
    double fac = 0.5 / sqrt( transfCoor.at(1) );
    double theta = transfCoor.at(2);
    double dPhidr = 0;
    double dPhidt = 0;
    double drdx = cos(alpha);
    double drdy = sin(alpha);
    double dtdx = ( -1 ) * sin(alpha);
    double dtdy = cos(alpha);

    switch ( number ) {
    case 1:
        dPhidr = ( -1 ) * fac * sin(0.5 * theta);
        dPhidt = fac * cos(0.5 * theta);
        break;
    case 2:
        dPhidr = fac * cos(0.5 * theta);
        dPhidt = fac * sin(0.5 * theta);
        break;
    case 3:
        dPhidr = fac * sin(0.5 * theta) * ( 2 * sin(theta) * sin(theta) - cos(theta) );
        dPhidt = fac * cos(theta) * cos(1.5 * theta);
        break;
    case 4:
        dPhidr = fac * cos(0.5 * theta) * ( cos(theta) + 2 * sin(theta) * sin(theta) );
        dPhidt = ( -1 ) * fac * cos(theta) * sin(1.5 * theta);
        break;
    }

    answer.at(1) = dPhidr * drdx + dPhidt * dtdx;
    answer.at(2) = dPhidr * drdy + dPhidt * dtdy;
#endif
}

void LinElBranchFunction :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const double &iR, const double &iTheta) const
{
    oEnrFunc.resize(4, 0.0);
    oEnrFunc [ 0 ] = sqrt(iR) * sin(0.5 * iTheta);
    oEnrFunc [ 1 ] = sqrt(iR) * sin(0.5 * iTheta) * sin(iTheta);
    oEnrFunc [ 2 ] = sqrt(iR) * cos(0.5 * iTheta);
    oEnrFunc [ 3 ] = sqrt(iR) * cos(0.5 * iTheta) * sin(iTheta);
}

void LinElBranchFunction :: evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const double &iR, const double &iTheta) const
{
    // Evaluate the enrichment function derivatives using the chain rule:
    // dPdx = dPdr*drdx + dPdt*dtdx
    // dPdy = dPdr*drdy + dPdt*dtdy

    oEnrFuncDeriv.resize(4);

    const double drdx =  cos(iTheta);
    const double drdy =  sin(iTheta);
    const double dtdx = -( 1.0 / iR ) * sin(iTheta);
    const double dtdy =  ( 1.0 / iR ) * cos(iTheta);
    /*
     *      double dtdy = 0.0;
     *      double eps = 1.0e-12;
     *      if( fabs(cos(iTheta)) < eps )
     *      {
     *              dtdy = 0.0;
     *      }
     *      else
     *      {
     *              dtdy = 1.0/(iR*iR);
     *      }
     */
    // Psi 1
    const double dP1dr = ( 1.0 / ( 2.0 * sqrt(iR) ) ) * sin(0.5 * iTheta);
    const double dP1dt = 0.5 * sqrt(iR) * cos(0.5 * iTheta);
    oEnrFuncDeriv [ 0 ].setValues(2, dP1dr * drdx + dP1dt * dtdx, dP1dr * drdy + dP1dt * dtdy);

    // Psi 2
    const double dP2dr = ( 1.0 / ( 2.0 * sqrt(iR) ) ) * sin(0.5 * iTheta) * sin(iTheta);
    const double dP2dt = 0.5 * sqrt(iR) * cos(0.5 * iTheta) * sin(iTheta) + sqrt(iR) * sin(0.5 * iTheta) * cos(iTheta);
    oEnrFuncDeriv [ 1 ].setValues(2, dP2dr * drdx + dP2dt * dtdx, dP2dr * drdy + dP2dt * dtdy);

    // Psi 3
    const double dP3dr = ( 1.0 / ( 2.0 * sqrt(iR) ) ) * cos(0.5 * iTheta);
    const double dP3dt = -0.5 * sqrt(iR) * sin(0.5 * iTheta);
    oEnrFuncDeriv [ 2 ].setValues(2, dP3dr * drdx + dP3dt * dtdx, dP3dr * drdy + dP3dt * dtdy);

    // Psi 4
    const double dP4dr = ( 1.0 / ( 2.0 * sqrt(iR) ) ) * cos(0.5 * iTheta) * sin(iTheta);
    const double dP4dt = -0.5 * sqrt(iR) * sin(0.5 * iTheta) * sin(iTheta) + sqrt(iR) * cos(0.5 * iTheta) * cos(iTheta);
    oEnrFuncDeriv [ 3 ].setValues(2, dP4dr * drdx + dP4dt * dtdx, dP4dr * drdy + dP4dt * dtdy);
}
} // end namespace oofem
