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

#include "combinedzzsiee.h"
#include "domain.h"
#include "element.h"
#include "conTable.h"
#include "mathfem.h"

namespace oofem {
#define CZZSI_ZERO_INDICATOR_TOL 1.e-3

void
CombinedZZSIErrorEstimator :: setDomain(Domain *d) {
    FEMComponent :: setDomain(d);
    zzee.setDomain(d);
    siee.setDomain(d);
    this->giveRemeshingCrit()->setDomain(d);
}


int
CombinedZZSIErrorEstimator :: estimateError(EE_ErrorMode mode, TimeStep *tStep)
{
    int result = zzee.estimateError(mode, tStep);
    result += siee.estimateError(mode, tStep);
    return result;
}

double
CombinedZZSIErrorEstimator :: giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep)
{
    this->estimateError(equilibratedEM, tStep);
    if ( type == indicatorET ) {
        return siee.giveElementError(type, elem, tStep);
    } else {
        return zzee.giveElementError(type, elem, tStep);
    }
}

double
CombinedZZSIErrorEstimator :: giveValue(EE_ValueType type, TimeStep *tStep)
{
    return zzee.giveValue(type, tStep);
}

RemeshingCriteria *
CombinedZZSIErrorEstimator :: giveRemeshingCrit()
{
    if ( this->rc ) {
        return this->rc;
    }

    return ( this->rc = new CombinedZZSIRemeshingCriteria(1, this) );
}

IRResultType
CombinedZZSIErrorEstimator :: initializeFrom(InputRecord *ir)
{
    zzee.initializeFrom(ir);
    siee.initializeFrom(ir);

    return this->giveRemeshingCrit()->initializeFrom(ir);
}




CombinedZZSIRemeshingCriteria :: CombinedZZSIRemeshingCriteria(int n, ErrorEstimator *e) : RemeshingCriteria(n, e),
    zzrc(n, e), dirc(n, e)
{ }

double
CombinedZZSIRemeshingCriteria :: giveRequiredDofManDensity(int num, TimeStep *tStep, int relative)
{
    double indicatorVal, currDensity;
    double proposedDensity;
    this->estimateMeshDensities(tStep);

    dirc.giveNodeChar(num, tStep, indicatorVal, currDensity);
    if ( indicatorVal > dirc.giveMinIndicatorLimit() ) {
        return dirc.giveRequiredDofManDensity(num, tStep, relative);
    } else if ( fabs(indicatorVal) > CZZSI_ZERO_INDICATOR_TOL ) {
        //return zzrc.giveRequiredDofManDensity (num, tStep, relative);
        // transition between dirc and zzrc
        proposedDensity = zzrc.giveRequiredDofManDensity(num, tStep, relative);
        proposedDensity = min(proposedDensity, currDensity);
        proposedDensity = max( proposedDensity, dirc.giveMinIndicatorDensity() );
        if ( relative ) {
            return proposedDensity / currDensity;
        } else {
            return proposedDensity;
        }
    } else {
        return zzrc.giveRequiredDofManDensity(num, tStep, relative);
    }

    return 0.0; // to make compiler happy
}


RemeshingStrategy
CombinedZZSIRemeshingCriteria :: giveRemeshingStrategy(TimeStep *tStep)
{
    RemeshingStrategy s1, s2;
    this->estimateMeshDensities(tStep);

    s1 = zzrc.giveRemeshingStrategy(tStep);
    s2 = dirc.giveRemeshingStrategy(tStep);

    //if ((s1 == RemeshingFromPreviousState_RS) || (s2 == RemeshingFromPreviousState_RS)) return RemeshingFromPreviousState_RS;
    if ( ( s1 == RemeshingFromPreviousState_RS ) || ( s2 == RemeshingFromPreviousState_RS ) ) {
        return RemeshingFromCurrentState_RS;
    } else if ( ( s1 == RemeshingFromCurrentState_RS ) || ( s2 == RemeshingFromCurrentState_RS ) ) {
        return RemeshingFromCurrentState_RS;
    } else {
        return NoRemeshing_RS;
    }
}

int
CombinedZZSIRemeshingCriteria :: estimateMeshDensities(TimeStep *tStep)
{
    zzrc.estimateMeshDensities(tStep);
    dirc.estimateMeshDensities(tStep);
    return 1;
}

IRResultType
CombinedZZSIRemeshingCriteria :: initializeFrom(InputRecord *ir)
{
    zzrc.initializeFrom(ir);
    dirc.initializeFrom(ir);
    return IRRT_OK;
}


double
CombinedZZSIRemeshingCriteria :: giveDofManDensity(int num)
{
    int i, isize;
    ConnectivityTable *ct = domain->giveConnectivityTable();
    const IntArray *con;
    ZZRemeshingCriteriaInterface *interface;
    double density = 0.0;

    con = ct->giveDofManConnectivityArray(num);
    isize = con->giveSize();

    for ( i = 1; i <= isize; i++ ) {
        interface = ( ZZRemeshingCriteriaInterface * )
                    domain->giveElement( con->at(i) )->giveInterface(DirectErrorIndicatorRCInterfaceType);
        if ( !interface ) {
            _error("giveDofManDensity: element does not support ZZRemeshingCriteriaInterface");
        }

        if ( i == 1 ) {
            density = interface->ZZRemeshingCriteriaI_giveCharacteristicSize();
        } else {
            density = min( density, interface->ZZRemeshingCriteriaI_giveCharacteristicSize() );
        }
    }

    return density;
}


void
CombinedZZSIRemeshingCriteria :: setDomain(Domain *d) {
    FEMComponent :: setDomain(d);
    zzrc.setDomain(d);
    dirc.setDomain(d);
}
} // end namespace oofem
