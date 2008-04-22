/* $Header: /home/cvs/bp/oofem/sm/src/scalarerrorindicator.C,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


#include "scalarerrorindicator.h"
#include "directerrorindicatorrc.h"
#include "element.h"
#include "integrationrule.h"
#include "mathfem.h"


int
ScalarErrorIndicator :: estimateError(EE_ErrorMode mode, TimeStep *tStep) {
    if ( indicatorType == 1 ) {
        if ( mode == equilibratedEM ) {
            varType = IST_PrincipalDamageTensor;
        } else {
            varType = IST_PrincipalDamageTempTensor;
        }
    }

    return 1;
}

double
ScalarErrorIndicator :: giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep) {
    FloatArray val;
    IntegrationRule *iRule = elem->giveDefaultIntegrationRulePtr();
    int i, result = 1, nip = iRule->getNumberOfIntegrationPoints();
    double sval, maxVal = 0.0;

    if ( type != indicatorET ) {
        return 0.0;
    }

    if ( this->skipRegion( elem->giveRegionNumber() ) ) {
        return 0.0;
    }

    for ( i = 0; i < nip; i++ ) {
        result = elem->giveIPValue(val, iRule->getIntegrationPoint(i), varType, tStep);
        if ( result ) {
            sval = sqrt( dotProduct( val, val, val.giveSize() ) );
            if ( i == 0 ) {
                maxVal = sval;
            } else {
                maxVal = max(maxVal, sval);
            }
        }
    }

    return maxVal;
}


IRResultType
ScalarErrorIndicator :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    ErrorEstimator :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, indicatorType, IFT_ScalarErrorIndicator_vartype, "vartype"); // Macro
    if ( indicatorType != 1 ) {
        _error("instanciateFrom: usupported varType");
    }

    return this->giveRemeshingCrit()->initializeFrom(ir);
}

RemeshingCriteria *
ScalarErrorIndicator :: giveRemeshingCrit() {
    if ( this->rc ) {
        return this->rc;
    }

    return ( this->rc = new DirectErrorIndicatorRC(1, this) );
}
