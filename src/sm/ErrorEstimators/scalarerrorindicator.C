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

#include "../sm/ErrorEstimators/scalarerrorindicator.h"
#include "../sm/ErrorEstimators/directerrorindicatorrc.h"
#include "element.h"
#include "integrationrule.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "errorestimatortype.h"
#include "classfactory.h"

namespace oofem {
REGISTER_ErrorEstimator(ScalarErrorIndicator, EET_SEI);

int
ScalarErrorIndicator :: estimateError(EE_ErrorMode mode, TimeStep *tStep)
{
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
ScalarErrorIndicator :: giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep)
{
    FloatArray val;
    int result = 1;
    double sval, maxVal = 0.0;

    if ( type != indicatorET ) {
        return 0.0;
    }

    if ( this->skipRegion( elem->giveRegionNumber() ) ) {
        return 0.0;
    }

    for ( GaussPoint *gp: *elem->giveDefaultIntegrationRulePtr() ) {
        result = elem->giveIPValue(val, gp, varType, tStep);
        if ( result ) {
            sval = val.computeNorm();
            if ( gp->giveNumber() == 1 ) {
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
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    ErrorEstimator :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, indicatorType, _IFT_ScalarErrorIndicator_vartype);
    if ( indicatorType != 1 ) {
        OOFEM_ERROR("usupported varType");
    }

    return this->giveRemeshingCrit()->initializeFrom(ir);
}

RemeshingCriteria *
ScalarErrorIndicator :: giveRemeshingCrit()
{
    if ( !this->rc ) {
        this->rc.reset( new DirectErrorIndicatorRC(1, this) );
    }

    return this->rc.get();
}
} // end namespace oofem
