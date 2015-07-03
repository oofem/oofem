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

#include "errorestimator.h"
#include "remeshingcrit.h"
#include "inputrecord.h"

namespace oofem {
ErrorEstimator :: ErrorEstimator(int n, Domain *d) : FEMComponent(n, d)
{
    skippedNelems = 0;
    regionSkipMap.clear();
}

ErrorEstimator :: ~ErrorEstimator()
{
}

void
ErrorEstimator :: setDomain(Domain *d)
{
    FEMComponent :: setDomain(d);
    this->giveRemeshingCrit()->setDomain(d);
}



IRResultType
ErrorEstimator :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    regionSkipMap.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, regionSkipMap, _IFT_ErrorEstimator_regionskipmap);
    this->IStype = IST_StressTensor;
    int val = ( int ) this->IStype;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_ErrorEstimator_IStype);
    this->IStype = ( InternalStateType ) val;

    return IRRT_OK;
}

void ErrorEstimator :: reinitialize()
{
    this->rc->reinitialize();
}

bool ErrorEstimator :: skipRegion(int reg)
{
    if ( reg <= regionSkipMap.giveSize() ) {
        return regionSkipMap.at(reg) > 0;
    } else {
        return false;
    }
}
} // end namespace oofem
