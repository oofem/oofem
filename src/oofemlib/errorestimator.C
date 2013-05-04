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

#include "errorestimator.h"
#include "remeshingcrit.h"

namespace oofem {

ErrorEstimator :: ErrorEstimator ( int n, Domain* d ) : FEMComponent ( n, d )
{
    rc = NULL;
    skippedNelems = 0;
    regionSkipMap.resize ( 0 );
}

ErrorEstimator :: ~ErrorEstimator()
{
    delete rc;
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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    regionSkipMap.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, regionSkipMap, _IFT_ErrorEstimator_regionskipmap);
    this->IStype = IST_StressTensor;
    int val = (int) this->IStype;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_ErrorEstimator_IStype);
    this->IStype = (InternalStateType) val;

    return IRRT_OK;
}

void ErrorEstimator :: reinitialize()
{
    this->rc->reinitialize();
}

bool ErrorEstimator :: skipRegion ( int reg )
{
    if ( reg <= regionSkipMap.giveSize() ) {
        return regionSkipMap.at ( reg ) > 0;
    } else {
        return false;
    }
}

} // end namespace oofem
