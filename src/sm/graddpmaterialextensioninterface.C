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



#include "domain.h"
#include "nonlocalbarrier.h"
#include "graddpmaterialextensioninterface.h"

#ifdef __PARALLEL_MODE
 #include "parallel.h"
#endif

#include <list>


namespace oofem {
// flag forcing the inclusion of all elements with volume inside support of weight function.
// This forces inclusion of all integration points of these elements, even if weight is zero
// If not defined (default) only integration points with nonzero weight are included.
// #define NMEI_USE_ALL_ELEMENTS_IN_SUPPORT


// constructor
GradDpMaterialExtensionInterface :: GradDpMaterialExtensionInterface(Domain *d)  : Interface()
{
    dom = d;

    cl = 0.;
    cl0 = 0.;
    averType = 0;
    beta = 0.;
    zeta = 0.;
}




IRResultType
GradDpMaterialExtensionInterface :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro


    // read the characteristic length
    IR_GIVE_FIELD(ir, cl, _IFT_GradDpMaterialExtensionInterface_cl);
    if ( cl < 0.0 ) {
        cl = 0.0;
    }

    cl0 = cl;
    // special averaging
    // averType = 0 ... classical
    // averType = 1 ... distance-based
    // averType = 1 ... stress-based
    averType = 0;

    IR_GIVE_OPTIONAL_FIELD(ir, averType, _IFT_GradDpMaterialExtensionInterface_averagingtype);
    if ( averType == 1 ) {
        IR_GIVE_FIELD(ir, beta, _IFT_GradDpMaterialExtensionInterface_beta);
        IR_GIVE_FIELD(ir, zeta, _IFT_GradDpMaterialExtensionInterface_zeta);
    } else if ( averType == 2 ) {
        IR_GIVE_FIELD(ir, beta, _IFT_GradDpMaterialExtensionInterface_beta);
    }

    return IRRT_OK;
}





void
GradDpMaterialExtensionInterface :: giveDistanceBasedCharacteristicLength(const FloatArray &gpCoords)
{
    double distance = zeta * cl0; // Initially distance from the boundary is set to the maximum value
    double temp;
    int ib, nbarrier = dom->giveNumberOfNonlocalBarriers();
    for ( ib = 1; ib <= nbarrier; ib++ ) { //Loop over all the nonlocal barriers to find minimum distance from the boundary
        temp = dom->giveNonlocalBarrier(ib)->calculateMinimumDistanceFromBoundary(gpCoords, zeta * cl0);
        if ( distance > temp ) { //Check to find minimum distance from boundary from all nonlocal boundaries
            distance = temp;
        }
    }

    //Calculate interaction radius based on the minimum distance from the nonlocal boundaries
    cl = ( ( 1 - beta ) / ( zeta * cl0 ) * distance + beta ) * cl0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
} // end namespace oofem
