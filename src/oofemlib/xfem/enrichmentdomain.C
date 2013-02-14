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

#include "mathfem.h"
#include "alist.h"
#include "enrichmentdomain.h"
#include "element.h"
#include "dofmanager.h"
#include <algorithm>

namespace oofem {
    /*
    IRResultType EnrichmentDomain :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    
    //IR_GIVE_FIELD(ir, geometry, IFT_EnrichmentItem_geometryItemNr, "geometryitem"); // Macro
    
    return IRRT_OK;
}
*/




IRResultType DofManList :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IntArray idList;
    IR_GIVE_FIELD(ir, idList, IFT_Identification, "list"); // Macro
    for ( int i = 1; i<=idList.giveSize(); i++) {
        this->dofManList.push_back( idList.at(i) );
    }
    return IRRT_OK;
    
}


bool DofManList :: isDofManagerEnriched(DofManager *dMan)
{
    int dManNumber = dMan->giveNumber();
    std::list< int > :: iterator p;
    p = std::find(this->dofManList.begin( ), this->dofManList.end( ), dManNumber);
    
    if ( p == this->dofManList.end( ) ) {
        return false;
    } else {
        return true;
    }
}




/// WholeDomain


/*
IRResultType WholeDomain :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    //IntArray idList;
    //IR_GIVE_FIELD(ir, idList, IFT_Identification, "all"); // Macro
    //for ( int i = 1; i<=idList.giveSize(); i++) {
    //    this->dofManList.push_back( idList.at(i) );
    //}
    return IRRT_OK;
    
}
*/
} // end namespace oofem
