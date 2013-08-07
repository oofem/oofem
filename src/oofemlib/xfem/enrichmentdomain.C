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
#include "enrichmentitem.h"
#include "element.h"
#include "dofmanager.h"
#include "connectivitytable.h"
#include "classfactory.h"
#include "enrichmentfunction.h"
#include "xfemmanager.h"

#include <algorithm>

#include <cmath>

namespace oofem {

REGISTER_EnrichmentDomain( DofManList )
REGISTER_EnrichmentDomain( WholeDomain )
REGISTER_EnrichmentDomain( EDBGCircle )

REGISTER_EnrichmentDomain( EDCrack )

//REGISTER_EnrichmentDomain( BasicGeometryDomain<Line> )

// General 

EnrichmentDomain::EnrichmentDomain()
{

}

void EnrichmentDomain_BG :: CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan)
{
	iEnrItem.updateNodeEnrMarker(ixFemMan, *this);
}

void DofManList :: CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan)
{
	iEnrItem.updateNodeEnrMarker(ixFemMan, *this);
}


IRResultType DofManList :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IntArray idList;
    IR_GIVE_FIELD(ir, idList, _IFT_DofManList_list);
    for ( int i = 1; i<=idList.giveSize(); i++) {
        this->dofManList.push_back( idList.at(i) );
    }
    return IRRT_OK;
    
}

void WholeDomain :: CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan)
{
	iEnrItem.updateNodeEnrMarker(ixFemMan, *this);
}

} // end namespace oofem
