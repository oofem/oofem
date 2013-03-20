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
#include "conTable.h"
#include <algorithm>

namespace oofem {

// General 

bool 
EnrichmentDomain :: isElementEnriched(Element *element) 
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        if ( this->isDofManagerEnriched( element->giveDofManager(i) ) ) {
            return true;
        }
    }
    return false;
}

void
EnrichmentDomain :: updateEnrichmentDomain()
{
    if ( DofManList *ded = dynamic_cast< DofManList * > (this) )  {
//        ded->updateEnrichmentDomain();    
    }
}



// DofMan list

IRResultType DofManList :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IntArray idList;
    IR_GIVE_FIELD(ir, idList, IFT_RecordIDField, "list"); // Macro
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


void
DofManList :: addDofManagers(IntArray &dofManNumbers)
{
    for ( int i = 1; i <= dofManNumbers.giveSize(); i++) {
        std::list< int > :: iterator p;
        p = std::find(this->dofManList.begin( ), this->dofManList.end( ), dofManNumbers.at(i));
        if ( p == this->dofManList.end( ) ) { // if new node
            this->dofManList.push_back( dofManNumbers.at(i) );
        }       
    }
}


void
DofManList :: updateEnrichmentDomain(IntArray &dofManNumbers)
{
    this->addDofManagers(dofManNumbers);

}


// Circle
bool 
EDBGCircle :: isDofManagerEnriched(DofManager *dMan)
{ 
#if 1
    FloatArray coords; 
    coords = *(dMan->giveCoordinates());
    return this->bg->isInside(coords);
#else
    
    int node = dMan->giveGlobalNumber();
    // gets neighbouring elements of a node
    Domain *d = dMan->giveDomain();
    const IntArray *neighbours = d->giveConnectivityTable()->giveDofManConnectivityArray(node);
    for ( int i = 1; i <= neighbours->giveSize(); i++ ) {
        // for each of the neighbouring elements finds out whether it interacts with this EnrichmentItem
        if ( isElementEnriched( d->giveElement( neighbours->at(i) ) ) ) {
            return true;
        }
    }

    return false;
#endif
};


bool
EDBGCircle :: isElementEnriched(const Element *element) 
{
#if 1
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        if ( this->isDofManagerEnriched( element->giveDofManager(i) ) ) {
            return true;
        }
    }
    return false;
#else
    Circle *c = static_cast < Circle * > ( this->bg );
    int numIntersections = c->computeNumberOfIntersectionPoints(element);
    //int numIntersections = this->bg->computeNumberOfIntersectionPoints(element);
    if ( numIntersections > 0 ) {
        return true;
    } else {
        return false;
    }
#endif
};

} // end namespace oofem
