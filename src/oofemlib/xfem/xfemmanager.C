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

#include "xfemmanager.h"
#include "inputrecord.h"
#include "intarray.h"
#include "conTable.h"
#include "flotarry.h"
#include "alist.h"
#include "domain.h"
#include "enrichmentitem.h"
#include "enrichmentdomain.h"
#include "element.h"
#include "dofmanager.h"
#include "cltypes.h"
#include "xfemelementinterface.h"
#include "usrdefsub.h"
#include "masterdof.h"
#include "datareader.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
XfemManager :: XfemManager(Domain *domain)
{
    this->domain = domain;
    this->enrichmentItemList = new AList< EnrichmentItem >(0);
    numberOfEnrichmentItems = -1;
}

XfemManager :: ~XfemManager()
{
    delete enrichmentItemList;
}

void
XfemManager :: clear()
{
    delete enrichmentItemList;
    enrichmentItemList = NULL;
    numberOfEnrichmentItems = -1;
}


void XfemManager :: giveActiveEIsFor(IntArray &answer, const Element *elem)
{
    for ( int i = 1; i <= this->giveNumberOfEnrichmentItems(); i++ ) {
        if ( this->giveEnrichmentItem(i)->isElementEnriched(elem) ) {
            answer.followedBy( enrichmentItemList->at(i)->giveNumber() );
        }
    }
}

bool XfemManager :: isElementEnriched(const Element *elem)
{
    // Loop over all EI which asks if el is enriched. 
    for ( int i = 1; i <= this->giveNumberOfEnrichmentItems(); i++ ){
        if ( this->giveEnrichmentItem(i)->isElementEnriched(elem) ){ 
            return true; 
        };
    } 
    return false;
}


EnrichmentItem *XfemManager :: giveEnrichmentItem(int n)
{
    // Returns the n-th enrichment item.
    if ( enrichmentItemList->includes(n) ) {
        return enrichmentItemList->at(n);
    } else {
        OOFEM_ERROR2("giveEnrichmentItem: undefined enrichmentItem (%d)", n);
    }
    return NULL;
}


// Old method: strange workflow in this method
///@todo: Broken but not in use and should be rewritten anyway
XfemManager :: XfemType XfemManager :: computeNodeEnrichmentType(int nodeNumber)
{
    XfemType ret;
    int intersectionCount = 0;

    // elements surrounding one node
    const IntArray *neighbours = this->giveDomain()->giveConnectivityTable()->giveDofManConnectivityArray(nodeNumber);
    for ( int i = 1; i <= neighbours->giveSize(); i++ ) {
        IntArray interactedEnrEl;
        // list of the EI's that are active in the element
        giveActiveEIsFor( interactedEnrEl, this->giveDomain()->giveElement( neighbours->at(i) ) );
        for ( int j = 1; j <= interactedEnrEl.giveSize(); j++ ) {
            // sums up number of intersections the geometry have with a given element
            //intersectionCount += this->giveEnrichmentItem( interactedEnrEl.at(j) )->computeNumberOfIntersectionPoints( this->giveDomain()->giveElement( neighbours->at(i) ) );
        }

        interactedEnrEl.zero();
    }
    // very specialized
    // only for 2d. Won't work if several ei are active in the neighboring element to a node.
    // one node could also have several TYPEs, Tip + inclusion etc.
    if ( intersectionCount == 0 ) {
        ret = STANDARD;
    } else if ( intersectionCount == 1 ) {
        ret = TIP;
    } else /*if ( intersectionCount > 1 )*/ {
        ret = SPLIT;
    }

    return ret;
}



void 
XfemManager :: createEnrichedDofs()
{   
    // Creates new dofs due to enrichment and appends them to the dof managers
    ///@todo: need to add check if dof already exists in the dofmanager
    int nrDofMan = this->giveDomain()->giveNumberOfDofManagers();
    IntArray dofIdArray;
 
    for (int j = 1; j <= this->giveNumberOfEnrichmentItems(); j++ ) {
        EnrichmentItem *ei = this->giveEnrichmentItem(j);
        for ( int k = 1; k <= ei->giveNumberOfEnrichmentDomains(); k++ ) {
            for ( int i = 1; i <= nrDofMan; i++ ) {
                DofManager *dMan = this->giveDomain()->giveDofManager(i); 
                if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,k) ) {
                    ei->computeDofManDofIdArray(dofIdArray, dMan, k);
                    int nDofs = dMan->giveNumberOfDofs();
                    for ( int m = 1; m<= dofIdArray.giveSize(); m++ ) {                      
                        dMan->appendDof( new MasterDof( nDofs + m, dMan, ( DofIDItem ) ( dofIdArray.at(m) ) ) );   
                    }
                }
            }        
        }
    }

}



IRResultType XfemManager :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, numberOfEnrichmentItems, _IFT_XfemManager_numberOfEnrichmentItems);
    return IRRT_OK;
}


int XfemManager :: instanciateYourself(DataReader *dr)
{
    const char *__proc = "instanciateYourself"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    std :: string name;

    enrichmentItemList->growTo(numberOfEnrichmentItems);
    for ( int i = 1; i <= numberOfEnrichmentItems; i++ ) {
        InputRecord *mir = dr->giveInputRecord(DataReader :: IR_enrichItemRec, i);
        result = mir->giveRecordKeywordField(name);

        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, "", mir, result);
        }

        EnrichmentItem *ei = CreateUsrDefEnrichmentItem( name.c_str(), i, this, this->giveDomain() );
        if ( ei == NULL ) {
            OOFEM_ERROR2( "XfemManager::instanciateYourself: unknown enrichment item (%s)", name.c_str() );
        }

        ei->initializeFrom(mir);
        ei->instanciateYourself(dr);
        this->enrichmentItemList->put(i, ei);
        this->createEnrichedDofs();
    }
    return 1;
}



void XfemManager :: updateIntegrationRule()
{
    for ( int i = 1; i <= this->giveDomain()->giveNumberOfElements(); i++ ) {
        Element *el = this->giveDomain()->giveElement(i);
        XfemElementInterface *xei = static_cast< XfemElementInterface * >( el->giveInterface(XfemElementInterfaceType) );
        xei->XfemElementInterface_updateIntegrationRule();
    }
}



#define XFEMMAN_STATE_SIZE 5

#define SAVE_COMPONENTS(size,type,giveMethod)   \
  {                                             \
    type* obj;                                  \
    for ( i = 1; i <= size; i++ ) {             \
        obj = giveMethod(i);                    \
        if ( ( mode & CM_Definition ) ) {       \
            ct =  (int) obj->giveClassID();     \
            if ( !stream->write(& ct, 1) ) {    \
                THROW_CIOERR(CIO_IOERR);        \
            }                                   \
        }                                       \
        if ( ( iores = obj->saveContext(stream, mode) ) != CIO_OK ) { \
            THROW_CIOERR(iores);                \
        }                                       \
    }                                           \
  }

#define RESTORE_COMPONENTS(size,type,resizeMethod,creator,giveMethod,setMethod) \
  {                                         \
    type *obj;                              \
    if ( mode & CM_Definition ) {           \
        resizeMethod(size);                 \
    }                                       \
    for ( i = 1; i <= size; i++ ) {         \
        if ( mode & CM_Definition ) {       \
            if ( !stream->read(& ct, 1) ) { \
                THROW_CIOERR(CIO_IOERR);    \
            }                               \
            compId = ( classType ) ct;      \
            obj = creator(compId, 0, this); \
        } else {                            \
            obj = giveMethod(i);            \
        }                                   \
        if ( ( iores = obj->restoreContext(stream, mode) ) != CIO_OK ) { \
            THROW_CIOERR(iores);            \
        }                                   \
        if ( mode & CM_Definition ) {       \
            setMethod(i, obj);              \
        }                                   \
    }                                       \
  }

///@todo: not fixed yet
#if 0
contextIOResultType XfemManager :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    int i,ct;

    if (mode & CM_Definition) {
        int _state[XFEMMAN_STATE_SIZE];
        _state[0]=this->domainIndex;
        _state[1]=this->dofIdPos;
        _state[2]=this->numberOfEnrichmentItems;
        _state[3]=this->numberOfEnrichmentFunctions;
        _state[4]=this->numberOfGeometryItems;
        if ( !stream->write(_state, XFEMMAN_STATE_SIZE) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
    // store enrichment items
    SAVE_COMPONENTS(this->numberOfEnrichmentItems,EnrichmentItem,this->giveEnrichmentItem);
    // store enrichment functions
    SAVE_COMPONENTS(this->numberOfEnrichmentFunctions,EnrichmentFunction,this->giveEnrichmentFunction);
    // store geometry items
    SAVE_COMPONENTS(this->numberOfGeometryItems,BasicGeometry,this->giveGeometry);
    // store fictPosition map
    int _size = fictPosition->giveSize();
    if ( !stream->write(& _size, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    for (i=1; i<=_size; i++) {
        if ((iores = this->fictPosition->at(i)->storeYourself(stream, mode)) != CIO_OK) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType XfemManager:: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    int i, ct;
    classType compId;

    if (mode & CM_Definition) {
        int _state[XFEMMAN_STATE_SIZE];
        if ( !stream->read(_state, XFEMMAN_STATE_SIZE) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        this->domainIndex = _state[0];
        this->dofIdPos = _state[1];
        this->numberOfEnrichmentItems = _state[2];
        this->numberOfEnrichmentFunctions = _state[3];
        this->numberOfGeometryItems = _state[4];
    }

    {  // restore enrichment items
        EnrichmentItem *obj;
        if ( mode & CM_Definition ) {
            enrichmentItemList->clear();
            enrichmentItemList->growTo(this->numberOfEnrichmentItems);
        }
        for ( i = 1; i <= this->numberOfEnrichmentItems; i++ ) {
            if ( mode & CM_Definition ) {
                if ( !stream->read(& ct, 1) ) {
                    THROW_CIOERR(CIO_IOERR);
                }
                compId = ( classType ) ct;
                obj = CreateUsrDefEnrichmentItem(compId, i, this, this->giveDomain());
            } else {
                obj = this->giveEnrichmentItem(i);
            }
            if ( ( iores = obj->restoreContext(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
            if ( mode & CM_Definition ) {
                enrichmentItemList->put(i, obj);
            }
        }
    }

    {  // restore enrichment functions
        EnrichmentFunction *obj;
        if ( mode & CM_Definition ) {
            enrichmentFunctionList->clear();
            enrichmentFunctionList->growTo(this->numberOfEnrichmentFunctions);
        }
        for ( i = 1; i <= this->numberOfEnrichmentFunctions; i++ ) {
            if ( mode & CM_Definition ) {
                if ( !stream->read(& ct, 1) ) {
                    THROW_CIOERR(CIO_IOERR);
                }
                compId = ( classType ) ct;
                obj = CreateUsrDefEnrichmentFunction(compId, i, this->giveDomain());
            } else {
                obj = this->giveEnrichmentFunction(i);
            }
            if ( ( iores = obj->restoreContext(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
            if ( mode & CM_Definition ) {
                enrichmentFunctionList->put(i, obj);
            }
        }
    }

    {  // restore geometry items
        BasicGeometry *obj;
        if ( mode & CM_Definition ) {
            geometryList->clear();
            geometryList->growTo(this->numberOfGeometryItems);
        }
        for ( i = 1; i <= this->numberOfGeometryItems; i++ ) {
            if ( mode & CM_Definition ) {
                if ( !stream->read(& ct, 1) ) {
                    THROW_CIOERR(CIO_IOERR);
                }
                compId = ( classType ) ct;
                obj = CreateUsrDefGeometry(compId);
            } else {
                obj = this->giveGeometry(i);
            }
            if ( ( iores = obj->restoreContext(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
            if ( mode & CM_Definition ) {
                geometryList->put(i, obj);
            }
        }
    }
    // restore fictPosition map
    int _size;
    IntArray *p;
    if ( !stream->read(& _size, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    this->fictPosition->clear();
    this->fictPosition->growTo(_size);
    for (i=1; i<=_size; i++) {
        p = new IntArray();
        if ((iores = p->restoreYourself(stream, mode)) != CIO_OK) {
            THROW_CIOERR(iores);
        }
        this->fictPosition->put(i, p);
    }

    return CIO_OK;
}

#endif

} // end namespace oofem
