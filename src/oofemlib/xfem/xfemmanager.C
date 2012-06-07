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

#include "inputrecord.h"
#include "intarray.h"
#include "femcmpnn.h"
#include "engngm.h"
#include "geometry.h"
#include "conTable.h"
#include "flotarry.h"
#include "alist.h"
#include "domain.h"
#include "fei2dtrlin.h"
#include "enrichmentitem.h"
#include "integrationrule.h"
#include "xfemmanager.h"
#include "element.h"
#include "dofmanager.h"
#include "delaunay.h"
#include "cltypes.h"
#include "patch.h"
#include "enrichmentfunction.h"
#include "gaussintegrationrule.h"
#include "fei2dquadlin.h"
#include "xfemelementinterface.h"
#include "oofem_limits.h"
#include "usrdefsub.h"
#include "masterdof.h"
#include "xfem/xfemmanager.h"
#include "patchintegrationrule.h"
#include "datareader.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
XfemManager :: XfemManager(EngngModel *emodel, int domainIndex)
{
    this->emodel = emodel;
    this->domainIndex = domainIndex;
    this->enrichmentFunctionList = new AList< EnrichmentFunction >(0);
    this->enrichmentItemList = new AList< EnrichmentItem >(0);
    this->geometryList = new AList< BasicGeometry >(0);
    this->fictPosition = new AList< IntArray >();
    numberOfEnrichmentItems = 0;
    numberOfEnrichmentFunctions = 0;
    numberOfGeometryItems = 0;
    dofIdPos = 0;
}

XfemManager :: ~XfemManager()
{
    delete enrichmentItemList;
    delete geometryList;
    delete enrichmentFunctionList;
    delete fictPosition;
}

void
XfemManager :: clear()
{
    delete enrichmentItemList;
    enrichmentItemList = NULL;

    delete geometryList;
    geometryList = NULL;

    delete enrichmentFunctionList;
    enrichmentFunctionList = NULL;

    delete fictPosition;
    fictPosition = NULL;

    numberOfEnrichmentItems = 0;
    numberOfEnrichmentFunctions = 0;
    numberOfGeometryItems = 0;
    dofIdPos = 0;
}



Domain *XfemManager :: giveDomain() { return emodel->giveDomain(domainIndex); }

void XfemManager :: getInteractedEI(IntArray &answer, Element *elem)
{
    int count = 0;
    for ( int i = 1; i <= this->giveNumberOfEnrichmentItems(); i++ ) {
        if ( this->giveEnrichmentItem(i)->interacts(elem) ) {
            count++;
            answer.resize(count);
            answer.at(count) = enrichmentItemList->at(i)->giveNumber();
        }
    }
}

bool XfemManager :: isInteracted(Element *elem)
{
    IntArray interactedEnrItems;
    getInteractedEI(interactedEnrItems, elem);
    if ( interactedEnrItems.giveSize() > 0 ) {
        return true;
    } else {
        return false;
    }
}

EnrichmentItem *XfemManager :: giveEnrichmentItem(int n)
// Returns the n-th enrichment item.
{
    if ( enrichmentItemList->includes(n) ) {
        return enrichmentItemList->at(n);
    } else {
        OOFEM_ERROR2("giveEnrichmentItem: undefined enrichmentItem (%d)", n);
    }

    return NULL;
}

BasicGeometry *XfemManager :: giveGeometry(int n)
// Returns the n-th geometry.
{
    if ( geometryList->includes(n) ) {
        return geometryList->at(n);
    } else {
        OOFEM_ERROR2("giveGeometry: undefined geometry (%d)", n);
    }

    return NULL;
}

EnrichmentFunction *XfemManager :: giveEnrichmentFunction(int n)
// Returns the n-th enrichment function.
{
    if ( enrichmentFunctionList->includes(n) ) {
        return enrichmentFunctionList->at(n);
    } else {
        OOFEM_ERROR2("giveEnrichmentFunction: undefined enrichmentFunction (%d)", n);
    }

    return NULL;
}

int XfemManager :: computeFictPosition()
{
    // gives for a particular node position of its fictitious node
    // it is supposed that the fictitious nodes are at the very end
    // this is supposed to be used for creation of locationArray for a dofmanager
    // the function returns simultaneously the last dof
    int nrNodes = emodel->giveDomain(1)->giveNumberOfDofManagers();
    int count = this->giveDomain()->giveEngngModel()->giveNumberOfEquations(EID_MomentumBalance);
    IntArray edofs;
    for ( int j = 1; j <= nrNodes; j++ ) {
        IntArray *dofs = new IntArray();

        for ( int i = 1; i <= this->enrichmentItemList->giveSize(); i++ ) {
            int dofSize = enrichmentItemList->at(i)->getDofIdArray()->giveSize();
            if ( enrichmentItemList->at(i)->isDofManEnriched(j) ) {
                edofs.resize(dofSize);
                for ( int k = 1; k <= dofSize; k++ ) {
                    count++;
                    edofs.at(k) = count;
                }

                dofs->followedBy(edofs);
            }
        }

        fictPosition->put(j, dofs);
    }

    return count;
}

bool XfemManager :: isEnriched(int nodeNumber)
{
    XfemManager :: XfemType nodeEnrType = this->computeNodeEnrichmentType(nodeNumber);
    if ( nodeEnrType != 0 ) {
        return true;
    } else {
        return false;
    }
}

XfemManager :: XfemType XfemManager :: computeNodeEnrichmentType(int nodeNumber)
{
    XfemType ret;
    int intersectionCount = 0;
    const IntArray *neighbours = emodel->giveDomain(domainIndex)->giveConnectivityTable()->giveDofManConnectivityArray(nodeNumber);
    for ( int i = 1; i <= neighbours->giveSize(); i++ ) {
        IntArray interactedEnrEl;
        getInteractedEI( interactedEnrEl, emodel->giveDomain(domainIndex)->giveElement( neighbours->at(i) ) );
        for ( int j = 1; j <= interactedEnrEl.giveSize(); j++ ) {
            intersectionCount += this->giveEnrichmentItem( interactedEnrEl.at(j) )->computeNumberOfIntersectionPoints( emodel->giveDomain(domainIndex)->giveElement( neighbours->at(i) ) );
        }

        interactedEnrEl.zero();
    }

    if ( intersectionCount == 0 ) {
        ret = STANDARD;
    } else if ( intersectionCount == 1 ) {
        ret = TIP;
    } else if ( intersectionCount > 1 ) {
        ret = SPLIT;
    }

    return ret;
}

IRResultType XfemManager :: initializeFrom(InputRecord *ir) {
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    if ( !ir->hasField(IFT_XfemManager_name, "xfemmanager") ) {
        OOFEM_ERROR("XfemManager::instanciateFrom: bad record");
    }

    IR_GIVE_FIELD(ir, numberOfEnrichmentFunctions, IFT_XfemManager_numberOfEnrichmentFunctions, "numberofenrichmentfunctions"); // Macro// Macro
    IR_GIVE_FIELD(ir, numberOfEnrichmentItems, IFT_XfemManager_numberOfEnrichmentItems, "numberofenrichmentitems"); // Macro
    IR_GIVE_FIELD(ir, numberOfGeometryItems, IFT_XfemManager_numberOfGeometryItems, "numberofgeometryitems");
    return IRRT_OK;
}

int XfemManager :: instanciateYourself(DataReader *dr)
{
    const char *__proc = "instanciateYourself"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    int i;
    std :: string name;
    EnrichmentItem *ei;
    EnrichmentFunction *ef;
    BasicGeometry *ge;
    InputRecord *mir;

    enrichmentFunctionList->growTo(numberOfEnrichmentFunctions);
    for ( i = 0; i < numberOfEnrichmentFunctions; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_enrichFuncRec, i + 1);
        result = mir->giveRecordKeywordField(name);

        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }

        ef = CreateUsrDefEnrichmentFunction( name.c_str(), i + 1, emodel->giveDomain(1) );
        if ( ef == NULL ) {
            OOFEM_ERROR2( "XfemManager::instanciateYourself: unknown enrichment function (%s)", name.c_str() );
        }

        enrichmentFunctionList->put(i + 1, ef);
        ef->initializeFrom(mir);
    }

    geometryList->growTo(numberOfGeometryItems);
    for ( i = 0; i < numberOfGeometryItems; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_geoRec, i + 1);
        result = mir->giveRecordKeywordField(name);
        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }

        ge = CreateUsrDefGeometry( name.c_str() );

        if ( ge == NULL ) {
            OOFEM_ERROR2( "XfemManager::instanciateYourself: unknown geometry (%s)", name.c_str() );
        }

        geometryList->put(i + 1, ge);
        ge->initializeFrom(mir);
    }

    enrichmentItemList->growTo(numberOfEnrichmentItems);
    for ( i = 0; i < numberOfEnrichmentItems; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_enrichItemRec, i + 1);
        result = mir->giveRecordKeywordField(name);

        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }

        ei = CreateUsrDefEnrichmentItem( name.c_str(), i + 1, this, emodel->giveDomain(1) );

        if ( ei == NULL ) {
            OOFEM_ERROR2( "XfemManager::instanciateYourself: unknown enrichment item (%s)", name.c_str() );
        }


        ei->initializeFrom(mir);
        int eindofs = ei->giveNumberOfDofs();
        IntArray dofIds(eindofs);
        for ( int j = 1; j <= eindofs; j++ ) {
            dofIds.at(j) = this->allocateNewDofID();
        }

        ei->setDofIdArray(dofIds);
        enrichmentItemList->put(i + 1, ei);
    }

#ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated enrichment items ", nnode)
#endif
    return 1;
}

DofIDItem XfemManager :: allocateNewDofID()
{
    int answer = 0;
    if ( dofIdPos < ( X_N - X_1 ) ) {
        answer = dofIdPos + X_1;
        dofIdPos++;
    }

    return ( DofIDItem ) answer;
}

void XfemManager :: updateIntegrationRule()
{
    for ( int i = 1; i <= this->giveDomain()->giveNumberOfElements(); i++ ) {
        Element *el = this->giveDomain()->giveElement(i);
        XfemElementInterface *xei = ( XfemElementInterface * ) el->giveInterface(XfemElementInterfaceType);
        xei->XfemElementInterface_updateIntegrationRule();
    }
}

void XfemManager :: updateGeometry(TimeStep *tStep)
{
    for ( int i = 1; i <= this->enrichmentItemList->giveSize(); i++ ) {
        enrichmentItemList->at(i)->updateGeometry(tStep);
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
                obj = CreateUsrDefEnrichmentItem(compId, i, this, emodel->giveDomain(domainIndex));
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
                obj = CreateUsrDefEnrichmentFunction(compId, i, emodel->giveDomain(domainIndex));
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

} // end namespace oofem
