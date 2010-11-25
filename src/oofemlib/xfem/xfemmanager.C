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

namespace oofem {
XfemManager :: XfemManager(EngngModel *emodel, int domainIndex) {
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

XfemManager :: ~XfemManager() {
    delete enrichmentItemList;
    delete geometryList;
    delete enrichmentFunctionList;
    delete fictPosition;
}
Domain *XfemManager :: giveDomain() { return emodel->giveDomain(domainIndex); }

void XfemManager :: getInteractedEI(IntArray &answer, Element *elem) {
    int count = 0;
    for ( int i = 1; i <= this->giveNumberOfEnrichmentItems(); i++ ) {
        if ( this->giveEnrichmentItem(i)->interacts(elem) ) {
            count++;
            answer.resize(count);
            answer.at(count) = enrichmentItemList->at(i)->giveNumber();
        }
    }
}

bool XfemManager :: isInteracted(Element *elem) {
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

int XfemManager :: computeFictPosition() {
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

bool XfemManager :: isEnriched(int nodeNumber) {
    XfemManager :: XfemType nodeEnrType = this->computeNodeEnrichmentType(nodeNumber);
    if ( nodeEnrType != 0 ) {
        return true;
    } else {
        return false;
    }
}

XfemManager :: XfemType XfemManager :: computeNodeEnrichmentType(int nodeNumber) {
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

IRResultType
XfemManager :: initializeFrom(InputRecord *ir) {
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

int
XfemManager :: instanciateYourself(DataReader *dr) {
    const char *__proc = "instanciateYourself"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    int i;
    char name [ MAX_NAME_LENGTH ];
    EnrichmentItem *ei;
    EnrichmentFunction *ef;
    BasicGeometry *ge;
    InputRecord *mir;

    enrichmentFunctionList->growTo(numberOfEnrichmentFunctions);
    for ( i = 0; i < numberOfEnrichmentFunctions; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_enrichFuncRec, i + 1);
        result = mir->giveRecordKeywordField(name, MAX_NAME_LENGTH);

        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }

        ef = CreateUsrDefEnrichmentFunction( name, i + 1, emodel->giveDomain(1) );
        if ( ef == NULL ) {
            OOFEM_ERROR2("XfemManager::instanciateYourself: unknown enrichment function (%s)", name);
        }

        enrichmentFunctionList->put(i + 1, ef);
        ef->initializeFrom(mir);
    }

    geometryList->growTo(numberOfGeometryItems);
    for ( i = 0; i < numberOfGeometryItems; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_geoRec, i + 1);
        result = mir->giveRecordKeywordField(name, MAX_NAME_LENGTH);
        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }

        ge = CreateUsrDefGeometry(name);

        if ( ge == NULL ) {
            OOFEM_ERROR2("XfemManager::instanciateYourself: unknown geometry (%s)", name);
        }

        geometryList->put(i + 1, ge);
        ge->initializeFrom(mir);
    }

    enrichmentItemList->growTo(numberOfEnrichmentItems);
    for ( i = 0; i < numberOfEnrichmentItems; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_enrichItemRec, i + 1);
        result = mir->giveRecordKeywordField(name, MAX_NAME_LENGTH);

        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }

        ei = CreateUsrDefEnrichmentItem( name, i + 1, this, emodel->giveDomain(1) );

        if ( ei == NULL ) {
            OOFEM_ERROR2("XfemManager::instanciateYourself: unknown enrichment item (%s)", name);
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

DofID XfemManager :: allocateNewDofID() {
    int answer = 0;
    if ( dofIdPos < ( X_N - X_1 ) ) {
        answer = dofIdPos + X_1;
        dofIdPos++;
    }

    return answer;
}


void XfemManager :: updateIntegrationRule() {
    for ( int i = 1; i <= this->giveDomain()->giveNumberOfElements(); i++ ) {
        Element *el = this->giveDomain()->giveElement(i);
        XfemElementInterface *xei = ( XfemElementInterface * ) el->giveInterface(XfemElementInterfaceType);
        xei->XfemElementInterface_updateIntegrationRule();
    }
}

void
XfemManager :: updateGeometry(TimeStep *tStep) {
    for ( int i = 1; i <= this->enrichmentItemList->giveSize(); i++ ) {
        enrichmentItemList->at(i)->updateGeometry(tStep);
    }
}
} // end namespace oofem
