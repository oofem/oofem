/*
 * File:   xfemmanager.h
 * Author: chamrova
 *
 * Created on October 30, 2008, 10:43 AM
 */

#ifndef _XFEMMANAGER_H
#define _XFEMMANAGER_H

#include "enrichmentitem.h"
#include "alist.h"
#include "domain.h"
#include "delaunay.h"
#include "gausspnt.h"
#include "engngm.h"
#include "patch.h"

namespace oofem {
/* this class manages the xfem part as well as takes over some functions which would appear
 * in the Domain and Node class */
class XfemManager
{
protected:
    /// associated EngineeringModel
    EngngModel *emodel;
    /// index of the associated domain
    int domainIndex;
    /// Enrichment item list
    AList< EnrichmentItem > *enrichmentItemList;
    /// Geometry list
    AList< BasicGeometry > *geometryList;
    /// Enrichment function list
    AList< EnrichmentFunction > *enrichmentFunctionList;
    /// map giving for a node a position of its fictitious node
    AList< IntArray > *fictPosition;
    /// index of next available dofId from pool
    int dofIdPos;
    int numberOfEnrichmentItems;
    int numberOfEnrichmentFunctions;
    int numberOfGeometryItems;

public:
    enum XfemType {
        SPLIT = 1, TIP = 4, STANDARD = 0
    };
    /// Constructor
    XfemManager(EngngModel *emodel, int index);
    /// destructor
    ~XfemManager();
    /* gets interacted enrichment items for a particular element, the enrichment items
     * are referenced by a number from the domain */
    void getInteractedEI(IntArray &answer, Element *elem);
    /// checks whether an element is interacted
    bool isInteracted(Element *elem);
    /// checks whether a node is interacted
    bool isEnriched(int nodeNumber);
    /// Accessor
    EnrichmentItem *giveEnrichmentItem(int n);
    /// Accessor
    BasicGeometry *giveGeometry(int n);
    /// Accessor
    EnrichmentFunction *giveEnrichmentFunction(int n);
    int giveNumberOfEnrichmentItems() { return enrichmentItemList->giveSize(); }
    /// computes for each node position of its fictitious node
    int computeFictPosition();
    /// computes the type of node enrichment, returns zero if the node is not enriched
    XfemType computeNodeEnrichmentType(int nodeNumber);
    /// Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);
    /// Instantiates the xfem components
    int instanciateYourself(DataReader *dr);
    const char *giveClassName() const { return "XfemManager"; }
    const char *giveInputRecordName() const { return "XfemManager"; }
    /// wrapper for updating the integration rule
    void updateIntegrationRule();
    /// gives Domain
    Domain *giveDomain();
    /// accessor
    IntArray *giveFictPosition(int nodeNumber) { return fictPosition->at(nodeNumber); }
    /// Stores receiver state to output stream.
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL) { return CIO_OK; }
    /// Restores the receiver state previously written in stream.
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL) { return CIO_OK; }
    /// geometry update; calls individual enrichment item updateGeometry method
    void updateGeometry(TimeStep *tStep);
protected:
    // changes dofIdPos to next index
    DofID allocateNewDofID();
};
} // end namespace oofem
#endif  /* _XFEMMANAGER_H */


