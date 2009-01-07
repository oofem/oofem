/* 
 * File:   xfemmanager.h
 * Author: chamrova
 *
 * Created on October 30, 2008, 10:43 AM
 */

#ifndef _XFEMMANAGER_H
#define	_XFEMMANAGER_H

#include "enrichmentitem.h"
#include "alist.h"
#include "domain.h"
#include "delaunay.h"
#include "gausspnt.h"

/* this class manages the xfem part as well as takes over some functions which would appear
 in the Element class */
class XfemManager {
    
protected:
    Domain *domain;
     /// Enrichment item list
    AList<EnrichmentItem>* enrichmentItemList ;
    /// Geometry list
    AList<Geometry>* geometryList ;
    /// Enrichment function list
    AList<EnrichmentFunction>* enrichmentFunctionList ;
    /// map giving for a node a position of its fictitious node
    IntArray fictPosition;

public:
    enum XfemType {
    SPLIT = 1,
    TIP = 4,
    STANDARD = 0
    };
    // constructor
    XfemManager(Domain* aDomain);
    // destructor
    ~XfemManager();
    /* gets interacted enrichment items for a particular element, the enrichment items
     are referenced by a number from the domain */
    void getInteractedEI(IntArray& answer, Element* elem);
    /* partitions an element into triangles,
    input are the vertices - together */
    void partitionElement(AList<Triangle>* triangles, AList<FloatArray>* together);
    // checks whether an element is interacted
    bool isInteracted(Element *elem);
    // updates integration points of an element, if the element is interacted
    void updateIntegrationRule(Element *elem);
    // wrap up for updateIntegrationRule for all the elemnts
    void updateIntegrationRules();
    // joins together the nodes of an element and intersection points
    void prepareNodesForDelaunay(AList<FloatArray>* answer, Element * elem);
    void setEnrichmentFunctionList(EnrichmentFunction *enrichmentFunction);
    void setEnrichmentItemList(EnrichmentItem* enrichmentitem);
    void setGeometryList(Geometry* geometry);
    EnrichmentItem* giveEnrichmentItem(int n);
    Geometry* giveGeometry(int n);
    EnrichmentFunction* giveEnrichmentFunction (int n);
    int giveNumberOfEnrichmentItems() {return enrichmentItemList->giveSize();}
    // creates enriched part of B matrix
    void createEnrBmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    // computes for each node position of its fictitious node
    void computeFictPosition();
    XfemType computeNodeEnrichmentType(int nodeNumber);
};

#endif	/* _XFEMMANAGER_H */


