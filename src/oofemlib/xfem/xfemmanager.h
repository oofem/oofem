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

/* this class manages the xfem part as well as takes over some functions which would appear
 in the Element class */
class XfemManager {
    
protected:
    Domain *domain;
public:
    enum XfemEleType {
    SPLIT = 2,
    TIP = 1,
    STANDARD
    };
    // constructor
    XfemManager(Domain* aDomain);
    /* gets interacted enrichment items for a particular element, the enrichment items
     are referenced by a number from the domain */
    void getInteractedEI(IntArray& answer, Element* elem);
    // checks whether to enrich a node
    bool enrichNode(int nodeNumber);
    // somehow does not seem to work
    //    XfemEleType getElementType(Element * elem);

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
   
};

#endif	/* _XFEMMANAGER_H */

