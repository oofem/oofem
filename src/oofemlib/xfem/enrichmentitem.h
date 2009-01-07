/*
 * File:   enrichmentitem.h
 * Author: chamrova
 *
 * Created on October 26, 2008, 11:42 AM
 */

#ifndef _ENRICHMENTITEM_H
#define	_ENRICHMENTITEM_H

#include "femcmpnn.h"
#include "domain.h"
#include "geometry.h"
#include "flotmtrx.h"
#include "enrichmentfunction.h"

class Geometry;

class EnrichmentItem : public FEMComponent
{
    public:
         // constructor; there is no destructor here, since Geometry is delete in Domain
         EnrichmentItem (int n,Domain* aDomain);
         // does nothing, but has to be here since it is pure virtual in the abstract class
         IRResultType initializeFrom (InputRecord* ir);
         // does nothing, but has to be here since it is pure virtual in the abstract class
         const char* giveClassName () const { return "EnrichmentItem" ; }
         // get method
         Geometry* giveGeometry();
         // assignes geometry to a particular enrichment item
         void setGeometry (Geometry *geometry) {this->geometry = geometry;}
         bool interacts(Element* element);
         void computeIntersectionPoints(AList<FloatArray>* intersectionPoints, Element *element);
         double computeNumberOfIntersectionPoints(Element *element);
         EnrichmentFunction* giveEnrichmentFunction() {return ef;}
         void setEnrichmentFunction(EnrichmentFunction *ef);
    protected:
         Geometry* geometry;
         EnrichmentFunction *ef;
};

class CrackTip : public EnrichmentItem {
    public:
        CrackTip (int n,Domain* aDomain) : EnrichmentItem(n, aDomain) {}
    protected:
};

class CrackInterior : public EnrichmentItem {
    public:
        CrackInterior (int n,Domain* aDomain) : EnrichmentItem(n, aDomain) {}

};
#endif	/* _ENRICHMENTITEM_H */


