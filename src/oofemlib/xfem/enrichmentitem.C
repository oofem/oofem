/*
 * File:   enrichmentitem.C
 * Author: chamrova
 *
 * Created on October 26, 2008, 11:44 AM
 */

#include "flotmtrx.h"
#include "enrichmentitem.h"
#include "geometry.h"
#include "element.h"
#include "enrichmentfunction.h"
#include "cltypes.h"
#include "cmath"

EnrichmentItem::EnrichmentItem(int n, Domain* aDomain): FEMComponent(n, aDomain){
    geometry = NULL;
    ef = NULL;
}

Geometry*  EnrichmentItem::giveGeometry ()
{
    return this->geometry;
}

IRResultType EnrichmentItem::initializeFrom (InputRecord* ir)
{
    return IRRT_OK;
}

// checks whether enrichment item interact the element
bool EnrichmentItem::interacts(Element* element){
    if(geometry->intersects(element)){
        return true;
    }
    else return false;
}

void EnrichmentItem::computeIntersectionPoints(AList<FloatArray>* intersectionPoints, Element *element){
    geometry->computeIntersectionPoints(element,intersectionPoints);
}

double EnrichmentItem::computeNumberOfIntersectionPoints(Element *element){
    return geometry->computeNumberOfIntersectionPoints(element);
}

void EnrichmentItem::setEnrichmentFunction(EnrichmentFunction *ef){
    this->ef = ef;
}



