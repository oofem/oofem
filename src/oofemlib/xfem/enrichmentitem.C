/* 
 * File:   enrichmentitem.C
 * Author: chamrova
 *
 * Created on October 26, 2008, 11:44 AM
 */

#include "enrichmentitem.h"
#include "geometry.h"
#include "element.h"
#include "enrichmentfunction.h"
#include "cltypes.h"

EnrichmentItem::EnrichmentItem(int n, Domain* aDomain): FEMComponent(n, aDomain){
    std::cout << "enrichmentitem constructor" << std::endl;
    geometry = NULL;
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
    if(geometry->interacts(element)){
        return true;
    }    
    else return false;
}

//XfemEleType EnrichmentItem::getInteractionType(Element *element){
//    XfemEleType xfeEleType = STANDARD;
//    if(this->interacts(element)){
//        if(geometry->giveNumberOfIntersectionPoints(element) == 2)
//            xfeEleType = SPLIT;
//        else
//            xfeEleType = TIP;
//    }
//    else xfeEleType = STANDARD;
//    return xfeEleType;
//}

void EnrichmentItem::getIntersectionPoints(AList<FloatArray>* intersectionPoints, Element *element){
    return geometry->giveIntersectionPoints(element,intersectionPoints);
}
