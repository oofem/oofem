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
#include "conTable.h"
#include "oofem_limits.h"
#include "usrdefsub.h"

EnrichmentItem::EnrichmentItem(int n, Domain* aDomain) : FEMComponent(n, aDomain) {
    geometry = NULL;
    ef = NULL;
}

BasicGeometry* EnrichmentItem::giveGeometry() {
    return this->geometry;
}

IRResultType EnrichmentItem::initializeFrom(InputRecord* ir) {
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    int geometryItemNr = 0;
    int enrichmentFunctionNr = 0;

    IR_GIVE_FIELD(ir, geometryItemNr, IFT_EnrichmentItem_geometryItemNr, "geometryitem"); // Macro
    IR_GIVE_FIELD(ir, enrichmentFunctionNr, IFT_EnrichmentItem_enrichmentFunctionNr, "enrichmentfunction"); // Macro
    BasicGeometry* enrItemGeometry = this->giveDomain()->giveEngngModel()->giveXfemManager(1)->giveGeometry(geometryItemNr);
    EnrichmentFunction* enrItemFunction = this->giveDomain()->giveEngngModel()->giveXfemManager(1)->giveEnrichmentFunction(enrichmentFunctionNr);
    this->geometry = enrItemGeometry;
    this->geometry->printYourself();
    this->setEnrichmentFunction(enrItemFunction);
    // this should go into enrichmentfunction probably
    enrItemFunction->insertEnrichmentItem(this);
    enrItemFunction->setActive(this);
    return IRRT_OK;
}

bool EnrichmentItem::interacts(Element* element) {
    return this->geometry->intersects(element);
}

bool EnrichmentItem::isOutside(BasicGeometry *bg) {
    return this->geometry->isOutside(bg);
}

void EnrichmentItem::computeIntersectionPoints(AList<FloatArray>* intersectionPoints, Element *element) {
    geometry->computeIntersectionPoints(element, intersectionPoints);
}

double EnrichmentItem::computeNumberOfIntersectionPoints(Element *element) {
    return geometry->computeNumberOfIntersectionPoints(element);
}

void EnrichmentItem::setEnrichmentFunction(EnrichmentFunction *ef) {
    this->ef = ef;
}

bool EnrichmentItem::isDofManEnriched(int nodeNumber) {
    bool ret = false;
    // gets neighbouring elements of a node
    const IntArray *neighbours = domain->giveConnectivityTable()->giveDofManConnectivityArray(nodeNumber);
    for (int i = 1; i <= neighbours->giveSize(); i++) {
        // for each of the neighbouring elements finds out whether it interacts with this EnrichmentItem
        if (this->interacts(domain->giveElement(neighbours->at(i)))) {
            ret = true;
            break;
        }
    }
    return ret;
}

IRResultType Inclusion::initializeFrom(InputRecord* ir) {
    EnrichmentItem::initializeFrom(ir);
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    int material = 0;
    IR_GIVE_FIELD(ir, material, IFT_EnrichmentItem_materialNr, "material"); // Macro
    this->mat = this->giveDomain()->giveMaterial(material);
    return IRRT_OK;
}
