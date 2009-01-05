
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

XfemManager::XfemManager(Domain* aDomain) {
    this->domain = aDomain;
}

void XfemManager::getInteractedEI(IntArray& answer, Element * elem) {
    int count = 0;
    for (int i = 1; i <= domain->giveNumberOfEnrichmentItems(); i++) {
        if (domain->giveEnrichmentItem(i)->interacts(elem)) {
            count++;
            answer.resize(count);
            answer.at(count) = domain->giveEnrichmentItem(i)->giveNumber();
        }
    }
}
// this will probably go into EnrichmentDetector

bool XfemManager::enrichNode(int nodeNumber) {
    // when do we enrich a node ? if its support is completely divided or only cut
    bool ret;
    IntArray interactedEnrEl;
    const IntArray *neighbours = domain->giveConnectivityTable()->giveDofManConnectivityArray(nodeNumber);
    for (int i = 1; i <= neighbours->giveSize(); i++) {
        getInteractedEI(interactedEnrEl, domain->giveElement(neighbours->at(i)));
    }
    if(interactedEnrEl.giveSize() > 0) ret = true;
    else ret = false;
    std::cout << ret << std::endl;
    return ret;
}

//XfemEleType XfemManager::getElementType(Element * elem) {
//    IntArray interactedEnrItems;
//    getInteractedEI(interactedEnrItems,elem);
//    return domain->giveEnrichmentItem(interactedEnrItems.at(1))->getInteractionType(elem);
//}

bool XfemManager::isInteracted(Element *elem) {
    IntArray interactedEnrItems;
    getInteractedEI(interactedEnrItems, elem);
    if (interactedEnrItems.giveSize() > 0) {
        return true;
    } else return false;
}

void XfemManager::prepareNodesForDelaunay(AList<FloatArray>* answer, Element * elem) {
    IntArray interactedEI;
    getInteractedEI(interactedEI, elem);
    // in intersecPoints the points of Element with interaction to EnrichmentItem will be stored
    AList<FloatArray>* intersecPoints = new AList<FloatArray > (0);
    for (int i = 1; i <= interactedEI.giveSize(); i++) {
        domain->giveEnrichmentItem(interactedEI.at(i))->getIntersectionPoints(intersecPoints, elem);
    }
    /* nodes of an element are copied to a different memory location
     so that the whole container of points for triangulation can be dealt with
     more easily (e.g. deleted)*/
    for (int i = 1; i <= elem->giveNumberOfNodes(); i++) {
        FloatArray *nodesCopy = new FloatArray();
        FloatArray *node = elem->giveDofManager(i)->giveCoordinates();
        *nodesCopy = *node;
        int sz = answer->giveSize();
        answer->put(sz + 1, nodesCopy);
    }
    for (int i = 1; i <= intersecPoints->giveSize(); i++) {
        int sz = answer->giveSize();
        answer->put(sz + 1, intersecPoints->at(i));
        intersecPoints->unlink(i);
    }
    delete intersecPoints;
}

void XfemManager::partitionElement(AList<Triangle>* answer, AList<FloatArray>* together) {
    Delaunay* dl = new Delaunay();
    dl->triangulate(together, answer);
    delete dl;
}

void XfemManager::updateIntegrationRules() {
    for (int i = 1; i <= this->domain->giveNumberOfElements(); i++) {
        this->updateIntegrationRule(this->domain->giveElement(i));
    }
}

void XfemManager::updateIntegrationRule(Element *elem) {
    if (this->isInteracted(elem)) {
        AList<Triangle>* triangles = new AList<Triangle > (0);
        // all the points coming into triangulation
        AList<FloatArray>* together = new AList<FloatArray > (0);
        this->prepareNodesForDelaunay(together, elem);
        this->partitionElement(triangles, together);
        /* now here I need to update the standard integrationRulesArray of the Parent Element
         * the problem is that there is no modifier function available in Element
         * the current points can be updated through elem->giveDefaultIntegrationRulePtr()
         * ->getIntegrationPoint->setCoordinates, but the new ones cannot be added
         * can I create a modifier for Element class?
        */
        for (int i = 1; i <= triangles->giveSize(); i++) {
            Patch *patch = new Patch(i, domain, triangles->at(i), elem);
            patch->setIntegrationRule();
            patch->convertIntoParental();
            delete patch;
        }
        delete triangles;
        delete together;
    }

}




