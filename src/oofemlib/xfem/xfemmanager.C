#include "intarray.h"
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

XfemManager::XfemManager(Domain* aDomain) {
    this->domain = aDomain;
    this->enrichmentFunctionList = new AList<EnrichmentFunction > (0);
    this->enrichmentItemList = new AList<EnrichmentItem > (0);
    this->geometryList = new AList<Geometry > (0);
}

XfemManager::~XfemManager() {
    delete enrichmentItemList;
    delete geometryList;
    delete enrichmentFunctionList;
}

void XfemManager::getInteractedEI(IntArray& answer, Element * elem) {
    int count = 0;
    for (int i = 1; i <= this->giveNumberOfEnrichmentItems(); i++) {
        if (this->giveEnrichmentItem(i)->interacts(elem)) {
            count++;
            answer.resize(count);
            answer.at(count) = this->giveEnrichmentItem(i)->giveNumber();
        }
    }
}

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
        this->giveEnrichmentItem(interactedEI.at(i))->computeIntersectionPoints(intersecPoints, elem);
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

// temporary function

void XfemManager::setEnrichmentItemList(EnrichmentItem * enrichmentitem) {
    int sz = this->enrichmentItemList->giveSize();
    enrichmentItemList->put(sz + 1, enrichmentitem);
}
// temporary function

void XfemManager::setGeometryList(Geometry *geoentity) {
    int sz = this->geometryList->giveSize();
    this->geometryList->put(sz + 1, geoentity);
}

void XfemManager::setEnrichmentFunctionList(EnrichmentFunction *enrichmentFunction) {
    int sz = this->enrichmentFunctionList->giveSize();
    this->enrichmentFunctionList->put(sz + 1, enrichmentFunction);
}

EnrichmentItem* XfemManager::giveEnrichmentItem(int n)
// Returns the n-th enrichment item.
{
    if (enrichmentItemList -> includes(n)) {
        return enrichmentItemList -> at(n);
    } else {
        OOFEM_ERROR2("giveEnrichmentItem: undefined enrichmentItem (%d)", n);
    }
    return NULL;
}

Geometry* XfemManager::giveGeometry(int n)
// Returns the n-th geometry.
{
    if (geometryList -> includes(n)) {
        return geometryList -> at(n);
    } else {
        OOFEM_ERROR2("giveGeometry: undefined geometry (%d)", n);
    }
    return NULL;
}

EnrichmentFunction* XfemManager::giveEnrichmentFunction(int n)
// Returns the n-th enrichment function.
{
    if (enrichmentFunctionList -> includes(n)) {
        return enrichmentFunctionList -> at(n);
    } else {
        OOFEM_ERROR2("giveEnrichmentFunction: undefined enrichmentFunction (%d)", n);
    }
    return NULL;
}

void XfemManager::createEnrBmatrixAt(GaussPoint *gp, FloatMatrix &answer) {
    // transformation of gauss point into global coordinates,
    // interpolation implementation-specific here, I cannot get from gp to interpolation
    FEI2dQuadLin interpolation(1, 2);
    IntArray nodes(4);
    for (int i = 1; i <= nodes.giveSize(); i++) {
        nodes.at(i) = gp->giveElement()->giveDofManagerNumber(i);
    }
    FloatArray gpGlob;
    interpolation.local2global(gpGlob, domain, nodes, *gp->giveCoordinates(), 0.0);
    // evaluation of N,dNdx
    FloatMatrix dNdx;
    interpolation.evaldNdx(dNdx, domain, nodes, *gp->giveCoordinates(), 0.0);
    FloatArray N;
    interpolation.edgeEvalN(N, *gp->giveCoordinates(), 0.0);
    // obtaining of Enrichment Items
    IntArray interactedEI;
    this->getInteractedEI(interactedEI, gp->giveElement());
    for (int i = 1; i <= interactedEI.giveSize(); i++) {
        EnrichmentItem *er = this->giveEnrichmentItem(interactedEI.at(i));
        // enrichment function at the global gauss point
        double efgp = er->giveEnrichmentFunction()->evaluateFunctionAt(&gpGlob);
        // derivative of enrichment function at the global gauss point
        FloatArray efgpD;
        er->giveEnrichmentFunction()->evaluateDerivativeAt(efgpD, &gpGlob);
        int enrCount = 0;
        answer.resize(1,1);
        for (int j = 1; j <= gp->giveElement()->giveNumberOfDofManagers(); j++) {
            if (this->computeNodeEnrichmentType(nodes.at(j)) != 0) {
                enrCount++;
                FloatArray *node = domain->giveDofManager(nodes.at(j))->giveCoordinates();
                double efnode = er->giveEnrichmentFunction()->evaluateFunctionAt(node);
                // matrix to be added anytime a node is enriched
                FloatMatrix toAdd(3, 2);
                toAdd.zero();
                double aa = dNdx.at(j, 1)*(efgp - efnode) + N.at(j) * efgpD.at(1);
                double bb = dNdx.at(j, 2)*(efgp - efnode) + N.at(j) * efgpD.at(2);
                toAdd.at(1, 1) = aa;
                toAdd.at(2, 2) = bb;
                toAdd.at(3, 1) = bb;
                toAdd.at(3, 2) = aa;
                answer.addSubMatrix(toAdd, 1, 2*enrCount - 1);     
            }    
        }
    }
}

void XfemManager::computeFictPosition() {
    // gives for a particular node position of its fictitious node
    // it is supposed that the fictitious nodes are at the very end
    // this is supposed to be used for creation of locationArray for a dofmanager
    int nrNodes = domain->giveNumberOfDofManagers();
    this->fictPosition.resize(nrNodes);
    int count = 0;
    for (int i = 1; i <= nrNodes; i++) {
        int enr = this->computeNodeEnrichmentType(i);
        if (enr > 0) {
            int sz = fictPosition.giveSize();
            fictPosition.resize(sz + 1);
            fictPosition.at(i) = nrNodes + count + 1;
            count = count + enr;
        }
    }
}

XfemManager::XfemType XfemManager::computeNodeEnrichmentType(int nodeNumber) {
    XfemType ret;
    int intersectionCount = 0;
    const IntArray *neighbours = domain->giveConnectivityTable()->giveDofManConnectivityArray(nodeNumber);
    for (int i = 1; i <= neighbours->giveSize(); i++) {
        IntArray interactedEnrEl;
        getInteractedEI(interactedEnrEl, domain->giveElement(neighbours->at(i)));
        for (int j = 1; j <= interactedEnrEl.giveSize(); j++) {
            intersectionCount += this->giveEnrichmentItem(interactedEnrEl.at(j))->computeNumberOfIntersectionPoints(domain->giveElement(neighbours->at(i)));
        }
        interactedEnrEl.zero();
    }
    if (intersectionCount == 0) ret = STANDARD;
    else if (intersectionCount == 1) ret = TIP;
    else if (intersectionCount > 1) ret = SPLIT;
    return ret;
}

