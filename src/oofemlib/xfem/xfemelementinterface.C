#include "xfemelementinterface.h"
#include "fei2dquadlin.h"
#include "patch.h"
#include "patchintegrationrule.h"

void XfemElementInterface :: XfemElementInterface_partitionElement(AList< Triangle > *answer, AList< FloatArray > *together) {
    Delaunay *dl = new Delaunay();
    dl->triangulate(together, answer);
    delete dl;
}

void XfemElementInterface :: XfemElementInterface_updateIntegrationRule() {
    XfemManager *xf = this->element->giveDomain()->giveEngngModel()->giveXfemManager(1);
    if ( xf->isInteracted(element) ) {
        AList< Triangle >triangles;
        // all the points coming into triangulation
        AList< FloatArray >together;
        this->XfemElementInterface_prepareNodesForDelaunay(& together);
        this->XfemElementInterface_partitionElement(& triangles, & together);
        PatchIntegrationRule *pir = new PatchIntegrationRule(1, element, &triangles);
        /*
         *      std::cout << triangles.at(0)->giveNrVertices();
         */
        int pointNr = 3 * triangles.giveSize();
        integrationDomain integDomain = element->giveIntegrationDomain();
        MaterialMode matMode = element->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0)->giveMaterialMode();
        pir->setUpIntegrationPoints(integDomain, pointNr, matMode);
        // here the old integration rule is deleted and patch integration rule is set
        element->setIntegrationRule(1, pir);
    }
}

void XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(AList< FloatArray > *answer) {
    XfemManager *xf = this->element->giveDomain()->giveEngngModel()->giveXfemManager(1);
    IntArray interactedEI;
    xf->getInteractedEI(interactedEI, element);
    // in intersecPoints the points of Element with interaction to EnrichmentItem will be stored
    AList< FloatArray > *intersecPoints = new AList< FloatArray >(0);
    for ( int i = 1; i <= interactedEI.giveSize(); i++ ) {
        xf->giveEnrichmentItem( interactedEI.at(i) )->computeIntersectionPoints(intersecPoints, element);
    }

    // nodes of an element are copied to a different memory location
    // so that the whole container of points for triangulation can be dealt with
    // more easily (e.g. deleted)
    for ( int i = 1; i <= element->giveNumberOfNodes(); i++ ) {
        FloatArray *node = element->giveDofManager(i)->giveCoordinates();
        FloatArray *nodesCopy = new FloatArray(*node);
        int sz = answer->giveSize();
        answer->put(sz + 1, nodesCopy);
    }

    for ( int i = 1; i <= intersecPoints->giveSize(); i++ ) {
        int sz = answer->giveSize();
        answer->put( sz + 1, intersecPoints->at(i) );
        intersecPoints->unlink(i);
    }

    delete intersecPoints;
}

