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
        AList< Triangle >triangles2;
        // all the points coming into triangulation
        AList< FloatArray >together1;
        AList< FloatArray >together2;
        this->XfemElementInterface_prepareNodesForDelaunay(& together1, & together2);
        this->XfemElementInterface_partitionElement(& triangles, & together1);
        this->XfemElementInterface_partitionElement(& triangles2, & together2);
        for(int i = 1; i <= triangles2.giveSize(); i++){
            int sz = triangles.giveSize();
            triangles.put(sz + 1, triangles2.at(i));
            triangles2.unlink(i);
        }
        PatchIntegrationRule *pir = new PatchIntegrationRule(1, element, &triangles);
        int pointNr = 3 * triangles.giveSize();
        integrationDomain integDomain = element->giveIntegrationDomain();
        MaterialMode matMode = element->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0)->giveMaterialMode();
        pir->setUpIntegrationPoints(integDomain, pointNr, matMode);
        // here the old integration rule is deleted and patch integration rule is set
        element->setIntegrationRule(1, pir);
    } 
}

void XfemElementInterface :: XfemElementInterface_prepareNodesForDelaunay(AList< FloatArray > *answer1, AList< FloatArray > *answer2) {
    XfemManager *xf = this->element->giveDomain()->giveEngngModel()->giveXfemManager(1);
    IntArray interactedEI;
    xf->getInteractedEI(interactedEI, element);
    // in intersecPoints the points of Element with interaction to EnrichmentItem will be stored
    AList< FloatArray > *intersecPoints = new AList< FloatArray >(0);
    for ( int i = 1; i <= interactedEI.giveSize(); i++ ) {
        xf->giveEnrichmentItem( interactedEI.at(i) )->computeIntersectionPoints(intersecPoints, element);
    }
    
    for ( int i = 1; i <= intersecPoints->giveSize(); i++ ) {
        int sz = answer1->giveSize();
        answer1->put(sz + 1, intersecPoints->at(i));
        FloatArray *node = intersecPoints->at(i);
        FloatArray *nodesCopy = new FloatArray(*node);
        int sz2 = answer2->giveSize();
        answer2->put(sz2 + 1, nodesCopy);
        
    }
    if(intersecPoints->giveSize() == 2) {
       
       double x1 = intersecPoints->at(1)->at(1);
       double x2 = intersecPoints->at(2)->at(1);
       double y1 = intersecPoints->at(1)->at(2);
       double y2 = intersecPoints->at(2)->at(2);
       for(int i = 1; i <= this->element->giveNumberOfDofManagers(); i++){
          double x = element->giveDofManager(i)->giveCoordinates()->at(1);
          double y = element->giveDofManager(i)->giveCoordinates()->at(2);
          double det = (x1 - x)*(y2 - y) - (x2 - x)*(y1 - y);
          FloatArray *node = element->giveDofManager(i)->giveCoordinates();
          FloatArray *nodesCopy = new FloatArray(*node);
          if(det > 0.00001) { 
             int sz = answer1->giveSize();
             answer1->put(sz + 1, nodesCopy);
           }
          else if (det < (-1)*0.00001){
             int sz = answer2->giveSize();
             answer2->put(sz + 1, nodesCopy);
          }
       }
    }
    
    // nodes of an element are copied to a different memory location
    // so that the whole container of points for triangulation can be dealt with
    // more easily (e.g. deleted)

    for ( int i = 1; i <= intersecPoints->giveSize(); i++ ) {
        intersecPoints->unlink(i);
    }
    delete intersecPoints; 
}

