#include "flotarry.h"
#include "gausspnt.h"
#include "integrationrule.h"
#include "patch.h"
#include "gaussintegrationrule.h"
#include "delaunay.h"
#include <math.h>
#include "fei2dtrlin.h"
#include "crosssection.h"
#include "enrichmentitem.h"

Patch :: Patch(Element *parent) : BasicGeometry() {
    this->parent = parent;
    this->gps = new AList<GaussPoint>();
}

Patch::Patch(Element *parent, AList<FloatArray> *vertices) : BasicGeometry() {
    this->parent = parent;
    this->vertices = vertices;
    this->gps = new AList<GaussPoint>();
}

Patch::~Patch(){
    for(int i = 1; i <= gps->giveSize(); i++){
       gps->unlink(i);
    }
    delete gps;
}

void Patch::computeMaterial(){
    XfemManager *xf = parent->giveDomain()->giveEngngModel()->giveXfemManager(1);
    EnrichmentItem *er = xf->giveEnrichmentItem(1);
    if(er->isOutside(this)) this->mat = parent->giveMaterial();
    else {
        this->mat = er->giveMaterial();
    }
}

void Patch::addGps(GaussPoint *gp){
    int sz = gps->giveSize();
    gps->put(sz + 1, gp);
}

bool Patch::hasGaussPoint(GaussPoint *gp){
    bool ret = false;
    for(int i = 1; i <= gps->giveSize(); i++){
       if(gp == gps->at(i)){
          ret = true;
          break;
       }
    }
    return ret;
}

FEI2dTrLin TrianglePatch :: interpolation(1, 2);

void TrianglePatch :: convertGPIntoParental(GaussPoint *gp) {
    FloatArray global;
    const FloatArray ** coords = new const FloatArray*[this->giveNrVertices()];
    // this we should put into the function before
    for(int i = 1; i <= this->giveNrVertices(); i++){
        coords[i-1] = new FloatArray(*this->giveVertex(i));
    }
    this->interpolation.local2global(global, coords, *gp->giveCoordinates(), 1.0);
    for(int i = 1; i <= this->giveNrVertices(); i++){
        delete coords[i-1];
    }
    delete [] coords;
    FloatArray local;
    parent->computeLocalCoordinates(local, global);
    gp->setCoordinates(local);
}

