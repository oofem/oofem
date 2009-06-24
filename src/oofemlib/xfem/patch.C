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

Patch :: Patch(Element *parent, int material) : BasicGeometry() {
    this->parent = parent;
    this->material = material;
}

Patch::Patch(Element *parent, AList<FloatArray> *vertices) : BasicGeometry() {
    this->parent = parent;
    this->vertices = vertices;
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




