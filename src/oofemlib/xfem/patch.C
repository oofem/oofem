#include "flotarry.h"
#include "gausspnt.h"
#include "integrationrule.h"
#include "patch.h"
#include "gaussintegrationrule.h"
#include "delaunay.h"
#include <math.h>
#include "fei2dtrlin.h"
#include "crosssection.h"

Patch :: Patch(Element *parent) : BasicGeometry() {
    this->parent = parent;
}

FEI2dTrLin TrianglePatch :: interpolation(1, 2);

void TrianglePatch :: convertGPIntoParental(GaussPoint *gp) {
    Triangle *tr = ( Triangle * ) this;
    int nrVertices = this->giveNrVertices();
    FloatArray N;
    this->interpolation.evalN(N, * gp->giveCoordinates(), 1.0);
    // stores the parental coordinates of one of the vertices
    FloatArray lCoor;
    FloatArray xcoor(nrVertices);
    FloatArray ycoor(nrVertices);

    for ( int i = 1; i <= nrVertices; i++ ) {
        FloatArray *point = tr->giveVertex(i);
        this->parent->computeLocalCoordinates(lCoor, * point);
        xcoor.at(i) = lCoor.at(1);
        ycoor.at(i) = lCoor.at(2);
    }

    // stores the parental coordinates of the Gauss Points
    FloatArray gpCoorInParental(2);
    gpCoorInParental.at(1) = dotProduct(N, xcoor, nrVertices);
    gpCoorInParental.at(2) = dotProduct(N, ycoor, nrVertices);
    gp->setCoordinates(gpCoorInParental);
    // this has to change, the computation works only in this particular case
    gp->setWeight( 2 * tr->getArea() * gp->giveWeight() );
}

