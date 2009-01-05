#include "flotarry.h"
#include "gausspnt.h"
#include "integrationrule.h"
#include "patch.h"
#include "gaussintegrationrule.h"
#include "delaunay.h"
#include <math.h>
#include "fei2dtrlin.h"
#include "crosssection.h"

FEI2dTrLin Patch::interpolation(1, 2);

Patch::Patch(int n, Domain *aDomain, Triangle* shape, Element *parent) : Element(n, aDomain) {
    this->shape = shape;
    this->parent = parent;
}

void Patch::setIntegrationRule() {
    this->numberOfIntegrationRules = 1;
    this->integrationRulesArray = new IntegrationRule*[1];
    this->integrationRulesArray[0] = new GaussIntegrationRule(1, this);
    int nPoints = this->shape->getVertices()->giveSize();
    this->integrationRulesArray[0]->setUpIntegrationPoints(_Triangle, nPoints, _PlaneStress);
}

double Patch::computeVolumeAround(GaussPoint* gp) {
    double determinant, weight, thickness, volume;
    // this is just a transform of what Delaunay outputs into what Transformation Jacobian takes
    const FloatArray ** coords = new const FloatArray*[3];
    for(int j = 0; j < this->getShape()->getVertices()->giveSize(); j++){
        coords[j] = this->getShape()->getVertices()->at(j + 1);
    }
    determinant = fabs(this->interpolation.giveTransformationJacobian(coords, *gp->giveCoordinates(), 0.0));
    delete [] coords;
    weight = gp->giveWeight();
    thickness = this->parent-> giveCrossSection() -> give('t');
    volume = determinant * weight * thickness;
    return volume;
}

double Patch::computeVolume() {
    double volume = 0;
    for (int j = 0; j < integrationRulesArray[0]->getNumberOfIntegrationPoints(); j++) {
          this->computeVolumeAround(integrationRulesArray[0]->getIntegrationPoint(j));
    }
    return volume;
}

void Patch::convertGPIntoParental(GaussPoint *gp){
    FloatArray N;
    this->interpolation.evalN(N, *gp->giveCoordinates(), 1.0);
    // stores the parental coordinates of one of the vertices
    FloatArray lCoor;
    int nrOfVertices = this->getShape()->getVertices()->giveSize();
    FloatArray xcoor(nrOfVertices);
    FloatArray ycoor(nrOfVertices);
    for(int i = 1; i <= nrOfVertices; i++){
         FloatArray *point = this->getShape()->getVertices()->at(i);
         this->parent->computeLocalCoordinates(lCoor, *point);
         xcoor.at(i) = lCoor.at(1);
         ycoor.at(i) = lCoor.at(2);
    }
    // stores the parental coordinates of the Gauss Points
    FloatArray gpCoorInParental(2);
    gpCoorInParental.at(1) = dotProduct(N,xcoor,nrOfVertices);
    gpCoorInParental.at(2) = dotProduct(N,ycoor,nrOfVertices);
    gp->setCoordinates(gpCoorInParental);
    gp->setWeight(2 * this->computeVolume() * gp->giveWeight());
}

void Patch::convertIntoParental() {
    for (int j = 0; j < integrationRulesArray[0]->getNumberOfIntegrationPoints(); j++) {
        GaussPoint *gp = integrationRulesArray[0]->getIntegrationPoint(j);
        this->convertGPIntoParental(gp);
    }
}
