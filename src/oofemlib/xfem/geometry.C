#include <math.h>
#include "geometry.h"
#include "element.h"
#include "dofmanager.h"
#include <iostream>

IRResultType Geometry::initializeFrom(InputRecord* ir) {
    return IRRT_OK;
}

bool Geometry::intersects(Element* element) {
    if (computeNumberOfIntersectionPoints(element) == 0)
        return false;
    else return true;
}

Line::Line(int n, Domain* aDomain, FloatArray* pointA, FloatArray* pointB) : Geometry(n, aDomain), pointA(pointA), pointB(pointB) {
}

Line::~Line() {
    delete pointA;
    delete pointB;
}

double Line::computeDistanceTo(FloatArray* point) {

    double a = pointA->at(2) - pointB->at(2);
    double b = pointB->at(1) - pointA->at(1);
    double c = pointA->at(1) * pointB->at(2) - pointB->at(1) * pointA->at(2);
    double l = pointA->distance(pointB);
    return (a * point->at(1) + b * point->at(2) + c) / l;

}

void Line::computeProjection(FloatArray& answer) {
    answer.resize(2);
    answer = *pointB - *pointA;
}

double Line::computeTangentialDistanceToEnd(FloatArray* point) {
    FloatArray projection;
    this->computeProjection(projection);
    FloatArray t(2);
    t = projection * (1.0 / projection.computeNorm());
    return dotProduct(*point - *pointB, t, 2);
}

int Line::computeNumberOfIntersectionPoints(Element* element) {
    int count = 0;
    int nrNodes = element->giveNumberOfDofManagers();
    FloatArray signedDist(nrNodes);
    FloatArray tanSignDist(nrNodes);
    // here I need to get max value and min value in the FloatArray
    // there is no function for that in FloatArray
    // below I set some values to start from, the values should be related to element size
    // rather than very big doubles
    double maxDist = -1000000.0;
    double maxTanDist = -1000000.0;
    double minDist = 1000000.0;
    double minTanDist = 1000000.0;
    for (int i = 1; i <= nrNodes; i++) {
        signedDist.at(i) = computeDistanceTo(element->giveDofManager(i)->giveCoordinates());
        tanSignDist.at(i) = computeTangentialDistanceToEnd(element->giveDofManager(i)->giveCoordinates());
    }
    for (int i = 1; i <= nrNodes; i++) {
        // finding out max and min values
        if (signedDist.at(i) > maxDist) maxDist = signedDist.at(i);
        if (tanSignDist.at(i) > maxTanDist) maxTanDist = tanSignDist.at(i);
        if (signedDist.at(i) < minDist) minDist = signedDist.at(i);
        if (tanSignDist.at(i) < minTanDist) minTanDist = tanSignDist.at(i);
    }
    if ((maxDist * minDist) < 0.0) {
        if (maxTanDist < 0.0) count = 2;
        else if ((maxTanDist * minTanDist) < 0) count = 1;
        else count = 0;
    }
    return count;
}

void Line::computeIntersectionPoints(Element* element, AList<FloatArray>* intersecPoints) {
    for (int i = 1; i <= element->giveNumberOfDofManagers(); i++) {
        int n1 = i;
        int n2 = 0;
        if (i < element->giveNumberOfDofManagers()) {
            n2 = i + 1;
        } else {
            n2 = 1;
        }
        double lsn1 = computeDistanceTo(element->giveDofManager(n1)->giveCoordinates());
        double lsn2 = computeDistanceTo(element->giveDofManager(n2)->giveCoordinates());
        double lst1 = computeTangentialDistanceToEnd(element->giveDofManager(n1)->giveCoordinates());
        double lst2 = computeTangentialDistanceToEnd(element->giveDofManager(n2)->giveCoordinates());
        if (lsn1 * lsn2 < 0 && lst1 > 0 && lst2 > 0) {
            double r = lsn1 / (lsn1 - lsn2);
            if (i <= element->giveNumberOfDofManagers()) {
                FloatArray* answer = new FloatArray(2);
                for (int j = 1; j <= answer->giveSize(); j++) {
                    answer->at(j) = (1 - r) * element->giveDofManager(n1)->giveCoordinate(j)
                            + r * element->giveDofManager(n2)->giveCoordinate(j);
                }
                int sz = intersecPoints->giveSize();
                intersecPoints->put(sz + 1, answer);
            }
        }
    }
}

double Line::computeInclinationAngle() {
    double y = pointB->at(2) - pointA->at(2);
    double x = pointB->at(1) - pointA->at(1);
    return atan2(y, x);
}

void Line::computeTransformationMatrix(FloatMatrix &answer){
    answer.resize(2,2);
    double alpha = this->computeInclinationAngle();
    answer.at(1,1) = cos(alpha);
    answer.at(1,2) = sin(alpha);
    answer.at(2,1) = (-1)*sin(alpha);
    answer.at(2,2) = cos(alpha);
}

 void Line::transformIntoPolar(FloatArray *point, FloatArray &answer){
     FloatArray* xp;
     FloatMatrix Qt;
     this->computeTransformationMatrix(Qt);
     FloatArray help = *point - *pointB;
     xp = Qt.Times(&help);
     answer.resize(2);
     answer.at(1) = xp->computeNorm();
     answer.at(2) = atan2(xp->at(2), xp->at(1));
     delete xp;
 }


