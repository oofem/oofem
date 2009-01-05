#include "geometry.h"
#include "element.h"
#include "dofmanager.h"

IRResultType Geometry::initializeFrom(InputRecord* ir) {
    return IRRT_OK;
}

bool Geometry::interacts(Element* element) {
    if(giveNumberOfIntersectionPoints(element) == 0)
        return false;
    else return true;
}

Line::Line(int n, Domain* aDomain, FloatArray* pointA, FloatArray* pointB) : Geometry(n, aDomain), pointA(pointA), pointB(pointB) {
}

Line::~Line() {
    delete pointA;
    delete pointB;
}

double Line::giveDistanceTo(FloatArray* point) {

    double a = pointA->at(2) - pointB->at(2);
    double b = pointB->at(1) - pointA->at(1);
    double c = pointA->at(1) * pointB->at(2) - pointB->at(1) * pointA->at(2);
    return a * point->at(1) + b * point->at(2) + c;

}

int Line::giveNumberOfIntersectionPoints(Element* element){
     int count = 0;
     // keeps the distance from the line to each of the element nodes
     FloatArray signedDist(element->giveNumberOfDofManagers());
     for (int i = 1; i <= signedDist.giveSize(); i++) {
        signedDist.at(i) = giveDistanceTo(element->giveDofManager(i)->giveCoordinates());
     }
     for (int i = 1; i <= signedDist.giveSize(); i++) {
        // the whole loop around element is done here
        int n2 = 0;
        if(i < signedDist.giveSize()) n2 = i+1;
        else n2 = 1;
        // if the sign of the points creating an edge is negative, then there is an intersection
        double mult = signedDist.at(i)*signedDist.at(n2);
        if(mult < 0) count++;
    } 
    return count;
}
void Line::giveIntersectionPoints(Element* element, FloatArray ** intersecPoints) {
    FloatArray signedDist(element->giveNumberOfDofManagers());
    for (int i = 1; i <= signedDist.giveSize(); i++) {
        signedDist.at(i) = giveDistanceTo(element->giveDofManager(i)->giveCoordinates());
    }
    int count = 0;
    for (int i = 1; i <= signedDist.giveSize(); i++) {
        int n2 = 0;
        if(i < signedDist.giveSize()) n2 = i+1;
        else n2 = 1;
        double mult = signedDist.at(i)*signedDist.at(n2);
        if(mult < 0) {
            count++;
            double r = signedDist.at(i) / (signedDist.at(i) - signedDist.at(n2));
            FloatArray* answer = new FloatArray(2);
            for (int j = 1; j <= answer->giveSize(); j++) {
                answer->at(j) = (1 - r) * element->giveDofManager(i)->giveCoordinate(j)
                        + r * element->giveDofManager(n2)->giveCoordinate(j);
            }
            intersecPoints[count] = answer;
        }
        
    }

}

void Line::giveIntersectionPoints(Element* element, AList<FloatArray>* intersecPoints) {

    for (int i = 1; i <= element->giveNumberOfDofManagers(); i++) {
        int n1 = i;
        int n2 = 0;
        if (i < element->giveNumberOfDofManagers()) {
            n2 = i + 1;
        } else {
            n2 = 1;
        }
        double lsn1 = giveDistanceTo(element->giveDofManager(n1)->giveCoordinates());
        double lsn2 = giveDistanceTo(element->giveDofManager(n2)->giveCoordinates());
        if (lsn1 * lsn2 < 0) {
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

