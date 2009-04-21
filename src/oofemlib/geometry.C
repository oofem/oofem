#include <math.h>
#include "alist.h"
#include "geometry.h"
#include "element.h"
#include "dofmanager.h"

BasicGeometry :: BasicGeometry() : Geometry() {
    this->vertices = new AList< FloatArray >(0);
}

BasicGeometry :: ~BasicGeometry() {
    delete vertices;
}

FloatArray *BasicGeometry :: giveVertex(int n) {
    return this->vertices->at(n);
}

void BasicGeometry :: setVertex(FloatArray *vertex) {
    int sz = this->vertices->giveSize();
    this->vertices->put(sz + 1, vertex);
}

bool BasicGeometry :: intersects(Element *element) {
    bool ret = 0;
    int ip = this->computeNumberOfIntersectionPoints(element);
    if ( ip > 0 ) {
        ret = true;
    } else {
        ret = false;
    }

    return ret;
}

Line :: Line(FloatArray *pointA, FloatArray *pointB) : BasicGeometry() {
    this->vertices->put(1, pointA);
    this->vertices->put(2, pointB);
}

double Line :: computeDistanceTo(FloatArray *point) {
    FloatArray *pointA = this->vertices->at(1);
    FloatArray *pointB = this->vertices->at(2);
    double a = pointA->at(2) - pointB->at(2);
    double b = pointB->at(1) - pointA->at(1);
    double c = pointA->at(1) * pointB->at(2) - pointB->at(1) * pointA->at(2);
    double l = pointA->distance(pointB);
    return ( a * point->at(1) + b * point->at(2) + c ) / l;
}

void Line :: computeProjection(FloatArray &answer) {
    answer.resize(2);
    answer = * vertices->at(2) - * vertices->at(1);
}

double Line :: computeTangentialDistanceToEnd(FloatArray *point) {
    FloatArray projection;
    this->computeProjection(projection);
    FloatArray t(2);
    t = projection * ( 1.0 / projection.computeNorm() );
    return dotProduct(* point - * vertices->at(2), t, 2);
}

int Line :: computeNumberOfIntersectionPoints(Element *element) {
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
    for ( int i = 1; i <= nrNodes; i++ ) {
        signedDist.at(i) = computeDistanceTo( element->giveDofManager(i)->giveCoordinates() );
        tanSignDist.at(i) = computeTangentialDistanceToEnd( element->giveDofManager(i)->giveCoordinates() );
    }

    for ( int i = 1; i <= nrNodes; i++ ) {
        // finding out max and min values
        if ( signedDist.at(i) > maxDist ) {
            maxDist = signedDist.at(i);
        }

        if ( tanSignDist.at(i) > maxTanDist ) {
            maxTanDist = tanSignDist.at(i);
        }

        if ( signedDist.at(i) < minDist ) {
            minDist = signedDist.at(i);
        }

        if ( tanSignDist.at(i) < minTanDist ) {
            minTanDist = tanSignDist.at(i);
        }
    }

    if ( ( maxDist * minDist ) < 0.0 ) {
        if ( maxTanDist < 0.0 ) {
            count = 2;
        } else if ( ( maxTanDist * minTanDist ) < 0 ) {
            count = 1;
        } else                                                     {
            count = 0;
        }
    }

    return count;
}

void Line :: computeIntersectionPoints(Element *element, AList< FloatArray > *intersecPoints) {
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        int n1 = i;
        int n2 = 0;
        if ( i < element->giveNumberOfDofManagers() ) {
            n2 = i + 1;
        } else {
            n2 = 1;
        }

        double lsn1 = computeDistanceTo( element->giveDofManager(n1)->giveCoordinates() );
        double lsn2 = computeDistanceTo( element->giveDofManager(n2)->giveCoordinates() );
        double lst1 = computeTangentialDistanceToEnd( element->giveDofManager(n1)->giveCoordinates() );
        double lst2 = computeTangentialDistanceToEnd( element->giveDofManager(n2)->giveCoordinates() );
        if ( lsn1 * lsn2 < 0 && lst1 < 0 && lst2 < 0 ) {
            double r = lsn1 / ( lsn1 - lsn2 );
            if ( i <= element->giveNumberOfDofManagers() ) {
                FloatArray *answer = new FloatArray(2);
                for ( int j = 1; j <= answer->giveSize(); j++ ) {
                    answer->at(j) = ( 1 - r ) * element->giveDofManager(n1)->giveCoordinate(j)
                    + r *element->giveDofManager(n2)->giveCoordinate(j);
                }

                int sz = intersecPoints->giveSize();
                intersecPoints->put(sz + 1, answer);
            }
        }
    }
}

double Line :: computeInclinationAngle() {
    FloatArray *pointA = this->vertices->at(1);
    FloatArray *pointB = this->vertices->at(2);
    double y = pointB->at(2) - pointA->at(2);
    double x = pointB->at(1) - pointA->at(1);
    return atan2(y, x);
}

void Line :: computeTransformationMatrix(FloatMatrix &answer) {
    answer.resize(2, 2);
    double alpha = this->computeInclinationAngle();
    answer.at(1, 1) = cos(alpha);
    answer.at(1, 2) = sin(alpha);
    answer.at(2, 1) = ( -1 ) * sin(alpha);
    answer.at(2, 2) = cos(alpha);
}

void Line :: transformIntoPolar(FloatArray *point, FloatArray &answer) {
    FloatArray *xp;
    FloatMatrix Qt;
    this->computeTransformationMatrix(Qt);
    FloatArray help = * point - * vertices->at(2);
    xp = Qt.Times(& help);
    answer.resize(2);
    answer.at(1) = xp->computeNorm();
    answer.at(2) = atan2( xp->at(2), xp->at(1) );
    delete xp;
}

IRResultType Line :: initializeFrom(InputRecord *ir) {
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    FloatArray *start = new FloatArray(2);
    FloatArray *end = new FloatArray(2);
    IR_GIVE_FIELD(ir, * start, IFT_Line_start, "start"); // Macro
    IR_GIVE_FIELD(ir, * end, IFT_Line_end, "end"); // Macro
    vertices->put(1, start);
    vertices->put(2, end);
    return IRRT_OK;
}

Triangle :: Triangle(FloatArray *p1, FloatArray *p2, FloatArray *p3) : BasicGeometry() {
    this->vertices->put(1, p1);
    this->vertices->put(2, p2);
    this->vertices->put(3, p3);
}

double Triangle :: getArea() {
    return fabs( 0.5 * ( vertices->at(1)->at(1) * ( vertices->at(2)->at(2) - vertices->at(3)->at(2) )
                        + vertices->at(2)->at(1) * ( vertices->at(3)->at(2) - vertices->at(1)->at(2) ) +
                        vertices->at(3)->at(1) * ( vertices->at(1)->at(2) - vertices->at(2)->at(2) ) ) );
}

double Triangle :: getRadiusOfCircumCircle() {
    return 0.25 * vertices->at(1)->distance( vertices->at(2) ) *
           vertices->at(2)->distance( vertices->at(3) ) *
           vertices->at(1)->distance( vertices->at(3) ) / this->getArea();
}

void Triangle :: computeBarycentrCoor(FloatArray &answer) {
    double c = vertices->at(1)->distance( vertices->at(2) );
    double a = vertices->at(2)->distance( vertices->at(3) );
    double b = vertices->at(1)->distance( vertices->at(3) );

    // just to avoid mutliple multiplication
    double aPow = a * a;
    double bPow = b * b;
    double cPow = c * c;

    answer.resize(3);
    answer.at(1) = aPow * ( ( -1 ) * aPow + bPow + cPow );
    answer.at(2) = bPow * ( aPow - bPow + cPow );
    answer.at(3) = cPow * ( aPow + bPow - cPow );
}

void Triangle :: computeCenterOfCircumCircle(FloatArray &answer) {
    FloatArray bar;
    this->computeBarycentrCoor(bar);
    double sum = bar.at(1) + bar.at(2) + bar.at(3);
    // center of the circumcircle
    answer.resize(2);
    for ( int i = 1; i <= answer.giveSize(); i++ ) {
        answer.at(i) = ( bar.at(1) * vertices->at(1)->at(i) + bar.at(2) * vertices->at(2)->at(i) + bar.at(3) * vertices->at(3)->at(i) ) / sum;
    }
}

void Triangle :: printYourself() {
    for ( int i = 1; i <= vertices->giveSize(); i++ ) {
        vertices->at(i)->printYourself();
    }
}

bool Triangle :: isOrientedAnticlockwise() {
    FloatMatrix fm(3, 2);
    for ( int i = 1; i <= fm.giveNumberOfRows(); i++ ) {
        for ( int j = 1; j <= fm.giveNumberOfColumns(); j++ ) {
            fm.at(i, j) = this->giveVertex(i)->at(j);
        }
    }

    FloatMatrix fm2(3, 1);
    for ( int i = 1; i <= fm.giveNumberOfRows(); i++ ) {
        fm2.at(i, 1) = 1;
    }

    fm.addSubMatrix(fm2, 1, 3);
    if ( fm.giveDeterminant() > 0.0001 ) {
        return true;
    } else {
        return false;
    }
}

void Triangle :: changeToAnticlockwise() {
    FloatArray *p2e = new FloatArray( *this->giveVertex(3) );
    FloatArray *p3e = new FloatArray( *this->giveVertex(2) );
    this->vertices->remove(2);
    this->vertices->remove(3);
    this->vertices->put(2, p2e);
    this->vertices->put(3, p3e);
}

Circle :: Circle(FloatArray *center, double radius) {
    this->vertices->put(1, center);
    this->radius = radius;
}

double Circle :: computeDistanceTo(FloatArray *point) {
    return vertices->at(1)->distance(point) - radius;
}

IRResultType Circle :: initializeFrom(InputRecord *ir) {
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    FloatArray *center = new FloatArray(2);
    IR_GIVE_FIELD(ir, * center, IFT_Circle_center, "center"); // Macro
    IR_GIVE_FIELD(ir, radius, IFT_Circle_radius, "radius"); // Macro
    vertices->put(1, center);
    return IRRT_OK;
}


