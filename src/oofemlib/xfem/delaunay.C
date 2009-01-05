#include "delaunay.h"
#include "flotarry.h"
#include "intarray.h"
#include "alist.h"
#include <math.h>

Triangle::Triangle(FloatArray* p1, FloatArray* p2, FloatArray* p3){
    vertices = new AList<FloatArray>(3);
    vertices->put(1,p1);
    vertices->put(2,p2);
    vertices->put(3,p3);
}

Triangle::~Triangle(){
    for (int i = 1; i <= vertices->giveSize(); i++){
        // this is done so that the vertices are not destroyed here, since they are shared
        vertices->unlink(i);
    }
    delete vertices;
}

double Triangle::getArea(){
    return fabs(0.5 * (vertices->at(1)->at(1)*(vertices->at(2)->at(2) - vertices->at(3)->at(2))
            + vertices->at(2)->at(1)*(vertices->at(3)->at(2) - vertices->at(1)->at(2)) +
            vertices->at(3)->at(1)*(vertices->at(1)->at(2) - vertices->at(2)->at(2))));
}

double Triangle::getRadiusOfCircumCircle(){
    return 0.25 * vertices->at(1)->distance(vertices->at(2)) * 
                  vertices->at(2)->distance(vertices->at(3)) *
                  vertices->at(1)->distance(vertices->at(3)) / this->getArea();
}

void Triangle::computeBarycentrCoor(FloatArray& answer){
    double c = vertices->at(1)->distance(vertices->at(2));
    double a = vertices->at(2)->distance(vertices->at(3));
    double b = vertices->at(1)->distance(vertices->at(3));
    
    // just to avoid mutliple multiplication
    double aPow = a*a;
    double bPow = b*b;
    double cPow = c*c;
    
    answer.resize(3);
    answer.at(1) = aPow * ((-1) * aPow + bPow + cPow);
    answer.at(2) = bPow * (aPow - bPow + cPow);
    answer.at(3) = cPow * (aPow + bPow - cPow);
}

void Triangle::computeCenterOfCircumCircle(FloatArray& answer){
    FloatArray bar;
    this->computeBarycentrCoor(bar);
    double sum = bar.at(1) + bar.at(2) + bar.at(3);
    // center of the circumcircle
    answer.resize(2);
    for(int i = 1; i <= answer.giveSize(); i++){
        answer.at(i) = (bar.at(1) * vertices->at(1)->at(i) + bar.at(2) * vertices->at(2)->at(i) + bar.at(3) * vertices->at(3)->at(i)) / sum;
    }
}

void Triangle::printYourself(){
    for(int i = 1; i <= vertices->giveSize(); i++){
        vertices->at(i)->printYourself();
    }
}
AList<FloatArray>* Triangle::getVertices(){
    return vertices;
}

bool Delaunay::colinear(FloatArray *p1, FloatArray *p2, FloatArray *p3) {
    double dist = p1->at(1)*(p2->at(2) - p3->at(2)) + p2->at(1)*(p3->at(2) - p1->at(2)) +
            p3->at(1)*(p1->at(2) - p2->at(2));
    // the tolerance probably needs a setter
    if ( dist < 0.0001 && dist > (-1)*0.0001) {
        return true;
    } else return false;
}

void Delaunay::printTriangles(AList<Triangle>* triangles) {
    for (int i = 1; i <= triangles->giveSize(); i++)
        triangles->at(i)->printYourself();
}

bool Delaunay::isInsideCC(FloatArray *p,  FloatArray* p1, FloatArray *p2, FloatArray *p3) {
    // he should not get destroyed, because he stores the coordinates
    Triangle tr(p1,p2,p3);
    double r = tr.getRadiusOfCircumCircle();
    FloatArray circumCenter;
    tr.computeCenterOfCircumCircle(circumCenter);
    double distance = circumCenter.distance(p);
    if (distance < r) return true;
    else return false;

}

void Delaunay::triangulate(AList<FloatArray>* vertices, AList<Triangle>* triangles){
    /// 4th order algorithm - four loops, only for testing puropses
    int n = vertices->giveSize();
    int count = 0;
    /// small shift of vertices
    for (int i = 1; i <= n; i++) {
        vertices->at(i)->at(1) += vertices->at(i)->at(1)*0.00001*double(rand())/RAND_MAX;
        vertices->at(i)->at(2) += vertices->at(i)->at(2)*0.00001*double(rand())/RAND_MAX;
    }
    for (int i = 1; i <= n; i++) {
        for (int j = i + 1; j <= n; j++) {
            for (int k = j + 1; k <= n; k++) {

                bool isTriangle = true;
                if (colinear(vertices->at(i), vertices->at(j),
                        vertices->at(k))) {
                    isTriangle = false;
                } else {
                    for (int a = 1; a <= n; a++) {
                        if (a != i && a != j && a != k) {
                            // checks whether a point a is inside a circumcircle of a triangle ijk
                            if (isInsideCC(vertices->at(a), vertices->at(i), vertices->at(j),
                                    vertices->at(k))) {
                                isTriangle = false;
                                break;
                            }
                        }
                    }
                }

                if (isTriangle) {
                    count++;
                    Triangle* triangle = new Triangle(vertices->at(i), vertices->at(j),
                            vertices->at(k));
                    triangles->put(count, triangle);
                }

            }
        }
    }
}

