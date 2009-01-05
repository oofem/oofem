/* 
 * File:   delaunay.h
 * Author: chamrova
 *
 * Created on October 28, 2008, 11:56 AM
 */

#ifndef _DELAUNAY_H
#define	_DELAUNAY_H

#include "flotarry.h"
#include "alist.h"

class Triangle {
    protected:
        AList<FloatArray>* vertices;
    public:
        Triangle(FloatArray* p1, FloatArray* p2, FloatArray* p3);
        ~Triangle();
        double getArea();
        void computeBarycentrCoor(FloatArray& answer);
        double getRadiusOfCircumCircle();
        void computeCenterOfCircumCircle(FloatArray& answer);
        void printYourself();
        AList<FloatArray>* getVertices();
};

// O(n4) algorithm, only for testing purposes

class Delaunay {
public:
    bool colinear(FloatArray *p1, FloatArray *p2, FloatArray *p3);
    void printTriangles(AList<Triangle>* triangles);
    bool isInsideCC(FloatArray *p,  FloatArray* p1, FloatArray *p2, FloatArray *p3);
    void triangulate(AList<FloatArray>* vertices, AList<Triangle>* triangles);
};


#endif	/* _DELAUNAY_H */

