/*
 * File:   geometry.h
 * Author: chamrova
 *
 * Created on October 26, 2008, 12:04 PM
 */

#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include "domain.h"
#include "flotarry.h"
#include "node.h"
#include <iostream>

/** Abstract representation for Geometry
 *  gives basic virtual methods for accessing the nodes
 */

class Geometry
{
public:
    /// Constructor
    Geometry() {};
    /// Destructor
    ~Geometry() {};
    /// returns number of nodes
    virtual int giveNumberOfNodes() { return 0; }
    /// returns coordinates at a position n
    virtual FloatArray *giveCoordinates(int n) { return NULL; }
    virtual void printYourself() {}
};

/** Concrete representation of Geometry
 * Element inherits from this class
 */
class ElementGeometry : public Geometry
{
protected:
    IntArray nodeNumbers;
public:
    ElementGeometry() : Geometry() {}
};

/** Concrete representation of Geometry
 * Patch inherits from this class
 */
class BasicGeometry : public Geometry
{
protected:
    /// List of geometry vertices
    AList< FloatArray > *vertices;
public:
    /// Constructor
    BasicGeometry();
    /// Destructor
    ~BasicGeometry();
    /// computes normal signed distance between this object and a point
    virtual double computeDistanceTo(FloatArray *point) { return 0; }
    /// checks whether an element is interacted, Element reference will be later replaced by Geometry
    virtual bool intersects(Element *element) {return NULL;};
    /// gives number of intersection points of Geometry entity with an element, Element reference will be later replaced by Geometry
    virtual int computeNumberOfIntersectionPoints(Element *element) { return 0; }
    // gives intersection points between this Geometry and Element
    virtual void computeIntersectionPoints(Element *element, AList< FloatArray > *intersecPoints) {}
    /// Accessor
    FloatArray *giveVertex(int n);
    /// Modifier
    void setVertex(FloatArray *vertex);
    /// Accessor
    AList< FloatArray > *giveVertices() { return this->vertices; }
    /// Initializes the Geometry from the InputRecord
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    /// gives class name
    virtual const char *giveClassName() const { return NULL; }
    /// returns number of Geometry vertices
    int giveNrVertices() { return this->vertices->giveSize(); }
    virtual bool isOutside(BasicGeometry *bg) {return NULL; }
};

class Line : public BasicGeometry
{
public:
    Line() : BasicGeometry() {}
    ~Line() {}
    Line(FloatArray *pointA, FloatArray *pointB);
    /// computes normal distance to a point
    double computeDistanceTo(FloatArray *point);
    /// computes tangential distance to a point
    double computeTangentialDistanceToEnd(FloatArray *point);
    void computeProjection(FloatArray &answer);
    int computeNumberOfIntersectionPoints(Element *element);
    void computeIntersectionPoints(Element *element, AList< FloatArray > *intersecPoints);
    double computeInclinationAngle();
    void computeTransformationMatrix(FloatMatrix &answer);
    void transformIntoPolar(FloatArray *point, FloatArray &answer);
    IRResultType initializeFrom(InputRecord *ir);
    bool isPointInside(FloatArray *point);
    bool intersects(Element *element);
    bool isOutside(BasicGeometry *bg);
};

class Triangle : public BasicGeometry
{
public:
    Triangle(FloatArray *p1, FloatArray *p2, FloatArray *p3);
    double getArea();
    void computeBarycentrCoor(FloatArray &answer);
    double getRadiusOfCircumCircle();
    void computeCenterOfCircumCircle(FloatArray &answer);
    void printYourself();
    int computeNumberOfIntersectionPoints(Element *element) { return 0; }
    bool isOrientedAnticlockwise();
    void changeToAnticlockwise();
};

class Circle : public BasicGeometry
{
protected:
    double radius;
public:
    Circle() : BasicGeometry() {}
    ~Circle() {}
    Circle(FloatArray *center, double radius);
    // normal distance to the surface not to the centre
    double computeDistanceTo(FloatArray *point);
    IRResultType initializeFrom(InputRecord *ir);
    const char *giveClassName() const { return "Circle"; }
    bool intersects(Element *element);
    void computeIntersectionPoints(Element *element, AList< FloatArray > *intersecPoints);
    void computeIntersectionPoints(Line *l, AList< FloatArray > *intersecPoints);
    bool isOutside(BasicGeometry *bg);
    void printYourself();
};
#endif  /* _GEOMETRY_H */




