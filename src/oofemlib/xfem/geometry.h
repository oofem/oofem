/* 
 * File:   geometry.h
 * Author: chamrova
 *
 * Created on October 26, 2008, 12:04 PM
 */

#ifndef _GEOMETRY_H
#define	_GEOMETRY_H

#include "femcmpnn.h"
#include "domain.h"
#include "flotarry.h"

class Geometry : public FEMComponent
{
public:
         // constructor
         Geometry (int n,Domain* aDomain) : FEMComponent (n, aDomain) {}
          // computes normal signed distance between a line and a point
         virtual double giveDistanceTo (FloatArray *point) = 0;
         // pure virtual
         IRResultType initializeFrom (InputRecord* ir);
         // pure virtual
         const char* giveClassName () const { return "Geometry" ; }
         // checks whether the element is interacted
         bool interacts(Element* element);
         // gives number of intersection points of a geometry entity with an element
         virtual int giveNumberOfIntersectionPoints(Element* element) = 0;
         // gives intersection points
         virtual void giveIntersectionPoints(Element* element, AList<FloatArray>* intersecPoints) = 0;
         // gives intersection points
         virtual void giveIntersectionPoints(Element* element, FloatArray ** intersecPoints) = 0;
};

class Line : public Geometry
{
protected:
         FloatArray* pointA;
         FloatArray* pointB;
public:
        /* temporary constructor, in the future the points will be set
        by a setter function */
         Line (int n, Domain* aDomain, FloatArray* pointA, FloatArray* pointB) ; 
        ~Line();
         double giveDistanceTo (FloatArray *point);
         int giveNumberOfIntersectionPoints(Element* element);
         void giveIntersectionPoints(Element* element, AList<FloatArray>* intersecPoints);
          // computes normal signed distance between a line and a point
         void giveIntersectionPoints(Element* element, FloatArray ** intersecPoints);
};
#endif	/* _GEOMETRY_H */

