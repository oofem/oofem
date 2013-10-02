/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef geometry_h
#define geometry_h

#include "domain.h"
#include "floatarray.h"
#include "node.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#ifdef __BOOST_MODULE
#include <BoostInterface.h>
#endif

namespace oofem {

///@name Input fields for geometries
//@{
#define _IFT_Circle_Name "circle"
#define _IFT_Circle_radius "radius"
#define _IFT_Circle_center "center"

#define _IFT_Line_Name "line"
#define _IFT_Line_start "start"
#define _IFT_Line_end "end"

#define _IFT_PointSwarm_Name "pointswarm" // just temporary
#define _IFT_PointSwarm_nodeID "nodeid"


#define _IFT_PolygonLine_Name "polygonline"
//#define _IFT_PolygonLine_start "start"
//#define _IFT_PolygonLine_end "end"
#define _IFT_PolygonLine_points "points"

//@}

/**
 * Abstract representation of Geometry
 * @author chamrova
 * @author Erik Svenning
 */
class BasicGeometry //: public Geometry
{
protected:
    /// List of geometry vertices.
	// AList does not provide elementary operations like insert,
	// therefore use std::vector instead.
//    AList< FloatArray > *vertices;
	std::vector< FloatArray > mVertices;
public:
    /// Constructor.
    BasicGeometry();

    /// Copy constructor: should be implemented when a class deals with pointers
    BasicGeometry(const BasicGeometry &iBasicGeometry);

    /// Destructor.
    virtual ~BasicGeometry();
    /// Computes normal signed distance between this object and a point.
    virtual double computeDistanceTo(const FloatArray *point) { return 0; }

    virtual double computeTangentialSignDist(FloatArray *point) { return 0; }

    /// Functions for computing signed distance in normal and tangential direction.
    /// Used by XFEM level set functions.
    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const = 0;
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const = 0;


    /// Checks whether an element is interacted, Element reference will be later replaced by Geometry.
    virtual bool intersects(Element *element) { return false; }
    /// Gives number of intersection points of Geometry entity with an element, Element reference will be later replaced by Geometry.
    virtual int computeNumberOfIntersectionPoints(Element *element) { return 0; }
    /// Gives intersection points between this Geometry and Element.
    virtual void computeIntersectionPoints(Element *element, std::vector< FloatArray > &oIntersectionPoints) { }

    /// Accessor.
//    FloatArray *giveVertex(int n);
    const FloatArray &giveVertex(int n) const {return mVertices[n-1];}

    /// Modifier.
    void setVertex(FloatArray *vertex);

    void insertVertexFront(const FloatArray &iP) {mVertices.insert( mVertices.begin(), iP );}
    void insertVertexBack(const FloatArray &iP) {mVertices.push_back(iP);}


    /// Accessor.
//    AList< FloatArray > *giveVertices() { return this->vertices; }
    /// Initializes the Geometry from the InputRecord.
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    /// Gives class name.
    virtual const char *giveClassName() const { return NULL; }
    /**
     * Returns classType id of receiver. Intended for run time
     * type checking. Every derived class have to overload this method.
     * @see classType.
     * @return Class type of receiver.
     */
    virtual classType giveClassID() const { return BasicGeometryClass; }
    /// Returns number of Geometry vertices.
    int giveNrVertices() const { return mVertices.size(); }
    virtual bool isOutside(BasicGeometry *bg) { return false; }
    virtual bool isInside(Element *el) { return false; }
    virtual bool isInside(FloatArray &point) { return false; }
    virtual void printYourself() { }
    /**
     * Stores the state of receiver to output stream.
     * @param stream Context stream.
     * @param mode Determines amount of info in stream.
     * @param obj Special parameter, used to pass optional parameters.
     * @return contextIOResultType.
     * @exception ContextIOERR If error encountered.
     */
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL) {return CIO_OK;}
    /**
     * Restores the state of receiver from output stream.
     * @param stream Context file.
     * @param mode Determines amount of info in stream.
     * @param obj Special parameter for sending extra information.
     * @return contextIOResultType.
     * @exception ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL) {return CIO_OK;}

#ifdef __OOFEG
    virtual void draw(oofegGraphicContext &gc) { }
#endif
};

class Line : public BasicGeometry
{
public:
    Line() : BasicGeometry() { }
    virtual ~Line() { }
    Line(FloatArray *pointA, FloatArray *pointB);

    virtual double computeDistanceTo(const FloatArray *point);
    /// Computes tangential distance to a point

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const {OOFEM_ERROR("Line::computeNormalSignDist -- not implemented");};
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const {OOFEM_ERROR("Line::computeTangentialSignDist -- not implemented");};


    double computeTangentialDistanceToEnd(FloatArray *point);
    void computeProjection(FloatArray &answer);
    virtual int computeNumberOfIntersectionPoints(Element *element);
    virtual void computeIntersectionPoints(Element *element, std::vector< FloatArray > &oIntersectionPoints);
    double computeInclinationAngle();
    void computeTransformationMatrix(FloatMatrix &answer);
    void transformIntoPolar(FloatArray *point, FloatArray &answer);
    virtual IRResultType initializeFrom(InputRecord *ir);
    bool isPointInside(FloatArray *point);
    virtual bool intersects(Element *element);
    virtual bool isOutside(BasicGeometry *bg);
};

class Triangle : public BasicGeometry
{
public:
    Triangle(FloatArray *p1, FloatArray *p2, FloatArray *p3);
    virtual ~Triangle() { }

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const {OOFEM_ERROR("Triangle::computeNormalSignDist -- not implemented");};
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const {OOFEM_ERROR("Triangle::computeTangentialSignDist -- not implemented");};

    double getArea();
    void computeBarycentrCoor(FloatArray &answer);
    double getRadiusOfCircumCircle();
    void computeCenterOfCircumCircle(FloatArray &answer);
    virtual void printYourself();
    virtual int computeNumberOfIntersectionPoints(Element *element) { return 0; }
    bool isOrientedAnticlockwise();
    void changeToAnticlockwise();
};

class Circle : public BasicGeometry
{
protected:
    double radius;
    const double mTangSignDist;
public:
    Circle() : BasicGeometry(), radius(0.0), mTangSignDist(1.0) { }
    virtual ~Circle() { }
    Circle(FloatArray *center, double radius);
    /// Computes the normal distance to the surface not to the center.
    virtual double computeDistanceTo(const FloatArray *point);

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const;

    // Irrelevant for a closed interface: we can always consider ourselves to be "inside" a closed interface in
    // tangential direction. Therefore, we may return any positive number.
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const {oDist = mTangSignDist;};

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "Circle"; }
    virtual bool intersects(Element *element);
    virtual void computeIntersectionPoints(Element *element, std::vector< FloatArray > &oIntersectionPoints);
    virtual void computeIntersectionPoints(Line *l, std::vector< FloatArray > &oIntersectionPoints);
    virtual int computeNumberOfIntersectionPoints(Element *element);
    virtual bool isOutside(BasicGeometry *bg);
    virtual bool isInside(Element *element);
    virtual bool isInside(FloatArray &point);
    virtual void printYourself();
};

class PolygonLine : public BasicGeometry
{
    static int nextLineIdNumber;
    int stepInd;
public:
	PolygonLine();
    virtual ~PolygonLine() { }
    /// Computes the normal distance to the surface not to the center.
    virtual double computeDistanceTo(const FloatArray *point);

    virtual double computeTangentialSignDist(FloatArray *point);

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const;
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "PolygonLine"; }

#ifdef __BOOST_MODULE
    virtual bool boundingBoxIntersects(Element *element);
#endif

    virtual bool intersects(Element *element);
    virtual void computeIntersectionPoints(Element *element, std::vector< FloatArray > &oIntersectionPoints);
    virtual void computeIntersectionPoints(Line *l, std::vector< FloatArray > &oIntersectionPoints);
    virtual int computeNumberOfIntersectionPoints(Element *element);
    virtual bool isOutside(BasicGeometry *bg);
    virtual bool isInside(Element *element);
    virtual bool isInside(FloatArray &point);

#ifdef __BOOST_MODULE
    virtual void calcBoundingBox(bPoint2 &oLC, bPoint2 &oUC);
#endif

    virtual void printYourself();

    // For debugging
    virtual void printVTK();

    // Id for writing VTK
    int lineIdNumber;

#ifdef __BOOST_MODULE
    // Upper and lower corner
    bPoint2 LC, UC;
#endif

};


class PointSwarm : public BasicGeometry
{
protected:
    std::list< int > idList;
public:
    PointSwarm() : BasicGeometry() { }
    virtual ~PointSwarm() { }
    PointSwarm(std::list<int> pointsID);

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const {OOFEM_ERROR("PointSwarm::computeNormalSignDist -- not implemented");};
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const {OOFEM_ERROR("PointSwarm::computeTangentialSignDist -- not implemented");};

    /// Computes the normal distance to the surface not to the center.
   // virtual double computeDistanceTo(FloatArray *point);
   virtual IRResultType initializeFrom(InputRecord *ir);
   // virtual const char *giveClassName() const { return "Circle"; }
   // virtual bool intersects(Element *element);
   // virtual void computeIntersectionPoints(Element *element, AList< FloatArray > *intersecPoints);
   // virtual void computeIntersectionPoints(Line *l, AList< FloatArray > *intersecPoints);
   // virtual bool isOutside(BasicGeometry *bg);
   // virtual bool isInside(Element *element);
   // virtual bool isInside(FloatArray &point);
    
};

} // end namespace oofem
#endif  // geometry_h




