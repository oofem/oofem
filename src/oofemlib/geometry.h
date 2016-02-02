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

#include "oofemcfg.h"
#include "error.h"
#include "floatarray.h"
#include "inputrecord.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#include <list>
#ifdef __BOOST_MODULE
 #include <BoostInterface.h>
#endif

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
#define _IFT_PolygonLine_points "points"

//@}

namespace oofem {

class Element;
class DynamicInputRecord;
class oofegGraphicContext;
class TipInfo;

/**
 * Abstract representation of Geometry
 * @author chamrova
 * @author Erik Svenning
 */
class OOFEM_EXPORT BasicGeometry //: public Geometry
{
protected:
    /// List of geometry vertices.
    std :: vector< FloatArray >mVertices;
public:
    /// Constructor.
    BasicGeometry();

    /// Copy constructor: should be implemented when a class deals with pointers
    BasicGeometry(const BasicGeometry & iBasicGeometry);

    /// Destructor.
    virtual ~BasicGeometry();

    virtual BasicGeometry *Clone() = 0;

    /// Computes normal signed distance between this object and a point.
    virtual double computeDistanceTo(const FloatArray *point) { return 0; }


    // For debugging
    virtual void printVTK(int iTStepIndex, int iIndex) {};

    /// Functions for computing signed distance in normal and tangential direction.
    /// Used by XFEM level set functions.
    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const = 0;
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const = 0;
    virtual void computeLocalCoordinates(FloatArray &oLocCoord, const FloatArray &iPoint) const { OOFEM_ERROR("not implemented."); }
    virtual void giveSubPolygon(std :: vector< FloatArray > &oPoints, const double &iXiStart, const double &iXiEnd) const { OOFEM_ERROR("?"); }
    virtual void giveGlobalCoordinates(FloatArray &oGlobalCoord, const double &iArcPos) const {OOFEM_ERROR("Not implemented.")};

    /// Computes tangential direction at given local coordinate (arcPos)
    virtual void giveTangent(FloatArray &oTangent, const double &iArcPosition) const {printf("BasicGeometry::giveTangent() not implemented.\n");}


    /// Checks whether an element is interacted, Element reference will be later replaced by Geometry.
    virtual bool intersects(Element *element) { return false; }
    /// Gives number of intersection points of Geometry entity with an element, Element reference will be later replaced by Geometry.
    virtual int computeNumberOfIntersectionPoints(Element *element) { return 0; }
    /// Gives intersection points between this Geometry and Element.
    virtual void computeIntersectionPoints(Element *element, std :: vector< FloatArray > &oIntersectionPoints) { }

    inline const FloatArray &giveVertex(int n) const { return mVertices [ n - 1 ]; }

    void setVertices(const std::vector<FloatArray> &iVertices) {mVertices = iVertices;}

    void removeDuplicatePoints(const double &iTolSquare);

    void insertVertexFront(const FloatArray &iP) { mVertices.insert(mVertices.begin(), iP); }
    void insertVertexBack(const FloatArray &iP) { mVertices.push_back(iP); }


    /// Initializes the Geometry from the InputRecord.
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    virtual void giveInputRecord(DynamicInputRecord &input) { OOFEM_ERROR("not implemented"); }
    /// Gives class name.
    virtual const char *giveClassName() const { return NULL; }
    std :: string errorInfo(const char *func) const { return std :: string(giveClassName()) + func; }
    /// Returns number of Geometry vertices.
    int giveNrVertices() const { return (int)mVertices.size(); }
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
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL) { return CIO_OK; }
    /**
     * Restores the state of receiver from output stream.
     * @param stream Context file.
     * @param mode Determines amount of info in stream.
     * @param obj Special parameter for sending extra information.
     * @return contextIOResultType.
     * @exception ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL) { return CIO_OK; }

#ifdef __OOFEG
    virtual void draw(oofegGraphicContext &gc) { }
#endif


    /**
     * Computes the distance between two lines.
     * Line 1 has start point iP1 and end point iP2.
     * Line 2 has start point iQ1 and end point iQ2.
     */
    static double computeLineDistance(const FloatArray &iP1, const FloatArray &iP2, const FloatArray &iQ1, const FloatArray &iQ2);

    /**
     * Returns start and end tip of the geometry, if applicable.
     */
    virtual bool giveTips(TipInfo &oStartTipInfo, TipInfo &oEndTipInfo) const {return false;}

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius) {OOFEM_ERROR("Not implemented.")};

};

class OOFEM_EXPORT Line : public BasicGeometry
{
public:
    Line() : BasicGeometry() { }
    virtual ~Line() { }
    Line(const FloatArray &iPointA, const FloatArray &iPointB);

    virtual BasicGeometry *Clone() { return new Line(*this); }

    virtual double computeDistanceTo(const FloatArray *point);
    /// Computes tangential distance to a point

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("not implemented"); }

    double computeTangentialDistanceToEnd(FloatArray *point);

    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const;

    void computeProjection(FloatArray &answer);
    virtual int computeNumberOfIntersectionPoints(Element *element);
    virtual void computeIntersectionPoints(Element *element, std :: vector< FloatArray > &oIntersectionPoints);
    double computeInclinationAngle();
    void computeTransformationMatrix(FloatMatrix &answer);
    void transformIntoPolar(FloatArray *point, FloatArray &answer);
    virtual IRResultType initializeFrom(InputRecord *ir);
    bool isPointInside(FloatArray *point);
    virtual bool intersects(Element *element);
    virtual bool isOutside(BasicGeometry *bg);

    double giveLength() const {return mVertices[0].distance( mVertices[1] );}
};

class OOFEM_EXPORT Triangle : public BasicGeometry
{
public:
    Triangle(const FloatArray & iP1, const FloatArray & iP2, const FloatArray & iP3);
    virtual ~Triangle() { }

    virtual BasicGeometry *Clone() { return new Triangle(*this); }

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("not implemented"); }
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const { OOFEM_ERROR("not implemented"); }

    double getArea();
    void computeBarycentrCoor(FloatArray &answer) const;
    double getRadiusOfCircumCircle();
    void computeCenterOfCircumCircle(FloatArray &answer) const;
    virtual void printYourself();
    virtual int computeNumberOfIntersectionPoints(Element *element) { return 0; }
    bool isOrientedAnticlockwise();
    void changeToAnticlockwise();

    /**
     * Checks if the projection of the the point iP onto the
     * triangle plane is inside the triangle.
     * @author Erik Svenning
     */
    bool pointIsInTriangle(const FloatArray &iP) const;

    /**
     * Split a triangle in four.
     */
    static void refineTriangle(std::vector<Triangle> &oRefinedTri, const Triangle &iTri);
};

class OOFEM_EXPORT Circle : public BasicGeometry
{
protected:
    double radius;
    const double mTangSignDist;
public:
    Circle() : BasicGeometry(), radius(0.0), mTangSignDist(1.0) { }
    virtual ~Circle() { }
    Circle(FloatArray &center, double radius);

    virtual BasicGeometry *Clone() { return new Circle(*this); }

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const;

    // Irrelevant for a closed interface: we can always consider ourselves to be "inside" a closed interface in
    // tangential direction. Therefore, we may return any positive number.
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const { oDist = mTangSignDist; }

    virtual void giveGlobalCoordinates(FloatArray &oGlobalCoord, const double &iArcPos) const;

    virtual void giveTangent(FloatArray &oTangent, const double &iArcPosition) const { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "Circle"; }
    virtual bool intersects(Element *element);
    virtual void computeIntersectionPoints(Element *element, std :: vector< FloatArray > &oIntersectionPoints);
    virtual void computeIntersectionPoints(Line *l, std :: vector< FloatArray > &oIntersectionPoints);
    virtual int computeNumberOfIntersectionPoints(Element *element);
    virtual bool isOutside(BasicGeometry *bg);
    virtual bool isInside(Element *element);
    virtual bool isInside(FloatArray &point);
    virtual void printYourself();

    double giveRadius() const {return radius;}

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius);

};

class OOFEM_EXPORT PolygonLine : public BasicGeometry
{
    bool mDebugVtk;
public:
    PolygonLine();
    virtual ~PolygonLine() { }

    virtual BasicGeometry *Clone() { return new PolygonLine(*this); }

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const;
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const;

    /// Computes arc length coordinate in the range [0,1]
    virtual void computeLocalCoordinates(FloatArray &oLocCoord, const FloatArray &iPoint) const;
    double computeLength() const;

    virtual void giveSubPolygon(std :: vector< FloatArray > &oPoints, const double &iXiStart, const double &iXiEnd) const;
    virtual void giveGlobalCoordinates(FloatArray &oGlobalCoord, const double &iArcPos) const;
    void giveNormal(FloatArray &oNormal, const double &iArcPosition) const;
    virtual void giveTangent(FloatArray &oTangent, const double &iArcPosition) const;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveClassName() const { return "PolygonLine"; }

#ifdef __BOOST_MODULE
    virtual bool boundingBoxIntersects(Element *element);
#endif

    virtual bool intersects(Element *element);
    virtual void computeIntersectionPoints(Element *element, std :: vector< FloatArray > &oIntersectionPoints);
    virtual void computeIntersectionPoints(Line *l, std :: vector< FloatArray > &oIntersectionPoints);
    void computeIntersectionPoints(const PolygonLine &iPolygonLine, std :: vector< FloatArray > &oIntersectionPoints) const;
    void computeIntersectionPoints(const FloatArray &iXStart, const FloatArray &iXEnd, std :: vector< FloatArray > &oIntersectionPoints) const;

    virtual int computeNumberOfIntersectionPoints(Element *element);
    virtual bool isOutside(BasicGeometry *bg);
    virtual bool isInside(Element *element);
    virtual bool isInside(FloatArray &point);

#ifdef __BOOST_MODULE
    virtual void calcBoundingBox(bPoint2 &oLC, bPoint2 &oUC);
#endif

    virtual void printYourself();

    // For debugging
    virtual void printVTK(int iTStepIndex, int iLineIndex);

#ifdef __BOOST_MODULE
    // Upper and lower corner
    bPoint2 LC, UC;
#endif

    virtual bool giveTips(TipInfo &oStartTipInfo, TipInfo &oEndTipInfo) const;

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius);

    /**
     * Keep only a part of the underlying geometry,
     * characterized by iArcPosStart and iArcPosEnd.
     */
    void cropPolygon(const double &iArcPosStart, const double &iArcPosEnd);

};


class OOFEM_EXPORT PointSwarm : public BasicGeometry
{
protected:
    std :: list< int >idList;
public:
    PointSwarm() : BasicGeometry() { }
    virtual ~PointSwarm() { }
    PointSwarm(std :: list< int >pointsID);

    virtual BasicGeometry *Clone() { return new PointSwarm(*this); }

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("not implemented"); }
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const { OOFEM_ERROR("not implemented"); }

    /// Computes the normal distance to the surface not to the center.
    // virtual double computeDistanceTo(FloatArray *point);
    virtual IRResultType initializeFrom(InputRecord *ir);
    // virtual const char *giveClassName() const { return "Circle"; }
    // virtual bool intersects(Element *element);
    // virtual void computeIntersectionPoints(Element *element, std :: vector< FloatArray > *intersecPoints);
    // virtual void computeIntersectionPoints(Line *l, std :: vector< FloatArray > *intersecPoints);
    // virtual bool isOutside(BasicGeometry *bg);
    // virtual bool isInside(Element *element);
    // virtual bool isInside(FloatArray &point);
};
} // end namespace oofem
#endif  // geometry_h
