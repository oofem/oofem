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

#include "mathfem.h"
#include "geometry.h"
#include "element.h"
#include "dofmanager.h"
#include "classfactory.h"
#include "inputrecord.h"
#include "dynamicinputrecord.h"
#include "xfem/xfemtolerances.h"
#include "feinterpol.h"
#include "xfem/tipinfo.h"


#include <fstream>
#include <limits>

namespace oofem {
REGISTER_Geometry(Line)
REGISTER_Geometry(Circle)
REGISTER_Geometry(PointSwarm)
REGISTER_Geometry(PolygonLine)

BasicGeometry :: BasicGeometry()
{ }

BasicGeometry :: BasicGeometry(const BasicGeometry &iBasicGeometry) :
    mVertices(iBasicGeometry.mVertices)
{ }

BasicGeometry :: ~BasicGeometry()
{ }


void  BasicGeometry :: removeDuplicatePoints(const double &iTolSquare)
{
    if ( mVertices.size() > 1 ) {
        for ( size_t i = mVertices.size()-1; i > 0; i--) {
            if ( distance_square(mVertices[i], mVertices[i-1]) < iTolSquare ) {
                mVertices.erase( mVertices.begin()+i );
            }
        }
    }
}

void BasicGeometry :: translate(const FloatArray &iTrans)
{
    for ( size_t i = 0; i < mVertices.size(); i++ ) {
        mVertices[i].add(iTrans);
    }
}


double BasicGeometry :: computeLineDistance(const FloatArray &iP1, const FloatArray &iP2, const FloatArray &iQ1, const FloatArray &iQ2)
{
    FloatArray u;
    u.beDifferenceOf(iP2, iP1);

    const double LengthP = u.computeNorm();

    FloatArray v;
    v.beDifferenceOf(iQ2, iQ1);
    const double LengthQ = v.computeNorm();


    // Regularization coefficients (to make it possible to solve when lines are parallel)
    const double c1 = (1.0e-14)*LengthP*LengthP;
    const double c2 = (1.0e-14)*LengthQ*LengthQ;

    const size_t minIter = 2;
    const size_t maxIter = 5;
    const double absTol = 1.0e-12;

    double xi = 0.0, eta = 0.0;

    FloatArray d;
    d = iP1;
    d.add(xi,u);
    d.add(-1.0, iQ1);
    d.add(-eta, v);

    FloatMatrix K(2,2), KInv;
    FloatArray dXi;

    bool lockXi = false, lockEta = false;

    for ( size_t iter = 0; iter < maxIter; iter++ ) {

        if ( xi < 0.0 ) {
            xi = 0.0;
            lockXi = true;
        }

        if ( xi > 1.0 ) {
            xi = 1.0;
            lockXi = true;
        }


        if ( eta < 0.0 ) {
            eta = 0.0;
            lockEta = true;
        }

        if ( eta > 1.0 ) {
            eta = 1.0;
            lockEta = true;
        }

        FloatArray R = {   d.dotProduct(u) + c1*xi,
                          -d.dotProduct(v) + c2*eta};

        if ( lockXi ) {
            R[0] = 0.0;
        }

        if ( lockEta ) {
            R[1] = 0.0;
        }

        const double res = R.computeNorm();
//        printf("iter: %lu res: %e\n", iter, res);

        if ( res < absTol && iter >= minIter ) {
//            printf("xi: %e eta: %e\n", xi, eta);
            break;
        }

        K(0,0) = -u.dotProduct(u)-c1;
        K(0,1) =  u.dotProduct(v);
        K(1,0) =  u.dotProduct(v);
        K(1,1) = -v.dotProduct(v)-c2;


        if ( lockXi ) {
            K(0,0) = -1.0;
            K(0,1) = K(1,0) = 0.0;
        }

        if ( lockEta ) {
            K(0,1) = K(1,0) = 0.0;
            K(1,1) = -1.0;
        }


        KInv.beInverseOf(K);

        dXi.beProductOf(KInv, R);

        xi  += dXi[0];
        eta += dXi[1];

        d = iP1;
        d.add(xi,u);
        d.add(-1.0, iQ1);
        d.add(-eta, v);
    }

    if ( xi < 0.0 ) {
        xi = 0.0;
    }

    if ( xi > 1.0 ) {
        xi = 1.0;
    }

    if ( eta < 0.0 ) {
        eta = 0.0;
    }

    if ( eta > 1.0 ) {
        eta = 1.0;
    }

    d = iP1;
    d.add(xi,u);
    d.add(-1.0, iQ1);
    d.add(-eta, v);

    const double dist = d.computeNorm();

    return dist;
}

bool Line :: intersects(Element *element)
{
    int ip = this->computeNumberOfIntersectionPoints(element);
    return ( ip > 0 );
}

Line :: Line(const FloatArray &iPointA, const FloatArray &iPointB) : BasicGeometry()
{
    mVertices.push_back(iPointA);
    mVertices.push_back(iPointB);
}

double Line :: computeDistanceTo(const FloatArray &point)
{
	// TODO: Is this function correct?! /ES
    const auto &pointA = mVertices [ 0 ];
    const auto &pointB = mVertices [ 1 ];
    double a = pointA.at(2) - pointB.at(2);
    double b = pointB.at(1) - pointA.at(1);
    double c = pointA.at(1) * pointB.at(2) - pointB.at(1) * pointA.at(2);
    double l = distance(pointA, pointB);
    return ( a * point.at(1) + b * point.at(2) + c ) / l;
}

void Line :: computeProjection(FloatArray &answer)
{
    answer.beDifferenceOf(mVertices [ 1 ], mVertices [ 0 ]);
}

double Line :: computeTangentialDistanceToEnd(const FloatArray &point)
{
    FloatArray projection;
    this->computeProjection(projection);
    FloatArray tmp;
    tmp.beDifferenceOf(point, mVertices [ 1 ]);
    return tmp.dotProduct(projection) / projection.computeNorm();
}

void Line :: computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinArcDist) const
{
    PolygonLine pl;
    pl.insertVertexBack( mVertices[0] );
    pl.insertVertexBack( mVertices[1] );

    pl.computeTangentialSignDist(oDist, iPoint, oMinArcDist);
}

int Line :: computeNumberOfIntersectionPoints(Element *element)
{

    int numIntersec = 0;
    const double relTol = 1.0e-3;
    const double LineLength = giveLength();
    const double absTol = relTol*std::max(LineLength, XfemTolerances::giveCharacteristicElementLength() );

    const int numEdges = element->giveInterpolation()->giveNumberOfEdges(element->giveGeometryType());

    for ( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ ) {
        auto bNodes = element->giveInterpolation()->boundaryGiveNodes(edgeIndex, element->giveGeometryType());

        const int nsLoc = bNodes.at(1);
        const int neLoc = bNodes.at( bNodes.giveSize() );

        FloatArray xS = element->giveNode(nsLoc)->giveCoordinates();
        xS.resizeWithValues(2);
        FloatArray xE = element->giveNode(neLoc)->giveCoordinates();
        xE.resizeWithValues(2);

        const double dist = BasicGeometry :: computeLineDistance(xS, xE, mVertices[0], mVertices[1]);

        if(dist < absTol) {
            numIntersec++;
        }
    }

    return numIntersec;
}

void Line :: computeIntersectionPoints(Element *element, std :: vector< FloatArray > &oIntersectionPoints)
{
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
        if ( lsn1 * lsn2 <= 0 && lst1 <= 0 && lst2 <= 0 ) {
            double r = lsn1 / ( lsn1 - lsn2 );
            if ( i <= element->giveNumberOfDofManagers() ) {
                FloatArray answer(2);
                for ( int j = 1; j <= answer.giveSize(); j++ ) {
                    answer.at(j) = ( 1 - r ) * element->giveDofManager(n1)->giveCoordinate(j)
                    + r *element->giveDofManager(n2)->giveCoordinate(j);
                }

                oIntersectionPoints.push_back(answer);
            }
        }
    }
}

double Line :: computeInclinationAngle()
{
    const FloatArray &pointA = mVertices [ 0 ];
    const FloatArray &pointB = mVertices [ 1 ];
    double y = pointB.at(2) - pointA.at(2);
    double x = pointB.at(1) - pointA.at(1);
    return atan2(y, x);
}

void Line :: computeTransformationMatrix(FloatMatrix &answer)
{
    answer.resize(2, 2);
    double alpha = this->computeInclinationAngle();
    answer.at(1, 1) =  cos(alpha);
    answer.at(1, 2) =  sin(alpha);
    answer.at(2, 1) = -sin(alpha);
    answer.at(2, 2) =  cos(alpha);
}

void Line :: transformIntoPolar(FloatArray *point, FloatArray &answer)
{
    FloatArray xp;
    FloatMatrix Qt;
    FloatArray help;
    this->computeTransformationMatrix(Qt);
    help.beDifferenceOf(* point, mVertices [ 1 ]);
    xp.beProductOf(Qt, help);
    answer.resize(2);
    answer.at(1) = xp.computeNorm();
    answer.at(2) = atan2( xp.at(2), xp.at(1) );
}

void Line :: initializeFrom(InputRecord &ir)
{
    mVertices.resize(2);
    IR_GIVE_FIELD(ir, mVertices [ 0 ], _IFT_Line_start);
    IR_GIVE_FIELD(ir, mVertices [ 1 ], _IFT_Line_end);
}

bool Line :: isPointInside(FloatArray *point)
{
    double maxX, minX, maxY, minY;
    if ( mVertices [ 0 ].at(1) > mVertices [ 1 ].at(1) ) {
        maxX = mVertices [ 0 ].at(1);
        minX = mVertices [ 1 ].at(1);
    } else {
        minX = mVertices [ 0 ].at(1);
        maxX = mVertices [ 1 ].at(1);
    }

    if ( mVertices [ 0 ].at(2) > mVertices [ 1 ].at(2) ) {
        maxY = mVertices [ 0 ].at(2);
        minY = mVertices [ 1 ].at(2);
    } else {
        minY = mVertices [ 0 ].at(2);
        maxY = mVertices [ 1 ].at(2);
    }

    if ( point->at(1) >= minX && point->at(1) <= maxX &&
         point->at(2) >= minY && point->at(2) <= maxY ) {
        return true;
    } else {
        return false;
    }
}

bool Line :: isOutside(BasicGeometry *bg)
{ // equivalent to up
    int count = 0;
    for ( int i = 1; i <= bg->giveNrVertices(); i++ ) {
        if ( this->computeDistanceTo( bg->giveVertex(i) ) > 0.1 ) {
            count++;
        }
    }

    return ( count != 0 );
}

Triangle :: Triangle(const FloatArray &iP1, const FloatArray &iP2, const FloatArray &iP3) : BasicGeometry()
{
    mVertices.push_back(iP1);
    mVertices.push_back(iP2);
    mVertices.push_back(iP3);
}

double Triangle :: getArea()
{
    return fabs( 0.5 * ( mVertices [ 0 ].at(1) * ( mVertices [ 1 ].at(2) - mVertices [ 2 ].at(2) )
                        + mVertices [ 1 ].at(1) * ( mVertices [ 2 ].at(2) - mVertices [ 0 ].at(2) ) +
                        mVertices [ 2 ].at(1) * ( mVertices [ 0 ].at(2) - mVertices [ 1 ].at(2) ) ) );
}

double Triangle :: getRadiusOfCircumCircle()
{
    return 0.25 * distance(mVertices [ 0 ], mVertices [ 1 ]) *
        distance(mVertices [ 1 ], mVertices [ 2 ]) *
        distance(mVertices [ 0 ], mVertices [ 2 ]) / this->getArea();
}

void Triangle :: computeBarycentrCoor(FloatArray &answer) const
{
    double c = distance(mVertices [ 0 ], mVertices [ 1 ]);
    double a = distance(mVertices [ 1 ], mVertices [ 2 ]);
    double b = distance(mVertices [ 0 ], mVertices [ 2 ]);

    // just to avoid multiple multiplication
    double aPow = a * a;
    double bPow = b * b;
    double cPow = c * c;

    answer.resize(3);
    answer.at(1) = aPow * ( -aPow + bPow + cPow );
    answer.at(2) = bPow * ( aPow - bPow + cPow );
    answer.at(3) = cPow * ( aPow + bPow - cPow );
}

void Triangle :: computeCenterOfCircumCircle(FloatArray &answer) const
{
    FloatArray bar;
    this->computeBarycentrCoor(bar);
    double sum = bar.at(1) + bar.at(2) + bar.at(3);
    // center of the circumcircle
    answer.resize(2);
    for ( int i = 1; i <= answer.giveSize(); i++ ) {
        answer.at(i) = ( bar.at(1) * mVertices [ 0 ].at(i) + bar.at(2) * mVertices [ 1 ].at(i) + bar.at(3) * mVertices [ 2 ].at(i) ) / sum;
    }
}

void Triangle :: printYourself()
{
    printf("Triangle: ");
    for ( size_t i = 0; i < mVertices.size(); i++ ) {
        mVertices [ i ].printYourself();
    }

    printf("\n");
}

bool Triangle :: isOrientedAnticlockwise()
{
    FloatMatrix fm(3, 3);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            fm.at(i, j) = this->giveVertex(i).at(j);
        }
        fm.at(i, 3) = 1.;
    }

    if ( fm.giveDeterminant() > 0.0001 ) {
        return true;
    } else {
        return false;
    }
}

void Triangle :: changeToAnticlockwise()
{
    std :: swap(mVertices [ 1 ], mVertices [ 2 ]);
}

bool Triangle :: pointIsInTriangle(const FloatArray &iP) const
{
    FloatArray P(iP);

    const double tol2 = 1.0e-18;

    // Compute triangle normal
    FloatArray p1p2;
    p1p2.beDifferenceOf(mVertices [ 1 ], mVertices [ 0 ]);

    FloatArray p1p3;
    p1p3.beDifferenceOf(mVertices [ 2 ], mVertices [ 0 ]);


    // Edge 1
    FloatArray t1;
    t1.beDifferenceOf(mVertices [ 1 ], mVertices [ 0 ]);
    if(t1.computeSquaredNorm() < tol2) {
        // The triangle is degenerated
        return false;
    }
    else {
        t1.normalize();
    }


    FloatArray a1;

    // Edge 2
    FloatArray t2;
    t2.beDifferenceOf(mVertices [ 2 ], mVertices [ 1 ]);
    if(t2.computeSquaredNorm() < tol2) {
        // The triangle is degenerated
        return false;
    }
    else {
        t2.normalize();
    }

    FloatArray a2;


    // Edge 3
    FloatArray t3;
    t3.beDifferenceOf(mVertices [ 0 ], mVertices [ 2 ]);
    if(t3.computeSquaredNorm() < tol2) {
        // The triangle is degenerated
        return false;
    }
    else {
        t3.normalize();
    }

    FloatArray a3;


    // Project point onto triangle plane
    FloatArray pProj = P;

    if( p1p2.giveSize() == 2 ) {
        // 2D
        a1 = {-t1[1], t1[0]};
        a2 = {-t2[1], t2[0]};
        a3 = {-t3[1], t3[0]};
    }
    else {
        // 3D
        FloatArray N;

        N.beVectorProductOf(p1p2, p1p3);

        if(N.computeSquaredNorm() < tol2) {
            // The triangle is degenerated
            return false;
        }
        else {
            N.normalize();
        }

        // Compute normal distance from triangle to point
        FloatArray p1p;
        p1p.beDifferenceOf(P, mVertices [ 0 ]);
        double d = p1p.dotProduct(N);

        pProj.add(-d, N);


        a1.beVectorProductOf(N, t1);
//        if(a1.computeSquaredNorm() < tol2) {
//            // The triangle is degenerated
//            return false;
//        }
//        else {
//            a1.normalize();
//        }

        a2.beVectorProductOf(N, t2);
//        if(a2.computeSquaredNorm() < tol2) {
//            // The triangle is degenerated
//            return false;
//        }
//        else {
//            a2.normalize();
//        }

        a3.beVectorProductOf(N, t3);
//        if(a3.computeSquaredNorm() < tol2) {
//            // The triangle is degenerated
//            return false;
//        }
//        else {
//            a3.normalize();
//        }

    }


    // Check if the point is on the correct side of all edges


    FloatArray p1pProj;
    p1pProj.beDifferenceOf(pProj, mVertices [ 0 ], mVertices[0].giveSize());
    if ( p1pProj.dotProduct(a1) < 0.0 ) {
        return false;
    }



    FloatArray p2pProj;
    p2pProj.beDifferenceOf(pProj, mVertices [ 1 ], mVertices[1].giveSize());
    if ( p2pProj.dotProduct(a2) < 0.0 ) {
        return false;
    }

    FloatArray p3pProj;
    p3pProj.beDifferenceOf(pProj, mVertices [ 2 ], mVertices[2].giveSize());
    if ( p3pProj.dotProduct(a3) < 0.0 ) {
        return false;
    }

    return true;
}

void Triangle :: refineTriangle(std::vector<Triangle> &oRefinedTri, const Triangle &iTri)
{
    const FloatArray &p1 = iTri.giveVertex(1);
    const FloatArray &p2 = iTri.giveVertex(2);
    const FloatArray &p3 = iTri.giveVertex(3);

    // Compute edge midpoints
    FloatArray q1, q2, q3;

    q1.beScaled(0.5, p1);
    q1.add(0.5, p2);

    q2.beScaled(0.5, p2);
    q2.add(0.5, p3);

    q3.beScaled(0.5, p3);
    q3.add(0.5, p1);

    oRefinedTri.push_back( Triangle(p1, q1, q3) );
    oRefinedTri.push_back( Triangle(q1, q2, q3) );
    oRefinedTri.push_back( Triangle(q1, p2, q2) );
    oRefinedTri.push_back( Triangle(q3, q2, p3) );
}

Circle :: Circle(FloatArray &center, double radius) :
    mTangSignDist(1.0)
{
    mVertices.push_back(center);

    this->radius = radius;
}

void Circle :: computeNormalSignDist(double &oDist, const FloatArray &iPoint) const
{
    oDist = distance(mVertices [ 0 ], iPoint) - radius;
}

void Circle :: giveGlobalCoordinates(FloatArray &oGlobalCoord, const double &iArcPos) const
{
    double angle = 2.0*M_PI*iArcPos;

    oGlobalCoord.resize(2);
    oGlobalCoord = { mVertices[0][0] + radius*cos(angle), mVertices[0][1] + radius*sin(angle) };
}

void Circle :: initializeFrom(InputRecord &ir)
{
    mVertices.resize(1);
    IR_GIVE_FIELD(ir, mVertices [ 0 ], _IFT_Circle_center);
    IR_GIVE_FIELD(ir, radius, _IFT_Circle_radius);
}

bool Circle :: intersects(Element *element)
{
    int count = 0;
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        const auto &nodeCoor = element->giveDofManager(i)->giveCoordinates();
        // distance from the node to the center of the circle
        double dist = distance(nodeCoor, mVertices [ 0 ]);
        if ( dist > this->radius ) {
            count++;
        }
    }

    if ( count == 0 || count == element->giveNumberOfDofManagers() ) {
        return false;
    } else {
        return true;
    }
}



bool
Circle :: isInside(const FloatArray &point)
{
    double dist = distance(this->giveVertex(1), point);
    if ( dist < this->radius ) {
        return true;
    }
    return false;
}


bool Circle :: isInside(Element *element)
{   // condition should maybe be that all nodes should be inside
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        const auto &nodeCoord = element->giveDofManager(i)->giveCoordinates();
        if ( isInside(nodeCoord) ) {
            return true;
        }
    }
    return false;
}




void Circle :: computeIntersectionPoints(Element *element, std :: vector< FloatArray > &oIntersectionPoints)
{
    if ( intersects(element) ) {
        for ( int i = 1; i <= element->giveNumberOfBoundarySides(); i++ ) {
            std :: vector< FloatArray >oneLineIntersects;
            ///@todo Move semantics or something would be useful here to avoid multiple copies.
            const auto &a = element->giveDofManager ( i )->giveCoordinates();
            const auto &b = ( i != element->giveNumberOfBoundarySides() ) ? 
                element->giveDofManager ( i + 1 )->giveCoordinates() :
                element->giveDofManager ( 1 )->giveCoordinates();

            Line l(a, b);
            computeIntersectionPoints(& l, oneLineIntersects);
            for ( int j = 1; j <= int ( oneLineIntersects.size() ); j++ ) {
                oIntersectionPoints.push_back(oneLineIntersects [ j - 1 ]);
            }
        }
    }
}

void Circle :: computeIntersectionPoints(Line *l, std :: vector< FloatArray > &oIntersectionPoints)
{
    double x1 = l->giveVertex(1).at(1);
    double y1 = l->giveVertex(1).at(2);
    double x2 = l->giveVertex(2).at(1);
    double y2 = l->giveVertex(2).at(2);
    double c1 = mVertices [ 0 ].at(1);
    double c2 = mVertices [ 0 ].at(2);
    double distX = x2 - x1;
    double distY = y2 - y1;
    double a = 0., b = 0., A, B, C;
    if ( distX != 0.0 ) {
        a = distY / distX;
        b = y1 - a * x1;
        A = 1 + a * a;
        B = 2 * ( ( -1 ) * c1 + b * a - a * c2 );
        C = c1 * c1 + b * b - 2 * c2 * b + c2 * c2 - radius * radius;
    } else {
        A = 1;
        B = ( -1 ) * 2 * c2;
        C = x1 * x1 - 2 * x1 * c1 + c1 * c1 + c2 * c2 - radius * radius;
    }

    double D = B * B - 4 * A * C;
    int sz = 0;
    if ( D < 0 ) {
        sz = 0;
    } else if ( D == 0 ) {
        sz = 1;
    } else {
        sz = 2;
    }

    for ( int i = 1; i <= sz; i++ ) {
        FloatArray point(2);
        double fn;
        if ( i == 1 ) {
            fn = sqrt(D);
        } else {
            fn = ( -1 ) * sqrt(D);
        }

        if ( distX != 0.0 ) {
            point.at(1) = ( ( -1 ) * B + fn ) / ( 2 * A );
            point.at(2) = a * point.at(1) + b;
        } else {
            point.at(1) = x1;
            point.at(2) = ( ( -1 ) * B + fn ) / ( 2 * A );
        }

        if ( l->isPointInside(& point) ) {
            oIntersectionPoints.push_back(point);
        }
    }
}

int
Circle :: computeNumberOfIntersectionPoints(Element *element)
{
    std :: vector< FloatArray >intersecPoints;

    this->computeIntersectionPoints(element, intersecPoints);
    return (int) intersecPoints.size();
}

bool Circle :: isOutside(BasicGeometry *bg)
{
    int count = 0;
    for ( int i = 1; i <= bg->giveNrVertices(); i++ ) {
        if ( 0.9999 * distance(bg->giveVertex(i), mVertices [ 0 ]) > this->radius ) {
            count++;
        }
    }

    if ( count != 0 ) {
        return true;
    } else {
        return false;
    }
}

void Circle :: printYourself()
{
    printf("Circle: radius = %e, center = ", this->radius);
    mVertices [ 0 ].printYourself();
    printf("\n");
}

void Circle :: giveBoundingSphere(FloatArray &oCenter, double &oRadius)
{
    oCenter = mVertices [ 0 ];
    oRadius = radius;
}


PolygonLine :: PolygonLine() : BasicGeometry()
{
    mDebugVtk = false;
#ifdef __BOOST_MODULE
    LC.x(0.0);
    LC.y(0.0);

    UC.x(0.0);
    UC.y(0.0);
#endif
}

void PolygonLine :: computeNormalSignDist(double &oDist, const FloatArray &iPoint) const
{
    FloatArray point = {iPoint[0], iPoint[1]};

    oDist = std :: numeric_limits< double > :: max();
    int numSeg = this->giveNrVertices() - 1;

    // TODO: This can probably be done in a nicer way.
    // Ensure that we work in 2d.
    const int dim = 2;

    for ( int segId = 1; segId <= numSeg; segId++ ) {
        // Crack segment
        const FloatArray &crackP1( this->giveVertex ( segId ) );

        const FloatArray &crackP2( this->giveVertex ( segId + 1 ) );

        double dist2 = 0.0;
        if ( segId == 1 ) {
            // Vector from start P1 to point X
            FloatArray u = {point.at(1) - crackP1.at(1), point.at(2) - crackP1.at(2)};

            // Line tangent vector
            FloatArray t = {crackP2.at(1) - crackP1.at(1), crackP2.at(2) - crackP1.at(2)};
            double l2 = t.computeSquaredNorm();

            if ( l2 > 0.0 ) {
                double l = t.normalize();
                double s = dot(u, t);

                if ( s > l ) {
                    // X is closest to P2
                    dist2 = distance_square(point, crackP2);
                } else {
                    double xi = s / l;
                    auto q = ( 1.0 - xi ) * crackP1 + xi * crackP2;
                    dist2 = distance_square(point, q);
                }
            } else {
                // If the points P1 and P2 coincide,
                // we can compute the distance to any
                // of these points.
                dist2 = distance_square(point, crackP1);
            }
        } else if ( segId == numSeg ) {
            // Vector from start P1 to point X
            FloatArray u = {point.at(1) - crackP1.at(1), point.at(2) - crackP1.at(2)};

            // Line tangent vector
            FloatArray t = {crackP2.at(1) - crackP1.at(1), crackP2.at(2) - crackP1.at(2)};
            double l2 = t.computeSquaredNorm();

            if ( l2 > 0.0 ) {
                double l = t.normalize();
                double s = dot(u, t);

                if ( s < 0.0 ) {
                    // X is closest to P1
                    dist2 = distance_square(point, crackP1);
                } else {
                    double xi = s / l;
                    auto q = ( 1.0 - xi ) * crackP1 + xi * crackP2;
                    dist2 = distance_square(point, q);
                }
            } else {
                // If the points P1 and P2 coincide,
                // we can compute the distance to any
                // of these points.
                dist2 = distance_square(point, crackP1);
            }
        } else {
            double arcPos = -1.0, dummy;
            dist2 = point.distance_square(crackP1, crackP2, arcPos, dummy);
        }

        if ( dist2 < oDist*oDist ) {
            FloatArray lineToP;
            lineToP.beDifferenceOf(point, crackP1, dim);

            FloatArray t;
            t.beDifferenceOf(crackP2, crackP1, dim);

            FloatArray n = {-t.at(2), t.at(1)};

            oDist = sgn( lineToP.dotProduct(n) ) * sqrt(dist2);
        }
    }
}

void PolygonLine :: computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinArcDist) const
{
    FloatArray point = iPoint;
    point.resizeWithValues(2);

    const int numSeg = this->giveNrVertices() - 1;

    double xi = 0.0, xiUnbounded = 0.0;

    if ( numSeg == 0 ) {
        FloatArray crackP1 = giveVertex ( 1 );
        oDist = distance(crackP1, iPoint);
        oMinArcDist = 0.0;
        return;
    }

    if ( numSeg == 1 ) {
        FloatArray crackP1 = giveVertex ( 1 );
        crackP1.resizeWithValues(2);
        FloatArray crackP2 = giveVertex ( 2 );
        crackP2.resizeWithValues(2);
        point.distance(crackP1, crackP2, xi, xiUnbounded);

        if( xiUnbounded < 0.0 ) {
            oDist = xiUnbounded * distance(crackP1, crackP2);
            oMinArcDist = 0.0;
            return;
        }

        if( xiUnbounded > 1.0 ) {
            oDist = -(xiUnbounded-1.0) * distance(crackP1, crackP2);
            oMinArcDist = 1.0;
            return;
        }

        const double L = computeLength();
        double distToStart  = xi*L;

        oDist = std::min(distToStart, (L - distToStart) );
        oMinArcDist = distToStart/L;
        return;
    }

    bool isBeforeStart = false, isAfterEnd = false;
    double distBeforeStart = 0.0, distAfterEnd = 0.0;

    ///////////////////////////////////////////////////////////////////
    // Check first segment
    FloatArray crackP1_start = giveVertex ( 1 );
    crackP1_start.resizeWithValues(2);
    FloatArray crackP2_start = giveVertex ( 2 );
    crackP2_start.resizeWithValues(2);
    const double distSeg_start = point.distance(crackP1_start, crackP2_start, xi, xiUnbounded);

    if( xiUnbounded < 0.0 ) {
        isBeforeStart = true;
        distBeforeStart = xiUnbounded * distance(crackP1_start, crackP2_start);
    }

    double arcPosPassed = distance(crackP1_start, crackP2_start);
    double distToStart  = xi * distance(crackP1_start, crackP2_start);

    double minGeomDist  = distSeg_start;




    ///////////////////////////////////////////////////////////////////
    // Check interior segments
    for ( int segId = 2; segId <= numSeg-1; segId++ ) {
        FloatArray crackP1 = giveVertex ( segId );
        crackP1.resizeWithValues(2);
        FloatArray crackP2 = giveVertex ( segId+1 );
        crackP2.resizeWithValues(2);

        const double distSeg = point.distance(crackP1, crackP2, xi, xiUnbounded);

        if(distSeg < minGeomDist) {
            isBeforeStart = false;
            minGeomDist = distSeg;
            distToStart = arcPosPassed + xi * distance(crackP1, crackP2);
        }

        arcPosPassed += distance(crackP1, crackP2);

    }



    ///////////////////////////////////////////////////////////////////
    // Check last segment
    FloatArray crackP1_end = giveVertex ( numSeg );
    crackP1_end.resizeWithValues(2);
    FloatArray crackP2_end = giveVertex ( numSeg+1 );
    crackP2_end.resizeWithValues(2);
    const double distSeg_end = point.distance(crackP1_end, crackP2_end, xi, xiUnbounded);

    if ( numSeg > 1 ) {
        if( xiUnbounded > 1.0 ) {
            arcPosPassed += xiUnbounded * distance(crackP1_end, crackP2_end);
        }
        else {
            arcPosPassed += xi * distance(crackP1_end, crackP2_end);
        }
    }

    if ( distSeg_end < minGeomDist ) {
        isBeforeStart = false;

        if( xiUnbounded > 1.0 ) {
            isAfterEnd = true;
            distAfterEnd = -(xiUnbounded-1.0) * distance(crackP1_end, crackP2_end);
        }

        distToStart = arcPosPassed;
    }

    ///////////////////////////////////////////////////////////////////
    // Return result

    if ( isBeforeStart ) {
        oDist = distBeforeStart;
        oMinArcDist = 0.0;
        return;
    }

    if ( isAfterEnd ) {
        oDist = distAfterEnd;
        oMinArcDist = 1.0;
        return;
    }

    const double L = computeLength();

    oDist = std::min(distToStart, (L - distToStart) );
    oMinArcDist = distToStart/L;
}

void PolygonLine :: computeLocalCoordinates(FloatArray &oLocCoord, const FloatArray &iPoint) const
{ }

double PolygonLine :: computeLength() const
{
    if ( mVertices.size() == 0 ) {
        return 0.0;
    }

    double L = 0.0;

    size_t numSeg = mVertices.size() - 1;

    if ( numSeg == 0 ) {
        return 0.0;
    }

    for ( size_t i = 0; i < numSeg; i++ ) {
        L += distance(mVertices [ i ], mVertices [ i + 1 ]);
    }

    return L;
}

void PolygonLine :: giveSubPolygon(std :: vector< FloatArray > &oPoints, const double &iXiStart, const double &iXiEnd) const
{
    double L = computeLength();
    double xSegStart = 0.0, xSegEnd = 0.0;
    double xiSegStart = 0.0, xiSegEnd = 0.0;
    size_t numSeg = mVertices.size() - 1;
    const double xiTol = 1.0e-9;
    if ( iXiStart < xiTol ) {
        // Add first point
        oPoints.push_back(mVertices [ 0 ]);
    }

    for ( size_t i = 0; i < numSeg; i++ ) {
        xSegEnd += distance(mVertices [ i ], mVertices [ i + 1 ]);

        xiSegStart = xSegStart / L;
        xiSegEnd = xSegEnd / L;

        if ( iXiStart > xiSegStart-xiTol && iXiStart < xiSegEnd+xiTol ) {
            // Start point is within the segment
            FloatArray p;
            double elXi = ( iXiStart - xiSegStart ) / ( xiSegEnd - xiSegStart );
            p.beScaled( ( 1.0 - elXi ), mVertices [ i ] );
            p.add(elXi, mVertices [ i + 1 ]);
            oPoints.push_back(p);
        }


        if ( iXiEnd > xiSegStart && iXiEnd < xiSegEnd ) {
            // End point is within the segment
            FloatArray p;
            double elXi = ( iXiEnd - xiSegStart ) / ( xiSegEnd - xiSegStart );
            p.beScaled( ( 1.0 - elXi ), mVertices [ i ] );
            p.add(elXi, mVertices [ i + 1 ]);
            oPoints.push_back(p);
        }

        if ( xiSegEnd > iXiStart && xiSegEnd < iXiEnd + xiTol ) {
            // End point of the segment is within range
            oPoints.push_back(mVertices [ i + 1 ]);
        }

        xSegStart = xSegEnd;
    }
}

void PolygonLine :: giveGlobalCoordinates(FloatArray &oGlobalCoord, const double &iArcPos) const
{
    double L = computeLength();
    double xSegStart = 0.0, xSegEnd = 0.0;
    double xiSegStart = 0.0, xiSegEnd = 0.0;
    size_t numSeg = mVertices.size() - 1;
    const double xiTol = 1.0e-9;
    if ( iArcPos < xiTol ) {
        oGlobalCoord = mVertices [ 0 ];
        return;
    }

    for ( size_t i = 0; i < numSeg; i++ ) {
        xSegEnd += distance(mVertices [ i ], mVertices [ i + 1 ]);

        xiSegStart = xSegStart / L;
        xiSegEnd = xSegEnd / L;

        if ( iArcPos > xiSegStart-xiTol && iArcPos < xiSegEnd+xiTol ) {
            // Point is within the segment
            FloatArray p;
            double elXi = ( iArcPos - xiSegStart ) / ( xiSegEnd - xiSegStart );
            p.beScaled( ( 1.0 - elXi ), mVertices [ i ] );
            p.add(elXi, mVertices [ i + 1 ]);
            oGlobalCoord = p;
            return;
        }

    }
}

void PolygonLine :: giveNormal(FloatArray &oNormal, const double &iArcPosition) const
{
    double L = computeLength();
    double xSegStart = 0.0, xSegEnd = 0.0;
    double xiSegStart = 0.0, xiSegEnd = 0.0;
    size_t numSeg = mVertices.size() - 1;
    const double xiTol = 1.0e-9;

    for ( size_t i = 0; i < numSeg; i++ ) {
        xSegEnd += distance(mVertices [ i ], mVertices [ i + 1 ]);

        xiSegStart = xSegStart / L;
        xiSegEnd = xSegEnd / L;

        if ( iArcPosition > xiSegStart-xiTol && iArcPosition < xiSegEnd+xiTol ) {
            // The given point is within the segment

            const FloatArray &p1 = mVertices [ i ];
            const FloatArray &p2 = mVertices [ i+1 ];

            FloatArray t = {p2(0) - p1(0), p2(1) - p1(1)};

            oNormal.resize(2);
            oNormal(0) = -t(1);
            oNormal(1) =  t(0);
            oNormal.normalize();

            return;
        }

    }

    OOFEM_ERROR("Arc position not found.")
}

void PolygonLine :: giveTangent(FloatArray &oTangent, const double &iArcPosition) const
{
    double L = computeLength();
    double xSegStart = 0.0, xSegEnd = 0.0;
    double xiSegStart = 0.0, xiSegEnd = 0.0;
    size_t numSeg = mVertices.size() - 1;
    const double xiTol = 1.0e-9;

    for ( size_t i = 0; i < numSeg; i++ ) {
        xSegEnd += distance(mVertices [ i ], mVertices [ i + 1 ]);

        xiSegStart = xSegStart / L;
        xiSegEnd = xSegEnd / L;

        if ( iArcPosition > xiSegStart-xiTol && iArcPosition < xiSegEnd+xiTol ) {
            // The given point is within the segment

            const FloatArray &p1 = mVertices [ i ];
            const FloatArray &p2 = mVertices [ i+1 ];

            oTangent = {p2(0) - p1(0), p2(1) - p1(1)};

            oTangent.normalize();

            return;
        }

    }

    OOFEM_ERROR("Arc position not found.")
}

void PolygonLine :: initializeFrom(InputRecord &ir)
{
    FloatArray points;
    IR_GIVE_FIELD(ir, points, _IFT_PolygonLine_points);

    int numPoints = points.giveSize() / 2;

    for ( int i = 1; i <= numPoints; i++ ) {
        mVertices.push_back({points.at(2 * ( i - 1 ) + 1), points.at( 2 * ( i   ) )});
    }

#ifdef __BOOST_MODULE
    // Precompute bounding box to speed up calculation of intersection points.
    calcBoundingBox(LC, UC);
#endif
}

void PolygonLine :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField( "PolygonLine", 1 );

    FloatArray points;
    int nVert = (int) mVertices.size();
    points.resize(nVert * 2);

    for ( int i = 0; i < nVert; i++ ) {
        points.at(2 * i + 1) = mVertices [ i ].at(1);
        points.at(2 * i + 2) = mVertices [ i ].at(2);
    }

    input.setField(points, _IFT_PolygonLine_points);
}

#ifdef __BOOST_MODULE
bool PolygonLine :: boundingBoxIntersects(Element *element)
{
    // Compute element bounding box
    bPoint2 eLC, eUC;
    eLC.x( element->giveNode(1)->giveCoordinate(1) );
    eLC.y( element->giveNode(1)->giveCoordinate(2) );

    eUC.x( element->giveNode(1)->giveCoordinate(1) );
    eUC.y( element->giveNode(1)->giveCoordinate(2) );

    int numNodes = element->giveNumberOfNodes();
    for ( int i = 2; i <= numNodes; i++ ) {
        eLC.x( min( eLC.x(), element->giveNode(i)->giveCoordinate(1) ) );
        eLC.y( min( eLC.y(), element->giveNode(i)->giveCoordinate(2) ) );

        eUC.x( max( eUC.x(), element->giveNode(i)->giveCoordinate(1) ) );
        eUC.y( max( eUC.y(), element->giveNode(i)->giveCoordinate(2) ) );
    }

    //printf("eLC: (%e, %e) eUC: (%e, %e) ", eLC.x(), eLC.y(), eUC.x(), eUC.y() );
    //printf(" LC: (%e, %e)  UC: (%e, %e)\n",  LC.x(),  LC.y(),  UC.x(),  UC.y() );


    // Check if there is any chance of overlap
    if ( !bOverlap(eLC, eUC, LC, UC) ) {
        // If the bounding boxes do not overlap,
        // there will certainly not be any intersections
        return false;
    }


    return true;
}
#endif

bool PolygonLine :: intersects(Element *element)
{
    printf("Warning: entering PolygonLine :: intersects(Element *element).\n");

#ifdef __BOOST_MODULE
    if ( !boundingBoxIntersects(element) ) {
        return false;
    }

    double distTol = 1.0e-9;

    int numSeg = this->giveNrVertices() - 1;


    // Loop over the crack segments and test each segment for
    // overlap with the element
    for ( int segId = 1; segId <= numSeg; segId++ ) {
        if ( element->giveGeometryType() == EGT_triangle_1 ) {
            // Crack segment
            bPoint2 crackP1( this->giveVertex ( segId )->at(1), this->giveVertex ( segId )->at(2) );
            bPoint2 crackP2( this->giveVertex ( segId + 1 )->at(1), this->giveVertex ( segId + 1 )->at(2) );
            bSeg2 crackSeg(crackP1, crackP2);


            // Triangle vertices
            bPoint2 x1( element->giveNode ( 1 )->giveCoordinate(1), element->giveNode ( 1 )->giveCoordinate(2) );
            bPoint2 x2( element->giveNode ( 2 )->giveCoordinate(1), element->giveNode ( 2 )->giveCoordinate(2) );
            bPoint2 x3( element->giveNode ( 3 )->giveCoordinate(1), element->giveNode ( 3 )->giveCoordinate(2) );

            bSeg2 edge1(x1, x2);
            bSeg2 edge2(x2, x3);
            bSeg2 edge3(x3, x1);

            double d1 = bDist(crackSeg, edge1);
            if ( d1 < distTol ) {
                return true;
            }

            double d2 = bDist(crackSeg, edge2);
            if ( d2 < distTol ) {
                return true;
            }

            double d3 = bDist(crackSeg, edge3);
            if ( d3 < distTol ) {
                return true;
            }
        }
    }



#endif

    return false;
}


bool PolygonLine :: isInside(Element *element)
{
    return false;
}

bool
PolygonLine :: isInside(const FloatArray &point)
{
    return false;
}



#ifdef __BOOST_MODULE
void PolygonLine :: calcBoundingBox(bPoint2 &oLC, bPoint2 &oUC)
{
    oLC.x( vertices->at(1)->at(1) );
    oLC.y( vertices->at(1)->at(2) );

    oUC.x( vertices->at(1)->at(1) );
    oUC.y( vertices->at(1)->at(2) );

    int numPoints = vertices->giveSize();
    for ( int i = 2; i <= numPoints; i++ ) {
        oLC.x( min( oLC.x(), vertices->at(i)->at(1) ) );
        oLC.y( min( oLC.y(), vertices->at(i)->at(2) ) );

        oUC.x( max( oUC.x(), vertices->at(i)->at(1) ) );
        oUC.y( max( oUC.y(), vertices->at(i)->at(2) ) );
    }
}
#endif


void PolygonLine :: computeIntersectionPoints(Element *element, std :: vector< FloatArray > &oIntersectionPoints)
{

    for ( int i = 1; i <= element->giveNumberOfBoundarySides(); i++ ) {
        std :: vector< FloatArray > oneLineIntersects;
        ///@todo Move semantics or something would be useful here to avoid multiple copies.
        const auto &xStart = element->giveDofManager ( i )->giveCoordinates();
        const auto &xEnd = ( i != element->giveNumberOfBoundarySides() ) ?
            element->giveDofManager ( i + 1 )->giveCoordinates() :
            element->giveDofManager ( 1 )->giveCoordinates();

        computeIntersectionPoints(xStart, xEnd, oneLineIntersects);

        for ( auto &interSect: oneLineIntersects ) {
            // Check that the intersection point has not already been identified.
            // This may happen if the crack intersects the element exactly at a node,
            // so that intersection is detected for both element edges in that node.

            // TODO: Set tolerance in a more transparent way.
            double distTol = 1.0e-9;

            bool alreadyFound = false;

            for ( auto &pInterSect: oIntersectionPoints ) {

                if ( distance(pInterSect, interSect) < distTol ) {
                    alreadyFound = true;
                    break;
                }
            }

            if ( !alreadyFound ) {
                oIntersectionPoints.push_back(interSect);
            }
        }


    }

#if 0
    printf("Warning: entering  PolygonLine :: computeIntersectionPoints(Element *element, std::vector< FloatArray > &oIntersectionPoints).\n");
#ifdef __BOOST_MODULE

    if ( !boundingBoxIntersects(element) ) {
        return;
    }


    for ( int i = 1; i <= element->giveNumberOfBoundarySides(); i++ ) {
        std :: vector< FloatArray >oneLineIntersects;
        ///@todo Move semantics or something would be useful here to avoid multiple copies.
        const auto &a = * element->giveDofManager ( i )->giveCoordinates();
        const auto &b = ( i != element->giveNumberOfBoundarySides() ) ?
            * element->giveDofManager ( i + 1 )->giveCoordinates() :
            * element->giveDofManager ( 1 )->giveCoordinates();

        Line l(a, b);

        computeIntersectionPoints(& l, oneLineIntersects);
        for ( FloatArray &interSect: oneLineIntersects ) {
            // Check that the intersection point has not already been identified.
            // This may happen if the crack intersects the element exactly at a node,
            // so that intersection is detected for both element edges in that node.

            // TODO: Set tolerance in a more transparent way.
            double distTol = 1.0e-9;

            bool alreadyFound = false;

            bPoint2 pNew( interSect.at(1), interSect.at(2) );

            for ( FloatArray &pInterSect: oIntersectionPoints ) {
                bPoint2 pOld( pInterSect.at(1), pInterSect.at(2) );

                if ( bDist(pOld, pNew) < distTol ) {
                    alreadyFound = true;
                    break;
                }
            }

            if ( !alreadyFound ) {
                oIntersectionPoints.push_back(interSect);
            }
        }
    }
#endif

#endif
}

void PolygonLine :: computeIntersectionPoints(Line *l, std :: vector< FloatArray > &oIntersectionPoints)
{
    printf("Warning: entering  PolygonLine :: computeIntersectionPoints(Line *l, std::vector< FloatArray > &oIntersectionPoints).\n");

#ifdef __BOOST_MODULE


    int numSeg = this->giveNrVertices() - 1;


    // Segment
    bPoint2 lineP1( l->giveVertex ( 1 )->at(1), l->giveVertex ( 1 )->at(2) );
    bPoint2 lineP2( l->giveVertex ( 2 )->at(1), l->giveVertex ( 2 )->at(2) );
    bSeg2 lineSeg(lineP1, lineP2);


    double distTol = 1.0e-9;

    bool foundOverlap = false;

    // Loop over the crack segments and test each segment for
    // overlap with the element
    for ( int segId = 1; segId <= numSeg; segId++ ) {
        // Crack segment
        bPoint2 crackP1( this->giveVertex ( segId )->at(1), this->giveVertex ( segId )->at(2) );
        bPoint2 crackP2( this->giveVertex ( segId + 1 )->at(1), this->giveVertex ( segId + 1 )->at(2) );
        bSeg2 crackSeg(crackP1, crackP2);

        bPoint2 intersectionPoint(0.0, 0.0);
        double d = bDist(crackSeg, lineSeg, & intersectionPoint);
        if ( d < distTol ) {
            if ( !foundOverlap ) {
                foundOverlap = true;
                oIntersectionPoints.emplace_back({intersectionPoint.x(), intersectionPoint.y()});
            }
        }
    }

#endif
}

void PolygonLine :: computeIntersectionPoints(const PolygonLine &iPolygonLine, std :: vector< FloatArray > &oIntersectionPoints) const
{
    int numSeg = this->giveNrVertices() - 1;
    for(int segIndex = 1; segIndex <= numSeg; segIndex++) {

        const FloatArray &xStart = this->giveVertex(segIndex);
        const FloatArray &xEnd = this->giveVertex(segIndex+1);

        iPolygonLine.computeIntersectionPoints(xStart, xEnd, oIntersectionPoints);
    }
}

void PolygonLine :: computeIntersectionPoints(const FloatArray &iXStart, const FloatArray &iXEnd, std :: vector< FloatArray > &oIntersectionPoints) const
{
    const double detTol = 1.0e-15;

    int numSeg = this->giveNrVertices() - 1;
    for(int segIndex = 1; segIndex <= numSeg; segIndex++) {

        const FloatArray &xStart = this->giveVertex(segIndex);
        const FloatArray &xEnd = this->giveVertex(segIndex+1);

        const FloatArray t1 = {xEnd(0) - xStart(0), xEnd(1) - xStart(1)};
        const FloatArray t2 = {iXEnd(0) - iXStart(0), iXEnd(1) - iXStart(1)};

        double xi1 = 0.0, xi2 = 0.0;
        int maxIter = 1;

        for(int iter = 0; iter < maxIter; iter++) {
            FloatArray temp = {iXStart(0) + xi2*t2(0) - xStart(0) - xi1*t1(0), iXStart(1) + xi2*t2(1) - xStart(1) - xi1*t1(1)};
            FloatArray res = {-t1.dotProduct(temp), t2.dotProduct(temp)};

            //printf("iter: %d res: %e\n", iter, res.computeNorm() );

            FloatMatrix K(2,2);
            K(0,0) = t1.dotProduct(t1);
            K(0,1) = -t1.dotProduct(t2);
            K(1,0) = -t1.dotProduct(t2);
            K(1,1) = t2.dotProduct(t2);

            double detK = K.giveDeterminant();

            if(detK < detTol) {
                return;
            }

            FloatMatrix KInv;
            KInv.beInverseOf(K);

            FloatArray dxi;
            dxi.beProductOf(KInv, res);

            xi1 -= dxi[0];
            xi2 -= dxi[1];
        }

//        printf("xi1: %e xi2: %e\n", xi1, xi2);


        if(xi1 >= 0.0 && xi1 <= 1.0 && xi2 >= 0.0 && xi2 <= 1.0) {
            FloatArray pos = xStart;
            pos.add(xi1, t1);
            oIntersectionPoints.push_back(pos);
        }

    }

}

int
PolygonLine :: computeNumberOfIntersectionPoints(Element *element)
{
    std :: vector< FloatArray >intersecPoints;
    this->computeIntersectionPoints(element, intersecPoints);
    return (int) intersecPoints.size();
}

bool PolygonLine :: isOutside(BasicGeometry *bg)
{
    return true;
}

void PolygonLine :: printYourself()
{
    printf("PolygonLine:\n");

    for(const FloatArray &x : mVertices) {
        x.printYourself();
    }

    printf("\n");
}

void PolygonLine :: printVTK(int iTStepIndex, int iLineIndex)
{
    // Debugging function: write crack geometry to vtk.

    // Write crack geometry to vtk
    std :: string vtkFileName;
    vtkFileName.append("crack");
    char lineIdNumberString [ 100 ];
    sprintf(lineIdNumberString, "%d", iLineIndex);
    vtkFileName.append(lineIdNumberString);

    vtkFileName.append("Step");
    char stepString [ 100 ];

    sprintf(stepString, "%d", iTStepIndex);
    vtkFileName.append(stepString);

    vtkFileName.append(".vtk");


    int numPoints = this->giveNrVertices();


    std :: ofstream file;
    file.open( vtkFileName.data() );

    // Write header
    file << "# vtk DataFile Version 2.0\n";
    file << "Geometry of a PolygonLine\n";
    file << "ASCII\n";

    file << "DATASET UNSTRUCTURED_GRID\n";


    // Write points
    file << "POINTS " << numPoints << "double\n";

    for ( int i = 1; i <= numPoints; i++ ) {
        file << this->giveVertex(i).at(1) << " " << this->giveVertex(i).at(2) << " 0.0\n";
    }


    // Write segments
    int numSeg = numPoints - 1;
    file << "CELLS " << numSeg << " " << numSeg * 3 << "\n";

    int numPointsPerSeg = 2;
    for ( int i = 0; i < numSeg; i++ ) {
        file << numPointsPerSeg << " " << i << " " << i + 1 << "\n";
    }


    // Write cell types
    file << "CELL_TYPES " << numSeg << "\n";
    int vtkCellType = 3;     // line segment
    for ( int i = 0; i < numSeg; i++ ) {
        file << vtkCellType << "\n";
    }

    file.close();
}


bool PolygonLine :: giveTips(TipInfo &oStartTipInfo, TipInfo &oEndTipInfo) const
{
    int nVert = giveNrVertices();
    if ( nVert > 1 ) {
        // Start tip
        TipInfo info1;
        const FloatArray &p1S = ( giveVertex(1) );
        const FloatArray &p2S = ( giveVertex(2) );

        // Tip position
        info1.mGlobalCoord = p1S;

        // Tip tangent
        info1.mTangDir.beDifferenceOf(p1S, p2S);
        info1.mTangDir.normalize();

        // Tip normal
        info1.mNormalDir = {
            -info1.mTangDir.at(2), info1.mTangDir.at(1)
        };

        info1.mTipIndex = 0;
        info1.mArcPos = 0.0;

        oStartTipInfo = info1;

        // End tip
        TipInfo info2;
        const FloatArray &p1E = ( giveVertex(nVert - 1) );
        const FloatArray &p2E = ( giveVertex(nVert) );

        // Tip position
        info2.mGlobalCoord = p2E;

        // Tip tangent
        info2.mTangDir.beDifferenceOf(p2E, p1E);
        info2.mTangDir.normalize();

        // Tip normal
        info2.mNormalDir = {
            -info2.mTangDir.at(2), info2.mTangDir.at(1)
        };

        info2.mTipIndex = 1;
        info2.mArcPos = 1.0;

        oEndTipInfo = info2;

        return true;
    }

    return false;
}

void PolygonLine :: giveBoundingSphere(FloatArray &oCenter, double &oRadius)
{
    int nVert = giveNrVertices();
    oCenter = {
        0.0, 0.0
    };
    oRadius = 0.0;

    if ( nVert > 0 ) {
        for ( int i = 1; i <= nVert; i++ ) {
            oCenter.add( giveVertex(i) );
        }

        oCenter.times( 1.0 / double( nVert ) );

        for ( int i = 1; i <= nVert; i++ ) {
            oRadius = std :: max( oRadius, distance(oCenter, giveVertex(i) ) );
        }
    }

}

void PolygonLine :: cropPolygon(const double &iArcPosStart, const double &iArcPosEnd)
{
    std::vector<FloatArray> points;
    giveSubPolygon(points, iArcPosStart, iArcPosEnd);
    setVertices(points);

    const double tol2 = 1.0e-18;
    removeDuplicatePoints(tol2);

}

void PointSwarm :: initializeFrom(InputRecord &ir)
{
    IntArray idList;

    IR_GIVE_FIELD(ir, idList, _IFT_PointSwarm_nodeID); // Macro

    for ( int i = 1; i <= idList.giveSize(); i++ ) {
        this->idList.push_back( idList.at(i) );
    }
}
} // end namespace oofem
