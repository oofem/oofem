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
#include "alist.h"
#include "geometry.h"
#include "element.h"
#include "dofmanager.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"


#include <fstream>
#include <limits>

namespace oofem {
REGISTER_Geometry(Line)
REGISTER_Geometry(Circle)
REGISTER_Geometry(PointSwarm)

BasicGeometry :: BasicGeometry()
{ }

BasicGeometry :: BasicGeometry(const BasicGeometry &iBasicGeometry) :
    mVertices(iBasicGeometry.mVertices)
{ }

BasicGeometry :: ~BasicGeometry()
{ }

// TODO: change to const FloatArray &iVertex
void BasicGeometry :: setVertex(FloatArray *vertex)
{
    mVertices.push_back(* vertex);
    delete vertex;
}

bool Line :: intersects(ElementGeometry *element)
{
    int ip = this->computeNumberOfIntersectionPoints(element);
    return ( ip > 0 );
}

Line :: Line(FloatArray *pointA, FloatArray *pointB) : BasicGeometry()
{
    mVertices.push_back(* pointA);
    delete pointA;
    mVertices.push_back(* pointB);
    delete pointB;
}

double Line :: computeDistanceTo(const FloatArray *point)
{
    const FloatArray &pointA = mVertices [ 0 ];
    const FloatArray &pointB = mVertices [ 1 ];
    double a = pointA.at(2) - pointB.at(2);
    double b = pointB.at(1) - pointA.at(1);
    double c = pointA.at(1) * pointB.at(2) - pointB.at(1) * pointA.at(2);
    double l = pointA.distance(pointB);
    return ( a * point->at(1) + b * point->at(2) + c ) / l;
}

void Line :: computeProjection(FloatArray &answer)
{
    answer.beDifferenceOf(mVertices [ 1 ], mVertices [ 0 ]);
}

double Line :: computeTangentialDistanceToEnd(FloatArray *point)
{
    FloatArray projection;
    this->computeProjection(projection);
    FloatArray tmp;
    tmp.beDifferenceOf(* point, mVertices [ 1 ]);
    return tmp.dotProduct(projection) / projection.computeNorm();
}

int Line :: computeNumberOfIntersectionPoints(ElementGeometry *element)
{
    int count = 0;
    int nrNodes = element->giveNumberOfDofManagers();
    FloatArray signedDist(nrNodes);
    FloatArray tanSignDist(nrNodes);

    for ( int i = 1; i <= nrNodes; i++ ) {
        signedDist.at(i) = computeDistanceTo( element->giveDofManager(i)->giveCoordinates() );
        tanSignDist.at(i) = computeTangentialDistanceToEnd( element->giveDofManager(i)->giveCoordinates() );
    }

    // here I need to get max value and min value in the FloatArray
    double maxDist = signedDist.at(1);
    double maxTanDist = tanSignDist.at(1);
    double minDist = signedDist.at(1);
    double minTanDist = tanSignDist.at(1);
    for ( int i = 2; i <= nrNodes; i++ ) {
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
        if ( maxTanDist <= 0.0 ) {
            count = 2;
        } else if ( ( maxTanDist * minTanDist ) <= 0 ) {
            count = 1;
        } else {
            count = 0;
        }
    }

    return count;
}

void Line :: computeIntersectionPoints(ElementGeometry *element, std :: vector< FloatArray > &oIntersectionPoints)
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

IRResultType Line :: initializeFrom(InputRecord *ir)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro

    FloatArray *start = new FloatArray(2);
    FloatArray *end = new FloatArray(2);
    IR_GIVE_FIELD(ir, * start, _IFT_Line_start);
    IR_GIVE_FIELD(ir, * end, _IFT_Line_end);

    mVertices.push_back(* start);
    delete start;
    mVertices.push_back(* end);
    delete end;
    return IRRT_OK;
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
        if ( this->computeDistanceTo( & ( bg->giveVertex(i) ) ) > 0.1 ) {
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
    return 0.25 * mVertices [ 0 ].distance(mVertices [ 1 ]) *
           mVertices [ 1 ].distance(mVertices [ 2 ]) *
           mVertices [ 0 ].distance(mVertices [ 2 ]) / this->getArea();
}

void Triangle :: computeBarycentrCoor(FloatArray &answer) const
{
    double c = mVertices [ 0 ].distance(mVertices [ 1 ]);
    double a = mVertices [ 1 ].distance(mVertices [ 2 ]);
    double b = mVertices [ 0 ].distance(mVertices [ 2 ]);

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
    // Compute triangle normal
    FloatArray p1p2;
    p1p2.beDifferenceOf(mVertices [ 1 ], mVertices [ 0 ]);

    FloatArray p1p3;
    p1p3.beDifferenceOf(mVertices [ 2 ], mVertices [ 0 ]);

    FloatArray N;
    N.beVectorProductOf(p1p2, p1p3);
    N.normalize();

    // Compute normal distance from triangle to point
    FloatArray p1p;
    p1p.beDifferenceOf(iP, mVertices [ 0 ]);
    double d = p1p.dotProduct(N);

    // Project point onto triangle plane
    FloatArray pProj = iP;
    pProj.add(-d, N);

    // Check if the point is on the correct side of all edges

    // Edge 1
    FloatArray t1;
    t1.beDifferenceOf(mVertices [ 1 ], mVertices [ 0 ]);
    t1.normalize();
    FloatArray a1;
    a1.beVectorProductOf(N, t1);
    a1.normalize();

    FloatArray p1pProj;
    p1pProj.beDifferenceOf(pProj, mVertices [ 0 ]);
    if ( p1pProj.dotProduct(a1) < 0.0 ) {
        return false;
    }


    // Edge 2
    FloatArray t2;
    t2.beDifferenceOf(mVertices [ 2 ], mVertices [ 1 ]);
    t2.normalize();
    FloatArray a2;
    a2.beVectorProductOf(N, t2);
    a2.normalize();

    FloatArray p2pProj;
    p2pProj.beDifferenceOf(pProj, mVertices [ 1 ]);
    if ( p2pProj.dotProduct(a2) < 0.0 ) {
        return false;
    }


    // Edge 3
    FloatArray t3;
    t3.beDifferenceOf(mVertices [ 0 ], mVertices [ 2 ]);
    t3.normalize();
    FloatArray a3;
    a3.beVectorProductOf(N, t3);
    a3.normalize();

    FloatArray p3pProj;
    p3pProj.beDifferenceOf(pProj, mVertices [ 2 ]);
    if ( p3pProj.dotProduct(a3) < 0.0 ) {
        return false;
    }

    return true;
}

Circle :: Circle(FloatArray *center, double radius) :
    mTangSignDist(1.0)
{
    mVertices.push_back(* center);
    delete center;

    this->radius = radius;
}

void Circle :: computeNormalSignDist(double &oDist, const FloatArray &iPoint) const
{
    oDist = mVertices [ 0 ].distance(iPoint) - radius;
}

IRResultType Circle :: initializeFrom(InputRecord *ir)
{

    IRResultType result; // Required by IR_GIVE_FIELD macro

    FloatArray *center = new FloatArray(2);
    IR_GIVE_FIELD(ir, * center, _IFT_Circle_center);
    IR_GIVE_FIELD(ir, radius, _IFT_Circle_radius);
    mVertices.push_back(* center);
    delete center;
    return IRRT_OK;
}

bool Circle :: intersects(ElementGeometry *element)
{
    int count = 0;
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        FloatArray *nodeCoor = element->giveDofManager(i)->giveCoordinates();
        // distance from the node to the center of the circle
        double dist = nodeCoor->distance(mVertices [ 0 ]);
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
Circle :: isInside(FloatArray &point)
{
    double dist = this->giveVertex(1).distance(point);
    if ( dist < this->radius ) {
        return true;
    }
    return false;
}


bool Circle :: isInside(ElementGeometry *element)
{   // condition should maybe be that all nodes should be inside
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        FloatArray nodeCoord = * element->giveDofManager(i)->giveCoordinates();
        if ( isInside(nodeCoord) ) {
            return true;
        }
    }
    return false;
}




void Circle :: computeIntersectionPoints(ElementGeometry *element, std :: vector< FloatArray > &oIntersectionPoints)
{
    if ( intersects(element) ) {
        for ( int i = 1; i <= element->giveNumberOfBoundarySides(); i++ ) {
            std :: vector< FloatArray >oneLineIntersects;
            FloatArray *a = new FloatArray( * ( element->giveDofManager ( i )->giveCoordinates() ) );
            FloatArray *b = NULL;
            if ( i != element->giveNumberOfBoundarySides() ) {
                b = new FloatArray( * ( element->giveDofManager ( i + 1 )->giveCoordinates() ) );
            } else {
                b = new FloatArray( * ( element->giveDofManager ( 1 )->giveCoordinates() ) );
            }

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
Circle :: computeNumberOfIntersectionPoints(ElementGeometry *element)
{
    std :: vector< FloatArray >intersecPoints;

    this->computeIntersectionPoints(element, intersecPoints);
    return intersecPoints.size();
}

bool Circle :: isOutside(BasicGeometry *bg)
{
    int count = 0;
    for ( int i = 1; i <= bg->giveNrVertices(); i++ ) {
        if ( 0.9999 * bg->giveVertex(i).distance(mVertices [ 0 ]) > this->radius ) {
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
    oDist = std :: numeric_limits< double > :: max();
    int numSeg = this->giveNrVertices() - 1;

    // TODO: This can probably be done in a nicer way.
    // Ensure that we work in 2d.
    FloatArray point;
    point.setValues( 2, iPoint.at(1), iPoint.at(2) );

    for ( int segId = 1; segId <= numSeg; segId++ ) {
        // Crack segment
        const FloatArray &crackP1( this->giveVertex ( segId ) );

        const FloatArray &crackP2( this->giveVertex ( segId + 1 ) );

        double dist = 0.0;
        if ( segId == 1 ) {
            // Vector from start P1 to point X
            FloatArray u;
            u.beDifferenceOf(point, crackP1);

            // Line tangent vector
            FloatArray t;
            t.beDifferenceOf(crackP2, crackP1);
            double l = norm(t);

            if ( l > 0.0 ) {
                t.normalize();
                double s = dot(u, t);

                //				if( s < 0.0 ) {
                //					// X is closest to P1
                //					dist = point.distance(crackP1);
                //				}
                //				else {
                if ( s > l ) {
                    // X is closest to P2
                    dist = point.distance(crackP2);
                } else {
                    double xi = s / l;
                    FloatArray q = ( 1.0 - xi ) * crackP1 + xi * crackP2;
                    dist = point.distance(q);
                }
                //				}
            } else {
                // If the points P1 and P2 coincide,
                // we can compute the distance to any
                // of these points.
                dist = point.distance(crackP1);
            }
        } else if ( segId == numSeg ) {
            // Vector from start P1 to point X
            FloatArray u;
            u.beDifferenceOf(point, crackP1);

            // Line tangent vector
            FloatArray t;
            t.beDifferenceOf(crackP2, crackP1);
            double l = norm(t);

            if ( l > 0.0 ) {
                t.normalize();
                double s = dot(u, t);

                if ( s < 0.0 ) {
                    // X is closest to P1
                    dist = point.distance(crackP1);
                } else {
                    //					if( s > l ) {
                    //						// X is closest to P2
                    //						dist = point.distance(crackP2);
                    //					}
                    //					else {
                    double xi = s / l;
                    FloatArray q = ( 1.0 - xi ) * crackP1 + xi * crackP2;
                    dist = point.distance(q);
                    //					}
                }
            } else {
                // If the points P1 and P2 coincide,
                // we can compute the distance to any
                // of these points.
                dist = point.distance(crackP1);
            }
        } else {
            double arcPos = -1.0;
            dist = point.distance(crackP1, crackP2, arcPos);
        }


        FloatArray t;
        t.beDifferenceOf(crackP2, crackP1);

        FloatArray n;
        n.setValues( 2, -t.at(2), t.at(1) );

        FloatArray lineToP;
        lineToP.beDifferenceOf(point, crackP1);

        double sign = sgn( lineToP.dotProduct(n) );

        if ( dist < fabs(oDist) ) {
            oDist = sign * dist;
        }
    }
}


void PolygonLine :: computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinArcDist) const
{
    double totalArcLength = computeLength();

    int numSeg = this->giveNrVertices() - 1;

    // TODO: This can probably be done in a nicer way.
    // Ensure that we work in 2d.
    FloatArray point;
    point.setValues( 2, iPoint.at(1), iPoint.at(2) );


    const FloatArray &crackPS( this->giveVertex(1) );
    const FloatArray &crackPE( this->giveVertex ( numSeg + 1 ) );

    double minDist = min( point.distance(crackPS), point.distance(crackPE) );
    double xiEl = 0.0;

    // Find the closest segment to determine the sign of the distance
    double minDistSeg = std :: numeric_limits< double > :: max();
    int minDistSegIndex = 0;
    double arcDistPassed = 0.0, minArcDist = 0.0;
    for ( int segId = 1; segId <= numSeg; segId++ ) {
        // Crack segment
        const FloatArray &crackP1( this->giveVertex ( segId ) );

        const FloatArray &crackP2( this->giveVertex ( segId + 1 ) );

        double distSeg = point.distance(crackP1, crackP2, xiEl);

        if ( distSeg < minDistSeg ) {
            minDistSeg = distSeg;
            minDistSegIndex = segId;

            if ( xiEl >= 0.0 && xiEl <= 1.0 ) {
                minArcDist = arcDistPassed + xiEl *crackP1.distance(crackP2) / totalArcLength;
            } else {
                if ( segId == 1 ) {
                    minArcDist = arcDistPassed + 0.0;
                } else if ( segId == numSeg ) {
                    minArcDist = arcDistPassed + 1.0;
                }
            }
        }

        arcDistPassed += crackP1.distance(crackP2) / totalArcLength;
    }

    oMinArcDist = minArcDist;

    if ( minDistSegIndex > 1 && minDistSegIndex < numSeg ) {
        // Interior segment is closest -> gamma is positive
        oDist = minDist;
        return;
    } else {
        if ( minDistSegIndex == 1 ) {
            const FloatArray &P1( this->giveVertex(1) );
            const FloatArray &P2( this->giveVertex(2) );

            FloatArray t;
            t.beDifferenceOf(P1, P2);

            t.normalize();

            FloatArray x_P1;
            x_P1.beDifferenceOf(point, P1);

            double sign = -x_P1.dotProduct(t);

            oDist = sign * minDist;
            return;
        } else if ( minDistSegIndex == numSeg ) {
            const FloatArray &P1( this->giveVertex ( minDistSegIndex ) );
            const FloatArray &P2( this->giveVertex ( minDistSegIndex + 1 ) );

            FloatArray t;
            t.beDifferenceOf(P2, P1);

            t.normalize();

            FloatArray x_P2;
            x_P2.beDifferenceOf(point, P2);

            double sign = -x_P2.dotProduct(t);

            oDist = sign * minDist;

            return;
        } else {
	  OOFEM_ERROR("PolygonLine :: computeTangentialSignDist: Index of minDistSegIndex not covered in loop: %d\n", minDistSegIndex);
        }
    }
}

void PolygonLine :: computeLocalCoordinates(FloatArray &oLocCoord, const FloatArray &iPoint) const
{ }

double PolygonLine :: computeLength() const
{
    double L = 0.0;

    size_t numSeg = mVertices.size() - 1;
    for ( size_t i = 0; i < numSeg; i++ ) {
        L += mVertices [ i ].distance(mVertices [ i + 1 ]);
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
        xSegEnd += mVertices [ i ].distance(mVertices [ i + 1 ]);

        xiSegStart = xSegStart / L;
        xiSegEnd        = xSegEnd / L;

        if ( iXiStart > xiSegStart && iXiStart < xiSegEnd ) {
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

        if ( xiSegEnd > iXiStart && xiSegEnd < iXiEnd ) {
            // End point of the segment is within range
            oPoints.push_back(mVertices [ i + 1 ]);
        }

        xSegStart = xSegEnd;
    }
}

IRResultType PolygonLine :: initializeFrom(InputRecord *ir)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro


    FloatArray *points = new FloatArray();
    IR_GIVE_FIELD(ir, * points, _IFT_PolygonLine_points);

    int numPoints = points->giveSize() / 2;

    for ( int i = 1; i <= numPoints; i++ ) {
        FloatArray *pos = new FloatArray(2);
        pos->at(1) = points->at(2 * ( i - 1 ) + 1);
        pos->at(2) = points->at( 2 * ( i   ) );
        mVertices.push_back(* pos);
        delete pos;
    }


#ifdef __BOOST_MODULE
    // Precompute bounding box to speed up calculation of intersection points.
    calcBoundingBox(LC, UC);
#endif

    delete points;

    return IRRT_OK;
}

void PolygonLine :: giveInputRecord(DynamicInputRecord &input)
{
    FloatArray points;
    int nVert = mVertices.size();
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

    //	printf("eLC: (%e, %e) eUC: (%e, %e) ", eLC.x(), eLC.y(), eUC.x(), eUC.y() );
    //	printf(" LC: (%e, %e)  UC: (%e, %e)\n",  LC.x(),  LC.y(),  UC.x(),  UC.y() );


    // Check if there is any chance of overlap
    if ( !bOverlap(eLC, eUC, LC, UC) ) {
        // If the bounding boxes do not overlap,
        // there will certainly not be any intersections
        return false;
    }


    return true;
}
#endif

bool PolygonLine :: intersects(ElementGeometry *element)
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


bool PolygonLine :: isInside(ElementGeometry *element)
{
    return false;
}

bool
PolygonLine :: isInside(FloatArray &point)
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


void PolygonLine :: computeIntersectionPoints(ElementGeometry *element, std :: vector< FloatArray > &oIntersectionPoints)
{
    printf("Warning: entering  PolygonLine :: computeIntersectionPoints(Element *element, std::vector< FloatArray > &oIntersectionPoints).\n");
#ifdef __BOOST_MODULE

    if ( !boundingBoxIntersects(element) ) {
        return;
    }


    for ( int i = 1; i <= element->giveNumberOfBoundarySides(); i++ ) {
        std :: vector< FloatArray >oneLineIntersects;
        FloatArray *a = new FloatArray( * ( element->giveDofManager ( i )->giveCoordinates() ) );
        FloatArray *b = NULL;
        if ( i != element->giveNumberOfBoundarySides() ) {
            b = new FloatArray( * ( element->giveDofManager ( i + 1 )->giveCoordinates() ) );
        } else {
            b = new FloatArray( * ( element->giveDofManager ( 1 )->giveCoordinates() ) );
        }

        Line l(a, b);

        computeIntersectionPoints(& l, oneLineIntersects);
        for ( int j = 1; j <= int ( oneLineIntersects.size() ); j++ ) {
            // Check that the intersection point has not already been identified.
            // This may happen if the crack intersects the element exactly at a node,
            // so that intersection is detected for both element edges in that node.

            // TODO: Set tolerance in a more transparent way.
            double distTol = 1.0e-9;

            bool alreadyFound = false;

            bPoint2 pNew( oneLineIntersects [ j - 1 ].at(1), oneLineIntersects [ j - 1 ].at(2) );

            int numPointsOld = oIntersectionPoints.size();
            for ( int k = 1; k <= numPointsOld; k++ ) {
                bPoint2 pOld( oIntersectionPoints [ k - 1 ].at(1), oIntersectionPoints [ k - 1 ].at(2) );

                if ( bDist(pOld, pNew) < distTol ) {
                    alreadyFound = true;
                    break;
                }
            }

            if ( !alreadyFound ) {
                //				int sz = intersecPoints->giveSize();
                //				intersecPoints->put( sz + 1, oneLineIntersects.at(j) );
                //				oneLineIntersects.unlink(j);
                oIntersectionPoints.push_back(oneLineIntersects [ j - 1 ]);
            }
        }
    }
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


                //				int sz = intersecPoints->giveSize();
                //				FloatArray *pos = new FloatArray(2);
                //				pos->at(1) = intersectionPoint.x();
                //				pos->at(2) = intersectionPoint.y();
                //				intersecPoints->put( sz + 1, pos );

                FloatArray pos(2);
                pos.at(1) = intersectionPoint.x();
                pos.at(2) = intersectionPoint.y();
                oIntersectionPoints.push_back(pos);
            }
        }
    }

#endif
}

int
PolygonLine :: computeNumberOfIntersectionPoints(ElementGeometry *element)
{
    std :: vector< FloatArray >intersecPoints;
    this->computeIntersectionPoints(element, intersecPoints);
    return intersecPoints.size();
}

bool PolygonLine :: isOutside(BasicGeometry *bg)
{
    return true;
}

void PolygonLine :: printYourself()
{
    printf("PolygonLine: start: ");
    mVertices [ 0 ].printYourself();
    printf(" end: ");
    mVertices [ 1 ].printYourself();
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



IRResultType PointSwarm :: initializeFrom(InputRecord *ir)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro
    IntArray idList;

    IR_GIVE_FIELD(ir, idList, _IFT_PointSwarm_nodeID); // Macro

    for ( int i = 1; i <= idList.giveSize(); i++ ) {
        this->idList.push_back( idList.at(i) );
    }
    return IRRT_OK;
}
} // end namespace oofem
