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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "delaunay.h"
#include "flotarry.h"
#include "intarray.h"
#include "alist.h"
#include "geometry.h"
#include "node.h"
#include "mathfem.h"

#include <map>

namespace oofem {
bool Delaunay :: colinear(FloatArray *p1, FloatArray *p2, FloatArray *p3)
{
    double dist = p1->at(1) * ( p2->at(2) - p3->at(2) ) + p2->at(1) * ( p3->at(2) - p1->at(2) ) +
                  p3->at(1) * ( p1->at(2) - p2->at(2) );
    // the tolerance probably needs a setter
    if ( dist < 0.0001 && dist > ( -1 ) * 0.0001 ) {
        return true;
    } else {
        return false;
    }
}

void Delaunay :: printTriangles(AList< Triangle > *triangles)
{
    for ( int i = 1; i <= triangles->giveSize(); i++ ) {
        triangles->at(i)->printYourself();
    }
}

bool Delaunay :: isInsideCC(FloatArray *p, FloatArray *p1,  FloatArray *p2,  FloatArray *p3)
{
    FloatArray *nodesCopy1 = new FloatArray(*p1);
    FloatArray *nodesCopy2 = new FloatArray(*p2);
    FloatArray *nodesCopy3 = new FloatArray(*p3);
    Triangle *tr = new Triangle(nodesCopy1, nodesCopy2, nodesCopy3);
    double r = tr->getRadiusOfCircumCircle();
    FloatArray circumCenter;
    tr->computeCenterOfCircumCircle(circumCenter);
    double distance = circumCenter.distance(p);
    delete tr;
    if ( distance < r ) {
        return true;
    } else {
        return false;
    }
}

void Delaunay :: triangulate(AList< FloatArray > *overtices, AList< Triangle > *triangles) {
    /// 4th order algorithm - four loops, only for testing puropses
    int n = overtices->giveSize();
    int count = 0;
    /// copy of vertices, since they will be shifted
    AList< FloatArray > *vertices = new AList< FloatArray >();
    std :: map< FloatArray *, FloatArray * >backToOld; // map for putting the shifted vertices into the old position
    for ( int i = 1; i <= overtices->giveSize(); i++ ) {
        FloatArray *ip = overtices->at(i);
        FloatArray *ipCopy = new FloatArray(*ip);
        vertices->put(i, ipCopy);
        backToOld [ ipCopy ] = ip;
    }

    // small shift of vertices
    for ( int i = 1; i <= n; i++ ) {
        vertices->at(i)->at(1) += vertices->at(i)->at(1) * 0.000001 * double ( rand() ) / RAND_MAX;
        vertices->at(i)->at(2) += vertices->at(i)->at(2) * 0.000001 * double ( rand() ) / RAND_MAX;
    }

    for ( int i = 1; i <= n; i++ ) {
        for ( int j = i + 1; j <= n; j++ ) {
            for ( int k = j + 1; k <= n; k++ ) {
                bool isTriangle = true;
                if ( colinear( vertices->at(i), vertices->at(j),
                              vertices->at(k) ) ) {
                    isTriangle = false;
                } else {
                    for ( int a = 1; a <= n; a++ ) {
                        if ( a != i && a != j && a != k ) {
                            // checks whether a point a is inside a circumcircle of a triangle ijk
                            if ( isInsideCC( vertices->at(a), vertices->at(i), vertices->at(j),
                                            vertices->at(k) ) ) {
                                isTriangle = false;
                                break;
                            }
                        }
                    }
                }

                if ( isTriangle ) {
                    count++;
                    // here we switch to old vertices
                    FloatArray *p1 = new FloatArray();
                    * p1 = * ( backToOld [ vertices->at(i) ] );
                    FloatArray *p2 = new FloatArray();
                    * p2 = * ( backToOld [ vertices->at(j) ] );
                    FloatArray *p3 = new FloatArray();
                    * p3 = * ( backToOld [ vertices->at(k) ] );
                    Triangle *triangle = new Triangle(p1, p2, p3);
                    if ( !triangle->isOrientedAnticlockwise() ) {
                        triangle->changeToAnticlockwise();
                    }

                    triangles->put(count, triangle);
                }
            }
        }
    }

    delete vertices;
}
} // end namespace oofem
