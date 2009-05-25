/* $Header: /home/cvs/bp/oofem/oofemlib/src/gaussintegrationrule.C,v 1.5.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


// file integrationRule.C

#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

// initialize class member

GaussIntegrationRule :: GaussIntegrationRule(int n, Element *e,
                                             int startIndx, int endIndx, bool dynamic) :
    IntegrationRule(n, e, startIndx, endIndx, dynamic) { }

GaussIntegrationRule :: GaussIntegrationRule(int n, Element *e) :
    IntegrationRule(n, e) { }

GaussIntegrationRule :: ~GaussIntegrationRule()
{ }



int
GaussIntegrationRule :: SetUpPointsOnLine(int nPoints, MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    int i;
    double weight;
    FloatArray *coord, *c, *w;

    switch ( nPoints ) {
    case 1:

        * arry = new GaussPoint * [ nPoints ];
        coord = new FloatArray(1);
        coord->at(1) = 0.0;
        weight = 2.0;
        ( * arry ) [ 0 ] = new GaussPoint(this, 1, coord, weight, mode);
        break;

    case 2:

        c = new FloatArray(2);
        w = new FloatArray(2);

        c->at(1) = -0.577350269189626;
        c->at(2) =  0.577350269189626;

        w->at(1) = 1.0;
        w->at(2) = 1.0;

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 2; i++ ) {
            coord = new FloatArray(1);
            coord->at(1) = c->at(i + 1);
            weight = w->at(i + 1);
            ( * arry ) [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
        }

        delete c;
        delete w;
        break;

    case 3:

        c = new FloatArray(3);
        w = new FloatArray(3);

        c->at(1) = -0.774596669241483;
        c->at(2) =  0.0;
        c->at(3) =  0.774596669241483;

        w->at(1) =  0.555555555555555;
        w->at(2) =  0.888888888888888;
        w->at(3) =  0.555555555555555;

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 3; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = c->at(i + 1);
            weight = w->at(i + 1);
            ( * arry ) [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
        }

        delete c;
        delete w;
        break;

    case 4:

        c = new FloatArray(4);
        w = new FloatArray(4);

        c->at(1) = -0.861136311594053;
        c->at(2) = -0.339981043584856;
        c->at(3) =  0.339981043584856;
        c->at(4) =  0.861136311594053;

        w->at(1) =  0.347854845137454;
        w->at(2) =  0.652145154862546;
        w->at(3) =  0.652145154862546;
        w->at(4) =  0.347854845137454;

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 4; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = c->at(i + 1);
            weight = w->at(i + 1);
            ( * arry ) [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
        }

        delete c;
        delete w;
        break;

    default:
        OOFEM_ERROR2("SetUpPointsOnLine: unsupported number of IPs (%d)", nPoints);
    }

    return nPoints;
}

int
GaussIntegrationRule :: SetUpPointsOnTriagle(int nPoints,
                                             MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    FloatArray *coord1;

    switch ( nPoints ) {
    case 1:

        * arry              = new GaussPoint * [ nPoints ];

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.333333333333;
        coord1->at(2)    = 0.333333333333;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, 0.5, mode);
        break;

    case 3:
        ( * arry )             = new GaussPoint * [ nPoints ];

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.166666666666666667;
        coord1->at(2)    = 0.166666666666666667;
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, 0.16666667, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.666666666666666667;
        coord1->at(2)    = 0.166666666666666667;
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, 0.16666667, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.166666666666666667;
        coord1->at(2)    = 0.666666666666666667;
        ( * arry ) [ 2 ]         = new GaussPoint(this, 3, coord1, 0.16666667, mode);
        break;

    case 4:

        * arry              = new GaussPoint * [ nPoints ];

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.2;
        coord1->at(2)    = 0.2;
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, 0.260416666666, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.6;
        coord1->at(2)    = 0.2;
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, 0.260416666666, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.2;
        coord1->at(2)    = 0.6;
        ( * arry ) [ 2 ]         = new GaussPoint(this, 3, coord1, 0.260416666666, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.333333333333333333;
        coord1->at(2)    = 0.333333333333333333;
        ( * arry ) [ 3 ]         = new GaussPoint(this, 4, coord1, -0.28125, mode);

        break;

    case 7:

        * arry              = new GaussPoint * [ nPoints ];

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.4701420641;
        coord1->at(2)    = 0.0597158717;
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, 0.06619705, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.4701420641;
        coord1->at(2)    = 0.4701420641;
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, 0.06619705, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.0597158717;
        coord1->at(2)    = 0.4701420641;
        ( * arry ) [ 2 ]         = new GaussPoint(this, 3, coord1, 0.06619705, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.1012865073;
        coord1->at(2)    = 0.1012865073;
        ( * arry ) [ 3 ]         = new GaussPoint(this, 4, coord1, 0.0629695902, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.7974269853;
        coord1->at(2)    = 0.1012865073;
        ( * arry ) [ 4 ]         = new GaussPoint(this, 5, coord1, 0.0629695902, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.1012865073;
        coord1->at(2)    = 0.7974269853;
        ( * arry ) [ 5 ]         = new GaussPoint(this, 6, coord1, 0.0629695902, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.333333333333333333;
        coord1->at(2)    = 0.333333333333333333;
        ( * arry ) [ 6 ]         = new GaussPoint(this, 7, coord1, 0.1125, mode);

        break;

    case 13:
        * arry              = new GaussPoint * [ nPoints ];

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.0651301029022;
        coord1->at(2)    = 0.0651301029022;
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, 0.0533472356088, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.8697397941956;
        coord1->at(2)    = 0.0651301029022;
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, 0.0533472356088, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.0651301029022;
        coord1->at(2)    = 0.8697397941956;
        ( * arry ) [ 2 ]         = new GaussPoint(this, 3, coord1, 0.0533472356088, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.3128654960049;
        coord1->at(2)    = 0.0486903154253;
        ( * arry ) [ 3 ]         = new GaussPoint(this, 4, coord1, 0.0771137608903, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.6384441885698;
        coord1->at(2)    = 0.3128654960049;
        ( * arry ) [ 4 ]         = new GaussPoint(this, 5, coord1, 0.0771137608903, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.0486903154253;
        coord1->at(2)    = 0.6384441885698;
        ( * arry ) [ 5 ]         = new GaussPoint(this, 6, coord1, 0.0771137608903, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.6384441885698;
        coord1->at(2)    = 0.0486903154253;
        ( * arry ) [ 6 ]         = new GaussPoint(this, 7, coord1, 0.0771137608903, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.3128654960049;
        coord1->at(2)    = 0.6384441885698;
        ( * arry ) [ 7 ]         = new GaussPoint(this, 8, coord1, 0.0771137608903, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.0486903154253;
        coord1->at(2)    = 0.3128654960049;
        ( * arry ) [ 8 ]         = new GaussPoint(this, 9, coord1, 0.0771137608903, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.2603459660790;
        coord1->at(2)    = 0.2603459660790;
        ( * arry ) [ 9 ]         = new GaussPoint(this, 10, coord1, 0.1756152576332, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.4793080678419;
        coord1->at(2)    = 0.2603459660790;
        ( * arry ) [ 10 ]         = new GaussPoint(this, 11, coord1, 0.1756152576332, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.2603459660790;
        coord1->at(2)    = 0.4793080678419;
        ( * arry ) [ 11 ]         = new GaussPoint(this, 12, coord1, 0.1756152576332, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.333333333333;
        coord1->at(2)    = 0.4793080678419;
        ( * arry ) [ 12 ]         = new GaussPoint(this, 13, coord1, -0.1495700444677, mode);

        break;

    default:
        OOFEM_ERROR2("SetUpPointsOnTriangle: unsupported number of IPs (%d)", nPoints);
    }

    return nPoints;
}

int
GaussIntegrationRule :: SetUpPointsOnSquare(int nPoints,
                                            MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    int i, j;
    double weight;
    FloatArray *coord, *c, *w;

    switch ( nPoints ) {
    case 1:

        * arry = new GaussPoint * [ nPoints ];
        coord = new FloatArray(2);
        coord->at(1) = 0.0;
        coord->at(2) = 0.0;
        weight = 4.0;
        ( * arry ) [ 0 ] = new GaussPoint(this, 1, coord, weight, mode);
        break;

    case 4:

        c = new FloatArray(2);
        w = new FloatArray(2);

        c->at(1) = -0.577350269189626;
        c->at(2) =  0.577350269189626;

        w->at(1) = 1.0;
        w->at(2) = 1.0;

        * arry  = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 2; i++ ) {
            for ( j = 0; j < 2; j++ ) {
                coord = new FloatArray(2);
                coord->at(1) = c->at(i + 1);
                coord->at(2) = c->at(j + 1);
                weight = w->at(i + 1) * w->at(j + 1);
                ( * arry ) [ 2 * i + j ] = new GaussPoint(this, 2 *i + j + 1, coord, weight, mode);
            }
        }

        delete c;
        delete w;
        break;

    case 9:

        c = new FloatArray(3);
        w = new FloatArray(3);

        c->at(1) = -0.774596669241483;
        c->at(2) =  0.0;
        c->at(3) =  0.774596669241483;

        w->at(1) =  0.555555555555555;
        w->at(2) =  0.888888888888888;
        w->at(3) =  0.555555555555555;

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 3; i++ ) {
            for ( j = 0; j < 3; j++ ) {
                coord  = new FloatArray(2);
                coord->at(1) = c->at(i + 1);
                coord->at(2) = c->at(j + 1);
                weight = w->at(i + 1) * w->at(j + 1);
                ( * arry ) [ 3 * i + j ] = new GaussPoint(this, 3 *i + j + 1, coord, weight, mode);
            }
        }

        delete c;
        delete w;
        break;

    case 16:

        c = new FloatArray(4);
        w = new FloatArray(4);

        c->at(1) = -0.861136311594053;
        c->at(2) = -0.339981043584856;
        c->at(3) =  0.339981043584856;
        c->at(4) =  0.861136311594053;

        w->at(1) =  0.347854845137454;
        w->at(2) =  0.652145154862546;
        w->at(3) =  0.652145154862546;
        w->at(4) =  0.347854845137454;

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 4; i++ ) {
            for ( j = 0; j < 4; j++ ) {
                coord  = new FloatArray(2);
                coord->at(1) = c->at(i + 1);
                coord->at(2) = c->at(j + 1);
                weight = w->at(i + 1) * w->at(j + 1);
                ( * arry ) [ 4 * i + j ] = new GaussPoint(this, 4 *i + j + 1, coord, weight, mode);
            }
        }

        delete c;
        delete w;
        break;

    case 64:

        c = new FloatArray(8);
        w = new FloatArray(8);

        c->at(1) = -0.960289856497536;
        c->at(2) = -0.796666477413627;
        c->at(3) = -0.525532409916329;
        c->at(4) = -0.183434642495650;
        c->at(5) = 0.183434642495650;
        c->at(6) = 0.525532409916329;
        c->at(7) = 0.796666477413627;
        c->at(8) = 0.96028985649753;

        w->at(1) = 0.101228536290375;
        w->at(2) = 0.222381034453374;
        w->at(3) = 0.313706645877887;
        w->at(4) = 0.362683783378362;
        w->at(5) = 0.362683783378362;
        w->at(6) = 0.313706645877887;
        w->at(7) = 0.222381034453374;
        w->at(8) = 0.10122853629037;

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 8; i++ ) {
            for ( j = 0; j < 8; j++ ) {
                coord  = new FloatArray(2);
                coord->at(1) = c->at(i + 1);
                coord->at(2) = c->at(j + 1);
                weight = w->at(i + 1) * w->at(j + 1);
                ( * arry ) [ 8 * i + j ] = new GaussPoint(this, 8 *i + j + 1, coord, weight, mode);
            }
        }

        delete c;
        delete w;
        break;

    default:
        OOFEM_ERROR2("SetUpPointsOnSquare: unsupported number of IPs (%d)", nPoints);
    }

    return nPoints;
}

int
GaussIntegrationRule :: SetUpPointsOnCube(int nPoints,
                                          MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    int i, j, k;
    double weight;
    FloatArray *coord, *c, *w;

    switch ( nPoints ) {
    case 1:

        * arry = new GaussPoint * [ nPoints ];
        coord = new FloatArray(3);
        coord->at(1) = 0.0;
        coord->at(2) = 0.0;
        coord->at(3) = 0.0;
        weight = 8.0;
        ( * arry ) [ 0 ] = new GaussPoint(this, 1, coord, weight, mode);
        break;

    case 8:

        c = new FloatArray(2);
        w = new FloatArray(2);

        c->at(1) = -0.577350269189626;
        c->at(2) =  0.577350269189626;

        w->at(1) = 1.0;
        w->at(2) = 1.0;

        * arry  = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 2; i++ ) {
            for ( j = 0; j < 2; j++ ) {
                for ( k = 0; k < 2; k++ ) {
                    coord = new FloatArray(3);
                    coord->at(1) = c->at(i + 1);
                    coord->at(2) = c->at(j + 1);
                    coord->at(3) = c->at(k + 1);
                    weight = w->at(i + 1) * w->at(j + 1) * w->at(k + 1);
                    ( * arry ) [ 4 * i + 2 * j + k ] = new GaussPoint(this, 4 *i + 2 *j + k + 1, coord, weight, mode);
                }
            }
        }

        delete c;
        delete w;
        break;

    case 27:

        c = new FloatArray(3);
        w = new FloatArray(3);

        c->at(1) = -0.774596669241483;
        c->at(2) =  0.0;
        c->at(3) =  0.774596669241483;

        w->at(1) =  0.555555555555555;
        w->at(2) =  0.888888888888888;
        w->at(3) =  0.555555555555555;

        * arry       = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 3; i++ ) {
            for ( j = 0; j < 3; j++ ) {
                for ( k = 0; k < 3; k++ ) {
                    coord  = new FloatArray(3);
                    coord->at(1) = c->at(i + 1);
                    coord->at(2) = c->at(j + 1);
                    coord->at(3) = c->at(k + 1);
                    weight = w->at(i + 1) * w->at(j + 1) * w->at(k + 1);
                    ( * arry ) [ 9 * i + 3 * j + k ] = new GaussPoint(this, 9 *i + 3 *j + k + 1, coord, weight, mode);
                }
            }
        }

        delete c;
        delete w;
        break;

    case 64:

        c = new FloatArray(4);
        w = new FloatArray(4);

        c->at(1) = -0.861136311594053;
        c->at(2) = -0.339981043584856;
        c->at(3) =  0.339981043584856;
        c->at(4) =  0.861136311594053;

        w->at(1) =  0.347854845137454;
        w->at(2) =  0.652145154862546;
        w->at(3) =  0.652145154862546;
        w->at(4) =  0.347854845137454;

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 4; i++ ) {
            for ( j = 0; j < 4; j++ ) {
                for ( k = 0; k < 4; k++ ) {
                    coord  = new FloatArray(3);
                    coord->at(1) = c->at(i + 1);
                    coord->at(2) = c->at(j + 1);
                    coord->at(3) = c->at(k + 1);
                    weight = w->at(i + 1) * w->at(j + 1) * w->at(k + 1);
                    ( * arry ) [ 16 * i + 4 * j + k ] = new GaussPoint(this, 16 *i + 4 *j + k + 1, coord, weight, mode);
                }
            }
        }

        delete c;
        delete w;
        break;

    default:
        OOFEM_ERROR2("SetUpPointsOnCube: unsupported number of IPs (%d)", nPoints);
    }

    return nPoints;
}


int
GaussIntegrationRule :: SetUpPointsOnTetrahedra(int nPoints,
                                                MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    // int         i,j,k;
    // double      weight;
    FloatArray *coord1;

    switch ( nPoints ) {
    case 1:

        * arry  = new GaussPoint * [ nPoints ];
        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.25;
        coord1->at(2)    = 0.25;
        coord1->at(3)    = 0.25;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, 1., mode);
        break;

    case 4:
        // quadratic formulae

        * arry  = new GaussPoint * [ nPoints ];

        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.58541020;
        coord1->at(2)    = 0.13819660;
        coord1->at(3)    = 0.13819660;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, 1. / 4., mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.13819660;
        coord1->at(2)    = 0.58541020;
        coord1->at(3)    = 0.13819660;
        ( * arry ) [ 1 ]        = new GaussPoint(this, 1, coord1, 1. / 4., mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.13819660;
        coord1->at(2)    = 0.13819660;
        coord1->at(3)    = 0.58541020;
        ( * arry ) [ 2 ]        = new GaussPoint(this, 1, coord1, 1. / 4., mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.13819660;
        coord1->at(2)    = 0.13819660;
        coord1->at(3)    = 0.13819660;
        ( * arry ) [ 3 ]        = new GaussPoint(this, 1, coord1, 1. / 4., mode);

        break;

    case 5:
        // cubic formulae

        * arry  = new GaussPoint * [ nPoints ];

        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.25;
        coord1->at(2)    = 0.25;
        coord1->at(3)    = 0.25;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, -4. / 5., mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.5;
        coord1->at(2)    = 1. / 6.;
        coord1->at(3)    = 1. / 6.;
        ( * arry ) [ 1 ]        = new GaussPoint(this, 1, coord1, 9. / 20., mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = 1. / 6.;
        coord1->at(2)    = 0.5;
        coord1->at(3)    = 1. / 6.;
        ( * arry ) [ 2 ]        = new GaussPoint(this, 1, coord1, 9. / 20., mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = 1. / 6.;
        coord1->at(2)    = 1. / 6.;
        coord1->at(3)    = 0.5;
        ( * arry ) [ 3 ]        = new GaussPoint(this, 1, coord1, 9. / 20., mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = 1. / 6.;
        coord1->at(2)    = 1. / 6.;
        coord1->at(3)    = 1. / 6.;
        ( * arry ) [ 4 ]        = new GaussPoint(this, 1, coord1, 9. / 20., mode);

        break;

    default:
        OOFEM_ERROR2("SetUpPointsOnTetrahedra: unsupported number of IPs (%d)", nPoints);
    }

    return nPoints;
}

int
GaussIntegrationRule :: getRequiredNumberOfIntegrationPoints(integrationDomain dType,
                                                             int approxOrder)
{
    int requiredNIP;
    if ( approxOrder < 0 ) {
        return 0;
    }

    switch ( dType ) {
    case _Line:
        requiredNIP = ( approxOrder + 1 ) / 2;
        if ( requiredNIP > 4 ) {
            return -1;
        }

        return requiredNIP;

    case _Triangle:
        if ( approxOrder <= 1 ) {
            return 1;
        }

        if ( approxOrder <= 3 ) {
            return 4;
        }

        if ( approxOrder <= 5 ) {
            return 7;
        }

        return -1;

    case _Square:
        requiredNIP = max( ( approxOrder + 1 ) / 2, 2 );
        requiredNIP *= requiredNIP;
        if ( requiredNIP > 16 ) {
            return -1;
        }

        return requiredNIP;

    case _Cube:
        requiredNIP = max( ( approxOrder + 1 ) / 2, 2 );
        requiredNIP *= requiredNIP * requiredNIP;
        if ( requiredNIP > 64 ) {
            return -1;
        }

        return requiredNIP;

    case _Tetrahedra:
        if ( approxOrder <= 1 ) {
            return 1;
        }

        if ( approxOrder <= 2 ) {
            return 4;
        }

        if ( approxOrder <= 3 ) {
            return 5;
        }

        return -1;

    default:
        OOFEM_ERROR("GaussIntegrationRule::setUpIntegrationPoints - unknown integrationDomain");
    }

    return -1;
}


int
GaussIntegrationRule :: SetUpPointsOn2DEmbeddedLine(int nPoints, MaterialMode mode, GaussPoint ***arry,
                                                    const FloatArray **coords)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    double weight, l;
    FloatArray *coord1;


    switch ( nPoints ) {
    case 1:
        * arry = new GaussPoint * [ nPoints ];
        coord1 = new FloatArray(2);
        l = 0.0;       //local coordinate on the line
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        weight = 2.0;
        ( * arry ) [ 0 ] = new GaussPoint(this, 1, coord1, weight, mode);
        break;

    case 2:

        * arry              = new GaussPoint * [ nPoints ];
        coord1             = new FloatArray(2);

        l = -0.577350269189626; //local coordinate on the line
        weight = 1.0;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l = 0.577350269189626;        //local coordinate on the line
        weight = 1.0;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, weight, mode);


        break;

    case 3:

        * arry              = new GaussPoint * [ nPoints ];
        coord1             = new FloatArray(2);

        l = -0.774596669241483; //local coordinate on the line
        weight =  0.555555555555555;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l = 0.0;        //local coordinate on the line
        weight = 0.888888888888888;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l = 0.774596669241483; //local coordinate on the line
        weight =  0.555555555555555;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( * arry ) [ 2 ]         = new GaussPoint(this, 3, coord1, weight, mode);
        break;

    case 4:


        * arry              = new GaussPoint * [ nPoints ];
        coord1             = new FloatArray(2);

        l = -0.861136311594053; //local coordinate on the line
        weight =  0.347854845137454;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l = -0.339981043584856;        //local coordinate on the line
        weight = 0.652145154862546;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l =  0.339981043584856; //local coordinate on the line
        weight =   0.652145154862546;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( * arry ) [ 2 ]         = new GaussPoint(this, 3, coord1, weight, mode);

        l =  0.861136311594053;        //local coordinate on the line
        weight =   0.347854845137454;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( * arry ) [ 3 ]         = new GaussPoint(this, 4, coord1, weight, mode);
        break;


    default:
        printf("SetUpPointsOnLine: such order of integration is not suported\n");
        exit(1);
    }

    // Initialize the material status at each Gauss point
    //for (int i=0; i<numberOfIntegrationPoints; i++)
    //  mat->giveStatus(gaussPointArray[i]);

    return nPoints;
}
