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

#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "mathfem.h"

namespace oofem {
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

    case 8:
    {
        double coord_rule8 [ 8 ] = {
            -0.960289856497536,
            -0.796666477413627,
            -0.525532409916329,
            -0.183434642495650,
            0.183434642495650,
            0.525532409916329,
            0.796666477413627,
            0.960289856497536
        };

        double weight_rule8 [ 8 ] = {
            0.101228536290375,
            0.222381034453374,
            0.313706645877887,
            0.362683783378362,
            0.362683783378362,
            0.313706645877887,
            0.222381034453374,
            0.101228536290375
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 8; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = coord_rule8 [ i ];
            weight = weight_rule8 [ i ];
            ( * arry ) [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
        }
    }

    break;

    case 16:
    {
        double coord_rule16 [ 16 ] = {
            -0.989400934991650,
            -0.944575023073233,
            -0.865631202387832,
            -0.755404408355003,
            -0.617876244402644,
            -0.458016777657227,
            -0.281603550779259,
            -0.095012509837637,
            0.095012509837637,
            0.281603550779259,
            0.458016777657227,
            0.617876244402644,
            0.755404408355003,
            0.865631202387832,
            0.944575023073233,
            0.989400934991650
        };

        double weight_rule16 [ 16 ] = {
            0.027152459411753,
            0.062253523938647,
            0.095158511682492,
            0.124628971255534,
            0.149595988816577,
            0.169156519395003,
            0.182603415044924,
            0.189450610455068,
            0.189450610455068,
            0.182603415044924,
            0.169156519395003,
            0.149595988816577,
            0.124628971255534,
            0.095158511682492,
            0.062253523938647,
            0.027152459411753
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 16; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = coord_rule16 [ i ];
            weight = weight_rule16 [ i ];
            ( * arry ) [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
        }
    }

    break;

    case 24:
    {
        double coord_rule24 [ 24 ] = {
            -0.995187219997021,
            -0.974728555971309,
            -0.938274552002733,
            -0.886415527004401,
            -0.820001985973903,
            -0.740124191578554,
            -0.648093651936976,
            -0.545421471388840,
            -0.433793507626045,
            -0.315042679696163,
            -0.191118867473616,
            -0.064056892862606,
            0.064056892862606,
            0.191118867473616,
            0.315042679696163,
            0.433793507626045,
            0.545421471388840,
            0.648093651936976,
            0.740124191578554,
            0.820001985973903,
            0.886415527004401,
            0.938274552002733,
            0.974728555971309,
            0.995187219997021
        };

        double weight_rule24 [ 24 ] = {
            0.012341229799991,
            0.028531388628934,
            0.044277438817419,
            0.059298584915436,
            0.073346481411080,
            0.086190161531953,
            0.097618652104114,
            0.107444270115966,
            0.115505668053726,
            0.121670472927803,
            0.125837456346828,
            0.127938195346752,
            0.127938195346752,
            0.125837456346828,
            0.121670472927803,
            0.115505668053726,
            0.107444270115966,
            0.097618652104114,
            0.086190161531953,
            0.073346481411080,
            0.059298584915436,
            0.044277438817419,
            0.028531388628934,
            0.012341229799991
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 24; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = coord_rule24 [ i ];
            weight = weight_rule24 [ i ];
            ( * arry ) [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
        }
    }

    break;

    case 32:
    {
        double coord_rule32 [ 32 ] = {
            -0.997263861849482,
            -0.985611511545268,
            -0.964762255587506,
            -0.934906075937740,
            -0.896321155766052,
            -0.849367613732570,
            -0.794483795967942,
            -0.732182118740290,
            -0.663044266930215,
            -0.587715757240762,
            -0.506899908932229,
            -0.421351276130635,
            -0.331868602282128,
            -0.239287362252137,
            -0.144471961582796,
            -0.048307665687738,
            0.048307665687738,
            0.144471961582796,
            0.239287362252137,
            0.331868602282128,
            0.421351276130635,
            0.506899908932229,
            0.587715757240762,
            0.663044266930215,
            0.732182118740290,
            0.794483795967942,
            0.849367613732570,
            0.896321155766052,
            0.934906075937740,
            0.964762255587506,
            0.985611511545268,
            0.997263861849482
        };

        double weight_rule32 [ 32 ] = {
            0.007018610009469,
            0.016274394730905,
            0.025392065309263,
            0.034273862913022,
            0.042835898022227,
            0.050998059262376,
            0.058684093478536,
            0.065822222776362,
            0.072345794108848,
            0.078193895787070,
            0.083311924226947,
            0.087652093004404,
            0.091173878695764,
            0.093844399080805,
            0.095638720079275,
            0.096540088514728,
            0.096540088514728,
            0.095638720079275,
            0.093844399080805,
            0.091173878695764,
            0.087652093004404,
            0.083311924226947,
            0.078193895787070,
            0.072345794108848,
            0.065822222776362,
            0.058684093478536,
            0.050998059262376,
            0.042835898022227,
            0.034273862913022,
            0.025392065309263,
            0.016274394730905,
            0.007018610009469
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 32; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = coord_rule32 [ i ];
            weight = weight_rule32 [ i ];
            ( * arry ) [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
        }
    }

    break;

    case 64:
    {
        double coord_rule64 [ 64 ] = {
            -0.999305041735772,
            -0.996340116771955,
            -0.991013371476744,
            -0.983336253884626,
            -0.973326827789911,
            -0.961008799652054,
            -0.946411374858403,
            -0.929569172131940,
            -0.910522137078503,
            -0.889315445995114,
            -0.865999398154093,
            -0.840629296252580,
            -0.813265315122798,
            -0.783972358943341,
            -0.752819907260532,
            -0.719881850171611,
            -0.685236313054233,
            -0.648965471254657,
            -0.611155355172393,
            -0.571895646202634,
            -0.531279464019895,
            -0.489403145707053,
            -0.446366017253464,
            -0.402270157963992,
            -0.357220158337668,
            -0.311322871990211,
            -0.264687162208767,
            -0.217423643740007,
            -0.169644420423993,
            -0.121462819296121,
            -0.072993121787799,
            -0.024350292663424,
            0.024350292663424,
            0.072993121787799,
            0.121462819296121,
            0.169644420423993,
            0.217423643740007,
            0.264687162208767,
            0.311322871990211,
            0.357220158337668,
            0.402270157963992,
            0.446366017253464,
            0.489403145707053,
            0.531279464019895,
            0.571895646202634,
            0.611155355172393,
            0.648965471254657,
            0.685236313054233,
            0.719881850171611,
            0.752819907260532,
            0.783972358943341,
            0.813265315122798,
            0.840629296252580,
            0.865999398154093,
            0.889315445995114,
            0.910522137078503,
            0.929569172131940,
            0.946411374858403,
            0.961008799652054,
            0.973326827789911,
            0.983336253884626,
            0.991013371476744,
            0.996340116771955,
            0.999305041735772
        };

        double weight_rule64 [ 64 ] = {
            0.001783280721693,
            0.004147033260559,
            0.006504457968980,
            0.008846759826363,
            0.011168139460130,
            0.013463047896718,
            0.015726030476025,
            0.017951715775697,
            0.020134823153530,
            0.022270173808383,
            0.024352702568711,
            0.026377469715055,
            0.028339672614260,
            0.030234657072403,
            0.032057928354851,
            0.033805161837142,
            0.035472213256882,
            0.037055128540240,
            0.038550153178615,
            0.039953741132721,
            0.041262563242623,
            0.042473515123654,
            0.043583724529323,
            0.044590558163757,
            0.045491627927418,
            0.046284796581314,
            0.046968182816210,
            0.047540165714830,
            0.047999388596458,
            0.048344762234803,
            0.048575467441503,
            0.048690957009140,
            0.048690957009140,
            0.048575467441503,
            0.048344762234803,
            0.047999388596458,
            0.047540165714830,
            0.046968182816210,
            0.046284796581314,
            0.045491627927418,
            0.044590558163757,
            0.043583724529323,
            0.042473515123654,
            0.041262563242623,
            0.039953741132721,
            0.038550153178615,
            0.037055128540240,
            0.035472213256882,
            0.033805161837142,
            0.032057928354851,
            0.030234657072403,
            0.028339672614260,
            0.026377469715055,
            0.024352702568711,
            0.022270173808383,
            0.020134823153530,
            0.017951715775697,
            0.015726030476025,
            0.013463047896718,
            0.011168139460130,
            0.008846759826363,
            0.006504457968980,
            0.004147033260559,
            0.001783280721693
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 64; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = coord_rule64 [ i ];
            weight = weight_rule64 [ i ];
            ( * arry ) [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
        }
    }

    break;

    default:
        OOFEM_ERROR2("SetUpPointsOnLine: unsupported number of IPs (%d)", nPoints);
    }

    return nPoints;
}

int
GaussIntegrationRule :: SetUpPointsOnTriangle(int nPoints,
                                             MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store coordinates and weights)
{
    FloatArray *coord1;

    switch ( nPoints ) {
    case 1:

        * arry              = new GaussPoint * [ nPoints ];

        coord1             = new FloatArray(2);
        coord1->at(1)    = 1. / 3.;
        coord1->at(2)    = 1. / 3.;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, 0.5, mode);
        break;

    case 3:
        ( * arry )             = new GaussPoint * [ nPoints ];

        coord1             = new FloatArray(2);
        coord1->at(1)    = 1. / 6.;
        coord1->at(2)    = 1. / 6.;
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, 1. / 6., mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 2. / 3.;
        coord1->at(2)    = 1. / 6.;
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, 1. / 6., mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 1. / 6.;
        coord1->at(2)    = 2. / 3.;
        ( * arry ) [ 2 ]         = new GaussPoint(this, 3, coord1, 1. / 6., mode);
        break;

    case 4:

        * arry              = new GaussPoint * [ nPoints ];

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.2;
        coord1->at(2)    = 0.2;
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, 25. / 96., mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.6;
        coord1->at(2)    = 0.2;
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, 25. / 96., mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.2;
        coord1->at(2)    = 0.6;
        ( * arry ) [ 2 ]         = new GaussPoint(this, 3, coord1, 25. / 96., mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 1. / 3.;
        coord1->at(2)    = 1. / 3.;
        ( * arry ) [ 3 ]         = new GaussPoint(this, 4, coord1, - 9. / 32., mode);

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
        coord1->at(1)    = 1. / 3.;
        coord1->at(2)    = 1. / 3.;
        ( * arry ) [ 6 ]         = new GaussPoint(this, 7, coord1, 9. / 80., mode);

        break;

    case 13:
        * arry              = new GaussPoint * [ nPoints ];

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.0651301029022;
        coord1->at(2)    = 0.0651301029022;
        ( * arry ) [ 0 ]         = new GaussPoint(this, 1, coord1, 0.0266736178044, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.8697397941956;
        coord1->at(2)    = 0.0651301029022;
        ( * arry ) [ 1 ]         = new GaussPoint(this, 2, coord1, 0.0266736178044, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.0651301029022;
        coord1->at(2)    = 0.8697397941956;
        ( * arry ) [ 2 ]         = new GaussPoint(this, 3, coord1, 0.0266736178044, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.3128654960049;
        coord1->at(2)    = 0.0486903154253;
        ( * arry ) [ 3 ]         = new GaussPoint(this, 4, coord1, 0.0385568804452, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.6384441885698;
        coord1->at(2)    = 0.3128654960049;
        ( * arry ) [ 4 ]         = new GaussPoint(this, 5, coord1, 0.0385568804452, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.0486903154253;
        coord1->at(2)    = 0.6384441885698;
        ( * arry ) [ 5 ]         = new GaussPoint(this, 6, coord1, 0.0385568804452, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.6384441885698;
        coord1->at(2)    = 0.0486903154253;
        ( * arry ) [ 6 ]         = new GaussPoint(this, 7, coord1, 0.0385568804452, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.3128654960049;
        coord1->at(2)    = 0.6384441885698;
        ( * arry ) [ 7 ]         = new GaussPoint(this, 8, coord1, 0.0385568804452, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.0486903154253;
        coord1->at(2)    = 0.3128654960049;
        ( * arry ) [ 8 ]         = new GaussPoint(this, 9, coord1, 0.0385568804452, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.2603459660790;
        coord1->at(2)    = 0.2603459660790;
        ( * arry ) [ 9 ]         = new GaussPoint(this, 10, coord1, 0.0878076288166, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.4793080678419;
        coord1->at(2)    = 0.2603459660790;
        ( * arry ) [ 10 ]         = new GaussPoint(this, 11, coord1, 0.0878076288166, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.2603459660790;
        coord1->at(2)    = 0.4793080678419;
        ( * arry ) [ 11 ]         = new GaussPoint(this, 12, coord1, 0.0878076288166, mode);

        coord1             = new FloatArray(2);
        coord1->at(1)    = 0.333333333333;
        coord1->at(2)    = 0.333333333333;
        ( * arry ) [ 12 ]         = new GaussPoint(this, 13, coord1, -0.0747850222339, mode);

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


    case 256:

        c = new FloatArray(16);
        w = new FloatArray(16);

        c->at(1) =  -0.989400934991650;
        c->at(2) =  -0.944575023073233;
        c->at(3) =  -0.865631202387832;
        c->at(4) =  -0.755404408355003;
        c->at(5) =  -0.617876244402644;
        c->at(6) =  -0.458016777657227;
        c->at(7) =  -0.281603550779259;
        c->at(8) =  -0.095012509837637;
        c->at(9) =   0.095012509837637;
        c->at(10) =   0.281603550779259;
        c->at(11) =   0.458016777657227;
        c->at(12) =   0.617876244402644;
        c->at(13) =   0.755404408355003;
        c->at(14) =   0.865631202387832;
        c->at(15) =   0.944575023073233;
        c->at(16) =   0.989400934991650;

        w->at(1) = 0.027152459411753;
        w->at(2) = 0.062253523938647;
        w->at(3) = 0.095158511682492;
        w->at(4) = 0.124628971255534;
        w->at(5) = 0.149595988816577;
        w->at(6) = 0.169156519395003;
        w->at(7) = 0.182603415044924;
        w->at(8) = 0.189450610455068;
        w->at(9) = 0.189450610455068;
        w->at(10) = 0.182603415044924;
        w->at(11) = 0.169156519395003;
        w->at(12) = 0.149595988816577;
        w->at(13) = 0.124628971255534;
        w->at(14) = 0.095158511682492;
        w->at(15) = 0.062253523938647;
        w->at(16) = 0.027152459411753;

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 16; i++ ) {
            for ( j = 0; j < 16; j++ ) {
                coord  = new FloatArray(2);
                coord->at(1) = c->at(i + 1);
                coord->at(2) = c->at(j + 1);
                weight = w->at(i + 1) * w->at(j + 1);
                ( * arry ) [ 16 * i + j ] = new GaussPoint(this, 16 *i + j + 1, coord, weight, mode);
            }
        }

        delete c;
        delete w;
        break;

    case 576:
    {
        double coord_rule24 [ 24 ] = {
            -0.995187219997021,
            -0.974728555971309,
            -0.938274552002733,
            -0.886415527004401,
            -0.820001985973903,
            -0.740124191578554,
            -0.648093651936976,
            -0.545421471388840,
            -0.433793507626045,
            -0.315042679696163,
            -0.191118867473616,
            -0.064056892862606,
            0.064056892862606,
            0.191118867473616,
            0.315042679696163,
            0.433793507626045,
            0.545421471388840,
            0.648093651936976,
            0.740124191578554,
            0.820001985973903,
            0.886415527004401,
            0.938274552002733,
            0.974728555971309,
            0.995187219997021
        };

        double weight_rule24 [ 24 ] = {
            0.012341229799991,
            0.028531388628934,
            0.044277438817419,
            0.059298584915436,
            0.073346481411080,
            0.086190161531953,
            0.097618652104114,
            0.107444270115966,
            0.115505668053726,
            0.121670472927803,
            0.125837456346828,
            0.127938195346752,
            0.127938195346752,
            0.125837456346828,
            0.121670472927803,
            0.115505668053726,
            0.107444270115966,
            0.097618652104114,
            0.086190161531953,
            0.073346481411080,
            0.059298584915436,
            0.044277438817419,
            0.028531388628934,
            0.012341229799991
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 24; i++ ) {
            for ( j = 0; j < 24; j++ ) {
                coord  = new FloatArray(2);
                coord->at(1) = coord_rule24 [ i ];
                coord->at(2) = coord_rule24 [ j ];
                weight = weight_rule24 [ i ] * weight_rule24 [ j ];
                ( * arry ) [ 24 * i + j ] = new GaussPoint(this, 24 *i + j + 1, coord, weight, mode);
            }
        }
    }

    break;

    case 1024:
    {
        double coord_rule32 [ 32 ] = {
            -0.997263861849482,
            -0.985611511545268,
            -0.964762255587506,
            -0.934906075937740,
            -0.896321155766052,
            -0.849367613732570,
            -0.794483795967942,
            -0.732182118740290,
            -0.663044266930215,
            -0.587715757240762,
            -0.506899908932229,
            -0.421351276130635,
            -0.331868602282128,
            -0.239287362252137,
            -0.144471961582796,
            -0.048307665687738,
            0.048307665687738,
            0.144471961582796,
            0.239287362252137,
            0.331868602282128,
            0.421351276130635,
            0.506899908932229,
            0.587715757240762,
            0.663044266930215,
            0.732182118740290,
            0.794483795967942,
            0.849367613732570,
            0.896321155766052,
            0.934906075937740,
            0.964762255587506,
            0.985611511545268,
            0.997263861849482
        };

        double weight_rule32 [ 32 ] = {
            0.007018610009469,
            0.016274394730905,
            0.025392065309263,
            0.034273862913022,
            0.042835898022227,
            0.050998059262376,
            0.058684093478536,
            0.065822222776362,
            0.072345794108848,
            0.078193895787070,
            0.083311924226947,
            0.087652093004404,
            0.091173878695764,
            0.093844399080805,
            0.095638720079275,
            0.096540088514728,
            0.096540088514728,
            0.095638720079275,
            0.093844399080805,
            0.091173878695764,
            0.087652093004404,
            0.083311924226947,
            0.078193895787070,
            0.072345794108848,
            0.065822222776362,
            0.058684093478536,
            0.050998059262376,
            0.042835898022227,
            0.034273862913022,
            0.025392065309263,
            0.016274394730905,
            0.007018610009469
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 32; i++ ) {
            for ( j = 0; j < 32; j++ ) {
                coord  = new FloatArray(2);
                coord->at(1) = coord_rule32 [ i ];
                coord->at(2) = coord_rule32 [ j ];
                weight = weight_rule32 [ i ] * weight_rule32 [ j ];
                ( * arry ) [ 32 * i + j ] = new GaussPoint(this, 32 *i + j + 1, coord, weight, mode);
            }
        }
    }

    break;

    case 4096:
    {
        double coord_rule64 [ 64 ] = {
            -0.999305041735772,
            -0.996340116771955,
            -0.991013371476744,
            -0.983336253884626,
            -0.973326827789911,
            -0.961008799652054,
            -0.946411374858403,
            -0.929569172131940,
            -0.910522137078503,
            -0.889315445995114,
            -0.865999398154093,
            -0.840629296252580,
            -0.813265315122798,
            -0.783972358943341,
            -0.752819907260532,
            -0.719881850171611,
            -0.685236313054233,
            -0.648965471254657,
            -0.611155355172393,
            -0.571895646202634,
            -0.531279464019895,
            -0.489403145707053,
            -0.446366017253464,
            -0.402270157963992,
            -0.357220158337668,
            -0.311322871990211,
            -0.264687162208767,
            -0.217423643740007,
            -0.169644420423993,
            -0.121462819296121,
            -0.072993121787799,
            -0.024350292663424,
            0.024350292663424,
            0.072993121787799,
            0.121462819296121,
            0.169644420423993,
            0.217423643740007,
            0.264687162208767,
            0.311322871990211,
            0.357220158337668,
            0.402270157963992,
            0.446366017253464,
            0.489403145707053,
            0.531279464019895,
            0.571895646202634,
            0.611155355172393,
            0.648965471254657,
            0.685236313054233,
            0.719881850171611,
            0.752819907260532,
            0.783972358943341,
            0.813265315122798,
            0.840629296252580,
            0.865999398154093,
            0.889315445995114,
            0.910522137078503,
            0.929569172131940,
            0.946411374858403,
            0.961008799652054,
            0.973326827789911,
            0.983336253884626,
            0.991013371476744,
            0.996340116771955,
            0.999305041735772
        };

        double weight_rule64 [ 64 ] = {
            0.001783280721693,
            0.004147033260559,
            0.006504457968980,
            0.008846759826363,
            0.011168139460130,
            0.013463047896718,
            0.015726030476025,
            0.017951715775697,
            0.020134823153530,
            0.022270173808383,
            0.024352702568711,
            0.026377469715055,
            0.028339672614260,
            0.030234657072403,
            0.032057928354851,
            0.033805161837142,
            0.035472213256882,
            0.037055128540240,
            0.038550153178615,
            0.039953741132721,
            0.041262563242623,
            0.042473515123654,
            0.043583724529323,
            0.044590558163757,
            0.045491627927418,
            0.046284796581314,
            0.046968182816210,
            0.047540165714830,
            0.047999388596458,
            0.048344762234803,
            0.048575467441503,
            0.048690957009140,
            0.048690957009140,
            0.048575467441503,
            0.048344762234803,
            0.047999388596458,
            0.047540165714830,
            0.046968182816210,
            0.046284796581314,
            0.045491627927418,
            0.044590558163757,
            0.043583724529323,
            0.042473515123654,
            0.041262563242623,
            0.039953741132721,
            0.038550153178615,
            0.037055128540240,
            0.035472213256882,
            0.033805161837142,
            0.032057928354851,
            0.030234657072403,
            0.028339672614260,
            0.026377469715055,
            0.024352702568711,
            0.022270173808383,
            0.020134823153530,
            0.017951715775697,
            0.015726030476025,
            0.013463047896718,
            0.011168139460130,
            0.008846759826363,
            0.006504457968980,
            0.004147033260559,
            0.001783280721693
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 64; i++ ) {
            for ( j = 0; j < 64; j++ ) {
                coord  = new FloatArray(2);
                coord->at(1) = coord_rule64 [ i ];
                coord->at(2) = coord_rule64 [ j ];
                weight = weight_rule64 [ i ] * weight_rule64 [ j ];
                ( * arry ) [ 64 * i + j ] = new GaussPoint(this, 64 *i + j + 1, coord, weight, mode);
            }
        }
    }
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

    case 512:
    {
        double coord_rule8 [ 8 ] = {
            -0.960289856497536,
            -0.796666477413627,
            -0.525532409916329,
            -0.183434642495650,
            0.183434642495650,
            0.525532409916329,
            0.796666477413627,
            0.960289856497536
        };

        double weight_rule8 [ 8 ] = {
            0.101228536290375,
            0.222381034453374,
            0.313706645877887,
            0.362683783378362,
            0.362683783378362,
            0.313706645877887,
            0.222381034453374,
            0.101228536290375
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 8; i++ ) {
            for ( j = 0; j < 8; j++ ) {
                for ( k = 0; k < 8; k++ ) {
                    coord  = new FloatArray(3);
                    coord->at(1) = coord_rule8 [ i ];
                    coord->at(2) = coord_rule8 [ j ];
                    coord->at(3) = coord_rule8 [ k ];
                    weight = weight_rule8 [ i ] * weight_rule8 [ j ] * weight_rule8 [ k ];
                    ( * arry ) [ 64 * i + 8 * j + k ] = new GaussPoint(this, 64 *i + 8 *j + k + 1, coord, weight, mode);
                }
            }
        }
    }

    break;

    case 4096:
    {
        double coord_rule16 [ 16 ] = {
            -0.989400934991650,
            -0.944575023073233,
            -0.865631202387832,
            -0.755404408355003,
            -0.617876244402644,
            -0.458016777657227,
            -0.281603550779259,
            -0.095012509837637,
            0.095012509837637,
            0.281603550779259,
            0.458016777657227,
            0.617876244402644,
            0.755404408355003,
            0.865631202387832,
            0.944575023073233,
            0.989400934991650
        };

        double weight_rule16 [ 16 ] = {
            0.027152459411753,
            0.062253523938647,
            0.095158511682492,
            0.124628971255534,
            0.149595988816577,
            0.169156519395003,
            0.182603415044924,
            0.189450610455068,
            0.189450610455068,
            0.182603415044924,
            0.169156519395003,
            0.149595988816577,
            0.124628971255534,
            0.095158511682492,
            0.062253523938647,
            0.027152459411753
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 16; i++ ) {
            for ( j = 0; j < 16; j++ ) {
                for ( k = 0; k < 16; k++ ) {
                    coord  = new FloatArray(3);
                    coord->at(1) = coord_rule16 [ i ];
                    coord->at(2) = coord_rule16 [ j ];
                    coord->at(3) = coord_rule16 [ k ];
                    weight = weight_rule16 [ i ] * weight_rule16 [ j ] * weight_rule16 [ k ];
                    ( * arry ) [ 256 * i + 16 * j + k ] = new GaussPoint(this, 256 *i + 16 *j + k + 1, coord, weight, mode);
                }
            }
        }
    }

    break;

    case 13824:
    {
        double coord_rule24 [ 24 ] = {
            -0.995187219997021,
            -0.974728555971309,
            -0.938274552002733,
            -0.886415527004401,
            -0.820001985973903,
            -0.740124191578554,
            -0.648093651936976,
            -0.545421471388840,
            -0.433793507626045,
            -0.315042679696163,
            -0.191118867473616,
            -0.064056892862606,
            0.064056892862606,
            0.191118867473616,
            0.315042679696163,
            0.433793507626045,
            0.545421471388840,
            0.648093651936976,
            0.740124191578554,
            0.820001985973903,
            0.886415527004401,
            0.938274552002733,
            0.974728555971309,
            0.995187219997021
        };

        double weight_rule24 [ 24 ] = {
            0.012341229799991,
            0.028531388628934,
            0.044277438817419,
            0.059298584915436,
            0.073346481411080,
            0.086190161531953,
            0.097618652104114,
            0.107444270115966,
            0.115505668053726,
            0.121670472927803,
            0.125837456346828,
            0.127938195346752,
            0.127938195346752,
            0.125837456346828,
            0.121670472927803,
            0.115505668053726,
            0.107444270115966,
            0.097618652104114,
            0.086190161531953,
            0.073346481411080,
            0.059298584915436,
            0.044277438817419,
            0.028531388628934,
            0.012341229799991
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 24; i++ ) {
            for ( j = 0; j < 24; j++ ) {
                for ( k = 0; k < 24; k++ ) {
                    coord  = new FloatArray(3);
                    coord->at(1) = coord_rule24 [ i ];
                    coord->at(2) = coord_rule24 [ j ];
                    coord->at(3) = coord_rule24 [ k ];
                    weight = weight_rule24 [ i ] * weight_rule24 [ j ] * weight_rule24 [ k ];
                    ( * arry ) [ 576 * i + 24 * j + k ] = new GaussPoint(this, 576 *i + 24 *j + k + 1, coord, weight, mode);
                }
            }
        }
    }

    break;

    case 32768:
    {
        double coord_rule32 [ 32 ] = {
            -0.997263861849482,
            -0.985611511545268,
            -0.964762255587506,
            -0.934906075937740,
            -0.896321155766052,
            -0.849367613732570,
            -0.794483795967942,
            -0.732182118740290,
            -0.663044266930215,
            -0.587715757240762,
            -0.506899908932229,
            -0.421351276130635,
            -0.331868602282128,
            -0.239287362252137,
            -0.144471961582796,
            -0.048307665687738,
            0.048307665687738,
            0.144471961582796,
            0.239287362252137,
            0.331868602282128,
            0.421351276130635,
            0.506899908932229,
            0.587715757240762,
            0.663044266930215,
            0.732182118740290,
            0.794483795967942,
            0.849367613732570,
            0.896321155766052,
            0.934906075937740,
            0.964762255587506,
            0.985611511545268,
            0.997263861849482
        };

        double weight_rule32 [ 32 ] = {
            0.007018610009469,
            0.016274394730905,
            0.025392065309263,
            0.034273862913022,
            0.042835898022227,
            0.050998059262376,
            0.058684093478536,
            0.065822222776362,
            0.072345794108848,
            0.078193895787070,
            0.083311924226947,
            0.087652093004404,
            0.091173878695764,
            0.093844399080805,
            0.095638720079275,
            0.096540088514728,
            0.096540088514728,
            0.095638720079275,
            0.093844399080805,
            0.091173878695764,
            0.087652093004404,
            0.083311924226947,
            0.078193895787070,
            0.072345794108848,
            0.065822222776362,
            0.058684093478536,
            0.050998059262376,
            0.042835898022227,
            0.034273862913022,
            0.025392065309263,
            0.016274394730905,
            0.007018610009469
        };

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 32; i++ ) {
            for ( j = 0; j < 32; j++ ) {
                for ( k = 0; k < 32; k++ ) {
                    coord  = new FloatArray(3);
                    coord->at(1) = coord_rule32 [ i ];
                    coord->at(2) = coord_rule32 [ j ];
                    coord->at(3) = coord_rule32 [ k ];
                    weight = weight_rule32 [ i ] * weight_rule32 [ j ] * weight_rule32 [ k ];
                    ( * arry ) [ 1024 * i + 32 * j + k ] = new GaussPoint(this, 1024 *i + 32 *j + k + 1, coord, weight, mode);
                }
            }
        }
    }

    break;

    default:
        OOFEM_ERROR2("SetUpPointsOnCube: unsupported number of IPs (%d)", nPoints);
    }

    return nPoints;
}


int
GaussIntegrationRule :: SetUpPointsOnTetrahedra(int nPoints,
                                                MaterialMode mode, GaussPoint ***arry)
{
    double a, b, c, w;
    FloatArray *coord1;

    switch ( nPoints ) {
    case 1:

        * arry  = new GaussPoint * [ nPoints ];
        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.25;
        coord1->at(2)    = 0.25;
        coord1->at(3)    = 0.25;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, 1./6., mode);
        break;

    case 4:
        // quadratic formulae

        * arry  = new GaussPoint * [ nPoints ];

        a = (5. + 3.*sqrt(5.))/20.;
        b = (5. - sqrt(5.))/20.;
        w = 1. / 24.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 1 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 2 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 3 ]        = new GaussPoint(this, 1, coord1, w, mode);

        break;

    case 5:
        // cubic formulae

        * arry  = new GaussPoint * [ nPoints ];

        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.25;
        coord1->at(2)    = 0.25;
        coord1->at(3)    = 0.25;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, -2. / 15., mode);

        a = 1. / 2.;
        b = 1. / 6.;
        w = 3. / 40.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 1 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 2 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 3 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 4 ]        = new GaussPoint(this, 1, coord1, w, mode);

        break;

    case 11:
        // exact x^4
        * arry  = new GaussPoint * [ nPoints ];

        a = 1. / 4.;
        w = -74. / 5625.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, w, mode);

        a = 5./70.;
        b = 11./14.;
        w = 343. / 45000.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( * arry ) [ 1 ]        = new GaussPoint(this, 2, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 2 ]        = new GaussPoint(this, 3, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 3 ]        = new GaussPoint(this, 4, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( * arry ) [ 4 ]        = new GaussPoint(this, 5, coord1, w, mode);

        a = (1. + sqrt(5./14.)) / 4.;
        b = (1. - sqrt(5./14.)) / 4.;
        w = 28. / 1125.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 5 ]        = new GaussPoint(this, 6, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 6 ]        = new GaussPoint(this, 7, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 7 ]        = new GaussPoint(this, 8, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( * arry ) [ 8 ]        = new GaussPoint(this, 9, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 9 ]        = new GaussPoint(this, 10, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 10 ]        = new GaussPoint(this, 11, coord1, w, mode);

        break;

    case 15:
        // Derived by Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
        // exact for x^5
        * arry  = new GaussPoint * [ nPoints ];

        a = 1. / 4.;
        w = 0.302836780970891856e-1;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( * arry ) [ 0 ]        = new GaussPoint(this, 1, coord1, w, mode);

        a = 0.0;
        b = 1. / 3.;
        w = 0.602678571428571597e-2;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 1 ]        = new GaussPoint(this, 2, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 2 ]        = new GaussPoint(this, 3, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 3 ]        = new GaussPoint(this, 4, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 4 ]        = new GaussPoint(this, 5, coord1, w, mode);

        ////////////// 
        a = 8. / 11.;
        b = 1. / 11.;
        w = 0.116452490860289742e-1;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 5 ]        = new GaussPoint(this, 6, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 6 ]        = new GaussPoint(this, 7, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 7 ]        = new GaussPoint(this, 8, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 8 ]        = new GaussPoint(this, 9, coord1, w, mode);

        ////////////// 
        a = 0.665501535736642813e-1;
        b = 0.433449846426335728;
        w = 0.109491415613864534e-1;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 9 ]        = new GaussPoint(this, 10, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 10 ]        = new GaussPoint(this, 11, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 11 ]        = new GaussPoint(this, 12, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( * arry ) [ 12 ]        = new GaussPoint(this, 13, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 13 ]        = new GaussPoint(this, 14, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 14 ]        = new GaussPoint(this, 15, coord1, w, mode);

        break;

    case 24:
        // Derived by Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
        // See also: Monomial cubature rules since “Stroud”: a compilation
        // Exact for x^6
        * arry  = new GaussPoint * [ nPoints ];
        a = 0.356191386222544953;
        b = 0.214602871259151684;
        w = 0.665379170969464506e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 0 ] = new GaussPoint(this, 1, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( * arry ) [ 1 ] = new GaussPoint(this, 2, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( * arry ) [ 2 ] = new GaussPoint(this, 3, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 3 ] = new GaussPoint(this, 4, coord1, w, mode);

        a = 0.877978124396165982;
        b = 0.406739585346113397e-1;
        w = 0.167953517588677620e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 4 ] = new GaussPoint(this, 5, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( * arry ) [ 5 ] = new GaussPoint(this, 6, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( * arry ) [ 6 ] = new GaussPoint(this, 7, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 7 ] = new GaussPoint(this, 8, coord1, w, mode);

        a = 0.329863295731730594e-1;
        b = 0.322337890142275646;
        w = 0.922619692394239843e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 8 ] = new GaussPoint(this, 9, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( * arry ) [ 9 ] = new GaussPoint(this, 10, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( * arry ) [ 10 ] = new GaussPoint(this, 11, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 11 ] = new GaussPoint(this, 12, coord1, w, mode);

        a = 0.603005664791649076;
        b = 0.269672331458315867;
        c = 0.636610018750175299e-1;
        w = 0.803571428571428248e-2;
        // 12 permutations follows here, a, b, c
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( * arry ) [ 12 ] = new GaussPoint(this, 13, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( * arry ) [ 13 ] = new GaussPoint(this, 14, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( * arry ) [ 14 ] = new GaussPoint(this, 15, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( * arry ) [ 15 ] = new GaussPoint(this, 16, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( * arry ) [ 16 ] = new GaussPoint(this, 17, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( * arry ) [ 17 ] = new GaussPoint(this, 18, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( * arry ) [ 18 ] = new GaussPoint(this, 19, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( * arry ) [ 19 ] = new GaussPoint(this, 20, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( * arry ) [ 20 ] = new GaussPoint(this, 21, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( * arry ) [ 21 ] = new GaussPoint(this, 22, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( * arry ) [ 22 ] = new GaussPoint(this, 23, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( * arry ) [ 23 ] = new GaussPoint(this, 24, coord1, w, mode);

        break;

    case 45:        
        // Derived by Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
        // See also: Monomial cubature rules since “Stroud”: a compilation
        // Exact for x^8
        * arry  = new GaussPoint * [ nPoints ];

        a = 1. / 4.;
        w = -0.393270066412926145e-1;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = a;
        coord1->at(3) = a;
        ( * arry ) [ 0 ] = new GaussPoint(this, 1, coord1, w, mode);

        a = 0.617587190300082967;
        b = 0.127470936566639015;
        w = 0.408131605934270525e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 1 ] = new GaussPoint(this, 2, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( * arry ) [ 2 ] = new GaussPoint(this, 3, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( * arry ) [ 3 ] = new GaussPoint(this, 4, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 4 ] = new GaussPoint(this, 5, coord1, w, mode);

        a = 0.903763508822103123;
        b = 0.320788303926322960e-1;
        w = 0.658086773304341943e-3;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 5 ] = new GaussPoint(this, 6, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( * arry ) [ 6 ] = new GaussPoint(this, 7, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( * arry ) [ 7 ] = new GaussPoint(this, 8, coord1, w, mode);
        
        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( * arry ) [ 8 ] = new GaussPoint(this, 9, coord1, w, mode);

        a = 0.450222904356718978;
        b = 0.497770956432810185e-1;
        w = 0.438425882512284693e-2;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 9 ]        = new GaussPoint(this, 10, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 10 ]        = new GaussPoint(this, 11, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 11 ]        = new GaussPoint(this, 12, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( * arry ) [ 12 ]        = new GaussPoint(this, 13, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 13 ]        = new GaussPoint(this, 14, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 14 ]        = new GaussPoint(this, 15, coord1, w, mode);

        a = 0.316269552601450060;
        b = 0.183730447398549945;
        w = 0.138300638425098166e-1;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 15 ] = new GaussPoint(this, 16, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 16 ] = new GaussPoint(this, 17, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( * arry ) [ 17 ] = new GaussPoint(this, 18, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( * arry ) [ 18 ] = new GaussPoint(this, 19, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( * arry ) [ 19 ] = new GaussPoint(this, 20, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( * arry ) [ 20 ] = new GaussPoint(this, 21, coord1, w, mode);

        a = 0.513280033360881072;
        b = 0.229177878448171174e-1;
        c = 0.231901089397150906;
        w = 0.424043742468372453e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( * arry ) [ 21 ] = new GaussPoint(this, 22, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( * arry ) [ 22 ] = new GaussPoint(this, 23, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( * arry ) [ 23 ] = new GaussPoint(this, 24, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( * arry ) [ 24 ] = new GaussPoint(this, 25, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( * arry ) [ 25 ] = new GaussPoint(this, 26, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( * arry ) [ 26 ] = new GaussPoint(this, 27, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( * arry ) [ 27 ] = new GaussPoint(this, 28, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( * arry ) [ 28 ] = new GaussPoint(this, 29, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( * arry ) [ 29 ] = new GaussPoint(this, 30, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( * arry ) [ 30 ] = new GaussPoint(this, 31, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( * arry ) [ 31 ] = new GaussPoint(this, 32, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( * arry ) [ 32 ] = new GaussPoint(this, 33, coord1, w, mode);

        a = 0.193746475248804382;
        b = 0.730313427807538396;
        c = 0.379700484718286102e-1;
        w = 0.223873973961420164e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( * arry ) [ 33 ] = new GaussPoint(this, 34, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( * arry ) [ 34 ] = new GaussPoint(this, 35, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( * arry ) [ 35 ] = new GaussPoint(this, 36, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( * arry ) [ 36 ] = new GaussPoint(this, 37, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( * arry ) [ 37 ] = new GaussPoint(this, 38, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( * arry ) [ 38 ] = new GaussPoint(this, 39, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( * arry ) [ 39 ] = new GaussPoint(this, 40, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( * arry ) [ 40 ] = new GaussPoint(this, 41, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( * arry ) [ 41 ] = new GaussPoint(this, 42, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( * arry ) [ 42 ] = new GaussPoint(this, 43, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( * arry ) [ 43 ] = new GaussPoint(this, 44, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( * arry ) [ 44 ] = new GaussPoint(this, 45, coord1, w, mode);

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
        if ( requiredNIP > 64 ) {
            return -1;
        }

        if ( requiredNIP <= 1) {
        	return 1;
        }

        if ( requiredNIP <= 2) {
        	return 2;
        }

        if ( requiredNIP <= 3) {
        	return 3;
        }

        if ( requiredNIP <= 4) {
        	return 4;
        }

        if ( requiredNIP <= 8) {
        	return 8;
        }

        if ( requiredNIP <= 16) {
        	return 16;
        }

        if ( requiredNIP <= 24) {
        	return 24;
        }

        if ( requiredNIP <= 32) {
        	return 32;
        }

        if ( requiredNIP <= 64) {
        	return 64;
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

        if ( approxOrder <= 4 ) {
            return 11; // unsure
        }

        if ( approxOrder <= 5 ) {
            return 15;
        }

        if ( approxOrder <= 6 ) {
            return 24;
        }

        if ( approxOrder <= 8 ) {
            return 45;
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
//   store coordinates and weights)
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
        OOFEM_ERROR2("GaussIntegrationRule :: SetUpPointsOn2DEmbeddedLine - order %d of integration is not suported\n",nPoints);
    }

    // Initialize the material status at each Gauss point
    //for (int i=0; i<numberOfIntegrationPoints; i++)
    //  mat->giveStatus(gaussPointArray[i]);

    return nPoints;
}
} // end namespace oofem
