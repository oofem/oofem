/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/lobattoir.C,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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


#include "lobattoir.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

namespace oofem {
// initialize class member

LobattoIntegrationRule :: LobattoIntegrationRule(int n, Element *e,
                                                 int startIndx, int endIndx, bool dynamic) :
    IntegrationRule(n, e, startIndx, endIndx, dynamic) { }


LobattoIntegrationRule :: ~LobattoIntegrationRule()
{ }



int
LobattoIntegrationRule :: SetUpPointsOnLine(int nPoints, Element *elem,
                                            MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Lobatto Integration Points
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

        c->at(1) = -1.0;
        c->at(2) =  1.0;

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

        c->at(1) = -1.0;
        c->at(2) =  0.0;
        c->at(3) =  1.0;

        w->at(1) =  0.333333333333333;
        w->at(2) =  1.333333333333333;
        w->at(3) =  0.333333333333333;

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

        c->at(1) = -1.0;
        c->at(2) = -0.447213595499958;
        c->at(3) =  0.447213595499958;
        c->at(4) =  1.0;

        w->at(1) =  0.166666666666667;
        w->at(2) =  0.833333333333333;
        w->at(3) =  0.833333333333333;
        w->at(4) =  0.166666666666667;

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

    case 5:

        c = new FloatArray(5);
        w = new FloatArray(5);

        c->at(1) = -1.0;
        c->at(2) = -0.654653670707977;
        c->at(3) =  0.0;
        c->at(4) =  0.654653670707977;
        c->at(5) =  1.0;

        w->at(1) =  0.1;
        w->at(2) =  0.544444444444444;
        w->at(3) =  0.711111111111111;
        w->at(4) =  0.544444444444444;
        w->at(5) =  0.1;


        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 5; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = c->at(i + 1);
            weight = w->at(i + 1);
            ( * arry ) [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
        }

        delete c;
        delete w;
        break;

    case 6:

        c = new FloatArray(6);
        w = new FloatArray(6);

        c->at(1) = -1.0;
        c->at(2) = -0.765055323929465;
        c->at(3) = -0.285231516480645;
        c->at(4) =  0.285231516480645;
        c->at(5) =  0.765055323929465;
        c->at(6) =  1.0;

        w->at(1) =  0.066666666666667;
        w->at(2) =  0.378474956297847;
        w->at(3) =  0.554858377035486;
        w->at(4) =  0.554858377035486;
        w->at(5) =  0.378474956297847;
        w->at(6) =  0.066666666666667;

        * arry = new GaussPoint * [ nPoints ];

        for ( i = 0; i < 6; i++ ) {
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
LobattoIntegrationRule :: SetUpPointsOnTriagle(int nPoints, Element *elem,
                                               MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    OOFEM_ERROR2("SetUpPointsOnTriangle: unsupported number of IPs (%d)", nPoints);
    return nPoints;
}

int
LobattoIntegrationRule :: SetUpPointsOnSquare(int nPoints, Element *elem,
                                              MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    OOFEM_ERROR2("SetUpPointsOnSquare: unsupported number of IPs (%d)", nPoints);
    return nPoints;
}

int
LobattoIntegrationRule :: SetUpPointsOnCube(int nPoints, Element *elem,
                                            MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    OOFEM_ERROR2("SetUpPointsOnCube: unsupported number of IPs (%d)", nPoints);
    return nPoints;
}


int
LobattoIntegrationRule :: SetUpPointsOnTetrahedra(int nPoints, Element *elem,
                                                  MaterialMode mode, GaussPoint ***arry)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    OOFEM_ERROR2("SetUpPointsOnTetrahedra: unsupported number of IPs (%d)", nPoints);
    return nPoints;
}

int
LobattoIntegrationRule :: getRequiredNumberOfIntegrationPoints(integrationDomain dType,
                                                               int approxOrder)
{
    int requiredNIP;
    if ( approxOrder < 0 ) {
        return 0;
    }

    switch ( dType ) {
    case _Line:
        if ( approxOrder <= 1 ) {
            return 1;
        }

        requiredNIP = ( int ) ceil( ( ( double ) approxOrder + 3.0 ) / 2. );
        if ( requiredNIP > 6 ) {
            return -1;
        }

        return requiredNIP;

    default:
        OOFEM_ERROR("LobattoIntegrationRule::setUpIntegrationPoints - unknown integrationDomain");
    }

    return -1;
}
} // end namespace oofem
