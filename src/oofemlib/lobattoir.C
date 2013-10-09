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

#include "lobattoir.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "mathfem.h"

namespace oofem {

LobattoIntegrationRule :: LobattoIntegrationRule(int n, Element *e,
                                                 int startIndx, int endIndx, bool dynamic) :
    IntegrationRule(n, e, startIndx, endIndx, dynamic) { }

LobattoIntegrationRule :: LobattoIntegrationRule(int n, Element *e) :
    IntegrationRule(n, e) { }

LobattoIntegrationRule :: ~LobattoIntegrationRule()
{ }


int
LobattoIntegrationRule :: SetUpPointsOnLine(int nPoints, MaterialMode mode)
{
    double weight;
    FloatArray *coord, *c, *w;

    switch ( nPoints ) {
    case 1:

        gaussPointArray = new GaussPoint * [ nPoints ];
        coord = new FloatArray(1);
        coord->at(1) = 0.0;
        weight = 2.0;
        gaussPointArray [ 0 ] = new GaussPoint(this, 1, coord, weight, mode);
        break;

    case 2:

        c = new FloatArray(2);
        w = new FloatArray(2);

        c->at(1) = -1.0;
        c->at(2) =  1.0;

        w->at(1) = 1.0;
        w->at(2) = 1.0;

        gaussPointArray = new GaussPoint * [ nPoints ];

        for ( int i = 0; i < 2; i++ ) {
            coord = new FloatArray(1);
            coord->at(1) = c->at(i + 1);
            weight = w->at(i + 1);
            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
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

        gaussPointArray = new GaussPoint * [ nPoints ];

        for ( int i = 0; i < 3; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = c->at(i + 1);
            weight = w->at(i + 1);
            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
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

        gaussPointArray = new GaussPoint * [ nPoints ];

        for ( int i = 0; i < 4; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = c->at(i + 1);
            weight = w->at(i + 1);
            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
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


        gaussPointArray = new GaussPoint * [ nPoints ];

        for ( int i = 0; i < 5; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = c->at(i + 1);
            weight = w->at(i + 1);
            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
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

        gaussPointArray = new GaussPoint * [ nPoints ];

        for ( int i = 0; i < 6; i++ ) {
            coord  = new FloatArray(1);
            coord->at(1) = c->at(i + 1);
            weight = w->at(i + 1);
            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
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
LobattoIntegrationRule :: SetUpPointsOnTriangle(int nPoints, MaterialMode mode)
{
    OOFEM_ERROR2("SetUpPointsOnTriangle: unsupported number of IPs (%d)", nPoints);
    return nPoints;
}

int
LobattoIntegrationRule :: SetUpPointsOnSquare(int nPoints, MaterialMode mode)
{
    OOFEM_ERROR2("SetUpPointsOnSquare: unsupported number of IPs (%d)", nPoints);
    return nPoints;
}

int
LobattoIntegrationRule :: SetUpPointsOnCube(int nPoints, MaterialMode mode)
{
    OOFEM_ERROR2("SetUpPointsOnCube: unsupported number of IPs (%d)", nPoints);
    return nPoints;
}


int
LobattoIntegrationRule :: SetUpPointsOnTetrahedra(int nPoints, MaterialMode mode)
{
    OOFEM_ERROR2("SetUpPointsOnTetrahedra: unsupported number of IPs (%d)", nPoints);
    return nPoints;
}

int
LobattoIntegrationRule :: getRequiredNumberOfIntegrationPoints(integrationDomain dType,
                                                               int approxOrder)
{
    int requiredNIP;

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
