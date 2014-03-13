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
    FloatArray coords_xi, weights;
    this->giveLineCoordsAndWeights(nPoints, coords_xi, weights);
    this->numberOfIntegrationPoints = nPoints;
    this->gaussPointArray  = new GaussPoint * [ nPoints ];

    for ( int i = 1; i <= nPoints; i++ ) {
        FloatArray *coord = new FloatArray(1);
        coord->at(1) = coords_xi.at(i);
        this->gaussPointArray [ i - 1 ] = new GaussPoint(this, i, coord, weights.at ( i ), mode);
    }

    this->intdomain = _Line;
    return numberOfIntegrationPoints;
}


int
LobattoIntegrationRule :: SetUpPointsOnSquare(int nPoints, MaterialMode mode)
//GaussIntegrationRule :: SetUpPointsOnSquare(int nPoints_xi1, int nPoints_xi2, MaterialMode mode)
{
    int nPoints_xi1 = ( int ) floor( sqrt( double ( nPoints ) ) );
    int nPoints_xi2 = nPoints_xi1;
    FloatArray coords_xi1, weights1, coords_xi2, weights2;
    this->giveLineCoordsAndWeights(nPoints_xi1, coords_xi1, weights1);
    this->giveLineCoordsAndWeights(nPoints_xi2, coords_xi2, weights2);
    this->numberOfIntegrationPoints = nPoints_xi1 * nPoints_xi2;
    this->gaussPointArray  = new GaussPoint * [ this->numberOfIntegrationPoints ];
    int count = 0;
    for ( int i = 1; i <= nPoints_xi1; i++ ) {
        for ( int j = 1; j <= nPoints_xi2; j++ ) {
            count++;
            FloatArray *coord = new FloatArray(2);
            coord->at(1) = coords_xi1.at(i);
            coord->at(2) = coords_xi2.at(j);
            this->gaussPointArray [ count - 1 ] = new GaussPoint(this, count, coord, weights1.at ( i ) *weights2.at ( j ), mode);
        }
    }

    this->intdomain = _Square;
    return this->numberOfIntegrationPoints;
}


int
LobattoIntegrationRule :: SetUpPointsOnCube(int nPoints, MaterialMode mode)
//GaussIntegrationRule :: SetUpPointsOnCube(int nPoints_xi1, int nPoints_xi2, int nPoints_xi3, MaterialMode mode)
{
    int nPoints_xi1 = ( int ) floor(cbrt( double ( nPoints ) ) + 0.5);
    int nPoints_xi2 = nPoints_xi1;
    int nPoints_xi3 = nPoints_xi1;
    FloatArray coords_xi1, weights1, coords_xi2, weights2, coords_xi3, weights3;
    this->giveLineCoordsAndWeights(nPoints_xi1, coords_xi1, weights1);
    this->giveLineCoordsAndWeights(nPoints_xi2, coords_xi2, weights2);
    this->giveLineCoordsAndWeights(nPoints_xi3, coords_xi3, weights3);
    this->numberOfIntegrationPoints = nPoints_xi1 * nPoints_xi2 * nPoints_xi3;
    this->gaussPointArray  = new GaussPoint * [ this->numberOfIntegrationPoints ];
    int count = 0;
    for ( int i = 1; i <= nPoints_xi1; i++ ) {
        for ( int j = 1; j <= nPoints_xi2; j++ ) {
            for ( int k = 1; k <= nPoints_xi3; k++ ) {
                count++;
                FloatArray *coord = new FloatArray(3);
                coord->at(1) = coords_xi1.at(i);
                coord->at(2) = coords_xi2.at(j);
                coord->at(3) = coords_xi3.at(k);
                this->gaussPointArray [ count - 1 ] = new GaussPoint(this, count, coord, weights1.at ( i ) *weights2.at ( j ) *weights3.at ( k ), mode);
            }
        }
    }

    this->intdomain = _Cube;
    return this->numberOfIntegrationPoints;
}


//int
//{
//    double weight;
//    FloatArray *coord, *c, *w;
//
//    switch ( nPoints ) {
//    case 1:
//
//        gaussPointArray = new GaussPoint * [ nPoints ];
//        coord = new FloatArray(1);
//        coord->at(1) = 0.0;
//        weight = 2.0;
//        gaussPointArray [ 0 ] = new GaussPoint(this, 1, coord, weight, mode);
//        break;
//
//    case 2:
//
//        c = new FloatArray(2);
//        w = new FloatArray(2);
//
//        c->at(1) = -1.0;
//        c->at(2) =  1.0;
//
//        w->at(1) = 1.0;
//        w->at(2) = 1.0;
//
//        gaussPointArray = new GaussPoint * [ nPoints ];
//
//        for ( int i = 0; i < 2; i++ ) {
//            coord = new FloatArray(1);
//            coord->at(1) = c->at(i + 1);
//            weight = w->at(i + 1);
//            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
//        }
//
//        delete c;
//        delete w;
//        break;
//
//    case 3:
//
//        c = new FloatArray(3);
//        w = new FloatArray(3);
//
//        c->at(1) = -1.0;
//        c->at(2) =  0.0;
//        c->at(3) =  1.0;
//
//        w->at(1) =  0.333333333333333;
//        w->at(2) =  1.333333333333333;
//        w->at(3) =  0.333333333333333;
//
//        gaussPointArray = new GaussPoint * [ nPoints ];
//
//        for ( int i = 0; i < 3; i++ ) {
//            coord  = new FloatArray(1);
//            coord->at(1) = c->at(i + 1);
//            weight = w->at(i + 1);
//            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
//        }
//
//        delete c;
//        delete w;
//        break;
//
//    case 4:
//
//        c = new FloatArray(4);
//        w = new FloatArray(4);
//
//        c->at(1) = -1.0;
//        c->at(2) = -0.447213595499958;
//        c->at(3) =  0.447213595499958;
//        c->at(4) =  1.0;
//
//        w->at(1) =  0.166666666666667;
//        w->at(2) =  0.833333333333333;
//        w->at(3) =  0.833333333333333;
//        w->at(4) =  0.166666666666667;
//
//        gaussPointArray = new GaussPoint * [ nPoints ];
//
//        for ( int i = 0; i < 4; i++ ) {
//            coord  = new FloatArray(1);
//            coord->at(1) = c->at(i + 1);
//            weight = w->at(i + 1);
//            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
//        }
//
//        delete c;
//        delete w;
//        break;
//
//    case 5:
//
//        c = new FloatArray(5);
//        w = new FloatArray(5);
//
//        c->at(1) = -1.0;
//        c->at(2) = -0.654653670707977;
//        c->at(3) =  0.0;
//        c->at(4) =  0.654653670707977;
//        c->at(5) =  1.0;
//
//        w->at(1) =  0.1;
//        w->at(2) =  0.544444444444444;
//        w->at(3) =  0.711111111111111;
//        w->at(4) =  0.544444444444444;
//        w->at(5) =  0.1;
//
//
//        gaussPointArray = new GaussPoint * [ nPoints ];
//
//        for ( int i = 0; i < 5; i++ ) {
//            coord  = new FloatArray(1);
//            coord->at(1) = c->at(i + 1);
//            weight = w->at(i + 1);
//            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
//        }
//
//        delete c;
//        delete w;
//        break;
//
//    case 6:
//
//        c = new FloatArray(6);
//        w = new FloatArray(6);
//
//        c->at(1) = -1.0;
//        c->at(2) = -0.765055323929465;
//        c->at(3) = -0.285231516480645;
//        c->at(4) =  0.285231516480645;
//        c->at(5) =  0.765055323929465;
//        c->at(6) =  1.0;
//
//        w->at(1) =  0.066666666666667;
//        w->at(2) =  0.378474956297847;
//        w->at(3) =  0.554858377035486;
//        w->at(4) =  0.554858377035486;
//        w->at(5) =  0.378474956297847;
//        w->at(6) =  0.066666666666667;
//
//        gaussPointArray = new GaussPoint * [ nPoints ];
//
//        for ( int i = 0; i < 6; i++ ) {
//            coord  = new FloatArray(1);
//            coord->at(1) = c->at(i + 1);
//            weight = w->at(i + 1);
//            gaussPointArray [ i ] = new GaussPoint(this, i + 1, coord, weight, mode);
//        }
//
//        delete c;
//        delete w;
//        break;
//
//    default:
//        OOFEM_ERROR("unsupported number of IPs (%d)", nPoints);
//    }
//
//    return nPoints;
//}

int
LobattoIntegrationRule :: SetUpPointsOnTriangle(int nPoints, MaterialMode mode)
{
    OOFEM_ERROR("unsupported number of IPs (%d)", nPoints);
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
        OOFEM_ERROR("unknown integrationDomain");
    }

    return -1;
}


void
LobattoIntegrationRule :: giveLineCoordsAndWeights(int nPoints, FloatArray &coords_xi, FloatArray &weights)
// Create arrays of coordinates and weights for Lobatto Integration Points of a line with 'nPoints' integrationpoints
{
    coords_xi.resize(nPoints);
    weights.resize(nPoints);

    switch ( nPoints ) {
    case 1:
        coords_xi = {0.0};
        weights = {2.0};
        break;

    case 2:
        coords_xi = {-1.0, 1.0};
        weights = {1.0, 1.0};
        break;

    case 3:
        coords_xi = {-1.0, 0.0, 1.0};
        weights = {0.333333333333333, 1.333333333333333, 0.333333333333333};
        break;

    case 4:
        coords_xi = {-1.0,
                     -0.447213595499958,
                      0.447213595499958,
                      1.0};
        weights = { 0.166666666666667,
                    0.833333333333333,
                    0.833333333333333,
                    0.166666666666667};
        break;

    case 5:
        coords_xi = {-1.0,
                     -0.654653670707977,
                      0.0,
                      0.654653670707977,
                      1.0};
        weights = { 0.1,
                    0.544444444444444,
                    0.711111111111111,
                    0.544444444444444,
                    0.1};
        break;

    case 6:
        coords_xi = {-1.0,
                     -0.765055323929465,
                     -0.285231516480645,
                      0.285231516480645,
                      0.765055323929465,
                      1.0};
        weights = { 0.066666666666667,
                    0.378474956297847,
                    0.554858377035486,
                    0.554858377035486,
                    0.378474956297847,
                    0.066666666666667};
        break;

    default:
        OOFEM_SIMPLE_ERROR("unsupported number of IPs (%d)", nPoints);
    }
}
} // end namespace oofem
