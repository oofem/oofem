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
    this->gaussPoints.resize( nPoints );

    for ( int i = 1; i <= nPoints; i++ ) {
        this->gaussPoints [ i - 1 ] = new GaussPoint(this, i, {coords_xi.at(i)}, weights.at ( i ), mode);
    }

    this->intdomain = _Line;
    return this->giveNumberOfIntegrationPoints();
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
    this->gaussPoints.resize( nPoints_xi1 * nPoints_xi2 );
    int count = 0;
    for ( int i = 1; i <= nPoints_xi1; i++ ) {
        for ( int j = 1; j <= nPoints_xi2; j++ ) {
            count++;
            this->gaussPoints [ count - 1 ] = new GaussPoint(this, count, {coords_xi1.at(i), coords_xi2.at(j)}, 
                                                             weights1.at ( i ) *weights2.at ( j ), mode);
        }
    }

    this->intdomain = _Square;
    return this->giveNumberOfIntegrationPoints();
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
    this->gaussPoints.resize( nPoints_xi1 * nPoints_xi2 * nPoints_xi3 );
    int count = 0;
    for ( int i = 1; i <= nPoints_xi1; i++ ) {
        for ( int j = 1; j <= nPoints_xi2; j++ ) {
            for ( int k = 1; k <= nPoints_xi3; k++ ) {
                count++;
                this->gaussPoints [ count - 1 ] = new GaussPoint(this, count, {coords_xi1.at(i), coords_xi2.at(j), coords_xi3.at(k)},
                                                                 weights1.at ( i ) *weights2.at ( j ) *weights3.at ( k ), mode);
            }
        }
    }

    this->intdomain = _Cube;
    return this->giveNumberOfIntegrationPoints();
}


int
LobattoIntegrationRule :: SetUpPointsOnTriangle(int nPoints, MaterialMode mode)
{
    OOFEM_ERROR("unsupported number of IPs (%d)", nPoints);
    return 0;
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
    switch ( nPoints ) {
    case 1:
        coords_xi = FloatArray{0.0};
        weights = FloatArray{2.0};
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

    default:
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
        if ( nPoints > 6 ) {
            OOFEM_LOG_WARNING("Unsupported number of IPs (%d) for LobattoIR, using 6 ips instead.", nPoints);
        }
        break;
    }
}
} // end namespace oofem
