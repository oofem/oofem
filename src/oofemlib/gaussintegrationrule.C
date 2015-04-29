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

#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "element.h"
#include "floatarray.h"
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
GaussIntegrationRule :: SetUpPointsOnLine(int nPoints, MaterialMode mode)
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
GaussIntegrationRule :: SetUpPointsOn2DEmbeddedLine(int nPoints, MaterialMode mode, const FloatArray &coord0, const FloatArray &coord1)
{
    FloatArray coords_xi, weights;
    this->giveLineCoordsAndWeights(nPoints, coords_xi, weights);
    this->gaussPoints.resize( nPoints );

    for ( int i = 1; i <= nPoints; i++ ) {
        double x = ( coords_xi.at(i) + 1.0 ) * 0.5;
        FloatArray subpatchCoord = {x};

        this->gaussPoints [ i - 1 ] = new GaussPoint(this, i, weights.at ( i ), mode);

        this->gaussPoints [ i - 1 ]->setSubPatchCoordinates(subpatchCoord);

        FloatArray globalCoord = { 
            ( 1. - x ) * coord0.at(1) + x * coord1.at(1),
            ( 1. - x ) * coord0.at(2) + x * coord1.at(2) };
        this->gaussPoints [ i - 1 ]->setGlobalCoordinates(globalCoord);

        FloatArray naturalCoord;
        this->giveElement()->computeLocalCoordinates(naturalCoord, globalCoord);
        this->gaussPoints [ i - 1 ]->setNaturalCoordinates(naturalCoord);
    }

    this->intdomain = _Embedded2dLine;
    return this->giveNumberOfIntegrationPoints();
}


int
GaussIntegrationRule :: SetUpPointsOnSquare(int nPoints, MaterialMode mode)
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
            this->gaussPoints [ count ] = new GaussPoint(this, count + 1, {coords_xi1.at(i), coords_xi2.at(j)},
                                                         weights1.at ( i ) *weights2.at ( j ), mode);
            count++;
        }
    }

    this->intdomain = _Square;
    return this->giveNumberOfIntegrationPoints();
}


int
GaussIntegrationRule :: SetUpPointsOnCube(int nPoints, MaterialMode mode)
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
                this->gaussPoints [ count ] = new GaussPoint(this, count + 1, {coords_xi1.at(i), coords_xi2.at(j), coords_xi3.at(k)},
                                                             weights1.at ( i ) *weights2.at ( j ) *weights3.at ( k ), mode);
                count++;
            }
        }
    }

    this->intdomain = _Cube;
    return this->giveNumberOfIntegrationPoints();
}


int GaussIntegrationRule :: SetUpPointsOnCubeLayers(int nPoints1, int nPoints2, int nPointsDepth, MaterialMode mode, const FloatArray &layerThickness)
{
    FloatArray coords_xi1, weights1, coords_xi2, weights2, coords_xi3, weights3;
    this->giveLineCoordsAndWeights(nPoints1, coords_xi1, weights1);
    this->giveLineCoordsAndWeights(nPoints2, coords_xi2, weights2);
    this->giveLineCoordsAndWeights(nPointsDepth, coords_xi3, weights3);
    int pointsPerLayer = nPoints1 * nPoints2 * nPointsDepth;
    this->gaussPoints.resize( pointsPerLayer * layerThickness.giveSize() );

    int count = 0;
    double totalThickness = layerThickness.sum();

    // bottom is scaled so that it goes from -1 to 1.
    double bottom = -1.;
    double scaledThickness;
    for ( int t = 1; t <= layerThickness.giveSize(); t++ ) {
        scaledThickness = layerThickness.at(t) / totalThickness;
        for ( int i = 1; i <= nPoints1; i++ ) {
            for ( int j = 1; j <= nPoints2; j++ ) {
                for ( int k = 1; k <= nPointsDepth; k++ ) {
                    this->gaussPoints [ count ] = new GaussPoint(this, count + 1, 
                                        {coords_xi1.at(i), coords_xi2.at(j), ( coords_xi3.at(k) + 1. ) * scaledThickness + bottom},
                                        weights1.at ( i ) *weights2.at ( j ) * ( weights3.at ( k ) *scaledThickness ), mode);
                    count++;
                }
            }
        }
        bottom += 2.0 * scaledThickness;
    }

    this->intdomain = _Cube;
    return this->giveNumberOfIntegrationPoints();
}


int
GaussIntegrationRule :: SetUpPointsOnTriangle(int nPoints, MaterialMode mode)
{
    FloatArray coords_xi1, coords_xi2, weights;
    this->giveTriCoordsAndWeights(nPoints, coords_xi1, coords_xi2, weights);
    this->gaussPoints.resize( nPoints );

    for ( int i = 1; i <= nPoints; i++ ) {
        this->gaussPoints [ i - 1 ] = new GaussPoint(this, i, {coords_xi1.at(i), coords_xi2.at(i)}, weights.at ( i ), mode);
    }

    this->intdomain = _Triangle;
    return this->giveNumberOfIntegrationPoints();
}


int
GaussIntegrationRule :: SetUpPointsOnTetrahedra(int nPoints, MaterialMode mode)
{
    FloatArray coords_xi1, coords_xi2, coords_xi3, weights;
    this->giveTetCoordsAndWeights(nPoints, coords_xi1, coords_xi2, coords_xi3, weights);
    this->gaussPoints.resize( nPoints );

    for ( int i = 1; i <= nPoints; i++ ) {
        this->gaussPoints [ i - 1 ] = new GaussPoint(this, i, 
                                {coords_xi1.at(i), coords_xi2.at(i), coords_xi3.at(i)}, weights.at ( i ), mode);
    }

    this->intdomain = _Tetrahedra;
    return this->giveNumberOfIntegrationPoints();
}


int
GaussIntegrationRule :: SetUpPointsOnWedge(int nPointsTri, int nPointsDepth, MaterialMode mode)
{
    FloatArray coords_xi1, coords_xi2, coords_xi3, weightsTri, weightsDepth;
    this->giveTriCoordsAndWeights(nPointsTri, coords_xi1, coords_xi2, weightsTri);
    this->giveLineCoordsAndWeights(nPointsDepth, coords_xi3, weightsDepth);
    this->gaussPoints.resize( nPointsTri * nPointsDepth );

    int count = 0;
    for ( int i = 1; i <= nPointsTri; i++ ) {
        for ( int j = 1; j <= nPointsDepth; j++ ) {
            this->gaussPoints [ count ] = new GaussPoint(this, count + 1, {coords_xi1.at(i), coords_xi2.at(i), coords_xi3.at(j)},
                                                         weightsTri.at ( i ) *weightsDepth.at ( j ), mode);
            count++;
        }
    }

    this->intdomain = _Wedge;
    return this->giveNumberOfIntegrationPoints();
}

int
GaussIntegrationRule :: SetUpPointsOnWedgeLayers(int nPointsTri, int nPointsDepth, MaterialMode mode, const FloatArray &layerThickness)
{
    FloatArray coords_xi1, coords_xi2, coords_xi3, weightsTri, weightsDepth;
    this->giveTriCoordsAndWeights(nPointsTri, coords_xi1, coords_xi2, weightsTri);
    this->giveLineCoordsAndWeights(nPointsDepth, coords_xi3, weightsDepth);
    int pointsPerLayer = nPointsTri * nPointsDepth;
    this->gaussPoints.resize( pointsPerLayer * layerThickness.giveSize() );

    int count = 0;
    double totalThickness = layerThickness.sum();

    // bottom is scaled so that it goes from -1 to 1.
    double bottom = -1.;
    double scaledThickness;
    for ( int k = 1; k <= layerThickness.giveSize(); k++ ) {
        scaledThickness = layerThickness.at(k) / totalThickness;
        for ( int i = 1; i <= nPointsTri; i++ ) {
            for ( int j = 1; j <= nPointsDepth; j++ ) {
                this->gaussPoints [ count ] = new GaussPoint(this, count + 1, 
                                {coords_xi1.at(i), coords_xi2.at(i), ( coords_xi3.at(j) + 1. ) * scaledThickness + bottom},
                                weightsTri.at ( i ) * ( weightsDepth.at ( j ) *scaledThickness ), mode);
                count++;
            }
        }
        bottom += 2.0 * scaledThickness;
    }

    this->intdomain = _Wedge;
    return this->giveNumberOfIntegrationPoints();
}

int
GaussIntegrationRule :: getRequiredNumberOfIntegrationPoints(integrationDomain dType,
                                                             int approxOrder)
{
    // Returns the required number of gp's to exactly integrate a polynomial of order 'approxOrder'
    int requiredNIP;
    if ( approxOrder < 0 ) {
        return 0;
    }

    switch ( dType ) {
    case _Line:
        requiredNIP = (int)ceil(( approxOrder + 1 ) / 2);

        if ( requiredNIP <= 8 ) {
            return requiredNIP;
        }

        if ( requiredNIP <= 16 ) {
            return 16;
        }

        if ( requiredNIP <= 24 ) {
            return 24;
        }

        if ( requiredNIP <= 32 ) {
            return 32;
        }

        if ( requiredNIP <= 64 ) {
            return 64;
        }

        return -1;

    case _Triangle:
        if ( approxOrder <= 1 ) {
            return 1;
        }

        if ( approxOrder <= 2 ) {
            return 3;
        }

        if ( approxOrder <= 3 ) {
            return 4;
        }

        if ( approxOrder <= 4 ) {
            return 6;
        }

        if ( approxOrder <= 5 ) {
            return 7;
        }

        if ( approxOrder <= 6 ) {
            return 12;
        }

        if ( approxOrder <= 7 ) {
            return 13;
        }

        if ( approxOrder <= 8 ) {
            return 16;
        }

        if ( approxOrder <= 9 ) {
            return 19;
        }

        if ( approxOrder <= 10 ) {
            return 25;
        }

        return -1;

    case _Square:
        requiredNIP = max( (int)ceil( (double)( approxOrder + 1 ) / (double)2), 2 );
        requiredNIP *= requiredNIP;
        if ( requiredNIP > 64 * 64 ) {
            return -1;
        }

        return requiredNIP;

    case _Cube:
        requiredNIP = max( (int)ceil( (double)( approxOrder + 1 ) / (double)2), 2 );
        requiredNIP *= requiredNIP * requiredNIP;
        if ( requiredNIP > 64 * 64 * 64 ) { // = 262 144 gp's maybe overkill :)
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

        break;
    
    ///@todo Assuming same approximation order for triangle as line. Not totally sure about these /JB
    case _Wedge:
        if ( approxOrder <= 1 ) { // (lin tri) x (lin line) = 1 x 1 = 1
            return 1;        
        }

        if ( approxOrder <= 2 ) { // (quad tri) x (quad line) = 3 x 2 = 6
            return 6;
        }

        if ( approxOrder <= 3 ) { // (cubic tri) x (cubic line) = 4 x 2 = 8
            return 8;
        }                
        
        if ( approxOrder <= 4 ) { // (quartic tri) x (quartic line) = 6 x 3 = 18
            return 18;
        }        
    default:
        OOFEM_ERROR("unknown integrationDomain");
    }

    return -1;
}


void
GaussIntegrationRule :: giveTetCoordsAndWeights(int nPoints, FloatArray &coords_xi1, FloatArray &coords_xi2, FloatArray &coords_xi3, FloatArray &weights)
{
    double a, b, c, w;
    coords_xi1.resize(nPoints);
    coords_xi2.resize(nPoints);
    coords_xi3.resize(nPoints);
    weights.resize(nPoints);
    switch ( nPoints ) {
    case 1:
        a = 1. / 4.;
        w = 1. / 6.;
        coords_xi1(0) = a;
        coords_xi2(0) = a;
        coords_xi3(0) = a;
        weights(0) = w;
        break;

    case 4:     // quadratic formulae
        a = ( 5. + 3. * sqrt(5.) ) / 20.;
        b = ( 5. - sqrt(5.) ) / 20.;
        w = 1. / 24.;
        coords_xi1(0) = a;
        coords_xi2(0) = b;
        coords_xi3(0) = b;
        weights(0) = w;
        coords_xi1(1) = b;
        coords_xi2(1) = a;
        coords_xi3(1) = b;
        weights(1) = w;
        coords_xi1(2) = b;
        coords_xi2(2) = b;
        coords_xi3(2) = a;
        weights(2) = w;
        coords_xi1(3) = b;
        coords_xi2(3) = b;
        coords_xi3(3) = b;
        weights(3) = w;
        break;

    case 5:     // cubic formulae
        a = 1. / 4.;
        w = -2. / 15.;
        coords_xi1(0) = a;
        coords_xi2(0) = a;
        coords_xi3(0) = a;
        weights(0) = w;
        a = 1. / 2.;
        b = 1. / 6.;
        w = 3. / 40.;
        coords_xi1(1) = a;
        coords_xi2(1) = b;
        coords_xi3(1) = b;
        weights(1) = w;
        coords_xi1(2) = b;
        coords_xi2(2) = a;
        coords_xi3(2) = b;
        weights(2) = w;
        coords_xi1(3) = b;
        coords_xi2(3) = b;
        coords_xi3(3) = a;
        weights(3) = w;
        coords_xi1(4) = b;
        coords_xi2(4) = b;
        coords_xi3(4) = b;
        weights(4) = w;
        break;

    case 11:     // exact x^4
        a = 1. / 4.;
        w = -74. / 5625.;
        coords_xi1(0) = a;
        coords_xi2(0) = a;
        coords_xi3(0) = a;
        weights(0) = w;
        a = 5. / 70.;
        b = 11. / 14.;
        w = 343. / 45000.;
        coords_xi1(1) = b;
        coords_xi2(1) = a;
        coords_xi3(1) = a;
        weights(1) = w;
        coords_xi1(2) = a;
        coords_xi2(2) = b;
        coords_xi3(2) = a;
        weights(2) = w;
        coords_xi1(3) = a;
        coords_xi2(3) = a;
        coords_xi3(3) = b;
        weights(3) = w;
        coords_xi1(4) = a;
        coords_xi2(4) = a;
        coords_xi3(4) = a;
        weights(4) = w;
        a = ( 1. + sqrt(5. / 14.) ) / 4.;
        b = ( 1. - sqrt(5. / 14.) ) / 4.;
        w = 28. / 1125.;
        coords_xi1(5) = a;
        coords_xi2(5) = a;
        coords_xi3(5) = b;
        weights(5) = w;
        coords_xi1(6) = a;
        coords_xi2(6) = b;
        coords_xi3(6) = a;
        weights(6) = w;
        coords_xi1(7) = a;
        coords_xi2(7) = b;
        coords_xi3(7) = b;
        weights(7) = w;
        coords_xi1(8) = b;
        coords_xi2(8) = a;
        coords_xi3(8) = a;
        weights(8) = w;
        coords_xi1(9) = b;
        coords_xi2(9) = a;
        coords_xi3(9) = b;
        weights(9) = w;
        coords_xi1(10) = b;
        coords_xi2(10) = b;
        coords_xi3(10) = a;
        weights(10) = w;
        break;

    case 15:     // exact for x^5
        // Derived by Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
        a = 1. / 4.;
        w = 0.302836780970891856e-1;
        coords_xi1(0) = a;
        coords_xi2(0) = a;
        coords_xi3(0) = a;
        weights(0) = w;
        a = 0.0;
        b = 1. / 3.;
        w = 0.602678571428571597e-2;
        coords_xi1(1) = a;
        coords_xi2(1) = b;
        coords_xi3(1) = b;
        weights(1) = w;
        coords_xi1(2) = b;
        coords_xi2(2) = a;
        coords_xi3(2) = b;
        weights(2) = w;
        coords_xi1(3) = b;
        coords_xi2(3) = b;
        coords_xi3(3) = a;
        weights(3) = w;
        coords_xi1(4) = b;
        coords_xi2(4) = b;
        coords_xi3(4) = b;
        weights(4) = w;
        a = 8. / 11.;
        b = 1. / 11.;
        w = 0.116452490860289742e-1;
        coords_xi1(5) = a;
        coords_xi2(5) = b;
        coords_xi3(5) = b;
        weights(5) = w;
        coords_xi1(6) = b;
        coords_xi2(6) = a;
        coords_xi3(6) = b;
        weights(6) = w;
        coords_xi1(7) = b;
        coords_xi2(7) = b;
        coords_xi3(7) = a;
        weights(7) = w;
        coords_xi1(8) = b;
        coords_xi2(8) = b;
        coords_xi3(8) = b;
        weights(8) = w;
        a = 0.665501535736642813e-1;
        b = 0.433449846426335728;
        w = 0.109491415613864534e-1;
        coords_xi1(9) = a;
        coords_xi2(9) = a;
        coords_xi3(9) = b;
        weights(9) = w;
        coords_xi1(10) = a;
        coords_xi2(10) = b;
        coords_xi3(10) = a;
        weights(10) = w;
        coords_xi1(11) = a;
        coords_xi2(11) = b;
        coords_xi3(11) = b;
        weights(11) = w;
        coords_xi1(12) = b;
        coords_xi2(12) = a;
        coords_xi3(12) = a;
        weights(12) = w;
        coords_xi1(13) = b;
        coords_xi2(13) = a;
        coords_xi3(13) = b;
        weights(13) = w;
        coords_xi1(14) = b;
        coords_xi2(14) = b;
        coords_xi3(14) = a;
        weights(14) = w;
        break;

    case 24:     // Exact for x^6
        // Derived by Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
        // See also: Monomial cubature rules since “Stroud”: a compilation
        a = 0.356191386222544953;
        b = 0.214602871259151684;
        w = 0.665379170969464506e-2;
        coords_xi1(0) = a;
        coords_xi2(0) = b;
        coords_xi3(0) = b;
        weights(0) = w;
        coords_xi1(1) = b;
        coords_xi2(1) = a;
        coords_xi3(1) = b;
        weights(1) = w;
        coords_xi1(2) = b;
        coords_xi2(2) = b;
        coords_xi3(2) = a;
        weights(2) = w;
        coords_xi1(3) = b;
        coords_xi2(3) = b;
        coords_xi3(3) = b;
        weights(3) = w;
        a = 0.877978124396165982;
        b = 0.406739585346113397e-1;
        w = 0.167953517588677620e-2;
        coords_xi1(4) = a;
        coords_xi2(4) = b;
        coords_xi3(4) = b;
        weights(4) = w;
        coords_xi1(5) = b;
        coords_xi2(5) = a;
        coords_xi3(5) = b;
        weights(5) = w;
        coords_xi1(6) = b;
        coords_xi2(6) = b;
        coords_xi3(6) = a;
        weights(6) = w;
        coords_xi1(7) = b;
        coords_xi2(7) = b;
        coords_xi3(7) = b;
        weights(7) = w;
        a = 0.329863295731730594e-1;
        b = 0.322337890142275646;
        w = 0.922619692394239843e-2;
        coords_xi1(8) = a;
        coords_xi2(8) = b;
        coords_xi3(8) = b;
        weights(8) = w;
        coords_xi1(9) = b;
        coords_xi2(9) = a;
        coords_xi3(9) = b;
        weights(9) = w;
        coords_xi1(10) = b;
        coords_xi2(10) = b;
        coords_xi3(10) = a;
        weights(10) = w;
        coords_xi1(11) = b;
        coords_xi2(11) = b;
        coords_xi3(11) = b;
        weights(11) = w;
        a = 0.603005664791649076;
        b = 0.269672331458315867;
        c = 0.636610018750175299e-1;
        w = 0.803571428571428248e-2;
        // 12 permutations follows here, a, b, c
        coords_xi1(12) = a;
        coords_xi2(12) = b;
        coords_xi3(12) = c;
        weights(12) = w;
        coords_xi1(13) = b;
        coords_xi2(13) = a;
        coords_xi3(13) = c;
        weights(13) = w;
        coords_xi1(14) = a;
        coords_xi2(14) = c;
        coords_xi3(14) = b;
        weights(14) = w;
        coords_xi1(15) = b;
        coords_xi2(15) = c;
        coords_xi3(15) = a;
        weights(15) = w;
        coords_xi1(16) = c;
        coords_xi2(16) = a;
        coords_xi3(16) = b;
        weights(16) = w;
        coords_xi1(17) = c;
        coords_xi2(17) = b;
        coords_xi3(17) = a;
        weights(17) = w;
        coords_xi1(18) = a;
        coords_xi2(18) = c;
        coords_xi3(18) = c;
        weights(18) = w;
        coords_xi1(19) = b;
        coords_xi2(19) = c;
        coords_xi3(19) = c;
        weights(19) = w;
        coords_xi1(20) = c;
        coords_xi2(20) = a;
        coords_xi3(20) = c;
        weights(20) = w;
        coords_xi1(21) = c;
        coords_xi2(21) = b;
        coords_xi3(21) = c;
        weights(21) = w;
        coords_xi1(22) = c;
        coords_xi2(22) = c;
        coords_xi3(22) = a;
        weights(22) = w;
        coords_xi1(23) = c;
        coords_xi2(23) = c;
        coords_xi3(23) = b;
        weights(23) = w;
        break;

    case 45:     // Exact for x^8
        // Derived by Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
        // See also: Monomial cubature rules since “Stroud”: a compilation
        a = 1. / 4.;
        w = -0.393270066412926145e-1;
        coords_xi1(0) = a;
        coords_xi2(0) = a;
        coords_xi3(0) = a;
        weights(0) = w;
        a = 0.617587190300082967;
        b = 0.127470936566639015;
        w = 0.408131605934270525e-2;
        coords_xi1(1) = a;
        coords_xi2(1) = b;
        coords_xi3(1) = b;
        weights(1) = w;
        coords_xi1(2) = b;
        coords_xi2(2) = a;
        coords_xi3(2) = b;
        weights(2) = w;
        coords_xi1(3) = b;
        coords_xi2(3) = b;
        coords_xi3(3) = a;
        weights(3) = w;
        coords_xi1(4) = b;
        coords_xi2(4) = b;
        coords_xi3(4) = b;
        weights(4) = w;
        a = 0.903763508822103123;
        b = 0.320788303926322960e-1;
        w = 0.658086773304341943e-3;
        coords_xi1(5) = a;
        coords_xi2(5) = b;
        coords_xi3(5) = b;
        weights(5) = w;
        coords_xi1(6) = b;
        coords_xi2(6) = a;
        coords_xi3(6) = b;
        weights(6) = w;
        coords_xi1(7) = b;
        coords_xi2(7) = b;
        coords_xi3(7) = a;
        weights(7) = w;
        coords_xi1(8) = b;
        coords_xi2(8) = b;
        coords_xi3(8) = b;
        weights(8) = w;
        a = 0.450222904356718978;
        b = 0.497770956432810185e-1;
        w = 0.438425882512284693e-2;
        coords_xi1(9) = a;
        coords_xi2(9) = a;
        coords_xi3(9) = b;
        weights(9) = w;
        coords_xi1(10) = a;
        coords_xi2(10) = b;
        coords_xi3(10) = a;
        weights(10) = w;
        coords_xi1(11) = a;
        coords_xi2(11) = b;
        coords_xi3(11) = b;
        weights(11) = w;
        coords_xi1(12) = b;
        coords_xi2(12) = a;
        coords_xi3(12) = a;
        weights(12) = w;
        coords_xi1(13) = b;
        coords_xi2(13) = a;
        coords_xi3(13) = b;
        weights(13) = w;
        coords_xi1(14) = b;
        coords_xi2(14) = b;
        coords_xi3(14) = a;
        weights(14) = w;
        a = 0.316269552601450060;
        b = 0.183730447398549945;
        w = 0.138300638425098166e-1;
        coords_xi1(15) = a;
        coords_xi2(15) = a;
        coords_xi3(15) = b;
        weights(15) = w;
        coords_xi1(16) = a;
        coords_xi2(16) = b;
        coords_xi3(16) = a;
        weights(16) = w;
        coords_xi1(17) = a;
        coords_xi2(17) = b;
        coords_xi3(17) = b;
        weights(17) = w;
        coords_xi1(18) = b;
        coords_xi2(18) = a;
        coords_xi3(18) = a;
        weights(18) = w;
        coords_xi1(19) = b;
        coords_xi2(19) = a;
        coords_xi3(19) = b;
        weights(19) = w;
        coords_xi1(20) = b;
        coords_xi2(20) = b;
        coords_xi3(20) = a;
        weights(20) = w;
        a = 0.513280033360881072;
        b = 0.229177878448171174e-1;
        c = 0.231901089397150906;
        w = 0.424043742468372453e-2;
        coords_xi1(21) = a;
        coords_xi2(21) = b;
        coords_xi3(21) = c;
        weights(21) = w;
        coords_xi1(22) = b;
        coords_xi2(22) = a;
        coords_xi3(22) = c;
        weights(22) = w;
        coords_xi1(23) = a;
        coords_xi2(23) = c;
        coords_xi3(23) = b;
        weights(23) = w;
        coords_xi1(24) = b;
        coords_xi2(24) = c;
        coords_xi3(24) = a;
        weights(24) = w;
        coords_xi1(25) = c;
        coords_xi2(25) = a;
        coords_xi3(25) = b;
        weights(25) = w;
        coords_xi1(26) = c;
        coords_xi2(26) = b;
        coords_xi3(26) = a;
        weights(26) = w;
        coords_xi1(27) = a;
        coords_xi2(27) = c;
        coords_xi3(27) = c;
        weights(27) = w;
        coords_xi1(28) = b;
        coords_xi2(28) = c;
        coords_xi3(28) = c;
        weights(28) = w;
        coords_xi1(29) = c;
        coords_xi2(29) = a;
        coords_xi3(29) = c;
        weights(29) = w;
        coords_xi1(30) = c;
        coords_xi2(30) = b;
        coords_xi3(30) = c;
        weights(30) = w;
        coords_xi1(31) = c;
        coords_xi2(31) = c;
        coords_xi3(31) = a;
        weights(31) = w;
        coords_xi1(32) = c;
        coords_xi2(32) = c;
        coords_xi3(32) = b;
        weights(32) = w;
        a = 0.193746475248804382;
        b = 0.730313427807538396;
        c = 0.379700484718286102e-1;
        w = 0.223873973961420164e-2;
        coords_xi1(33) = a;
        coords_xi2(33) = b;
        coords_xi3(33) = c;
        weights(33) = w;
        coords_xi1(34) = b;
        coords_xi2(34) = a;
        coords_xi3(34) = c;
        weights(34) = w;
        coords_xi1(35) = a;
        coords_xi2(35) = c;
        coords_xi3(35) = b;
        weights(35) = w;
        coords_xi1(36) = b;
        coords_xi2(36) = c;
        coords_xi3(36) = a;
        weights(36) = w;
        coords_xi1(37) = c;
        coords_xi2(37) = a;
        coords_xi3(37) = b;
        weights(37) = w;
        coords_xi1(38) = c;
        coords_xi2(38) = b;
        coords_xi3(38) = a;
        weights(38) = w;
        coords_xi1(39) = a;
        coords_xi2(39) = c;
        coords_xi3(39) = c;
        weights(39) = w;
        coords_xi1(40) = b;
        coords_xi2(40) = c;
        coords_xi3(40) = c;
        weights(40) = w;
        coords_xi1(41) = c;
        coords_xi2(41) = a;
        coords_xi3(41) = c;
        weights(41) = w;
        coords_xi1(42) = c;
        coords_xi2(42) = b;
        coords_xi3(42) = c;
        weights(42) = w;
        coords_xi1(43) = c;
        coords_xi2(43) = c;
        coords_xi3(43) = a;
        weights(43) = w;
        coords_xi1(44) = c;
        coords_xi2(44) = c;
        coords_xi3(44) = b;
        weights(44) = w;

        break;

    default:
        OOFEM_SERROR("unsupported number of IPs (%d)", nPoints);
    }
}

void
GaussIntegrationRule :: giveTriCoordsAndWeights(int nPoints, FloatArray &coords_xi1, FloatArray &coords_xi2, FloatArray &weights)
// Create arrays of coordinates and weights for Gauss Integration Points of a triangle with 'nPoints' integrationpoints
// Coordinates xi1 and xi2 are the two first area coordinates of the triangle.
// Dunavant quadrature rules for triangles in area coordinates. Taken from http://www.mems.rice.edu/~akin/Elsevier/Chap_10.pdf
{
    switch ( nPoints ) {
    case 1:
        coords_xi1 = FloatArray{0.333333333333};
        coords_xi2 = FloatArray{0.333333333333};
        weights = FloatArray{0.5};
        break;

    case 3:
        coords_xi1 = {
            0.166666666666667,
            0.666666666666667,
            0.166666666666667
        };
        coords_xi2 = {
            0.166666666666667,
            0.166666666666667,
            0.666666666666667
        };
        weights = {
            0.166666666666666,
            0.166666666666666,
            0.166666666666666
        };
        break;

    case 4:

        coords_xi1 = {
            0.333333333333333,
            0.200000000000000,
            0.200000000000000,
            0.600000000000000
        };
        coords_xi2 = {
            0.333333333333333,
            0.600000000000000,
            0.200000000000000,
            0.200000000000000
        };
        weights = {
            -0.281250000000000,
            0.260416666666667,
            0.260416666666667,
            0.260416666666667
        };
        break;

    case 6:
        coords_xi1 = {
            0.445948490915965,
            0.445948490915965,
            0.108103018168070,
            0.091576213509771,
            0.091576213509771,
            0.816847572980459
        };
        coords_xi2 = {
            0.108103018168070,
            0.445948490915965,
            0.445948490915965,
            0.816847572980459,
            0.091576213509771,
            0.091576213509771
        };
        weights = {
            0.111690794839006,
            0.111690794839006,
            0.111690794839006,
            0.054975871827661,
            0.054975871827661,
            0.054975871827661
        };
        break;

    case 7:
        coords_xi1 = {
            0.333333333333333,
            0.470142064105115,
            0.470142064105115,
            0.059715871789770,
            0.101286507323456,
            0.101286507323456,
            0.797426985353087
        };
        coords_xi2 = {
            0.333333333333333,
            0.059715871789770,
            0.470142064105115,
            0.470142064105115,
            0.797426985353087,
            0.101286507323456,
            0.101286507323456
        };
        weights = {
            0.112500000000000,
            0.066197076394253,
            0.066197076394253,
            0.066197076394253,
            0.062969590272414,
            0.062969590272414,
            0.062969590272414
        };
        break;

    case 12:
        coords_xi1 = {
            0.249286745170910,
            0.249286745170910,
            0.501426509658179,
            0.063089014491502,
            0.063089014491502,
            0.873821971016996,
            0.636502499121399,
            0.310352451033784,
            0.053145049844817,
            0.636502499121399,
            0.310352451033784,
            0.053145049844817
        };
        coords_xi2 = {
            0.501426509658179,
            0.249286745170910,
            0.249286745170910,
            0.873821971016996,
            0.063089014491502,
            0.063089014491502,
            0.053145049844817,
            0.636502499121399,
            0.310352451033784,
            0.053145049844817,
            0.636502499121399,
            0.310352451033784
        };
        weights = {
            0.058393137863189,
            0.058393137863189,
            0.058393137863189,
            0.025422453185104,
            0.025422453185104,
            0.025422453185104,
            0.041425537809187,
            0.041425537809187,
            0.041425537809187,
            0.041425537809187,
            0.041425537809187,
            0.041425537809187
        };
        break;

    case 13:
        coords_xi1 = {
            0.333333333333333,
            0.260345966079040,
            0.260345966079040,
            0.479308067841920,
            0.065130102902216,
            0.065130102902216,
            0.869739794195568,
            0.638444188569810,
            0.312865496004874,
            0.048690315425316,
            0.638444188569810,
            0.312865496004874,
            0.048690315425316
        };
        coords_xi2 = {
            0.333333333333333,
            0.479308067841920,
            0.260345966079040,
            0.260345966079040,
            0.869739794195568,
            0.065130102902216,
            0.065130102902216,
            0.048690315425316,
            0.638444188569810,
            0.312865496004874,
            0.048690315425316,
            0.638444188569810,
            0.312865496004874
        };
        weights = {
            -0.074785022233841,
            0.087807628716604,
            0.087807628716604,
            0.087807628716604,
            0.026673617804419,
            0.026673617804419,
            0.026673617804419,
            0.038556880445128,
            0.038556880445128,
            0.038556880445128,
            0.038556880445128,
            0.038556880445128,
            0.038556880445128
        };
        break;

    case 16:
        coords_xi1 = {
            0.333333333333333,
            0.459292588292723,
            0.459292588292723,
            0.081414823414554,
            0.170569307751760,
            0.170569307751760,
            0.658861384496480,
            0.050547228317031,
            0.050547228317031,
            0.898905543365938,
            0.728492392955404,
            0.263112829634638,
            0.008394777409958,
            0.728492392955404,
            0.263112829634638,
            0.008394777409958
        };
        coords_xi2 = {
            0.333333333333333,
            0.081414823414554,
            0.459292588292723,
            0.459292588292723,
            0.658861384496480,
            0.170569307751760,
            0.170569307751760,
            0.898905543365938,
            0.050547228317031,
            0.050547228317031,
            0.008394777409958,
            0.728492392955404,
            0.263112829634638,
            0.008394777409958,
            0.728492392955404,
            0.263112829634638
        };
        weights = {
            0.072157803838894,
            0.047545817133642,
            0.047545817133642,
            0.047545817133642,
            0.051608685267359,
            0.051608685267359,
            0.051608685267359,
            0.016229248811599,
            0.016229248811599,
            0.016229248811599,
            0.013615157087218,
            0.013615157087218,
            0.013615157087218,
            0.013615157087218,
            0.013615157087218,
            0.013615157087218
        };
        break;


    case 19:
        coords_xi1 = {
            0.333333333333333,
            0.489682519198738,
            0.489682519198738,
            0.020634961602525,
            0.437089591492937,
            0.437089591492937,
            0.125820817014127,
            0.188203535619033,
            0.188203535619033,
            0.623592928761935,
            0.044729513394453,
            0.044729513394453,
            0.910540973211095,
            0.741198598784498,
            0.221962989160766,
            0.036838412054736,
            0.741198598784498,
            0.221962989160766,
            0.036838412054736
        };
        coords_xi2 = {
            0.333333333333333,
            0.020634961602525,
            0.489682519198738,
            0.489682519198738,
            0.125820817014127,
            0.437089591492937,
            0.437089591492937,
            0.623592928761935,
            0.188203535619033,
            0.188203535619033,
            0.910540973211095,
            0.044729513394453,
            0.044729513394453,
            0.036838412054736,
            0.741198598784498,
            0.221962989160766,
            0.036838412054736,
            0.741198598784498,
            0.221962989160766
        };
        weights = {
            0.048567898141400,
            0.015667350113570,
            0.015667350113570,
            0.015667350113570,
            0.038913770502387,
            0.038913770502387,
            0.038913770502387,
            0.039823869463605,
            0.039823869463605,
            0.039823869463605,
            0.012788837829349,
            0.012788837829349,
            0.012788837829349,
            0.021641769688645,
            0.021641769688645,
            0.021641769688645,
            0.021641769688645,
            0.021641769688645,
            0.021641769688645
        };
        break;

    case 25:
        coords_xi1 = {
            0.333333333333333,
            0.485577633383657,
            0.485577633383657,
            0.028844733232685,
            0.109481575485037,
            0.109481575485037,
            0.781036849029926,
            0.550352941820999,
            0.307939838764121,
            0.141707219414880,
            0.550352941820999,
            0.307939838764121,
            0.141707219414880,
            0.728323904597411,
            0.246672560639903,
            0.025003534762686,
            0.728323904597411,
            0.246672560639903,
            0.025003534762686,
            0.923655933587500,
            0.066803251012200,
            0.009540815400299,
            0.923655933587500,
            0.066803251012200,
            0.009540815400299
        };
        coords_xi2 = {
            0.333333333333333,
            0.028844733232685,
            0.485577633383657,
            0.485577633383657,
            0.781036849029926,
            0.109481575485037,
            0.109481575485037,
            0.141707219414880,
            0.550352941820999,
            0.307939838764121,
            0.141707219414880,
            0.550352941820999,
            0.307939838764121,
            0.025003534762686,
            0.728323904597411,
            0.246672560639903,
            0.025003534762686,
            0.728323904597411,
            0.246672560639903,
            0.009540815400299,
            0.923655933587500,
            0.066803251012200,
            0.009540815400299,
            0.923655933587500,
            0.066803251012200
        };
        weights = {
            0.045408995191377,
            0.018362978878233,
            0.018362978878233,
            0.018362978878233,
            0.022660529717764,
            0.022660529717764,
            0.022660529717764,
            0.036378958422710,
            0.036378958422710,
            0.036378958422710,
            0.036378958422710,
            0.036378958422710,
            0.036378958422710,
            0.014163621265529,
            0.014163621265529,
            0.014163621265529,
            0.014163621265529,
            0.014163621265529,
            0.014163621265529,
            0.004710833481867,
            0.004710833481867,
            0.004710833481867,
            0.004710833481867,
            0.004710833481867,
            0.004710833481867
        };
        break;

    default:
        OOFEM_SERROR("unsupported number of IPs (%d)", nPoints);
    }
}



void
GaussIntegrationRule :: giveLineCoordsAndWeights(int nPoints, FloatArray &coords_xi, FloatArray &weights)
// Create arrays of coordinates and weights for Gauss Integration Points of a line with 'nPoints' integrationpoints
{
    switch ( nPoints ) {
    case 1:
        coords_xi = FloatArray{0.0};
        weights = FloatArray{2.0};
        break;

    case 2:
        coords_xi = {-0.577350269189626, 0.577350269189626};
        weights = {1.0, 1.0};
        break;

    case 3:
        coords_xi = {
            -0.774596669241483,
            0.0,
            0.774596669241483
        };
        weights = {
            0.555555555555555,
            0.888888888888888,
            0.555555555555555
        };
        break;

    case 4:
        coords_xi = {
            -0.861136311594053,
            -0.339981043584856,
            0.339981043584856,
            0.861136311594053
        };
        weights = {
            0.347854845137454,
            0.652145154862546,
            0.652145154862546,
            0.347854845137454
        };
        break;

    case 5:
        coords_xi = {
            -0.9061798459386639927976269,
            -0.5384693101056830910363144,
            0.0,
            0.5384693101056830910363144,
            0.9061798459386639927976269
        };
        weights = {
            0.2369268850561890875142640,
            0.4786286704993664680412915,
            0.5688888888888888888888889,
            0.4786286704993664680412915,
            0.2369268850561890875142640
        };
        break;

    case 6:
        coords_xi = {
            -0.2386191860831969086305017,
            -0.6612093864662645136613996,
            -0.9324695142031520278123016,
            0.9324695142031520278123016,
            0.6612093864662645136613996,
            0.2386191860831969086305017
        };
        weights = {
            0.4679139345726910473898703,
            0.3607615730481386075698335,
            0.1713244923791703450402961,
            0.1713244923791703450402961,
            0.3607615730481386075698335,
            0.4679139345726910473898703
        };
        break;

    case 7:
        coords_xi = {
            -0.9491079123427585245261897,
            -0.7415311855993944398638648,
            -0.4058451513773971669066064,
            0.0,
            0.4058451513773971669066064,
            0.7415311855993944398638648,
            0.9491079123427585245261897
        };
        weights = {
            0.1294849661688696932706114,
            0.2797053914892766679014678,
            0.3818300505051189449503698,
            0.4179591836734693877551020,
            0.3818300505051189449503698,
            0.2797053914892766679014678,
            0.1294849661688696932706114
        };
        break;

    case 8:
        coords_xi = {
            -0.960289856497536,
            -0.796666477413627,
            -0.525532409916329,
            -0.183434642495650,
            0.183434642495650,
            0.525532409916329,
            0.796666477413627,
            0.960289856497536
        };
        weights = {
            0.101228536290375,
            0.222381034453374,
            0.313706645877887,
            0.362683783378362,
            0.362683783378362,
            0.313706645877887,
            0.222381034453374,
            0.101228536290375
        };
        break;


    case 16:
        coords_xi = {
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
        weights = {
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
        break;

    case 24:
        coords_xi = {
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
        weights = {
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
        break;

    case 32:
        coords_xi = {
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
        weights = {
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
        break;

    case 64:
        coords_xi = {
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
        weights = {
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
        break;

    default:
        OOFEM_SERROR("unsupported number of IPs (%d)", nPoints);
    }
}
} // end namespace oofem
