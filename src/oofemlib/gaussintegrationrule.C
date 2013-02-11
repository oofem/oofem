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
GaussIntegrationRule :: SetUpPointsOnLine(int nPoints, MaterialMode mode)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store corrdinates and weights)
{
    FloatArray coords_xi, weights;
    this->giveLineCoordsAndWeights(nPoints, coords_xi, weights);
    this->numberOfIntegrationPoints = nPoints;
    this->gaussPointArray  = new GaussPoint * [ nPoints ];

    for ( int i = 1; i <= nPoints; i++ ) {
        FloatArray *coord = new FloatArray(1);
        coord->at(1) = coords_xi.at(i);
        this->gaussPointArray [ i - 1 ] = new GaussPoint(this, i, coord, weights.at(i), mode);
    }

    return nPoints;
}

int
GaussIntegrationRule :: SetUpPointsOnSquare(int nPoints, MaterialMode mode)
/*GaussIntegrationRule :: SetUpPointsOnSquare(int nPoints_xi1, int nPoints_xi2, MaterialMode mode)
 * Creates an array of nPoints Gauss Integration Points. The points are set up as a product between two
 * 1D integration rules.
 * (Don't confuse with GaussPoint - elem is only the container where to store corrdinates and weights)
 */
{
    int nPoints_xi1 = floor( sqrt( double( nPoints ) ) );
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
            this->gaussPointArray [ count - 1 ] = new GaussPoint(this, count, coord, weights1.at(i) * weights2.at(j), mode);
        }
    }

    return this->numberOfIntegrationPoints;
}

int
GaussIntegrationRule :: SetUpPointsOnCube(int nPoints, MaterialMode mode)
/*GaussIntegrationRule :: SetUpPointsOnCube(int nPoints_xi1, int nPoints_xi2, int nPoints_xi3, MaterialMode mode)
 * Creates an array of nPoints Gauss Integration Points. The points are set up as a product between three
 * 1D integration rules.
 * (Don't confuse with GaussPoint - elem is only the container where to store corrdinates and weights)
 */
{
    int nPoints_xi1 = floor(cbrt( double( nPoints ) ) + 0.5);
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
                this->gaussPointArray [ count - 1 ] = new GaussPoint(this, count, coord, weights1.at(i) * weights2.at(j) * weights3.at(k), mode);
            }
        }
    }

    return this->numberOfIntegrationPoints;
}



int
GaussIntegrationRule :: SetUpPointsOnTriangle(int nPoints, MaterialMode mode)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store coordinates and weights)
{
    FloatArray coords_xi1, coords_xi2, weights;
    this->giveTriCoordsAndWeights(nPoints, coords_xi1, coords_xi2, weights);
    this->numberOfIntegrationPoints = nPoints;
    this->gaussPointArray  = new GaussPoint * [ nPoints ];

    for ( int i = 1; i <= nPoints; i++ ) {
        FloatArray *coord = new FloatArray(2);
        coord->at(1) = coords_xi1.at(i);
        coord->at(2) = coords_xi2.at(i);
        this->gaussPointArray [ i - 1 ] = new GaussPoint(this, i, coord, weights.at(i), mode);
    }

    return nPoints;
}

int
GaussIntegrationRule :: SetUpPointsOnTetrahedra(int nPoints, MaterialMode mode)
{
    double a, b, c, w;
    FloatArray *coord1;
    this->numberOfIntegrationPoints = nPoints;
    switch ( nPoints ) {
    case 1:

        this->gaussPointArray  = new GaussPoint * [ nPoints ];
        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.25;
        coord1->at(2)    = 0.25;
        coord1->at(3)    = 0.25;
        ( this->gaussPointArray ) [ 0 ]        = new GaussPoint(this, 1, coord1, 1. / 6., mode);
        break;

    case 4:
        // quadratic formulae

        this->gaussPointArray  = new GaussPoint * [ nPoints ];

        a = ( 5. + 3. * sqrt(5.) ) / 20.;
        b = ( 5. - sqrt(5.) ) / 20.;
        w = 1. / 24.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 0 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 1 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 2 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 3 ]        = new GaussPoint(this, 1, coord1, w, mode);

        break;

    case 5:
        // cubic formulae

        this->gaussPointArray  = new GaussPoint * [ nPoints ];

        coord1 = new FloatArray(3);
        coord1->at(1)    = 0.25;
        coord1->at(2)    = 0.25;
        coord1->at(3)    = 0.25;
        ( this->gaussPointArray ) [ 0 ]        = new GaussPoint(this, 1, coord1, -2. / 15., mode);

        a = 1. / 2.;
        b = 1. / 6.;
        w = 3. / 40.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 1 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 2 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 3 ]        = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 4 ]        = new GaussPoint(this, 1, coord1, w, mode);

        break;

    case 11:
        // exact x^4
        this->gaussPointArray  = new GaussPoint * [ nPoints ];

        a = 1. / 4.;
        w = -74. / 5625.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 0 ]        = new GaussPoint(this, 1, coord1, w, mode);

        a = 5. / 70.;
        b = 11. / 14.;
        w = 343. / 45000.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 1 ]        = new GaussPoint(this, 2, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 2 ]        = new GaussPoint(this, 3, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 3 ]        = new GaussPoint(this, 4, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 4 ]        = new GaussPoint(this, 5, coord1, w, mode);

        a = ( 1. + sqrt(5. / 14.) ) / 4.;
        b = ( 1. - sqrt(5. / 14.) ) / 4.;
        w = 28. / 1125.;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 5 ]        = new GaussPoint(this, 6, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 6 ]        = new GaussPoint(this, 7, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 7 ]        = new GaussPoint(this, 8, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 8 ]        = new GaussPoint(this, 9, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 9 ]        = new GaussPoint(this, 10, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 10 ]        = new GaussPoint(this, 11, coord1, w, mode);

        break;

    case 15:
        // Derived by Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
        // exact for x^5
        this->gaussPointArray  = new GaussPoint * [ nPoints ];

        a = 1. / 4.;
        w = 0.302836780970891856e-1;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 0 ]        = new GaussPoint(this, 1, coord1, w, mode);

        a = 0.0;
        b = 1. / 3.;
        w = 0.602678571428571597e-2;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 1 ]        = new GaussPoint(this, 2, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 2 ]        = new GaussPoint(this, 3, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 3 ]        = new GaussPoint(this, 4, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 4 ]        = new GaussPoint(this, 5, coord1, w, mode);

        //////////////
        a = 8. / 11.;
        b = 1. / 11.;
        w = 0.116452490860289742e-1;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 5 ]        = new GaussPoint(this, 6, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 6 ]        = new GaussPoint(this, 7, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 7 ]        = new GaussPoint(this, 8, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 8 ]        = new GaussPoint(this, 9, coord1, w, mode);

        //////////////
        a = 0.665501535736642813e-1;
        b = 0.433449846426335728;
        w = 0.109491415613864534e-1;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 9 ]        = new GaussPoint(this, 10, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 10 ]        = new GaussPoint(this, 11, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 11 ]        = new GaussPoint(this, 12, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 12 ]        = new GaussPoint(this, 13, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 13 ]        = new GaussPoint(this, 14, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 14 ]        = new GaussPoint(this, 15, coord1, w, mode);

        break;

    case 24:
        // Derived by Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
        // See also: Monomial cubature rules since “Stroud”: a compilation
        // Exact for x^6
        this->gaussPointArray  = new GaussPoint * [ nPoints ];
        a = 0.356191386222544953;
        b = 0.214602871259151684;
        w = 0.665379170969464506e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 0 ] = new GaussPoint(this, 1, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 1 ] = new GaussPoint(this, 2, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 2 ] = new GaussPoint(this, 3, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 3 ] = new GaussPoint(this, 4, coord1, w, mode);

        a = 0.877978124396165982;
        b = 0.406739585346113397e-1;
        w = 0.167953517588677620e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 4 ] = new GaussPoint(this, 5, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 5 ] = new GaussPoint(this, 6, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 6 ] = new GaussPoint(this, 7, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 7 ] = new GaussPoint(this, 8, coord1, w, mode);

        a = 0.329863295731730594e-1;
        b = 0.322337890142275646;
        w = 0.922619692394239843e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 8 ] = new GaussPoint(this, 9, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 9 ] = new GaussPoint(this, 10, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 10 ] = new GaussPoint(this, 11, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 11 ] = new GaussPoint(this, 12, coord1, w, mode);

        a = 0.603005664791649076;
        b = 0.269672331458315867;
        c = 0.636610018750175299e-1;
        w = 0.803571428571428248e-2;
        // 12 permutations follows here, a, b, c
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 12 ] = new GaussPoint(this, 13, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 13 ] = new GaussPoint(this, 14, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 14 ] = new GaussPoint(this, 15, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 15 ] = new GaussPoint(this, 16, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 16 ] = new GaussPoint(this, 17, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 17 ] = new GaussPoint(this, 18, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 18 ] = new GaussPoint(this, 19, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 19 ] = new GaussPoint(this, 20, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 20 ] = new GaussPoint(this, 21, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 21 ] = new GaussPoint(this, 22, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 22 ] = new GaussPoint(this, 23, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 23 ] = new GaussPoint(this, 24, coord1, w, mode);

        break;

    case 45:
        // Derived by Patrick Keast, MODERATE-DEGREE TETRAHEDRAL QUADRATURE FORMULAS
        // See also: Monomial cubature rules since “Stroud”: a compilation
        // Exact for x^8
        this->gaussPointArray  = new GaussPoint * [ nPoints ];

        a = 1. / 4.;
        w = -0.393270066412926145e-1;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = a;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 0 ] = new GaussPoint(this, 1, coord1, w, mode);

        a = 0.617587190300082967;
        b = 0.127470936566639015;
        w = 0.408131605934270525e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 1 ] = new GaussPoint(this, 2, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 2 ] = new GaussPoint(this, 3, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 3 ] = new GaussPoint(this, 4, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 4 ] = new GaussPoint(this, 5, coord1, w, mode);

        a = 0.903763508822103123;
        b = 0.320788303926322960e-1;
        w = 0.658086773304341943e-3;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 5 ] = new GaussPoint(this, 6, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 6 ] = new GaussPoint(this, 7, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 7 ] = new GaussPoint(this, 8, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = b;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 8 ] = new GaussPoint(this, 9, coord1, w, mode);

        a = 0.450222904356718978;
        b = 0.497770956432810185e-1;
        w = 0.438425882512284693e-2;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 9 ]        = new GaussPoint(this, 10, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 10 ]        = new GaussPoint(this, 11, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 11 ]        = new GaussPoint(this, 12, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 12 ]        = new GaussPoint(this, 13, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 13 ]        = new GaussPoint(this, 14, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 14 ]        = new GaussPoint(this, 15, coord1, w, mode);

        a = 0.316269552601450060;
        b = 0.183730447398549945;
        w = 0.138300638425098166e-1;
        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 15 ] = new GaussPoint(this, 16, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 16 ] = new GaussPoint(this, 17, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = a;
        coord1->at(2)    = b;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 17 ] = new GaussPoint(this, 18, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 18 ] = new GaussPoint(this, 19, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = a;
        coord1->at(3)    = b;
        ( this->gaussPointArray ) [ 19 ] = new GaussPoint(this, 20, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1)    = b;
        coord1->at(2)    = b;
        coord1->at(3)    = a;
        ( this->gaussPointArray ) [ 20 ] = new GaussPoint(this, 21, coord1, w, mode);

        a = 0.513280033360881072;
        b = 0.229177878448171174e-1;
        c = 0.231901089397150906;
        w = 0.424043742468372453e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 21 ] = new GaussPoint(this, 22, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 22 ] = new GaussPoint(this, 23, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 23 ] = new GaussPoint(this, 24, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 24 ] = new GaussPoint(this, 25, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 25 ] = new GaussPoint(this, 26, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 26 ] = new GaussPoint(this, 27, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 27 ] = new GaussPoint(this, 28, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 28 ] = new GaussPoint(this, 29, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 29 ] = new GaussPoint(this, 30, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 30 ] = new GaussPoint(this, 31, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 31 ] = new GaussPoint(this, 32, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 32 ] = new GaussPoint(this, 33, coord1, w, mode);

        a = 0.193746475248804382;
        b = 0.730313427807538396;
        c = 0.379700484718286102e-1;
        w = 0.223873973961420164e-2;
        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 33 ] = new GaussPoint(this, 34, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 34 ] = new GaussPoint(this, 35, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 35 ] = new GaussPoint(this, 36, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 36 ] = new GaussPoint(this, 37, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 37 ] = new GaussPoint(this, 38, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 38 ] = new GaussPoint(this, 39, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = a;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 39 ] = new GaussPoint(this, 40, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = b;
        coord1->at(2) = c;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 40 ] = new GaussPoint(this, 41, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = a;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 41 ] = new GaussPoint(this, 42, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = b;
        coord1->at(3) = c;
        ( this->gaussPointArray ) [ 42 ] = new GaussPoint(this, 43, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = a;
        ( this->gaussPointArray ) [ 43 ] = new GaussPoint(this, 44, coord1, w, mode);

        coord1 = new FloatArray(3);
        coord1->at(1) = c;
        coord1->at(2) = c;
        coord1->at(3) = b;
        ( this->gaussPointArray ) [ 44 ] = new GaussPoint(this, 45, coord1, w, mode);

        break;

    default:
        OOFEM_ERROR2("SetUpPointsOnTetrahedra: unsupported number of IPs (%d)", nPoints);
    }

    return nPoints;
}







int
GaussIntegrationRule :: SetUpPointsOnWedge(int nPointsTri, int nPointsDepth, MaterialMode mode)
/* Creates an array of (nPointsTri*nPointsDepth) Gauss Integration Points. The points are set up as a product between a
 * triangular integration rule in the plane and a 1D integration rule in the depth direction.
 * (Don't confuse with GaussPoint - elem is only the container where to store coordinates and weights)
 */
{
    FloatArray coords_xi1, coords_xi2, coords_xi3, weightsTri, weightsDepth;
    this->giveTriCoordsAndWeights(nPointsTri, coords_xi1, coords_xi2, weightsTri);
    this->giveLineCoordsAndWeights(nPointsDepth, coords_xi3, weightsDepth);
    this->numberOfIntegrationPoints = nPointsTri * nPointsDepth;
    this->gaussPointArray = new GaussPoint * [ this->numberOfIntegrationPoints ];

    for ( int i = 1, ind = 0; i <= nPointsTri; i++ ) {
        for ( int j = 1; j <= nPointsDepth; j++ ) {
            FloatArray *coord = new FloatArray(3);
            coord->at(1) = coords_xi1.at(i);
            coord->at(2) = coords_xi2.at(i);
            coord->at(3) = coords_xi3.at(j);
            this->gaussPointArray [ ind ] = new GaussPoint(this, 1, coord, weightsTri.at(i) * weightsDepth.at(j), mode);
            ind++;
        }
    }

    return 1;
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

        if ( requiredNIP <= 1 ) {
            return 1;
        }

        if ( requiredNIP <= 2 ) {
            return 2;
        }

        if ( requiredNIP <= 3 ) {
            return 3;
        }

        if ( requiredNIP <= 4 ) {
            return 4;
        }

        if ( requiredNIP <= 8 ) {
            return 8;
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

        return requiredNIP;

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
        requiredNIP = max( ( approxOrder + 1 ) / 2, 2 );
        requiredNIP *= requiredNIP;
        if ( requiredNIP > 64 * 64 ) {
            return -1;
        }

        return requiredNIP;

    case _Cube:
        requiredNIP = max( ( approxOrder + 1 ) / 2, 2 );
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

        return -1;

    default:
        OOFEM_ERROR("GaussIntegrationRule::setUpIntegrationPoints - unknown integrationDomain");
    }

    return -1;
}


int
GaussIntegrationRule :: SetUpPointsOn2DEmbeddedLine(int nPoints, MaterialMode mode, const FloatArray **coords)
// creates array of nPoints Gauss Integration Points
// ( don't confuse with GaussPoint - elem is only the container where to
//   store coordinates and weights)
{
    double weight, l;
    FloatArray *coord1;


    switch ( nPoints ) {
    case 1:
        this->gaussPointArray = new GaussPoint * [ nPoints ];
        coord1 = new FloatArray(2);
        l = 0.0;       //local coordinate on the line
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        weight = 2.0;
        ( this->gaussPointArray ) [ 0 ] = new GaussPoint(this, 1, coord1, weight, mode);
        break;

    case 2:

        this->gaussPointArray              = new GaussPoint * [ nPoints ];
        coord1             = new FloatArray(2);

        l = -0.577350269189626; //local coordinate on the line
        weight = 1.0;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( this->gaussPointArray ) [ 0 ]         = new GaussPoint(this, 1, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l = 0.577350269189626;        //local coordinate on the line
        weight = 1.0;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( this->gaussPointArray ) [ 1 ]         = new GaussPoint(this, 2, coord1, weight, mode);


        break;

    case 3:

        this->gaussPointArray              = new GaussPoint * [ nPoints ];
        coord1             = new FloatArray(2);

        l = -0.774596669241483; //local coordinate on the line
        weight =  0.555555555555555;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( this->gaussPointArray ) [ 0 ]         = new GaussPoint(this, 1, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l = 0.0;        //local coordinate on the line
        weight = 0.888888888888888;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( this->gaussPointArray ) [ 1 ]         = new GaussPoint(this, 2, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l = 0.774596669241483; //local coordinate on the line
        weight =  0.555555555555555;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( this->gaussPointArray ) [ 2 ]         = new GaussPoint(this, 3, coord1, weight, mode);
        break;

    case 4:


        this->gaussPointArray              = new GaussPoint * [ nPoints ];
        coord1             = new FloatArray(2);

        l = -0.861136311594053; //local coordinate on the line
        weight =  0.347854845137454;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( this->gaussPointArray ) [ 0 ]         = new GaussPoint(this, 1, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l = -0.339981043584856;        //local coordinate on the line
        weight = 0.652145154862546;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( this->gaussPointArray ) [ 1 ]         = new GaussPoint(this, 2, coord1, weight, mode);

        coord1             = new FloatArray(2);

        l =  0.339981043584856; //local coordinate on the line
        weight =   0.652145154862546;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( this->gaussPointArray ) [ 2 ]         = new GaussPoint(this, 3, coord1, weight, mode);

        l =  0.861136311594053;        //local coordinate on the line
        weight =   0.347854845137454;
        coord1->at(1) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(1) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(1);
        coord1->at(2) = ( 1. - ( l + 1 ) * 0.5 ) * coords [ 0 ]->at(2) + ( l + 1 ) * 0.5 * coords [ 1 ]->at(2);
        ( this->gaussPointArray ) [ 3 ]         = new GaussPoint(this, 4, coord1, weight, mode);
        break;


    default:
        OOFEM_ERROR2("GaussIntegrationRule :: SetUpPointsOn2DEmbeddedLine - order %d of integration is not suported\n", nPoints);
    }

    // Initialize the material status at each Gauss point
    //for (int i=0; i<numberOfIntegrationPoints; i++)
    //  mat->giveStatus(gaussPointArray[i]);

    return nPoints;
}


void
GaussIntegrationRule :: giveTriCoordsAndWeights(int nPoints, FloatArray &coords_xi1, FloatArray &coords_xi2, FloatArray &weights)
// Create arrays of coordinates and weights for Gauss Integration Points of a trinagle with 'nPoints' integrationpoints
// Coordinates xi1 and xi2 are the two firts area coordinates of the tringle.
// Dunavant quadrature rules for triangles in area coordinates. Taken from http://www.mems.rice.edu/~akin/Elsevier/Chap_10.pdf
{
    coords_xi1.resize(nPoints);
    coords_xi2.resize(nPoints);
    weights.resize(nPoints);
    switch ( nPoints ) {
    case 1:
        coords_xi1.setValues(1, 0.333333333333);
        coords_xi2.setValues(1, 0.333333333333);
        weights.setValues(1, 0.5);
        break;

    case 3:
        coords_xi1.setValues(3,
                             0.166666666666667,
                             0.666666666666667,
                             0.166666666666667
                             );

        coords_xi2.setValues(3,
                             0.166666666666667,
                             0.166666666666667,
                             0.666666666666667
                             );
        weights.setValues(3,
                          0.166666666666666,
                          0.166666666666666,
                          0.166666666666666
                          );

        break;

    case 4:

        coords_xi1.setValues(4,
                             0.333333333333333,
                             0.200000000000000,
                             0.200000000000000,
                             0.600000000000000
                             );

        coords_xi2.setValues(4,
                             0.333333333333333,
                             0.600000000000000,
                             0.200000000000000,
                             0.200000000000000
                             );

        weights.setValues(4,
                          -0.281250000000000,
                          0.260416666666667,
                          0.260416666666667,
                          0.260416666666667
                          );
        break;

    case 6:

        coords_xi1.setValues(6,
                             0.445948490915965,
                             0.445948490915965,
                             0.108103018168070,
                             0.091576213509771,
                             0.091576213509771,
                             0.816847572980459
                             );

        coords_xi2.setValues(6,
                             0.108103018168070,
                             0.445948490915965,
                             0.445948490915965,
                             0.816847572980459,
                             0.091576213509771,
                             0.091576213509771
                             );

        weights.setValues(6,
                          0.111690794839006,
                          0.111690794839006,
                          0.111690794839006,
                          0.054975871827661,
                          0.054975871827661,
                          0.054975871827661
                          );
        break;

    case 7:

        coords_xi1.setValues(7,
                             0.333333333333333,
                             0.470142064105115,
                             0.470142064105115,
                             0.059715871789770,
                             0.101286507323456,
                             0.101286507323456,
                             0.797426985353087
                             );
        coords_xi2.setValues(7,
                             0.333333333333333,
                             0.059715871789770,
                             0.470142064105115,
                             0.470142064105115,
                             0.797426985353087,
                             0.101286507323456,
                             0.101286507323456
                             );

        weights.setValues(7,
                          0.112500000000000,
                          0.066197076394253,
                          0.066197076394253,
                          0.066197076394253,
                          0.062969590272414,
                          0.062969590272414,
                          0.062969590272414
                          );

        break;

    case 12:

        coords_xi1.setValues(12,
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
                             );

        coords_xi2.setValues(12,
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
                             );

        weights.setValues(12,
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
                          );

        break;

    case 13:
        coords_xi1.setValues(13,
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
                             );

        coords_xi2.setValues(13,
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
                             );

        weights.setValues(13,
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
                          );
        break;

    case 16:

        coords_xi1.setValues(16,
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
                             );
        coords_xi2.setValues(16,
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
                             );

        weights.setValues(16,
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
                          );

        break;


    case 19:

        coords_xi1.setValues(19,
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
                             );
        coords_xi2.setValues(19,
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
                             );

        weights.setValues(19,
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
                          );

        break;

    case 25:

        coords_xi1.setValues(25,
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
                             );
        coords_xi2.setValues(25,
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
                             );

        weights.setValues(25,
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
                          );

        break;

    default:
        OOFEM_ERROR2("giveTriCoordsAndWeights: unsupported number of IPs (%d)", nPoints);
    }
}



void
GaussIntegrationRule :: giveLineCoordsAndWeights(int nPoints, FloatArray &coords_xi, FloatArray &weights)
// Create arrays of coordinates and weights for Gauss Integration Points of a line with 'nPoints' integrationpoints
{
    coords_xi.resize(nPoints);
    weights.resize(nPoints);

    switch ( nPoints ) {
    case 1:

        coords_xi.setValues(1, 0.0);
        weights.setValues(1, 2.0);

        break;

    case 2:

        coords_xi.setValues(2,
                            -0.577350269189626,
                            0.577350269189626
                            );

        weights.setValues(2,
                          1.0,
                          1.0
                          );

        break;

    case 3:

        coords_xi.setValues(3,
                            -0.774596669241483,
                            0.0,
                            0.774596669241483
                            );

        weights.setValues(3,
                          0.555555555555555,
                          0.888888888888888,
                          0.555555555555555
                          );

        break;

    case 4:

        coords_xi.setValues(4,
                            -0.861136311594053,
                            -0.339981043584856,
                            0.339981043584856,
                            0.861136311594053
                            );

        weights.setValues(4,
                          0.347854845137454,
                          0.652145154862546,
                          0.652145154862546,
                          0.347854845137454
                          );

        break;

    case 5:

        coords_xi.setValues(5,
                            -0.9061798459386639927976269,
                            -0.5384693101056830910363144,
                            0.0,
                            0.5384693101056830910363144,
                            0.9061798459386639927976269
                            );

        weights.setValues(5,
                          0.2369268850561890875142640,
                          0.4786286704993664680412915,
                          0.5688888888888888888888889,
                          0.4786286704993664680412915,
                          0.2369268850561890875142640
                          );

        break;

    case 6:

        coords_xi.setValues(6,
                            -0.2386191860831969086305017,
                            -0.6612093864662645136613996,
                            -0.9324695142031520278123016,
                            0.9324695142031520278123016,
                            0.6612093864662645136613996,
                            0.2386191860831969086305017
                            );

        weights.setValues(6,
                          0.4679139345726910473898703,
                          0.3607615730481386075698335,
                          0.1713244923791703450402961,
                          0.1713244923791703450402961,
                          0.3607615730481386075698335,
                          0.4679139345726910473898703
                          );

        break;

    case 7:

        coords_xi.setValues(7,
                            -0.9491079123427585245261897,
                            -0.7415311855993944398638648,
                            -0.4058451513773971669066064,
                            0.0,
                            0.4058451513773971669066064,
                            0.7415311855993944398638648,
                            0.9491079123427585245261897
                            );

        weights.setValues(7,
                          0.1294849661688696932706114,
                          0.2797053914892766679014678,
                          0.3818300505051189449503698,
                          0.4179591836734693877551020,
                          0.3818300505051189449503698,
                          0.2797053914892766679014678,
                          0.1294849661688696932706114
                          );

        break;

    case 8:

        coords_xi.setValues(8,
                            -0.960289856497536,
                            -0.796666477413627,
                            -0.525532409916329,
                            -0.183434642495650,
                            0.183434642495650,
                            0.525532409916329,
                            0.796666477413627,
                            0.960289856497536
                            );

        weights.setValues(8,
                          0.101228536290375,
                          0.222381034453374,
                          0.313706645877887,
                          0.362683783378362,
                          0.362683783378362,
                          0.313706645877887,
                          0.222381034453374,
                          0.101228536290375
                          );

        break;


    case 16:
        coords_xi.setValues(16,
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
                            );

        weights.setValues(16,
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
                          );

        break;


    case 24:

        coords_xi.setValues(24,
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
                            );

        weights.setValues(24,
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
                          );

        break;



    case 32:

        coords_xi.setValues(32,
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
                            );

        weights.setValues(32,
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
                          );

        break;

    case 64:

        coords_xi.setValues(64,
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
                            );

        weights.setValues(64,
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
                          );

        break;

    default:
        OOFEM_ERROR2("SetUpPointsOnLine: unsupported number of IPs (%d)", nPoints);
    }
}
} // end namespace oofem
