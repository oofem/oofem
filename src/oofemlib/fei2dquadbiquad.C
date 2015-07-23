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

#include "fei2dquadbiquad.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
void
FEI2dQuadBiQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double u, v;

    u = lcoords.at(1);
    v = lcoords.at(2);

    double a[] = {
        0.5 * ( u - 1.0 ) * u, 1.0 - u * u, 0.5 * ( u + 1.0 ) * u
    };
    double b[] = {
        0.5 * ( v - 1.0 ) * v, 1.0 - v * v, 0.5 * ( v + 1.0 ) * v
    };

    answer.resize(9);
    answer.at(1) = a [ 0 ] * b [ 0 ];
    answer.at(5) = a [ 1 ] * b [ 0 ];
    answer.at(2) = a [ 2 ] * b [ 0 ];

    answer.at(8) = a [ 0 ] * b [ 1 ];
    answer.at(9) = a [ 1 ] * b [ 1 ];
    answer.at(6) = a [ 2 ] * b [ 1 ];

    answer.at(4) = a [ 0 ] * b [ 2 ];
    answer.at(7) = a [ 1 ] * b [ 2 ];
    answer.at(3) = a [ 2 ] * b [ 2 ];
}


void
FEI2dQuadBiQuad :: giveDerivatives(FloatMatrix &dN, const FloatArray &lc)
{
    double u = lc.at(1);
    double v = lc.at(2);

    double a[] = {
        0.5 * ( u - 1.0 ) * u, 1.0 - u * u, 0.5 * ( u + 1.0 ) * u
    };
    double b[] = {
        0.5 * ( v - 1.0 ) * v, 1.0 - v * v, 0.5 * ( v + 1.0 ) * v
    };

    double da[] = {
        u - 0.5, -2.0 * u, u + 0.5
    };
    double db[] = {
        v - 0.5, -2.0 * v, v + 0.5
    };

    dN.resize(9, 2);

    dN.at(1, 1) = da [ 0 ] * b [ 0 ];
    dN.at(5, 1) = da [ 1 ] * b [ 0 ];
    dN.at(2, 1) = da [ 2 ] * b [ 0 ];
    dN.at(8, 1) = da [ 0 ] * b [ 1 ];
    dN.at(9, 1) = da [ 1 ] * b [ 1 ];
    dN.at(6, 1) = da [ 2 ] * b [ 1 ];
    dN.at(4, 1) = da [ 0 ] * b [ 2 ];
    dN.at(7, 1) = da [ 1 ] * b [ 2 ];
    dN.at(3, 1) = da [ 2 ] * b [ 2 ];

    dN.at(1, 2) = a [ 0 ] * db [ 0 ];
    dN.at(5, 2) = a [ 1 ] * db [ 0 ];
    dN.at(2, 2) = a [ 2 ] * db [ 0 ];
    dN.at(8, 2) = a [ 0 ] * db [ 1 ];
    dN.at(9, 2) = a [ 1 ] * db [ 1 ];
    dN.at(6, 2) = a [ 2 ] * db [ 1 ];
    dN.at(4, 2) = a [ 0 ] * db [ 2 ];
    dN.at(7, 2) = a [ 1 ] * db [ 2 ];
    dN.at(3, 2) = a [ 2 ] * db [ 2 ];
}


IntegrationRule *FEI2dQuadBiQuad :: giveIntegrationRule(int order)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Square, order + 6);
    iRule->SetUpPointsOnSquare(points, _Unknown);
    return iRule;
}
} // end namespace oofem
