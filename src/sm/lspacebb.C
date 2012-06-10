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

#include "lspacebb.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"

#ifdef __OOFEG
 #include "engngm.h"
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
#endif

namespace oofem {
LSpaceBB :: LSpaceBB(int n, Domain *aDomain) : LSpace(n, aDomain)
{ }

void
LSpaceBB :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [6x24] strain-displacement matrix {B} of the receiver, eva-
// luated at aGaussPoint.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    int i;
    FloatMatrix dnx, dnx0;
    FloatArray coord(3);

    answer.resize(6, 24);
    answer.zero();
    coord.zero();


    LSpace :: interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));
    LSpace :: interpolation.evaldNdx(dnx0, coord, FEIElementGeometryWrapper(this));

    // deviatoric part fully integrated, volumetric part in one point
    // here we follow BBar approach
    //
    // construct BB = B(gp) + Pv [B(0)-B(gp)]
    // where Pv is volumetric projection mtrx
    // B(gp) is original geometrical matrix evalueated at gp
    // B(0)  is geometrical matrix evalueated at centroid
    //
    // assemble Pv [B(0)-B(gp)]
    for ( i = 1; i <= 8; i++ ) {
        answer.at(1, 3 * i - 2) = answer.at(2, 3 * i - 2) = answer.at(3, 3 * i - 2) = ( dnx0.at(i, 1) - dnx.at(i, 1) ) / 3.0;
        answer.at(1, 3 * i - 1) = answer.at(2, 3 * i - 1) = answer.at(3, 3 * i - 1) = ( dnx0.at(i, 2) - dnx.at(i, 2) ) / 3.0;
        answer.at(1, 3 * i - 0) = answer.at(2, 3 * i - 0) = answer.at(3, 3 * i - 0) = ( dnx0.at(i, 3) - dnx.at(i, 3) ) / 3.0;
    }

    // add B(gp)
    for ( i = 1; i <= 8; i++ ) {
        answer.at(1, 3 * i - 2) += dnx.at(i, 1);
        answer.at(2, 3 * i - 1) += dnx.at(i, 2);
        answer.at(3, 3 * i - 0) += dnx.at(i, 3);
    }

    for ( i = 1; i <= 8; i++ ) {
        answer.at(4, 3 * i - 1) += dnx.at(i, 3);
        answer.at(4, 3 * i - 0) += dnx.at(i, 2);

        answer.at(5, 3 * i - 2) += dnx.at(i, 3);
        answer.at(5, 3 * i - 0) += dnx.at(i, 1);

        answer.at(6, 3 * i - 2) += dnx.at(i, 2);
        answer.at(6, 3 * i - 1) += dnx.at(i, 1);
    }
}
} // end namespace oofem
