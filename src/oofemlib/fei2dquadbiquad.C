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

#include "fei2dquadbiquad.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {

void
FEI2dQuadBiQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi, eta, n9;

    answer.resize(9);

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    n9 = ( 1. - ksi * ksi ) * ( 1. - eta * eta );

    answer.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25 * ( ksi + eta - 1. ) - 0.25*n9;
    answer.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25 * ( -ksi + eta - 1. ) - 0.25*n9;
    answer.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25 * ( -ksi - eta - 1. ) - 0.25*n9;
    answer.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25 * ( ksi - eta - 1. ) - 0.25*n9;
    answer.at(5) = 0.5 * ( 1. - ksi * ksi ) * ( 1. + eta ) - 0.5*n9;
    answer.at(6) = 0.5 * ( 1. - ksi ) * ( 1. - eta * eta ) - 0.5*n9;
    answer.at(7) = 0.5 * ( 1. - ksi * ksi ) * ( 1. - eta ) - 0.5*n9;
    answer.at(8) = 0.5 * ( 1. + ksi ) * ( 1. - eta * eta ) - 0.5*n9;
    answer.at(9) = n9;
}


void
FEI2dQuadBiQuad :: giveDerivatives(FloatMatrix &dn, const FloatArray &lc)
{
    double ksi, eta, dn9dxi, dn9deta;
    ksi = lc.at(1);
    eta = lc.at(2);
    dn.resize(9, 2);
    
    // dn/dxi
    dn9dxi = 2.0 * ksi * ( 1. - eta * eta );
    dn.at(1,1) =  0.25 * ( 1. + eta ) * ( 2.0 * ksi + eta ) - 0.25 * dn9dxi;
    dn.at(2,1) = -0.25 * ( 1. + eta ) * ( -2.0 * ksi + eta ) - 0.25 * dn9dxi;
    dn.at(3,1) = -0.25 * ( 1. - eta ) * ( -2.0 * ksi - eta ) - 0.25 * dn9dxi;
    dn.at(4,1) =  0.25 * ( 1. - eta ) * ( 2.0 * ksi - eta ) - 0.25 * dn9dxi;
    dn.at(5,1) = -ksi * ( 1. + eta ) - 0.5 * dn9dxi;
    dn.at(6,1) = -0.5 * ( 1. - eta * eta ) - 0.5 * dn9dxi;
    dn.at(7,1) = -ksi * ( 1. - eta ) - 0.5 * dn9dxi;
    dn.at(8,1) =  0.5 * ( 1. - eta * eta ) - 0.5 * dn9dxi;
    dn.at(9,1) =  dn9dxi;

    // dn/deta
    dn9deta = ( 1. - ksi * ksi ) * 2.0 * eta;
    dn.at(1,2) =  0.25 * ( 1. + ksi ) * ( 2.0 * eta + ksi ) - 0.25 * dn9deta;
    dn.at(2,2) =  0.25 * ( 1. - ksi ) * ( 2.0 * eta - ksi ) - 0.25 * dn9deta;
    dn.at(3,2) = -0.25 * ( 1. - ksi ) * ( -2.0 * eta - ksi ) - 0.25 * dn9deta;
    dn.at(4,2) = -0.25 * ( 1. + ksi ) * ( -2.0 * eta + ksi ) - 0.25 * dn9deta;
    dn.at(5,2) =  0.5 * ( 1. - ksi * ksi ) - 0.5 * dn9deta;
    dn.at(6,2) = -eta * ( 1. - ksi ) - 0.5 * dn9deta;
    dn.at(7,2) = -0.5 * ( 1. - ksi * ksi ) - 0.5 * dn9deta;
    dn.at(8,2) = -eta * ( 1. + ksi ) - 0.5 * dn9deta;
    dn.at(9,2) =  dn9deta;
}

} // end namespace oofem
