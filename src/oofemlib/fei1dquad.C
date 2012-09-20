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

#include "fei1dquad.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {
void
FEI1dQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer.resize(3);
    answer.zero();

    answer.at(1) = ksi* ( ksi - 1. ) * 0.5;
    answer.at(2) = ksi* ( 1. + ksi ) * 0.5;
    answer.at(3) = (1. - ksi * ksi);
}

void
FEI1dQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double J = this->giveTransformationJacobian(lcoords,cellgeo);
    double ksi = lcoords.at(1);
    answer.resize(1, 3);
    answer.zero();

    answer.at(1, 1) = (-1./2. + ksi)/J;
    answer.at(1, 2) =  (1./2. + ksi)/J;
    answer.at(1, 3) =  -2.*ksi/J;
}

void
FEI1dQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n(3);
    answer.resize(1);

    this->evalN(n, lcoords, cellgeo);
    answer.at(1) = ( n.at(1) * cellgeo.giveVertexCoordinates(1)->at(cindx) +
            n.at(2) * cellgeo.giveVertexCoordinates(2)->at(cindx) + n.at(3) * cellgeo.giveVertexCoordinates(3)->at(cindx) );
}

int
FEI1dQuad :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    double x1, x2, x3;
    double a,b,c;

    x1 = cellgeo.giveVertexCoordinates(1)->at(cindx);
    x2 = cellgeo.giveVertexCoordinates(2)->at(cindx);
    x3 = cellgeo.giveVertexCoordinates(3)->at(cindx);

    a = 0.5*(x1+x2)-x3;
    b = 0.5*(x2-x1);
    c = x3 - coords.at(1);

    answer.resize(1);
    if(fabs(a)<10.e-6)
    {
        double ksi = ( 2.0 * coords.at(1) - ( x1 + x2 ) ) / ( x2 - x1 );
        answer.at(1) = clamp(ksi, -1., 1.);
        return fabs(ksi) <= 1.0;
    }
    else
    {
        double ksi1 = (-b+sqrt(b*b-4.*a*c))/(2.*a);
        double ksi2 = (-b-sqrt(b*b-4.*a*c))/(2.*a);

        if( (fabs(ksi1) <= 1.) && (fabs(ksi2) <= 1.) ) { // Two roots, element must be bad
            answer.at(1) = 0.;
            return 0;
        } else if ( fabs(ksi1) <= 1. ) {
            answer.at(1) = ksi1;
            return 1;
        } else if( fabs(ksi2) <= 1. ) {
            answer.at(1) = ksi2;
            return 1;
        } else {
            answer.at(1) = 0.;
            return 0;
        }
    }
}

double
FEI1dQuad :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{

    double x1, x2, x3, J, ksi;


    x1 = cellgeo.giveVertexCoordinates(1)->at(cindx);
    x2 = cellgeo.giveVertexCoordinates(2)->at(cindx);
    x3 = cellgeo.giveVertexCoordinates(3)->at(cindx);
    ksi = lcoords.at(1);

    J = 1./2.*(2*ksi-1)*x1 + 1./2.*(2*ksi+1)*x2 - 2.*ksi*x3;
    return J;

}

double
FEI1dQuad :: giveLength(const FEICellGeometry &cellgeo) const
{
    return fabs( cellgeo.giveVertexCoordinates(2)->at(cindx) - cellgeo.giveVertexCoordinates(1)->at(cindx) );
}
} // end namespace oofem
