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

#include "fei1dlin.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {

double
FEI1dLin :: giveLength(const FEICellGeometry &cellgeo) const
{
    return fabs( cellgeo.giveVertexCoordinates(1)->at(cindx) - cellgeo.giveVertexCoordinates(2)->at(cindx) );
}

void
FEI1dLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

void
FEI1dLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l = this->computeLength(cellgeo);
    answer.resize(2, 1);

    answer.at(1, 1) = -1.0 / l;
    answer.at(2, 1) =  1.0 / l;
}

void
FEI1dLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n(2);
    answer.resize(1);

    this->evalN(n, lcoords, cellgeo);
    answer.at(1) = ( n.at(1) * cellgeo.giveVertexCoordinates(1)->at(cindx) +
                    n.at(2) * cellgeo.giveVertexCoordinates(2)->at(cindx) );
}

int
FEI1dLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    double ksi, x1, x2;
    answer.resize(1);


    x1 = cellgeo.giveVertexCoordinates(1)->at(cindx);
    x2 = cellgeo.giveVertexCoordinates(2)->at(cindx);

    answer.at(1) = ksi = ( 2.0 * coords.at(1) - ( x1 + x2 ) ) / ( x2 - x1 );
    return ( fabs(ksi) <= 1.0 ) ? 1 : 0;
}

double
FEI1dLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l = computeLength(cellgeo);
    return 0.5 * l;
}

double
FEI1dLin :: computeLength(const FEICellGeometry &cellgeo)
{
    return ( cellgeo.giveVertexCoordinates(2)->at(cindx) - cellgeo.giveVertexCoordinates(1)->at(cindx) );
}
} // end namespace oofem
