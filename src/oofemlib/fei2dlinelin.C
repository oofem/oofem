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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#include "fei2dlinelin.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "mathfem.h"

namespace oofem {
void FEI2dLineLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time)
{
    double xi = lcoords(0);
    answer.resize(2);
    answer.at(1) = ( 1. - xi ) * 0.5;
    answer.at(2) = ( 1. + xi ) * 0.5;
}

void FEI2dLineLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time)
{
    // Not meaningful to return anything.
    answer.resize(0,0);
}

void FEI2dLineLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time)
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo, time);
    answer.resize(max(xind,yind));
    answer.zero();
    answer.at(xind) = ( n(0) * cellgeo.giveVertexCoordinates(1)->at(xind) +
                        n(1) * cellgeo.giveVertexCoordinates(2)->at(xind) );
    answer.at(yind) = ( n(0) * cellgeo.giveVertexCoordinates(1)->at(yind) +
                        n(1) * cellgeo.giveVertexCoordinates(2)->at(yind) );
}

int FEI2dLineLin :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo, double time)
{
    double xi;
    double x2_x1, y2_y1;

    x2_x1 = cellgeo.giveVertexCoordinates(2)->at(xind) - cellgeo.giveVertexCoordinates(1)->at(xind);
    y2_y1 = cellgeo.giveVertexCoordinates(2)->at(yind) - cellgeo.giveVertexCoordinates(1)->at(yind);

    // Projection of the global coordinate gives the value interpolated in [0,1].
    xi = (x2_x1*gcoords(0) + y2_y1*gcoords(1))/(sqrt(x2_x1*x2_x1 + y2_y1*y2_y1));
    // Map to [-1,1] domain.
    xi = xi*2 - 1;

    answer.resize(1);
    answer(0) = xi;
    return true;
}

double FEI2dLineLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time)
{
    double x2_x1, y2_y1;
    x2_x1 = cellgeo.giveVertexCoordinates(2)->at(xind) - cellgeo.giveVertexCoordinates(1)->at(xind);
    y2_y1 = cellgeo.giveVertexCoordinates(2)->at(yind) - cellgeo.giveVertexCoordinates(1)->at(yind);
    return sqrt(x2_x1*x2_x1 + y2_y1*y2_y1)/2.0;
}

double FEI2dLineLin :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    double x2_x1, y2_y1;
    x2_x1 = cellgeo.giveVertexCoordinates(2)->at(xind) - cellgeo.giveVertexCoordinates(1)->at(xind);
    y2_y1 = cellgeo.giveVertexCoordinates(2)->at(yind) - cellgeo.giveVertexCoordinates(1)->at(yind);
    return sqrt(x2_x1*x2_x1 + y2_y1*y2_y1);
}

} // end namespace oofem
