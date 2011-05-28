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

#include "fei2dlinequad.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"

namespace oofem {
void FEI2dLineQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time)
{
    double xi = lcoords(0);
    answer.resize(3);
    answer(0) = 0.5*(xi-1.0)*xi;
    answer(1) = 0.5*(xi+1.0)*xi;
    answer(2) = 1.0-xi*xi;
}

void FEI2dLineQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time)
{
    // Not meaningful to return anything.
    answer.resize(0,0);
}

void FEI2dLineQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time)
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo, time);
    answer.resize(2);
    answer.at(1) = ( n(0) * cellgeo.giveVertexCoordinates(1)->at(xind) +
                     n(1) * cellgeo.giveVertexCoordinates(2)->at(xind) +
                     n(2) * cellgeo.giveVertexCoordinates(3)->at(xind) );
    answer.at(2) = ( n(0) * cellgeo.giveVertexCoordinates(1)->at(yind) +
                     n(1) * cellgeo.giveVertexCoordinates(2)->at(yind) +
                     n(2) * cellgeo.giveVertexCoordinates(3)->at(yind) );
}

int FEI2dLineQuad :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo, double time)
{
    double x1_x2, y1_y2, px_x3, py_y3, x3_x2_x1, y3_y2_y1;
    double b0, b1, b2, b3;

    x1_x2 = cellgeo.giveVertexCoordinates(1)->at(xind) - cellgeo.giveVertexCoordinates(2)->at(xind);
    y1_y2 = cellgeo.giveVertexCoordinates(1)->at(yind) - cellgeo.giveVertexCoordinates(2)->at(yind);
    px_x3 = gcoords.at(1) - cellgeo.giveVertexCoordinates(3)->at(xind);
    py_y3 = gcoords.at(2) - cellgeo.giveVertexCoordinates(3)->at(yind);
    x3_x2_x1 = 2*cellgeo.giveVertexCoordinates(3)->at(xind) - cellgeo.giveVertexCoordinates(2)->at(xind) - cellgeo.giveVertexCoordinates(1)->at(xind);
    y3_y2_y1 = 2*cellgeo.giveVertexCoordinates(3)->at(yind) - cellgeo.giveVertexCoordinates(2)->at(yind) - cellgeo.giveVertexCoordinates(1)->at(yind);

    b0 = 0.50*(x1_x2*px_x3 + y1_y2*py_y3);
    b1 = 0.25*(x1_x2*x1_x2 + y1_y2*y1_y2) + x3_x2_x1*px_x3 + y3_y2_y1*py_y3;
    b2 = 0.75*(x1_x2*x3_x2_x1 + y1_y2*y3_y2_y1);
    b3 = 0.50*(x3_x2_x1*x3_x2_x1 + y3_y2_y1*y3_y2_y1);

    double r[3];
    int roots;
    cubic(b3, b2, b1, b0, &r[0], &r[1], &r[2], &roots);

    // Copy them over, along with boundary cases.
    int points = 2;
    double p[5] = {-1.0, 1.0};
    for (int i = 0; i < roots; i++) {
        if (r[i] > -1.0 && r[i] < 1.0) {
            // The cubic solver has pretty bad accuracy, performing a single newton iteration
            // which typically improves the solution by many many orders of magnitude.
            // TODO: Move this to the cubic solver.
            r[i] -= (b0 + b1*r[i] + b2*r[i]*r[i] + b3*r[i]*r[i]*r[i])/(b1 + 2*b2*r[i] + 3*b3*r[i]*r[i]);
            p[points] = r[i];
            points++;
        }
    }

    double min_distance2 = 0.0, min_xi = 0, distance2;
    FloatArray f(2);
    answer.resize(1);

    for (int i = 0; i < points; i++) {
        answer(0) = p[i];
        this->local2global(f, answer, cellgeo, 0.0);
        distance2 = f.distance_square(gcoords);
        if ( i == 0 || distance2 < min_distance2 ) {
            min_distance2 = distance2;
            min_xi = answer(0);
        }
    }

    answer(0) = min_xi;
    return true;
}

double FEI2dLineQuad :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time)
{
    OOFEM_ERROR("FEI2DLineQuad :: edgeComputeLength - Not implemented");
    return 0.0;
}

double FEI2dLineQuad :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI2DLineQuad :: edgeComputeLength - Not implemented");
    return 0.0;
}

} // end namespace oofem
