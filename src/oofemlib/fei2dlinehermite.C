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

#include "fei2dlinehermite.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {

double FEI2dLineHermite :: giveLength(const FEICellGeometry &cellgeo) const
{
    double x2_x1, y2_y1;
    x2_x1 = cellgeo.giveVertexCoordinates(2)->at(xind) - cellgeo.giveVertexCoordinates(1)->at(xind);
    y2_y1 = cellgeo.giveVertexCoordinates(2)->at(yind) - cellgeo.giveVertexCoordinates(1)->at(yind);
    return sqrt(x2_x1*x2_x1 + y2_y1*y2_y1);
}

void FEI2dLineHermite :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    double l = this->giveLength(cellgeo);

    answer.resize(4);
    answer.zero();

    answer.at(1) = 0.25 * (1.0 - ksi) * ( 1.0 - ksi) * (2.0 + ksi);
    answer.at(2) =  0.125 * l * (1.0 - ksi) * (1.0 - ksi) * (1.0 + ksi);
    answer.at(3) = 0.25 * (1.0 + ksi) * ( 1.0 + ksi) * (2.0 - ksi);
    answer.at(4) = -0.125 * l * (1.0 + ksi) * (1.0 + ksi) * (1.0 - ksi);
}

void FEI2dLineHermite :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    // This is dNds projected on the direction if the linear edge. If any other interpolation is used for geometry, this can't be used anymore.
    FloatArray dNds;
    this->edgeEvaldNds(dNds, 1, lcoords, cellgeo);
    // Tangent line to project on
    FloatArray vec(2);
    vec.at(1) = cellgeo.giveVertexCoordinates(2)->at(xind) - cellgeo.giveVertexCoordinates(1)->at(xind);
    vec.at(2) = cellgeo.giveVertexCoordinates(2)->at(yind) - cellgeo.giveVertexCoordinates(1)->at(yind);
    vec.normalize();

    answer.beDyadicProductOf(dNds, vec);
}

void FEI2dLineHermite :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n;
    this->evalN(n, lcoords, cellgeo);
    answer.resize(max(xind,yind));
    answer.zero();
    answer.at(xind) = ( n.at(1) * cellgeo.giveVertexCoordinates(1)->at(xind) +
                        n.at(2) * cellgeo.giveVertexCoordinates(2)->at(xind) );
    answer.at(yind) = ( n.at(1) * cellgeo.giveVertexCoordinates(1)->at(yind) +
                        n.at(2) * cellgeo.giveVertexCoordinates(2)->at(yind) );
}

int FEI2dLineHermite :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    double xi;
    double x2_x1, y2_y1;

    x2_x1 = cellgeo.giveVertexCoordinates(2)->at(xind) - cellgeo.giveVertexCoordinates(1)->at(xind);
    y2_y1 = cellgeo.giveVertexCoordinates(2)->at(yind) - cellgeo.giveVertexCoordinates(1)->at(yind);

    // Projection of the global coordinate gives the value interpolated in [0,1].
    xi = (x2_x1*gcoords(0) + y2_y1*gcoords(1))/(sqrt(x2_x1*x2_x1 + y2_y1*y2_y1));
    // Map to [-1,1] domain.
    xi = xi*2.0 - 1.0;

    answer.resize(1);
    answer(0) = clamp(xi, -1., 1.);
    return false;
}

void FEI2dLineHermite :: edgeEvaldNds(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l_inv = this->giveLength(cellgeo);
    double ksi = lcoords.at(1);

    answer.resize(4);
    answer.zero();

    answer.at(1) =  1.5 * (ksi * ksi - 1.0) * l_inv;
    answer.at(2) =  0.25 * (ksi - 1.0) * (3.0 * ksi + 1.0);
    answer.at(3) = -1.5 * (ksi * ksi - 1.0) * l_inv;
    answer.at(4) =  0.25 * (ksi + 1.0) * (3.0 * ksi - 1.0);
}

void FEI2dLineHermite :: edgeEvald2Nds2(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l_inv = this->giveLength(cellgeo);
    double ksi = lcoords.at(1);

    answer.resize(4);
    answer.zero();

    answer.at(1) =  l_inv * 6.0 * ksi * l_inv;
    answer.at(2) =  l_inv * (3.0 * ksi - 1.0);
    answer.at(3) = -l_inv * 6.0 * ksi * l_inv;
    answer.at(4) =  l_inv * (3.0 * ksi + 1.0);
}

void FEI2dLineHermite :: edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double dx, dy;

    dx = cellgeo.giveVertexCoordinates(2)->at(xind) - cellgeo.giveVertexCoordinates(1)->at(xind);
    dy = cellgeo.giveVertexCoordinates(2)->at(yind) - cellgeo.giveVertexCoordinates(1)->at(yind);

    double L = sqrt(dx*dx + dy*dy);

    normal.resize(2);
    normal.at(1) = dy/L;
    normal.at(2) = -dx/L;
}

double FEI2dLineHermite :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x2_x1, y2_y1;
    x2_x1 = cellgeo.giveVertexCoordinates(2)->at(xind) - cellgeo.giveVertexCoordinates(1)->at(xind);
    y2_y1 = cellgeo.giveVertexCoordinates(2)->at(yind) - cellgeo.giveVertexCoordinates(1)->at(yind);
    return sqrt(x2_x1*x2_x1 + y2_y1*y2_y1)*0.5;
}

} // end namespace oofem
