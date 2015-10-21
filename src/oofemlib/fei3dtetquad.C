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

#include "fei3dtetquad.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
double
FEI3dTetQuad :: giveVolume(const FEICellGeometry &cellgeo) const
{
    // Use linear approximation.

    double x1, x2, x3, x4;     //, x5, x6, x7, x8, x9, x10;
    double y1, y2, y3, y4;     //, y5, y6, y7, y8, y9, y10;
    double z1, z2, z3, z4;     //, z5, z6, z7, z8, z9, z10;

    x1 = cellgeo.giveVertexCoordinates(1)->at(1);
    y1 = cellgeo.giveVertexCoordinates(1)->at(2);
    z1 = cellgeo.giveVertexCoordinates(1)->at(3);
    x2 = cellgeo.giveVertexCoordinates(2)->at(1);
    y2 = cellgeo.giveVertexCoordinates(2)->at(2);
    z2 = cellgeo.giveVertexCoordinates(2)->at(3);
    x3 = cellgeo.giveVertexCoordinates(3)->at(1);
    y3 = cellgeo.giveVertexCoordinates(3)->at(2);
    z3 = cellgeo.giveVertexCoordinates(3)->at(3);
    x4 = cellgeo.giveVertexCoordinates(4)->at(1);
    y4 = cellgeo.giveVertexCoordinates(4)->at(2);
    z4 = cellgeo.giveVertexCoordinates(4)->at(3);
#if 0
    x5=cellgeo.giveVertexCoordinates(5)->at(1);     y5=cellgeo.giveVertexCoordinates(5)->at(2);     z5=cellgeo.giveVertexCoordinates(5)->at(3);
    x6=cellgeo.giveVertexCoordinates(6)->at(1);     y6=cellgeo.giveVertexCoordinates(6)->at(2);     z6=cellgeo.giveVertexCoordinates(6)->at(3);
    x7=cellgeo.giveVertexCoordinates(7)->at(1);     y7=cellgeo.giveVertexCoordinates(7)->at(2);     z7=cellgeo.giveVertexCoordinates(7)->at(3);
    x8=cellgeo.giveVertexCoordinates(8)->at(1);     y8=cellgeo.giveVertexCoordinates(8)->at(2);     z8=cellgeo.giveVertexCoordinates(8)->at(3);
    x9=cellgeo.giveVertexCoordinates(9)->at(1);     y9=cellgeo.giveVertexCoordinates(9)->at(2);     z9=cellgeo.giveVertexCoordinates(9)->at(3);
    x10=cellgeo.giveVertexCoordinates(10)->at(1);   y10=cellgeo.giveVertexCoordinates(10)->at(2);   z10=cellgeo.giveVertexCoordinates(10)->at(3);

    double area = x1*y3*z2 - x1*y2*z3 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1 + x1*y2*z4 - x1*y4*z2 - x2*y1*z4 + x2*y4*z1 + x4*y1*z2 - x4*y2*z1 - x1*y3*z4 + x1*y4*z3
     + x3*y1*z4 - x3*y4*z1 - x4*y1*z3 + x4*y3*z1 - 2*x1*y2*z6 + 2*x1*y3*z5 - 2*x1*y5*z3 + 2*x1*y6*z2 + 2*x2*y1*z6 + x2*y3*z4 - x2*y4*z3 - 2*x2*y6*z1 - 2*x3*y1*z5 - x3*y2*z4
     + x3*y4*z2 + 2*x3*y5*z1 + x4*y2*z3 - x4*y3*z2 + 2*x5*y1*z3 - 2*x5*y3*z1 - 2*x6*y1*z2 + 2*x6*y2*z1 - 2*x1*y2*z7 + 2*x1*y3*z6 - 2*x1*y4*z5 + 2*x1*y5*z4 - 2*x1*y6*z3 + 2*x1*y7*z2
     + 2*x2*y1*z7 - 2*x2*y3*z5 + 2*x2*y5*z3 - 2*x2*y7*z1 - 2*x3*y1*z6 + 2*x3*y2*z5 - 2*x3*y5*z2 + 2*x3*y6*z1 + 2*x4*y1*z5 - 2*x4*y5*z1 - 2*x5*y1*z4 - 2*x5*y2*z3 + 2*x5*y3*z2
     + 2*x5*y4*z1 + 2*x6*y1*z3 - 2*x6*y3*z1 - 2*x7*y1*z2 + 2*x7*y2*z1 + 2*x1*y2*z8 - 2*x1*y8*z2 - 2*x2*y1*z8 + 2*x2*y4*z5
     + - 2*x2*y5*z4 + 2*x2*y8*z1 - 2*x4*y2*z5 + 2*x4*y5*z2 + 2*x5*y2*z4 - 2*x5*y4*z2 + 2*x8*y1*z2 - 2*x8*y2*z1 + 2*x1*y2*z9 - 2*x1*y3*z8 + 2*x1*y4*z7 - 4*x1*y5*z6 + 4*x1*y6*z5
     + - 2*x1*y7*z4 + 2*x1*y8*z3 - 2*x1*y9*z2 - 2*x2*y1*z9 - 2*x2*y3*z7 - 2*x2*y4*z6 + 2*x2*y6*z4 + 2*x2*y7*z3 + 2*x2*y9*z1 + 2*x3*y1*z8 + 2*x3*y2*z7 - 2*x3*y7*z2 - 2*x3*y8*z1
     + - 2*x4*y1*z7 + 2*x4*y2*z6 - 2*x4*y6*z2 + 2*x4*y7*z1 + 4*x5*y1*z6 - 4*x5*y6*z1 - 4*x6*y1*z5 - 2*x6*y2*z4 + 2*x6*y4*z2 + 4*x6*y5*z1 + 2*x7*y1*z4 - 2*x7*y2*z3 + 2*x7*y3*z2
     + - 2*x7*y4*z1 - 2*x8*y1*z3 + 2*x8*y3*z1 + 2*x9*y1*z2 - 2*x9*y2*z1 - 4*x1*y5*z7 + 4*x1*y7*z5 + 4*x2*y5*z6 - 4*x2*y6*z5 + 2*x3*y4*z6 - 2*x3*y6*z4 - 2*x4*y3*z6 + 2*x4*y6*z3
     + 4*x5*y1*z7 - 4*x5*y2*z6 + 4*x5*y6*z2 - 4*x5*y7*z1 + 4*x6*y2*z5 + 2*x6*y3*z4 - 2*x6*y4*z3 - 4*x6*y5*z2 - 4*x7*y1*z5 + 4*x7*y5*z1 - 2*x1*y3*z10 - 2*x1*y4*z9 + 4*x1*y5*z8
     + - 4*x1*y6*z7 + 4*x1*y7*z6 - 4*x1*y8*z5 + 2*x1*y9*z4 + 2*x1*y10*z3 + 2*x2*y3*z9 + 2*x2*y4*z8 + 4*x2*y5*z7 - 4*x2*y7*z5 - 2*x2*y8*z4 - 2*x2*y9*z3 + 2*x3*y1*z10 - 2*x3*y2*z9
     + - 2*x3*y4*z7 - 4*x3*y5*z6 + 4*x3*y6*z5 + 2*x3*y7*z4 + 2*x3*y9*z2 - 2*x3*y10*z1 + 2*x4*y1*z9
     + - 2*x4*y2*z8 + 2*x4*y3*z7 - 2*x4*y7*z3 + 2*x4*y8*z2 - 2*x4*y9*z1 - 4*x5*y1*z8 - 4*x5*y2*z7 + 4*x5*y3*z6 - 4*x5*y6*z3 + 4*x5*y7*z2 + 4*x5*y8*z1 + 4*x6*y1*z7 - 4*x6*y3*z5
     + 4*x6*y5*z3 - 4*x6*y7*z1 - 4*x7*y1*z6 + 4*x7*y2*z5 - 2*x7*y3*z4 + 2*x7*y4*z3 - 4*x7*y5*z2 + 4*x7*y6*z1 + 4*x8*y1*z5 + 2*x8*y2*z4 - 2*x8*y4*z2 - 4*x8*y5*z1 - 2*x9*y1*z4
     + 2*x9*y2*z3 - 2*x9*y3*z2 + 2*x9*y4*z1 - 2*x10*y1*z3 + 2*x10*y3*z1 + 2*x1*y4*z10 + 4*x1*y5*z9 - 4*x1*y9*z5 - 2*x1*y10*z4 + 2*x2*y3*z10 - 4*x2*y5*z8 - 4*x2*y6*z7 + 4*x2*y7*z6
     + 4*x2*y8*z5 - 2*x2*y10*z3 - 2*x3*y2*z10 - 2*x3*y4*z8 + 4*x3*y5*z7 - 4*x3*y7*z5 + 2*x3*y8*z4 + 2*x3*y10*z2 - 2*x4*y1*z10 + 2*x4*y3*z8 + 4*x4*y5*z6 - 4*x4*y6*z5 - 2*x4*y8*z3
     + 2*x4*y10*z1 - 4*x5*y1*z9 + 4*x5*y2*z8 - 4*x5*y3*z7 - 4*x5*y4*z6 + 4*x5*y6*z4 + 4*x5*y7*z3 - 4*x5*y8*z2 + 4*x5*y9*z1 + 4*x6*y2*z7 + 4*x6*y4*z5 - 4*x6*y5*z4 - 4*x6*y7*z2
     + - 4*x7*y2*z6 + 4*x7*y3*z5 - 4*x7*y5*z3 + 4*x7*y6*z2 - 4*x8*y2*z5 - 2*x8*y3*z4 + 2*x8*y4*z3 + 4*x8*y5*z2 + 4*x9*y1*z5 - 4*x9*y5*z1 + 2*x10*y1*z4 + 2*x10*y2*z3 - 2*x10*y3*z2 -
     + 2*x10*y4*z1 + 4*x1*y6*z9 - 4*x1*y7*z8 + 4*x1*y8*z7 - 4*x1*y9*z6
     + - 2*x2*y4*z10 - 4*x2*y5*z9 + 4*x2*y9*z5 + 2*x2*y10*z4 + 2*x3*y4*z9 + 4*x3*y5*z8 + 4*x3*y6*z7 - 4*x3*y7*z6 - 4*x3*y8*z5 - 2*x3*y9*z4 + 2*x4*y2*z10 - 2*x4*y3*z9 - 4*x4*y5*z7
     + 4*x4*y7*z5 + 2*x4*y9*z3 - 2*x4*y10*z2 + 4*x5*y2*z9 - 4*x5*y3*z8 + 4*x5*y4*z7 - 4*x5*y7*z4 + 4*x5*y8*z3 - 4*x5*y9*z2 - 4*x6*y1*z9 - 4*x6*y3*z7 + 4*x6*y7*z3 + 4*x6*y9*z1
     + 4*x7*y1*z8 + 4*x7*y3*z6 - 4*x7*y4*z5 + 4*x7*y5*z4 - 4*x7*y6*z3 - 4*x7*y8*z1 - 4*x8*y1*z7 + 4*x8*y3*z5 - 4*x8*y5*z3 + 4*x8*y7*z1 + 4*x9*y1*z6
     + - 4*x9*y2*z5 + 2*x9*y3*z4 - 2*x9*y4*z3 + 4*x9*y5*z2 - 4*x9*y6*z1 - 2*x10*y2*z4 + 2*x10*y4*z2 - 4*x1*y6*z10 + 4*x1*y10*z6 + 4*x2*y6*z9 - 4*x2*y7*z8 + 4*x2*y8*z7 - 4*x2*y9*z6
     + - 4*x3*y5*z9 + 4*x3*y9*z5 - 4*x4*y5*z8 + 4*x4*y6*z7 - 4*x4*y7*z6 + 4*x4*y8*z5 + 4*x5*y3*z9 + 4*x5*y4*z8 - 4*x5*y8*z4 - 4*x5*y9*z3 + 4*x6*y1*z10 - 4*x6*y2*z9 - 4*x6*y4*z7
     + 4*x6*y7*z4 + 4*x6*y9*z2 - 4*x6*y10*z1 + 4*x7*y2*z8 + 4*x7*y4*z6 - 4*x7*y6*z4 - 4*x7*y8*z2 - 4*x8*y2*z7 - 4*x8*y4*z5 + 4*x8*y5*z4 + 4*x8*y7*z2
     + 4*x9*y2*z6 - 4*x9*y3*z5 + 4*x9*y5*z3 - 4*x9*y6*z2 - 4*x10*y1*z6 + 4*x10*y6*z1 - 4*x1*y7*z10 - 4*x1*y8*z9 + 4*x1*y9*z8 + 4*x1*y10*z7 + 4*x2*y6*z10 - 4*x2*y10*z6 - 4*x3*y6*z9
     + 4*x3*y7*z8 - 4*x3*y8*z7 + 4*x3*y9*z6 + 4*x4*y5*z9 - 4*x4*y9*z5 - 4*x5*y4*z9 - 16*x5*y6*z7 + 16*x5*y7*z6 + 4*x5*y9*z4 - 4*x6*y2*z10 + 4*x6*y3*z9 + 16*x6*y5*z7 - 16*x6*y7*z5
     + - 4*x6*y9*z3 + 4*x6*y10*z2 + 4*x7*y1*z10 - 4*x7*y3*z8 - 16*x7*y5*z6 + 16*x7*y6*z5 + 4*x7*y8*z3 - 4*x7*y10*z1 + 4*x8*y1*z9 + 4*x8*y3*z7 - 4*x8*y7*z3 - 4*x8*y9*z1 - 4*x9*y1*z8
     + - 4*x9*y3*z6 + 4*x9*y4*z5 - 4*x9*y5*z4 + 4*x9*y6*z3 + 4*x9*y8*z1 - 4*x10*y1*z7 + 4*x10*y2*z6 - 4*x10*y6*z2 + 4*x10*y7*z1 + 4*x1*y8*z10 - 4*x1*y10*z8 + 4*x2*y7*z10 - 4*x2*y8*z9
     + 4*x2*y9*z8 - 4*x2*y10*z7 - 4*x3*y6*z10 + 4*x3*y10*z6 - 4*x4*y6*z9 + 4*x4*y7*z8 - 4*x4*y8*z7 + 4*x4*y9*z6 + 4*x6*y3*z10 + 4*x6*y4*z9 - 4*x6*y9*z4 - 4*x6*y10*z3 - 4*x7*y2*z10
     + - 4*x7*y4*z8 + 4*x7*y8*z4 + 4*x7*y10*z2 - 4*x8*y1*z10 + 4*x8*y2*z9 + 4*x8*y4*z7 - 4*x8*y7*z4 - 4*x8*y9*z2 + 4*x8*y10*z1 - 4*x9*y2*z8 - 4*x9*y4*z6 + 4*x9*y6*z4 + 4*x9*y8*z2
     + 4*x10*y1*z8 + 4*x10*y2*z7 - 4*x10*y3*z6 + 4*x10*y6*z3 - 4*x10*y7*z2 - 4*x10*y8*z1 + 4*x1*y9*z10 - 4*x1*y10*z9 - 4*x2*y8*z10 + 4*x2*y10*z8 + 4*x3*y7*z10 + 4*x3*y8*z9 - 4*x3*y9*z8
     + - 4*x3*y10*z7 + 4*x4*y6*z10 - 4*x4*y10*z6 + 16*x5*y6*z9 - 16*x5*y7*z8 + 16*x5*y8*z7 - 16*x5*y9*z6 - 4*x6*y4*z10 - 16*x6*y5*z9 + 16*x6*y9*z5 + 4*x6*y10*z4 - 4*x7*y3*z10 + 16*x7*y5*z8
     + - 16*x7*y8*z5 + 4*x7*y10*z3 + 4*x8*y2*z10 - 4*x8*y3*z9 - 16*x8*y5*z7 + 16*x8*y7*z5 + 4*x8*y9*z3 - 4*x8*y10*z2 - 4*x9*y1*z10 + 4*x9*y3*z8 + 16*x9*y5*z6 - 16*x9*y6*z5 - 4*x9*y8*z3
     + 4*x9*y10*z1 + 4*x10*y1*z9 - 4*x10*y2*z8 + 4*x10*y3*z7 + 4*x10*y4*z6 - 4*x10*y6*z4 - 4*x10*y7*z3 + 4*x10*y8*z2 - 4*x10*y9*z1 - 4*x2*y9*z10 + 4*x2*y10*z9 + 4*x3*y8*z10 - 4*x3*y10*z8
     + - 4*x4*y7*z10 + 4*x4*y8*z9 - 4*x4*y9*z8 + 4*x4*y10*z7 + 4*x7*y4*z10 - 4*x7*y10*z4 - 4*x8*y3*z10 - 4*x8*y4*z9 + 4*x8*y9*z4 + 4*x8*y10*z3 + 4*x9*y2*z10 + 4*x9*y4*z8 - 4*x9*y8*z4
     + - 4*x9*y10*z2 - 4*x10*y2*z9 + 4*x10*y3*z8 - 4*x10*y4*z7 + 4*x10*y7*z4 - 4*x10*y8*z3 + 4*x10*y9*z2 - 4*x3*y9*z10 + 4*x3*y10*z9 - 4*x4*y8*z10 + 4*x4*y10*z8 - 16*x5*y8*z9 + 16*x5*y9*z8
     + 4*x8*y4*z10 + 16*x8*y5*z9 - 16*x8*y9*z5 - 4*x8*y10*z4 + 4*x9*y3*z10 - 16*x9*y5*z8 + 16*x9*y8*z5 - 4*x9*y10*z3 - 4*x10*y3*z9 - 4*x10*y4*z8 + 4*x10*y8*z4 + 4*x10*y9*z3 + 4*x4*y9*z10
     + - 4*x4*y10*z9 + 16*x6*y7*z10 - 16*x6*y10*z7 - 16*x7*y6*z10 + 16*x7*y10*z6 - 4*x9*y4*z10 + 4*x9*y10*z4 + 4*x10*y4*z9 + 16*x10*y6*z7 - 16*x10*y7*z6 - 4*x10*y9*z4 - 16*x6*y9*z10
     + 16*x6*y10*z9 + 16*x7*y8*z10 - 16*x7*y10*z8 - 16*x8*y7*z10 + 16*x8*y10*z7 + 16*x9*y6*z10 - 16*x9*y10*z6 - 16*x10*y6*z9 + 16*x10*y7*z8 - 16*x10*y8*z7 + 16*x10*y9*z6 + 16*x8*y9*z10
     + - 16*x8*y10*z9 - 16*x9*y8*z10 + 16*x9*y10*z8 + 16*x10*y8*z9 - 16*x10*y9*z8;

     //printf("Q Area=%f\n", area);
#endif
    double area = x1 * y3 * z2 - x1 * y2 * z3 + x2 * y1 * z3 - x2 * y3 * z1 - x3 * y1 * z2 + x3 * y2 * z1 + x1 * y2 * z4 - x1 * y4 * z2 - x2 * y1 * z4 + x2 * y4 * z1 + x4 * y1 * z2 - x4 * y2 * z1 -
                  x1 * y3 * z4 + x1 * y4 * z3 + x3 * y1 * z4 - x3 * y4 * z1 - x4 * y1 * z3 + x4 * y3 * z1 + x2 * y3 * z4 - x2 * y4 * z3 - x3 * y2 * z4 + x3 * y4 * z2 + x4 * y2 * z3 - x4 * y3 * z2;

    //printf("L Area=%f\n", area);

    area = area / 6.0;
    return area;
}

void
FEI3dTetQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x1 = lcoords(0);
    double x2 = lcoords(1);
    double x3 = lcoords(2);
    double x4 = 1.0 - x1 - x2 - x3;

    answer.resize(10);
    answer(0) = x1 * ( 2 * x1 - 1 );
    answer(1) = x2 * ( 2 * x2 - 1 );
    answer(2) = x3 * ( 2 * x3 - 1 );
    answer(3) = x4 * ( 2 * x4 - 1 );

    answer(4) = 4 * x1 * x2;
    answer(5) = 4 * x2 * x3;
    answer(6) = 4 * x3 * x1;
    answer(7) = 4 * x1 * x4;
    answer(8) = 4 * x2 * x4;
    answer(9) = 4 * x3 * x4;
}

double
FEI3dTetQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix, inv, dNduvw, coords;
    this->evaldNdxi(dNduvw, lcoords, cellgeo);
    coords.resize( 3, dNduvw.giveNumberOfRows() );
    for ( int i = 1; i <= dNduvw.giveNumberOfRows(); i++ ) {
        coords.setColumn(* cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
    inv.beInverseOf(jacobianMatrix);

    answer.beProductOf(dNduvw, inv);
    return jacobianMatrix.giveDeterminant();
}

void
FEI3dTetQuad :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords , const FEICellGeometry &cellgeo)
{
    double x1 = lcoords(0);
    double x2 = lcoords(1);
    double x3 = lcoords(2);
    double x4 = 1.0 - x1 - x2 - x3;

    answer.resize(10, 3);

    // dNj/dx1
    answer(0, 0) = 4 * x1 - 1;
    answer(1, 0) = 0;
    answer(2, 0) = 0;
    answer(3, 0) = -4 * x4 + 1;
    answer(4, 0) = 4 * x2;
    answer(5, 0) = 0;
    answer(6, 0) = 4 * x3;
    answer(7, 0) = 4 * ( x4 - x1 );
    answer(8, 0) = -4 * x2;
    answer(9, 0) = -4 * x3;

    // dNj/dx2
    answer(0, 1) = 0;
    answer(1, 1) = 4 * x2 - 1;
    answer(2, 1) = 0;
    answer(3, 1) = -4 * x4 + 1;
    answer(4, 1) = 4 * x1;
    answer(5, 1) = 4 * x3;
    answer(6, 1) = 0;
    answer(7, 1) = -4 * x1;
    answer(8, 1) = 4 * ( x4 - x2 );
    answer(9, 1) = -4 * x3;

    // dNj/dx3
    answer(0, 2) = 0;
    answer(1, 2) = 0;
    answer(2, 2) = 4 * x3 - 1;
    answer(3, 2) = -4 * x4 + 1;
    answer(4, 2) = 0;
    answer(5, 2) = 4 * x2;
    answer(6, 2) = 4 * x1;
    answer(7, 2) = -4 * x1;
    answer(8, 2) = -4 * x2;
    answer(9, 2) = 4 * ( x4 - x3 );
}


void
FEI3dTetQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray N;
    this->evalN(N, lcoords, cellgeo);
    answer.clear();
    for ( int i = 1; i <= N.giveSize(); i++ ) {
        answer.add( N.at(i), * cellgeo.giveVertexCoordinates(i) );
    }
}

#define POINT_TOL 1e-6

int
FEI3dTetQuad :: global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo)
{
    FloatArray res, delta, guess, lcoords_guess;
    FloatMatrix jac;
    double convergence_limit, error = 0.0;

    // find a suitable convergence limit
    convergence_limit = 1e-6 * this->giveCharacteristicLength(cellgeo);

    // setup initial guess
    lcoords_guess.resize( gcoords.giveSize() );
    lcoords_guess.zero();

    // apply Newton-Raphson to solve the problem
    for ( int nite = 0; nite < 10; nite++ ) {
        // compute the residual
        this->local2global(guess, lcoords_guess, cellgeo);
        res.beDifferenceOf(gcoords, guess);

        // check for convergence
        error = res.computeNorm();
        if ( error < convergence_limit ) {
            break;
        }

        // compute the corrections
        this->giveJacobianMatrixAt(jac, lcoords_guess, cellgeo);
        jac.solveForRhs(res, delta, true);

        // update guess
        lcoords_guess.add(delta);
    }
    if ( error > convergence_limit ) { // Imperfect, could give false negatives.
        answer.resize(4);
        answer.zero();
        return false;
    }

    answer.resize(4);
    answer(0) = lcoords_guess(0);
    answer(1) = lcoords_guess(1);
    answer(2) = lcoords_guess(2);

    bool inside = true;
    for ( int i = 0; i < 3; i++ ) {
        if ( answer(i) < ( 0. - POINT_TOL ) ) {
            answer(i) = 0.;
            inside = false;
        } else if ( answer(i) > ( 1. + POINT_TOL ) ) {
            answer(i) = 1.;
            inside = false;
        }
    }

    answer(3) = 1.0 - answer(0) - answer(1) - answer(2); // Do this afterwards, since it might get clamped.
    if ( answer(3) < 0. - POINT_TOL ) {
        return false;
    }
    return inside;
}


double
FEI3dTetQuad :: giveCharacteristicLength(const FEICellGeometry &cellgeo) const
{
    return cellgeo.giveVertexCoordinates(1)->distance( cellgeo.giveVertexCoordinates(2) );
}


void
FEI3dTetQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix dNduvw, coords;
    this->evaldNdxi(dNduvw, lcoords, cellgeo);
    coords.resize( 3, dNduvw.giveNumberOfRows() );
    for ( int i = 1; i <= dNduvw.giveNumberOfRows(); i++ ) {
        coords.setColumn(* cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
}


void
FEI3dTetQuad :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double xi = lcoords.at(1);
    answer.resize(3);
    answer(0) = 0.5 * ( xi - 1.0 ) * xi;
    answer(1) = 0.5 * ( xi + 1.0 ) * xi;
    answer(2) = 1.0 - xi * xi;
}

void
FEI3dTetQuad :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                             const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    ///@todo Implement this
    OOFEM_ERROR("Not supported");
}

void
FEI3dTetQuad :: edgeLocal2global(FloatArray &answer, int iedge,
                                 const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    FloatArray N;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeEvalN(N, iedge, lcoords, cellgeo);

    answer.clear();
    for ( int i = 0; i < N.giveSize(); ++i ) {
        answer.add( N(i), * cellgeo.giveVertexCoordinates( edgeNodes(i) ) );
    }
}


double
FEI3dTetQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray edgeNodes;
    this->computeLocalEdgeMapping(edgeNodes, iedge);
    ///@todo Implement this
    OOFEM_ERROR("Not supported");
    return -1;
}


void
FEI3dTetQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    edgeNodes.resize(3);

    if ( iedge == 1 ) { // edge between nodes 1 2
        edgeNodes(0) = 1;
        edgeNodes(1) = 2;
        edgeNodes(2) = 5;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        edgeNodes(0) = 2;
        edgeNodes(1) = 3;
        edgeNodes(2) = 6;
    } else if ( iedge == 3 ) { // edge between nodes 3 1
        edgeNodes(0) = 3;
        edgeNodes(1) = 1;
        edgeNodes(2) = 7;
    } else if ( iedge == 4 ) { // edge between nodes 1 4
        edgeNodes(0) = 1;
        edgeNodes(1) = 4;
        edgeNodes(2) = 8;
    } else if ( iedge == 5 ) { // edge between nodes 2 4
        edgeNodes(0) = 2;
        edgeNodes(1) = 4;
        edgeNodes(2) = 9;
    } else if ( iedge == 6 ) { // edge between nodes 3 4
        edgeNodes(0) = 3;
        edgeNodes(1) = 4;
        edgeNodes(2) = 10;
    } else {
        OOFEM_ERROR("wrong edge number (%d)", iedge);
    }
}

double
FEI3dTetQuad :: edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo)
{
    ///@todo Implement this
    OOFEM_ERROR("Not supported");
    return -1;
}

void
FEI3dTetQuad :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1. - l1 - l2;

    answer.resize(6);

    answer.at(1) = ( 2. * l1 - 1. ) * l1;
    answer.at(2) = ( 2. * l2 - 1. ) * l2;
    answer.at(3) = ( 2. * l3 - 1. ) * l3;
    answer.at(4) = 4. * l1 * l2;
    answer.at(5) = 4. * l2 * l3;
    answer.at(6) = 4. * l3 * l1;
}

void
FEI3dTetQuad :: surfaceLocal2global(FloatArray &answer, int isurf,
                                    const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray nodes;
    FloatArray N;
    this->computeLocalSurfaceMapping(nodes, isurf);
    this->surfaceEvalN(N, isurf, lcoords, cellgeo);

    answer.clear();
    for ( int i = 0; i < N.giveSize(); ++i ) {
        answer.add( N(i), * cellgeo.giveVertexCoordinates( nodes(i) ) );
    }
}

void
FEI3dTetQuad :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray snodes;
    this->computeLocalSurfaceMapping(snodes, isurf);

    FloatArray lcoords_tet(4);
    lcoords_tet.at(snodes.at(1)) = lcoords.at(1);
    lcoords_tet.at(snodes.at(2)) = lcoords.at(2);
    lcoords_tet.at(snodes.at(3)) = 1. - lcoords.at(1) - lcoords.at(2);

    FloatMatrix fullB;
    this->evaldNdx(fullB, lcoords_tet, cellgeo);
    answer.resize(snodes.giveSize(), 3);
    for ( int i = 1; i <= snodes.giveSize(); ++i ) {
        for ( int j = 1; j <= 3; ++j ) {
            answer.at(i, j) = fullB.at(snodes.at(i), j);
        }
    }
}

double
FEI3dTetQuad :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray snodes(3);
    FloatArray a, b;
    this->computeLocalSurfaceMapping(snodes, isurf);

    double l1, l2, l3;
    l1 = lcoords.at(1);
    l2 = lcoords.at(2);
    l3 = 1.0 - l1 - l2;

    FloatArray dNdxi(6), dNdeta(6);

    dNdxi(0) = 4.0 * l1 - 1.0;
    dNdxi(1) = 0.0;
    dNdxi(2) = -1.0 * ( 4.0 * l3 - 1.0 );
    dNdxi(3) = 4.0 * l2;
    dNdxi(4) = -4.0 * l2;
    dNdxi(5) = 4.0 * l3 - 4.0 * l1;

    dNdeta(0) = 0.0;
    dNdeta(1) = 4.0 * l2 - 1.0;
    dNdeta(2) = -1.0 * ( 4.0 * l3 - 1.0 );
    dNdeta(3) = 4.0 * l1;
    dNdeta(4) = 4.0 * l3 - 4.0 * l2;
    dNdeta(5) = -4.0 * l1;

    for ( int i = 0; i < 6; ++i ) {
        a.add( dNdxi(i),  * cellgeo.giveVertexCoordinates( snodes(i) ) );
        b.add( dNdeta(i), * cellgeo.giveVertexCoordinates( snodes(i) ) );
    }
    answer.beVectorProductOf(a, b);
    return answer.normalize();
}

double
FEI3dTetQuad :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo)
{
    FloatArray normal;
    return this->surfaceEvalNormal(normal, isurf, lcoords, cellgeo);
}

void
FEI3dTetQuad :: computeLocalSurfaceMapping(IntArray &surfNodes, int isurf)
{
    int aNode = 0, bNode = 0, cNode = 0, dNode = 0, eNode = 0, fNode = 0;
    surfNodes.resize(6);

    if ( isurf == 1 ) {
        aNode = 1;
        bNode = 3;
        cNode = 2;
        dNode = 7;
        eNode = 6;
        fNode = 5;
    } else if ( isurf == 2 ) {
        aNode = 1;
        bNode = 2;
        cNode = 4;
        dNode = 5;
        eNode = 9;
        fNode = 8;
    } else if ( isurf == 3 ) {
        aNode = 2;
        bNode = 3;
        cNode = 4;
        dNode = 6;
        eNode = 10;
        fNode = 9;
    } else if ( isurf == 4 ) {
        aNode = 1;
        bNode = 4;
        cNode = 3;
        dNode = 8;
        eNode = 10;
        fNode = 7;
    } else {
        OOFEM_ERROR("wrong surface number (%d)", isurf);
    }

    surfNodes.at(1) = aNode;
    surfNodes.at(2) = bNode;
    surfNodes.at(3) = cNode;
    surfNodes.at(4) = dNode;
    surfNodes.at(5) = eNode;
    surfNodes.at(6) = fNode;
}

double FEI3dTetQuad :: evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo)
{
    IntArray fNodes;
    this->computeLocalSurfaceMapping(fNodes, iEdge);

    const FloatArray &c1 = * cellgeo.giveVertexCoordinates( fNodes.at(1) );
    const FloatArray &c2 = * cellgeo.giveVertexCoordinates( fNodes.at(2) );
    const FloatArray &c3 = * cellgeo.giveVertexCoordinates( fNodes.at(3) );
    const FloatArray &c4 = * cellgeo.giveVertexCoordinates( fNodes.at(4) );
    const FloatArray &c5 = * cellgeo.giveVertexCoordinates( fNodes.at(5) );
    const FloatArray &c6 = * cellgeo.giveVertexCoordinates( fNodes.at(6) );

    // Expression derived in Mathematica:
    return (
               c1(2) * ( c2(1) * ( -2 * c3(0) -  3 * c4(0) +  5 * c5(0) +  5 * c6(0) ) +
                        c3(1) * ( 2 * c2(0)            -  5 * c4(0) -  5 * c5(0) +  3 * c6(0) ) +
                        c4(1) * ( 3 * c2(0) +  5 * c3(0)            -  4 * c5(0) - 24 * c6(0) ) +
                        c5(1) * ( -5 * c2(0) +  5 * c3(0) +  4 * c4(0)            -  4 * c6(0) ) +
                        c6(1) * ( -5 * c2(0) -  3 * c3(0) + 24 * c4(0) +  4 * c5(0) ) ) +
               c2(2) * ( c1(1) * ( 2 * c3(0) +  3 * c4(0) -  5 * c5(0) -  5 * c6(0) ) +
                        c3(1) * ( -2 * c1(0)                       +  5 * c4(0) -  3 * c5(0) +  5 * c6(0) ) +
                        c4(1) * ( -3 * c1(0)            -  5 * c3(0)            + 24 * c5(0) +  4 * c6(0) ) +
                        c5(1) * ( 5 * c1(0)            +  3 * c3(0) - 24 * c4(0)            -  4 * c6(0) ) +
                        c6(1) * ( 5 * c1(0)            -  5 * c3(0) -  4 * c4(0) +  4 * c5(0) ) ) +
               c3(2) * ( c1(1) * ( -2 * c2(0)            +  5 * c4(0) +  5 * c5(0) -  3 * c6(0) ) +
                        c2(1) * ( 2 * c1(0)                       -  5 * c4(0) +  3 * c5(0) -  5 * c6(0) ) +
                        c4(1) * ( -5 * c1(0) +  5 * c2(0)                       -  4 * c5(0) +  4 * c6(0) ) +
                        c5(1) * ( -5 * c1(0) -  3 * c2(0)            +  4 * c4(0)            + 24 * c6(0) ) +
                        c6(1) * ( 3 * c1(0) +  5 * c2(0)            -  4 * c4(0) - 24 * c5(0) ) ) +
               c4(2) * ( c1(1) * ( -3 * c2(0) -  5 * c3(0)            +  4 * c5(0) + 24 * c6(0) ) +
                        c2(1) * ( 3 * c1(0)            +  5 * c3(0)            - 24 * c5(0) -  4 * c6(0) ) +
                        c3(1) * ( 5 * c1(0) -  5 * c2(0)                       +  4 * c5(0) -  4 * c6(0) ) +
                        c5(1) * ( -4 * c1(0) + 24 * c2(0) -  4 * c3(0)                       - 16 * c6(0) ) +
                        c6(1) * ( -24 * c1(0) +  4 * c2(0) +  4 * c3(0)            + 16 * c5(0) ) ) +
               c5(2) * ( c1(1) * ( 5 * c2(0) -  5 * c3(0) -  4 * c4(0)            +  4 * c6(0) ) +
                        c2(1) * ( -5 * c1(0)            -  3 * c3(0) + 24 * c4(0)            +  4 * c6(0) ) +
                        c3(1) * ( 5 * c1(0) +  3 * c2(0)            -  4 * c4(0)            - 24 * c6(0) ) +
                        c4(1) * ( 4 * c1(0) - 24 * c2(0) +  4 * c3(0)                       + 16 * c6(0) ) +
                        c6(1) * ( -4 * c1(0) -  4 * c2(0) + 24 * c3(0) - 16 * c4(0) ) ) +
               c6(2) * ( c1(1) * ( 5 * c2(0) +  3 * c3(0) - 24 * c4(0) -  4 * c5(0) ) +
                        c2(1) * ( -5 * c1(0)            +  5 * c3(0) +  4 * c4(0) -  4 * c5(0) ) +
                        c3(1) * ( -3 * c1(0) -  5 * c2(0)            +  4 * c4(0) + 24 * c5(0) ) +
                        c4(1) * ( 24 * c1(0) -  4 * c2(0) -  4 * c3(0)            - 16 * c5(0) ) +
                        c5(1) * ( 4 * c1(0) +  4 * c2(0) - 24 * c3(0) + 16 * c4(0) ) )
               ) / 30.;
}

IntegrationRule *
FEI3dTetQuad :: giveIntegrationRule(int order)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Tetrahedra, order + 3);
    iRule->SetUpPointsOnTetrahedra(points, _Unknown);
    return iRule;
}

IntegrationRule *
FEI3dTetQuad :: giveBoundaryIntegrationRule(int order, int boundary)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, NULL);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, order + 2);
    iRule->SetUpPointsOnTriangle(points, _Unknown);
    return iRule;
}
} // end namespace oofem
