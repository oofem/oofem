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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "sm/Elements/Beams/libeam3dboundarymembrane.h"
#include "sm/Materials/structuralms.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(LIBeam3dBoundaryMembrane);

LIBeam3dBoundaryMembrane :: LIBeam3dBoundaryMembrane(int n, Domain *aDomain) : LIBeam3dBoundary(n, aDomain)
    // Constructor.
{}


void
LIBeam3dBoundaryMembrane :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode == 3 ) {
        answer = { E_xx, E_xy, E_yx, E_yy };
    } else {
        answer = { D_u, D_v, D_w, R_u, R_v, R_w };
    }
}


void
LIBeam3dBoundaryMembrane :: computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FloatArray unitCellSize;
    unitCellSize.resize(2);
    unitCellSize.at(1) = this->giveNode(3)->giveCoordinate(1);
    unitCellSize.at(2) = this->giveNode(3)->giveCoordinate(2);

    IntArray switches1, switches2;
    this->giveSwitches(switches1, this->location.at(1) );
    this->giveSwitches(switches2, this->location.at(2) );

    FloatMatrix k1, k2;
    k1.resize(6, 4);
    k2.resize(6, 4);

    k1.at(1, 1) = unitCellSize.at(1) * switches1.at(1);
    k1.at(1, 2) = unitCellSize.at(2) * switches1.at(2);
    k1.at(2, 3) = unitCellSize.at(1) * switches1.at(1);
    k1.at(2, 4) = unitCellSize.at(2) * switches1.at(2);

    k2.at(1, 1) = unitCellSize.at(1) * switches2.at(1);
    k2.at(1, 2) = unitCellSize.at(2) * switches2.at(2);
    k2.at(2, 3) = unitCellSize.at(1) * switches2.at(1);
    k2.at(2, 4) = unitCellSize.at(2) * switches2.at(2);

    answer.resize(12, 12);
    answer.beUnitMatrix();
    answer.resizeWithData(12, 16);

    answer.assemble(k1, { 1, 2, 3, 4, 5, 6 }, { 13, 14, 15, 16 });
    answer.assemble(k2, { 7, 8, 9, 10, 11, 12 }, { 13, 14, 15, 16 });
}
} // end namespace oofem
