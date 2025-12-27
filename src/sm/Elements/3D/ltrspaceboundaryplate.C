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

#include "sm/Elements/3D/ltrspaceboundaryplate.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei3dtetlin.h"
#include "classfactory.h"
#include "Materials/structuralms.h"


namespace oofem {
REGISTER_Element(LTRSpaceBoundaryPlate);

LTRSpaceBoundaryPlate :: LTRSpaceBoundaryPlate(int n, Domain *aDomain) :
    LTRSpaceBoundary(n, aDomain)
{}

void
LTRSpaceBoundaryPlate :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode == 5 ) {
        answer = { E_xx, E_xy, E_yx, E_yy, E_zx, E_zy, K_xx, K_yy, K_xy, K_yx };
    } else {
        answer = { D_u, D_v, D_w };
    }
}


void
LTRSpaceBoundaryPlate :: computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FloatArray unitCellSize;
    unitCellSize.resize(2);
    unitCellSize.at(1) = this->giveNode(5)->giveCoordinate(1);
    unitCellSize.at(2) = this->giveNode(5)->giveCoordinate(2);

    IntArray switches1, switches2, switches3, switches4;
    this->giveSwitches(switches1, this->location.at(1) );
    this->giveSwitches(switches2, this->location.at(2) );
    this->giveSwitches(switches3, this->location.at(3) );
    this->giveSwitches(switches4, this->location.at(4) );

    FloatArray w(4);
    for ( int i = 1; i <= 4; i++ ) {
        w.at(i) = this->giveNode(i)->giveCoordinate(3);
    }

    FloatMatrix k1, k2, k3, k4;
    k1.resize(3, 10);
    k2.resize(3, 10);
    k3.resize(3, 10);
    k4.resize(3, 10);

    k1.at(1, 1) = unitCellSize.at(1) * switches1.at(1);
    k1.at(1, 2) = unitCellSize.at(2) * switches1.at(2);
    k1.at(2, 3) = unitCellSize.at(1) * switches1.at(1);
    k1.at(2, 4) = unitCellSize.at(2) * switches1.at(2);
    k1.at(3, 5) = unitCellSize.at(1) * switches1.at(1);
    k1.at(3, 6) = unitCellSize.at(2) * switches1.at(2);
    k1.at(1, 7) = -w.at(1) * unitCellSize.at(1) * switches1.at(1);
    k1.at(2, 8) = -w.at(1) * unitCellSize.at(2) * switches1.at(2);
    k1.at(1, 9) = -w.at(1) * unitCellSize.at(2) * switches1.at(2);
    k1.at(2, 10) = -w.at(1) * unitCellSize.at(1) * switches1.at(1);

    k2.at(1, 1) = unitCellSize.at(1) * switches2.at(1);
    k2.at(1, 2) = unitCellSize.at(2) * switches2.at(2);
    k2.at(2, 3) = unitCellSize.at(1) * switches2.at(1);
    k2.at(2, 4) = unitCellSize.at(2) * switches2.at(2);
    k2.at(3, 5) = unitCellSize.at(1) * switches2.at(1);
    k2.at(3, 6) = unitCellSize.at(2) * switches2.at(2);
    k2.at(1, 7) = -w.at(2) * unitCellSize.at(1) * switches2.at(1);
    k2.at(2, 8) = -w.at(2) * unitCellSize.at(2) * switches2.at(2);
    k2.at(1, 9) = -w.at(2) * unitCellSize.at(2) * switches2.at(2);
    k2.at(2, 10) = -w.at(2) * unitCellSize.at(1) * switches2.at(1);

    k3.at(1, 1) = unitCellSize.at(1) * switches3.at(1);
    k3.at(1, 2) = unitCellSize.at(2) * switches3.at(2);
    k3.at(2, 3) = unitCellSize.at(1) * switches3.at(1);
    k3.at(2, 4) = unitCellSize.at(2) * switches3.at(2);
    k3.at(3, 5) = unitCellSize.at(1) * switches3.at(1);
    k3.at(3, 6) = unitCellSize.at(2) * switches3.at(2);
    k3.at(1, 7) = -w.at(3) * unitCellSize.at(1) * switches3.at(1);
    k3.at(2, 8) = -w.at(3) * unitCellSize.at(2) * switches3.at(2);
    k3.at(1, 9) = -w.at(3) * unitCellSize.at(2) * switches3.at(2);
    k3.at(2, 10) = -w.at(3) * unitCellSize.at(1) * switches3.at(1);

    k4.at(1, 1) = unitCellSize.at(1) * switches4.at(1);
    k4.at(1, 2) = unitCellSize.at(2) * switches4.at(2);
    k4.at(2, 3) = unitCellSize.at(1) * switches4.at(1);
    k4.at(2, 4) = unitCellSize.at(2) * switches4.at(2);
    k4.at(3, 5) = unitCellSize.at(1) * switches4.at(1);
    k4.at(3, 6) = unitCellSize.at(2) * switches4.at(2);
    k4.at(1, 7) = -w.at(4) * unitCellSize.at(1) * switches4.at(1);
    k4.at(2, 8) = -w.at(4) * unitCellSize.at(2) * switches4.at(2);
    k4.at(1, 9) = -w.at(4) * unitCellSize.at(2) * switches4.at(2);
    k4.at(2, 10) = -w.at(4) * unitCellSize.at(1) * switches4.at(1);

    answer.resize(12, 12);
    answer.beUnitMatrix();
    answer.resizeWithData(12, 22);

    answer.assemble(k1, { 1, 2, 3 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 });
    answer.assemble(k2, { 4, 5, 6 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 });
    answer.assemble(k3, { 7, 8, 9 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 });
    answer.assemble(k4, { 10, 11, 12 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 });
}
} // end namespace oofem
