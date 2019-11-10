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

#include "sm/Elements/3D/ltrspaceboundaryvoigt.h"
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
REGISTER_Element(LTRSpaceBoundaryVoigt);

LTRSpaceBoundaryVoigt :: LTRSpaceBoundaryVoigt(int n, Domain *aDomain) :
    LTRSpaceBoundary(n, aDomain)
{}

void
LTRSpaceBoundaryVoigt :: initializeFrom(InputRecord &ir)
{
    LTRSpaceBoundary :: initializeFrom(ir);
}

void
LTRSpaceBoundaryVoigt :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode == 5 ) {
        answer = { E_xx, E_yy, E_zz, G_yz, G_xz, G_xy };
    } else {
        answer = { D_u, D_v, D_w };
    }
}


void
LTRSpaceBoundaryVoigt :: computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FloatArray unitCellSize;
    unitCellSize.resize(3);
    unitCellSize.at(1) = this->giveNode(5)->giveCoordinate(1);
    unitCellSize.at(2) = this->giveNode(5)->giveCoordinate(2);
    unitCellSize.at(3) = this->giveNode(5)->giveCoordinate(3);

    IntArray switches1, switches2, switches3, switches4;
    this->giveSwitches(switches1, this->location.at(1) );
    this->giveSwitches(switches2, this->location.at(2) );
    this->giveSwitches(switches3, this->location.at(3) );
    this->giveSwitches(switches4, this->location.at(4) );

    FloatMatrix k21_node1, k22_node1, k21_node2, k22_node2, k21_node3, k22_node3, k21_node4, k22_node4;
    k21_node1.resize(3, 3);
    k22_node1.resize(3, 3);
    k21_node2.resize(3, 3);
    k22_node2.resize(3, 3);
    k21_node3.resize(3, 3);
    k22_node3.resize(3, 3);
    k21_node4.resize(3, 3);
    k22_node4.resize(3, 3);

    k21_node1.at(1, 1) = unitCellSize.at(1) * switches1.at(1);
    k21_node1.at(2, 2) = unitCellSize.at(2) * switches1.at(2);
    k21_node1.at(3, 3) = unitCellSize.at(3) * switches1.at(3);

    k22_node1.at(1, 2) = unitCellSize.at(3) * switches1.at(3);
    k22_node1.at(1, 3) = unitCellSize.at(2) * switches1.at(2);
    k22_node1.at(2, 1) = unitCellSize.at(3) * switches1.at(3);

    k21_node2.at(1, 1) = unitCellSize.at(1) * switches2.at(1);
    k21_node2.at(2, 2) = unitCellSize.at(2) * switches2.at(2);
    k21_node2.at(3, 3) = unitCellSize.at(3) * switches2.at(3);

    k22_node2.at(1, 2) = unitCellSize.at(3) * switches2.at(3);
    k22_node2.at(1, 3) = unitCellSize.at(2) * switches2.at(2);
    k22_node2.at(2, 1) = unitCellSize.at(3) * switches2.at(3);

    k21_node3.at(1, 1) = unitCellSize.at(1) * switches3.at(1);
    k21_node3.at(2, 2) = unitCellSize.at(2) * switches3.at(2);
    k21_node3.at(3, 3) = unitCellSize.at(3) * switches3.at(3);

    k22_node3.at(1, 2) = unitCellSize.at(3) * switches3.at(3);
    k22_node3.at(1, 3) = unitCellSize.at(2) * switches3.at(2);
    k22_node3.at(2, 1) = unitCellSize.at(3) * switches3.at(3);

    k21_node4.at(1, 1) = unitCellSize.at(1) * switches4.at(1);
    k21_node4.at(2, 2) = unitCellSize.at(2) * switches4.at(2);
    k21_node4.at(3, 3) = unitCellSize.at(3) * switches4.at(3);

    k22_node4.at(1, 2) = unitCellSize.at(3) * switches4.at(3);
    k22_node4.at(1, 3) = unitCellSize.at(2) * switches4.at(2);
    k22_node4.at(2, 1) = unitCellSize.at(3) * switches4.at(3);

    answer.resize(12, 12);
    answer.beUnitMatrix();
    answer.resizeWithData(12, 18);

    answer.assemble(k21_node1, { 1, 2, 3 }, { 13, 14, 15 });
    answer.assemble(k22_node1, { 1, 2, 3 }, { 16, 17, 18 });
    answer.assemble(k21_node2, { 4, 5, 6 }, { 13, 14, 15 });
    answer.assemble(k22_node2, { 4, 5, 6 }, { 16, 17, 18 });
    answer.assemble(k21_node3, { 7, 8, 9 }, { 13, 14, 15 });
    answer.assemble(k22_node3, { 7, 8, 9 }, { 16, 17, 18 });
    answer.assemble(k21_node4, { 10, 11, 12 }, { 13, 14, 15 });
    answer.assemble(k22_node4, { 10, 11, 12 }, { 16, 17, 18 });
}
} // end namespace oofem
