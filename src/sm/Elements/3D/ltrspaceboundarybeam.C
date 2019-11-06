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

#include "sm/Elements/3D/ltrspaceboundarybeam.h"
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
REGISTER_Element(LTRSpaceBoundaryBeam);

LTRSpaceBoundaryBeam :: LTRSpaceBoundaryBeam(int n, Domain *aDomain) :
    LTRSpaceBoundary(n, aDomain)
{}

void
LTRSpaceBoundaryBeam :: initializeFrom(InputRecord &ir)
{
    LTRSpaceBoundary :: initializeFrom(ir);
}

void
LTRSpaceBoundaryBeam :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode == 5 ) {
        answer = { E_xx, E_zx, K_xx };
    } else {
        answer = { D_u, D_v, D_w };
    }
}


void
LTRSpaceBoundaryBeam :: computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    double unitCellSize = this->giveNode(5)->giveCoordinate(1);

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
    k1.resize(3, 3);
    k2.resize(3, 3);
    k3.resize(3, 3);
    k4.resize(3, 3);

    k1.at(1, 1) = unitCellSize * switches1.at(1);
    k1.at(1, 3) = -w.at(1) * unitCellSize * switches1.at(1);
    k1.at(3, 2) = unitCellSize * switches1.at(1);

    k2.at(1, 1) = unitCellSize * switches2.at(1);
    k2.at(1, 3) = -w.at(2) * unitCellSize * switches2.at(1);
    k2.at(3, 2) = unitCellSize * switches2.at(1);

    k3.at(1, 1) = unitCellSize * switches3.at(1);
    k3.at(1, 3) = -w.at(3) * unitCellSize * switches3.at(1);
    k3.at(3, 2) = unitCellSize * switches3.at(1);

    k4.at(1, 1) = unitCellSize * switches4.at(1);
    k4.at(1, 3) = -w.at(4) * unitCellSize * switches4.at(1);
    k4.at(3, 2) = unitCellSize * switches4.at(1);

    answer.resize(12, 12);
    answer.beUnitMatrix();
    answer.resizeWithData(12, 15);

    answer.assemble(k1, { 1, 2, 3 }, { 13, 14, 15 });
    answer.assemble(k2, { 4, 5, 6 }, { 13, 14, 15 });
    answer.assemble(k3, { 7, 8, 9 }, { 13, 14, 15 });
    answer.assemble(k4, { 10, 11, 12 }, { 13, 14, 15 });
}
} // end namespace oofem
