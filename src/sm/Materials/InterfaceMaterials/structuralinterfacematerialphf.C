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

#include "structuralinterfacematerialphf.h"
#include "structuralinterfacematerialstatus.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"

namespace oofem {
StructuralInterfaceMaterialPhF :: StructuralInterfaceMaterialPhF(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{
    this->useNumericalTangent = false;
}



void
StructuralInterfaceMaterialPhF :: giveEngTraction_2d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const double damage, TimeStep *tStep)
{
    FloatArray jump3D(3), traction3D;
    jump3D = { jump.at(1), 0.0, jump.at(2) };
    this->giveEngTraction_3d(traction3D, gp, jump3D, damage, tStep);
    answer = { traction3D.at(1), traction3D.at(3) };

#ifdef DEBUG
    if ( fabs( traction3D.at(2) ) > 1.0e-3 ) {
        OOFEM_ERROR("Traction vector obtained from 3D state contains a nonzero thickness stress component")
    }
#endif
}





} // end namespace oofem
