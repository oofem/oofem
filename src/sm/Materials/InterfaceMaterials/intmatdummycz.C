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

#include "intmatdummycz.h"
#include "classfactory.h"
#include "gausspoint.h"

namespace oofem {
REGISTER_Material(IntMatDummyCZ);

IntMatDummyCZ :: IntMatDummyCZ(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{

}

IntMatDummyCZ :: ~IntMatDummyCZ()
{

}

void IntMatDummyCZ :: giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                                const FloatMatrix &F, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus(gp) );

    status->letTempJumpBe(jump);

    answer.resize(3);
    answer.zero();
}

void IntMatDummyCZ :: give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(3,3);
    answer.zero();
}

IRResultType IntMatDummyCZ :: initializeFrom(InputRecord *ir)
{
    return StructuralInterfaceMaterial :: initializeFrom(ir);
}

void IntMatDummyCZ :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);
}

void IntMatDummyCZ :: printYourself()
{
    printf("I am a IntMatDummyCZ.\n");
}

} /* namespace oofem */
