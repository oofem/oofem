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

#include "cohint.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {

REGISTER_Material(CohesiveInterfaceMaterial);

CohesiveInterfaceMaterial :: CohesiveInterfaceMaterial(int n, Domain *d) : StructuralInterfaceMaterial(n, d) { }

void
CohesiveInterfaceMaterial :: giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus(gp) );

    answer.resize(3);

    // normal part of elastic stress-strain law
    answer.at(1) = kn * jump.at(1);

    // shear part of elastic stress-strain law
    answer.at(2) = ks * jump.at(2);
    answer.at(3) = ks * jump.at(3);

    // update gp
    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);
}


void
CohesiveInterfaceMaterial :: give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = kn;
    answer.at(2, 2) = answer.at(3, 3) = ks;
}

int
CohesiveInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
}

IRResultType
CohesiveInterfaceMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    // elastic parameters
    IR_GIVE_FIELD(ir, kn, _IFT_CohesiveInterfaceMaterial_kn);
    IR_GIVE_FIELD(ir, ks, _IFT_CohesiveInterfaceMaterial_ks);

    return StructuralInterfaceMaterial :: initializeFrom(ir);
}

void
CohesiveInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);
    input.setField(this->kn, _IFT_CohesiveInterfaceMaterial_kn);
    input.setField(this->ks, _IFT_CohesiveInterfaceMaterial_ks);
}

} // namespace oofem
