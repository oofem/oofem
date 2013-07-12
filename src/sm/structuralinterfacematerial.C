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

#include "StructuralInterfaceMaterial.h"
#include "structuralinterfacematerialstatus.h"
#include "dynamicinputrecord.h"

namespace oofem {
int
StructuralInterfaceMaterial :: hasMaterialModeCapability(MaterialMode mode)
//  
// returns whether receiver supports given mode
//
{
    return mode == _1dInterface  ||  mode == _2dInterface ||
           mode == _3dInterface;
}


int
StructuralInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_InterfaceJump ) {
        answer = status->giveJump();
        return 1;

    } else if ( type == IST_InterfaceTraction ) {
        answer = status->giveTraction();
        return 1;

    } else if ( type == IST_InterfaceFirstPKTraction ) {
        answer = status->giveFirstPKTraction();
        return 1;

    } else if ( type == IST_DeformationGradientTensor ) {
        answer.beVectorForm( status->giveF() );
        return 1;

    } else {
        return Material :: giveIPValue(answer, gp, type, atTime);
    }
}


InternalStateValueType
StructuralInterfaceMaterial :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_InterfaceJump ) || ( type == IST_InterfaceTraction ) ||
        ( type == IST_InterfaceFirstPKTraction ) ) {
        return ISVT_VECTOR; 
    } else if ( type == IST_DeformationGradientTensor ) {
        return ISVT_TENSOR_G;
    } else {
        return Material :: giveIPValueType(type);
    }
}


int
StructuralInterfaceMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( ( type == IST_InterfaceJump ) || ( type == IST_InterfaceTraction ) ||
        ( type == IST_InterfaceFirstPKTraction ) ) {
            answer.setValues(3, 1, 2, 3);

    } else if ( type == IST_DeformationGradientTensor ) {
        answer.setValues(9, 1, 2, 3, 4, 5, 6, 7, 8, 9);
        return 1;
    } else {
        return Material :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}

// Currently not needed
IRResultType
StructuralInterfaceMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    this->Material :: initializeFrom(ir);
    return IRRT_OK;
}


void
StructuralInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    Material :: giveInputRecord(input);
}


} // end namespace oofem
