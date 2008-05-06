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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//
// FILE: internalstatetype.h
//

#ifndef internalstatetype_h
#define internalstatetype_h

#include "enumitem.h"

//
// following type determine the mode of some value.
// which can be requested from various specialized methods.
// particular specialized methods (for example method for computing the load vector)
// are general, i.e., they are able to compute response for
// both totalLoadVector and incrementalLoadVector charTypes.
// The particular type of response is then requested using parameter of ValueModeType type.
//

/**
 * Type  representing the physical meaning of element or constitutive model internal variable.
 * Values of this type are used, when these internal variables are requested.
 */
#define InternalStateType_DEF \
    ENUM_ITEM_WITH_VALUE(IST_Undedfined, 0) \
    ENUM_ITEM_WITH_VALUE(IST_StressTensor, 1) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalStressTensor, 2) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalStressTempTensor, 3) \
    ENUM_ITEM_WITH_VALUE(IST_StrainTensor, 4) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalStrainTensor, 5) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalStrainTempTensor, 6) \
    ENUM_ITEM_WITH_VALUE(IST_BeamForceMomentumTensor, 7) \
    ENUM_ITEM_WITH_VALUE(IST_BeamStrainCurvatureTensor, 8) \
    ENUM_ITEM_WITH_VALUE(IST_ShellForceMomentumTensor, 9) \
    ENUM_ITEM_WITH_VALUE(IST_ShellStrainCurvatureTensor, 10) \
    ENUM_ITEM_WITH_VALUE(IST_CurvatureTensor, 11) \
    ENUM_ITEM_WITH_VALUE(IST_DisplacementVector, 12) \
    ENUM_ITEM_WITH_VALUE(IST_DamageTensor, 13) \
    ENUM_ITEM_WITH_VALUE(IST_DamageInvTensor, 14) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalDamageTensor, 15) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalDamageTempTensor, 16) \
    ENUM_ITEM_WITH_VALUE(IST_CrackState, 17) \
    ENUM_ITEM_WITH_VALUE(IST_StressTensorTemp, 18) \
    ENUM_ITEM_WITH_VALUE(IST_StrainTensorTemp, 19) \
    ENUM_ITEM_WITH_VALUE(IST_ForceTensorTemp, 20) \
    ENUM_ITEM_WITH_VALUE(IST_MomentumTensorTemp, 21) \
    ENUM_ITEM_WITH_VALUE(IST_CurvatureTensorTemp, 22) \
    ENUM_ITEM_WITH_VALUE(IST_DisplacementVectorTemp, 23) \
    ENUM_ITEM_WITH_VALUE(IST_DamageTensorTemp, 24) \
    ENUM_ITEM_WITH_VALUE(IST_DamageInvTensorTemp, 25) \
    ENUM_ITEM_WITH_VALUE(IST_CrackStateTemp, 26) \
    ENUM_ITEM_WITH_VALUE(IST_PlasticStrainTensor, 27) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalPlasticStrainTensor, 28) \
    ENUM_ITEM_WITH_VALUE(IST_CylindricalStressTensor, 29) \
    ENUM_ITEM_WITH_VALUE(IST_CylindricalStrainTensor, 30) \
    ENUM_ITEM_WITH_VALUE(IST_MaxEquivalentStrainLevel, 31) \
    ENUM_ITEM_WITH_VALUE(IST_ErrorIndicatorLevel, 32) \
    ENUM_ITEM_WITH_VALUE(IST_InternalStressError, 33) \
    ENUM_ITEM_WITH_VALUE(IST_PrimaryUnknownError, 34) \
    ENUM_ITEM_WITH_VALUE(IST_RelMeshDensity, 35) \
    ENUM_ITEM_WITH_VALUE(IST_MicroplaneDamageValues, 36) \
    ENUM_ITEM_WITH_VALUE(IST_Temperature, 37) \
    ENUM_ITEM_WITH_VALUE(IST_MassConcentration_1, 38) \
    ENUM_ITEM_WITH_VALUE(IST_HydrationDegree, 39) \
    ENUM_ITEM_WITH_VALUE(IST_Humidity, 40) \
    ENUM_ITEM_WITH_VALUE(IST_Velocity, 41) \
    ENUM_ITEM_WITH_VALUE(IST_Pressure, 42) \
    ENUM_ITEM_WITH_VALUE(IST_VOFFraction, 43) \
    ENUM_ITEM_WITH_VALUE(IST_Density, 44) \
    ENUM_ITEM_WITH_VALUE(IST_MaterialInterfaceVal, 45) \
  \
    ENUM_ITEM(CrackStatuses) \
    ENUM_ITEM(CrackedFlag) \
    ENUM_ITEM(CrackDirs) \

enum InternalStateType {
    InternalStateType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


char *__InternalStateTypeToString(InternalStateType _value);


#endif // internalstatetype_h

