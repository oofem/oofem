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

#ifndef internalstatetype_h
#define internalstatetype_h

#include "enumitem.h"

namespace oofem {
//
// following type determine the mode of some value.
// which can be requested from various specialized methods.
// particular specialized methods (for example method for computing the load vector)
// are general, i.e., they are able to compute response for
// both totalLoadVector and incrementalLoadVector charTypes.
// The particular type of response is then requested using parameter of ValueModeType type.
//
#define InternalStateType_DEF \
    ENUM_ITEM_WITH_VALUE(IST_Undefined, 0) \
    ENUM_ITEM_WITH_VALUE(IST_StressTensor, 1) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalStressTensor, 2) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalStressTempTensor, 3) \
    ENUM_ITEM_WITH_VALUE(IST_StrainTensor, 4) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalStrainTensor, 5) \
    ENUM_ITEM_WITH_VALUE(IST_PrincipalStrainTempTensor, 6) \
    ENUM_ITEM_WITH_VALUE(IST_BeamForceMomentumTensor, 7) \
    ENUM_ITEM_WITH_VALUE(IST_BeamStrainCurvatureTensor, 8) \
    ENUM_ITEM_WITH_VALUE(IST_ShellMomentumTensor, 9)  \
    ENUM_ITEM_WITH_VALUE(IST_ShellForceTensor, 10) \
    ENUM_ITEM_WITH_VALUE(IST_ShellCurvatureTensor, 11) \
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
    ENUM_ITEM_WITH_VALUE(IST_MaterialNumber, 46) \
    ENUM_ITEM_WITH_VALUE(IST_ElementNumber, 47) \
    ENUM_ITEM_WITH_VALUE(IST_BoneVolumeFraction, 48) \
    ENUM_ITEM_WITH_VALUE(IST_PlasStrainEnerDens, 49) \
    ENUM_ITEM_WITH_VALUE(IST_ElasStrainEnerDens, 50) \
    ENUM_ITEM_WITH_VALUE(IST_TotalStrainEnerDens, 51) \
    ENUM_ITEM_WITH_VALUE(IST_DamageScalar, 52) \
    ENUM_ITEM_WITH_VALUE(IST_MaterialOrientation_x, 53) \
    ENUM_ITEM_WITH_VALUE(IST_MaterialOrientation_y, 54) \
    ENUM_ITEM_WITH_VALUE(IST_MaterialOrientation_z, 55) \
    ENUM_ITEM_WITH_VALUE(IST_TemperatureFlow, 56) \
    ENUM_ITEM_WITH_VALUE(IST_MassConcentrationFlow_1, 57) \
    ENUM_ITEM_WITH_VALUE(IST_HumidityFlow, 58) \
    ENUM_ITEM_WITH_VALUE(IST_CrackStatuses, 59) \
    ENUM_ITEM_WITH_VALUE(IST_CrackedFlag, 60) \
    ENUM_ITEM_WITH_VALUE(IST_CrackDirs, 61)     \
    ENUM_ITEM_WITH_VALUE(IST_CumPlasticStrain, 62) \
    ENUM_ITEM_WITH_VALUE(IST_CumPlasticStrain_2, 63) \
    ENUM_ITEM_WITH_VALUE(IST_StressWorkDensity, 64) \
    ENUM_ITEM_WITH_VALUE(IST_DissWorkDensity, 65) \
    ENUM_ITEM_WITH_VALUE(IST_FreeEnergyDensity, 66) \
    ENUM_ITEM_WITH_VALUE(IST_ThermalConductivityIsotropic, 67) \
    ENUM_ITEM_WITH_VALUE(IST_HeatCapacity, 68) \
    ENUM_ITEM_WITH_VALUE(IST_AverageTemperature, 69) \
    ENUM_ITEM_WITH_VALUE(IST_YoungModulusVirginPaste, 70) \
    ENUM_ITEM_WITH_VALUE(IST_PoissonRatioVirginPaste, 71) \
    ENUM_ITEM_WITH_VALUE(IST_YoungModulusConcrete, 72) \
    ENUM_ITEM_WITH_VALUE(IST_PoissonRatioConcrete, 73) \
    ENUM_ITEM_WITH_VALUE(IST_VolumetricPlasticStrain, 74) \
    ENUM_ITEM_WITH_VALUE(IST_DeviatoricStrain, 75) \
    ENUM_ITEM_WITH_VALUE(IST_DeviatoricStress, 76) \
    ENUM_ITEM_WITH_VALUE(IST_Viscosity, 77)                     \
    ENUM_ITEM_WITH_VALUE(IST_CharacteristicLength, 78)  \
    ENUM_ITEM_WITH_VALUE(IST_DeviatoricStrainMeasure, 79) \
    ENUM_ITEM_WITH_VALUE(IST_DeviatoricStressMeasure, 80) \
    ENUM_ITEM_WITH_VALUE(IST_vonMisesStress, 81) \
    ENUM_ITEM_WITH_VALUE(IST_CrackVector, 82) \
    ENUM_ITEM_WITH_VALUE(IST_PressureGradient, 83) \
    ENUM_ITEM_WITH_VALUE(IST_DissWork, 84) \
    ENUM_ITEM_WITH_VALUE(IST_DeltaDissWork, 85) \
    ENUM_ITEM_WITH_VALUE(IST_StressCapPos, 86) \
    ENUM_ITEM_WITH_VALUE(IST_TangentNorm, 87) \
    ENUM_ITEM_WITH_VALUE(IST_Tangent, 88) \
    ENUM_ITEM_WITH_VALUE(IST_DirectorField, 89) \
    ENUM_ITEM_WITH_VALUE(IST_CrackWidth, 90) \
    ENUM_ITEM_WITH_VALUE(IST_DeformationGradientTensor, 91) \
    ENUM_ITEM_WITH_VALUE(IST_FirstPKStressTensor, 92) \
    ENUM_ITEM_WITH_VALUE(IST_XFEMEnrichment, 93) \
    ENUM_ITEM_WITH_VALUE(IST_XFEMNumIntersecPoints, 94) \
    ENUM_ITEM_WITH_VALUE(IST_XFEMLevelSetPhi, 95) \
    ENUM_ITEM_WITH_VALUE(IST_Maturity, 96) \
    ENUM_ITEM_WITH_VALUE(IST_CauchyStressTensor, 97) \
    ENUM_ITEM_WITH_VALUE(IST_InterfaceJump, 98) \
    ENUM_ITEM_WITH_VALUE(IST_InterfaceTraction, 99) \
    ENUM_ITEM_WITH_VALUE(IST_InterfaceFirstPKTraction, 100) \
    ENUM_ITEM_WITH_VALUE(IST_StressTensor_Reduced, 101) \
    ENUM_ITEM_WITH_VALUE(IST_StrainTensor_Reduced, 102) \
    ENUM_ITEM_WITH_VALUE(IST_CrossSectionNumber, 103 ) \
    ENUM_ITEM_WITH_VALUE(IST_ShellStrainTensor, 104) \
    ENUM_ITEM_WITH_VALUE(IST_AbaqusStateVector, 105)


/**
 * Type  representing the physical meaning of element or constitutive model internal variable.
 * Values of this type are used, when these internal variables are requested.
 */
enum InternalStateType {
    InternalStateType_DEF
};

enum ElementCharSizeMethod {
    ECSM_Unknown,
    ECSM_SquareRootOfArea,
    ECSM_Projection,
    ECSM_ProjectionCentered,
    ECSM_Oliver1,
    ECSM_Oliver1modified,
    ECSM_Oliver2
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__InternalStateTypeToString(InternalStateType _value);
} // end namespace oofem
#endif // internalstatetype_h
