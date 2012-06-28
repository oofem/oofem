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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "error.h"
#include "chartype.h"
#include "internalstatetype.h"
#include "internalstatevaluetype.h"
#include "elementgeometrytype.h"
#include "unknowntype.h"
#include "materialmode.h"
#include "matresponsemode.h"
#include "valuemodetype.h"
#include "materialmappingalgorithmtype.h"
#include "meshpackagetype.h"
#include "domaintype.h"
#include "doftype.h"
#include "dofiditem.h"
#include "contextioerr.h"
#include "field.h"
#include "equationid.h"

#ifndef __MAKEDEPEND
 #include <cstring>
#endif

namespace oofem {

char cltypesGiveUnknownTypeModeKey(ValueModeType mode)
{
    switch ( mode ) {
    case VM_Unknown:      return 0;

    case VM_Total:        return 'u';

    case VM_Velocity:     return 'v';

    case VM_Acceleration: return 'a';

    default: OOFEM_ERROR("cltypesGiveUnknownTypeModeKey : unsupported ValueModeType");
    }

    return 0;
}


InternalStateValueType giveInternalStateValueType(InternalStateType type)
{
    switch ( type ) {
    case IST_StressTensor:
    case IST_StrainTensor:
    case IST_CurvatureTensor:
    case IST_DamageTensor:
    case IST_DamageInvTensor:
    case IST_PrincipalDamageTensor:
    case IST_PrincipalDamageTempTensor:
    case IST_StressTensorTemp:
    case IST_StrainTensorTemp:
    case IST_ForceTensorTemp:
    case IST_MomentumTensorTemp:
    case IST_CurvatureTensorTemp:
    case IST_DamageTensorTemp:
    case IST_DamageInvTensorTemp:
    case IST_PlasticStrainTensor:
    case IST_CylindricalStressTensor:
    case IST_CylindricalStrainTensor:
    case IST_MaxEquivalentStrainLevel:
    case IST_CumPlasticStrain:
    case IST_CumPlasticStrain_2:
    case IST_DeviatoricStrain:
    case IST_DeviatoricStress:
        return ISVT_TENSOR_S3;

    case IST_BeamForceMomentumTensor:
    case IST_BeamStrainCurvatureTensor:
    case IST_ShellForceMomentumTensor:
    case IST_ShellStrainCurvatureTensor:
        return ISVT_TENSOR_G;

    case IST_PrincipalStressTensor:
    case IST_PrincipalStressTempTensor:
    case IST_PrincipalStrainTensor:
    case IST_PrincipalStrainTempTensor:
    case IST_PrincipalPlasticStrainTensor:
    case IST_DisplacementVector:
    case IST_DisplacementVectorTemp:
    case IST_CrackState:
    case IST_CrackStateTemp:
    case IST_MicroplaneDamageValues:
    case IST_Velocity:
    case IST_MaterialOrientation_x:
    case IST_MaterialOrientation_y:
    case IST_MaterialOrientation_z:
    case IST_TemperatureFlow:
    case IST_MassConcentrationFlow_1:
    case IST_HumidityFlow:
    case IST_CrackDirs:
    case IST_CrackStatuses:
    case IST_CrackVector:
        return ISVT_VECTOR;

    case IST_ErrorIndicatorLevel:
    case IST_InternalStressError:
    case IST_PrimaryUnknownError:
    case IST_RelMeshDensity:
    case IST_Temperature:
    case IST_MassConcentration_1:
    case IST_HydrationDegree:
    case IST_Humidity:
    case IST_Pressure:
    case IST_VOFFraction:
    case IST_Density:
    case IST_MaterialInterfaceVal:
    case IST_MaterialNumber:
    case IST_ElementNumber:
    case IST_BoneVolumeFraction:
    case IST_PlasStrainEnerDens:
    case IST_ElasStrainEnerDens:
    case IST_TotalStrainEnerDens:
    case IST_DamageScalar:
    case IST_CrackedFlag:
    case IST_StressWorkDensity:
    case IST_DissWorkDensity:
    case IST_FreeEnergyDensity:
    case IST_ThermalConductivityIsotropic:
    case IST_HeatCapacity:
    case IST_AverageTemperature:
    case IST_YoungModulusVirginPaste:
    case IST_PoissonRatioVirginPaste:
    case IST_YoungModulusConcrete:
    case IST_PoissonRatioConcrete:
    case IST_Viscosity:
    case IST_vonMisesStress:
        return ISVT_SCALAR;

    default:
        return ISVT_UNDEFINED;
    }
}



ContextIOERR :: ContextIOERR(contextIOResultType e, const char *file, int line)
{
    error = e;
    this->file = file;
    this->line = line;
    this->msg = NULL;
}

ContextIOERR :: ContextIOERR(contextIOResultType e, const char *msg, const char *file, int line)
{
    error = e;
    this->file = file;
    this->line = line;
    this->msg  = msg;
}

ContextIOERR :: ~ContextIOERR() { }


void
ContextIOERR :: print()
{
    if ( msg ) {
        __OOFEM_ERROR3(file, line, "ContextIOERR encountered, error code: %d\n%s", error, msg);
    }  else {
        __OOFEM_ERROR2(file, line, "ContextIOERR encountered, error code: %d", error);
    }
}

/*
 * The enum declaration and to_string conversion
 * inspired by X-Macros technique, as described in Wikipedia entry on C preprocessor
 * (http://en.wikipedia.org/wiki/C_preprocessor)
 */

#define ENUM_ITEM(element) case element: return # element;

#define ENUM_ITEM_WITH_VALUE(element, val) case element: return # element;

#define TO_STRING_BODY(enum_def)                        \
    switch ( _value ) { \
        enum_def \
    default: return "Unknown"; \
    }

const char *__InternalStateTypeToString(InternalStateType _value) {
    TO_STRING_BODY(InternalStateType_DEF)
}

const char *__UnknownTypeToString(UnknownType _value) {
    TO_STRING_BODY(UnknownType_DEF)
}

const char *__dofTypeToString(dofType _value) {
    TO_STRING_BODY(dofType_DEF)
}

const char *__domainTypeToString(domainType _value) {
    TO_STRING_BODY(domainType_DEF)
}

const char *__MaterialModeToString(MaterialMode _value) {
    TO_STRING_BODY(MaterialMode_DEF)
}

const char *__Element_Geometry_TypeToString(Element_Geometry_Type _value) {
    TO_STRING_BODY(Element_Geometry_Type_DEF)
}

const char *__ValueModeTypeToString(ValueModeType _value) {
    TO_STRING_BODY(ValueModeType_DEF)
}

const char *__MatResponseModeToString(MatResponseMode _value) {
    TO_STRING_BODY(MatResponseMode_DEF)
}

const char *__DofIDItemToString(DofIDItem _value) {
    TO_STRING_BODY(DofIDItem_DEF)
}

const char *__CharTypeToString(CharType _value) {
    TO_STRING_BODY(CharType_DEF)
}

const char *__MaterialMappingAlgorithmTypeToString(MaterialMappingAlgorithmType _value) {
    TO_STRING_BODY(MaterialMappingAlgorithmType_DEF)
}

const char *__MeshPackageTypeToString(MeshPackageType _value) {
    TO_STRING_BODY(MeshPackageType_DEF)
}

const char *__EquationIDToString(EquationID _value) {
    TO_STRING_BODY(EquationID_DEF)
}

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef TO_STRING_BODY



#define ENUM_ITEM(element) \
    if (!strncasecmp(_value, # element, strlen(# element) ) ) { \
        return element; \
    }

#define ENUM_ITEM_WITH_VALUE(element, val) \
    if (!strncasecmp(_value, # element, strlen(# element) ) ) { \
        return element; \
    }

#define TO_STRING_BODY(enum_def) \
        enum_def \
        OOFEM_ERROR2("Can not convert string %s to value", _value);\
        return FieldType(0);

FieldType __StringToFieldType(const char * _value) {
    TO_STRING_BODY(FieldType_DEF)
}




} // end namespace oofem
