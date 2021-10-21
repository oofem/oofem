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
#include "xfem/xfemmanager.h"

#include <cstring>
#include <string>

namespace oofem {
char cltypesGiveUnknownTypeModeKey(ValueModeType mode)
{
    switch ( mode ) {
    case VM_Unknown:        return 0;

    case VM_Total:          return 'u';

    case VM_Velocity:       return 'v';

    case VM_Acceleration:   return 'a';

    case VM_TotalIntrinsic: return 'i';

    default: OOFEM_ERROR("unsupported ValueModeType");
    }

    return 0;
}


InternalStateValueType giveInternalStateValueType(InternalStateType type)
{
    switch ( type ) {
    case IST_StrainTensor:
    case IST_StrainTensorTemp:
    case IST_DeviatoricStrain:
    case IST_PlasticStrainTensor:
    case IST_ThermalStrainTensor:
    case IST_ElasticStrainTensor:
    case IST_CylindricalStrainTensor:
    case IST_CreepStrainTensor:
    case IST_ShellStrainTensor:
    case IST_CurvatureTensor:
    case IST_CurvatureTensorTemp:
    case IST_EigenStrainTensor:
    case IST_CrackStrainTensor:
        return ISVT_TENSOR_S3E;

    case IST_StressTensor:
    case IST_StressTensorTemp:
    case IST_CylindricalStressTensor:
    case IST_DeviatoricStress:
    case IST_CauchyStressTensor:

    case IST_ShellForceTensor:
    case IST_ShellForceTensorTemp:
    case IST_ShellMomentTensor:
    case IST_MomentTensor:
    case IST_MomentTensorTemp:

    ///@todo Should be have these are S3E or just S3?
    case IST_AutogenousShrinkageTensor:
    case IST_DryingShrinkageTensor:
    case IST_TotalShrinkageTensor:

    // Damage tensors
    case IST_DamageTensor:
    case IST_DamageTensorTemp:
    case IST_DamageInvTensor:
    case IST_DamageInvTensorTemp:
    case IST_PrincipalDamageTensor:
    case IST_PrincipalDamageTempTensor:
        return ISVT_TENSOR_S3;

    case IST_DeformationGradientTensor:
    case IST_FirstPKStressTensor:
    case IST_MacroSlipGradient:
    case IST_ReinfMembraneStress:
    //case IST_MaterialOrientation:
        return ISVT_TENSOR_G;

    case IST_BeamForceMomentTensor:
    case IST_BeamStrainCurvatureTensor:
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
    case IST_CrackStatusesTemp:
    case IST_CrackVector:
    case IST_2ndCrackVector:
    case IST_3rdCrackVector:      
    case IST_InterfaceFirstPKTraction:
    case IST_InterfaceTraction:
    case IST_InterfaceJump:
    case IST_InterfaceNormal:
    case IST_PrincStressVector1:
    case IST_PrincStressVector2:
    case IST_PrincStressVector3:
    case IST_MacroSlipVector:
    case IST_TransferStress:
        return ISVT_VECTOR;

    case IST_MaxEquivalentStrainLevel:
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
    case IST_CumPlasticStrain:
    case IST_CumPlasticStrain_2:
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
    case IST_VolumetricPlasticStrain:
    case IST_Viscosity:
    case IST_DeviatoricStrainMeasure:
    case IST_DeviatoricStressMeasure:
    case IST_vonMisesStress:
    case IST_XFEMEnrichment:
    case IST_XFEMNumIntersecPoints:
    case IST_XFEMLevelSetPhi:
    case IST_Maturity:
    case IST_CrossSectionNumber:
    case IST_CrackWidth:
    case IST_2ndCrackWidth:
    case IST_3rdCrackWidth: 
    case IST_TensileStrength:
    case IST_ResidualTensileStrength:
    case IST_CrackIndex:
    case IST_FiberStressNL:
    case IST_FiberStressLocal:
    case IST_CrackSlip:
    case IST_EquivalentTime:
    case IST_MoistureContent:
    case IST_IncrementCreepModulus:
    case IST_InternalSource:
        return ISVT_SCALAR;

    default:
        return ISVT_UNDEFINED;
    }
}


int giveInternalStateTypeSize(InternalStateValueType valType)
{
    switch ( valType ) {
    case ISVT_TENSOR_S3:
    case ISVT_TENSOR_S3E:
    case ISVT_TENSOR_G:
        return 9;

    case ISVT_VECTOR:
        return 3;

    case ISVT_SCALAR:
        return 1;

    default:
        return 0;
    }
}


InternalStateValueType giveInternalStateValueType(UnknownType type)
{
    if ( type == DisplacementVector || type == EigenVector || type == VelocityVector || type == DirectorField || type == MacroSlipVector || type == ResidualForce ) {
        return ISVT_VECTOR;
    } else if ( type == FluxVector || type == PressureVector || type == Temperature || type == Humidity || type == DeplanationFunction ) {
        return ISVT_SCALAR;
    } else {
        OOFEM_ERROR( "unsupported UnknownType %s", __UnknownTypeToString(type) );
        return ISVT_SCALAR; // To make compiler happy.
    }
}


ContextIOERR :: ContextIOERR(contextIOResultType e, const char *file, int line) :
    error(e),
    msg(nullptr),
    file(file),
    line(line)
{
    this->full_message = "ContextIOERR " + std::to_string(error) + " at line " + std::to_string(line) + " in file \"" + file + "\"";
}

ContextIOERR :: ContextIOERR(contextIOResultType e, const char *msg, const char *file, int line) :
    error(e),
    msg(msg),
    file(file),
    line(line)
{
    this->full_message = "ContextIOERR " + std::to_string(error) + " at line " + std::to_string(line) + " in file \"" + file + "\":" + this->msg;
}


void
ContextIOERR :: print()
{
    if ( msg ) {
        oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, NULL, file, line, 
                                  "ContextIOERR encountered, error code: %d\n%s", error, msg);
    }  else {
        oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, NULL, file, line, 
                                  "ContextIOERR encountered, error code: %d", error);
    }
    OOFEM_EXIT(1);
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

std :: string __DofIDItemToString(DofIDItem _value) {
    if ( _value >= MaxDofID ) {
        char tmp [ 1024 ];
        sprintf(tmp, "X_%d", _value - MaxDofID + 1);
        return tmp;
    }
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

const char *__XFEMStateTypeToString(XFEMStateType _value) {
    TO_STRING_BODY(XFEMStateType_DEF)
}
} // end namespace oofem
