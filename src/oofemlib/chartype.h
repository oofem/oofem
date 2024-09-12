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

#ifndef chartype_h
#define chartype_h

#include "enumitem.h"

namespace oofem {
/**
 * Type representing kind of characteristic value (of scalar, vector or tensorial character) or
 * unknown, which is required, requested, returned, or passed to/from various general services.
 * It typically describes the physical meaning of corresponding component.
 * Typically, many top base classes declare general services for requesting or computing
 * some "characteristic" values of given type. Then only one service for all values of sane type
 * (like vector, scalar) is declared, passing the type of required value (of CharType type) as parameter.
 * Particular implementation based on passed CharType value, usually invokes corresponding specialized services
 * and returns result. If passed CharType value is of unsupported value, error is generated.
 * @see ValueModeType type.
 */
#define CharType_DEF                              \
    ENUM_ITEM_WITH_VALUE(UnknownCharType, 0)                    \
    ENUM_ITEM_WITH_VALUE(StiffnessMatrix, 1)                    \
    ENUM_ITEM_WITH_VALUE(TangentStiffnessMatrix, 2)             \
    ENUM_ITEM_WITH_VALUE(SecantStiffnessMatrix, 3)              \
    ENUM_ITEM_WITH_VALUE(ElasticStiffnessMatrix, 4)             \
    ENUM_ITEM_WITH_VALUE(MassMatrix, 5)                         \
    ENUM_ITEM_WITH_VALUE(LumpedMassMatrix, 6)                   \
    ENUM_ITEM_WITH_VALUE(ConductivityMatrix, 9)                 \
    ENUM_ITEM_WITH_VALUE(CapacityMatrix, 10)                    \
    ENUM_ITEM_WITH_VALUE(InitialStressMatrix, 11)               \
    ENUM_ITEM_WITH_VALUE(LumpedInitialStressMatrix, 12)         \
    /* characteristic vectors */ \
    ENUM_ITEM_WITH_VALUE(ExternalForcesVector, 150)                                       \
    ENUM_ITEM_WITH_VALUE(InternalForcesVector, 151)                                       \
    ENUM_ITEM_WITH_VALUE(LastEquilibratedInternalForcesVector, 152)                       \
    ENUM_ITEM_WITH_VALUE(InertiaForcesVector, 160) \
    /* PFEM */	\
    ENUM_ITEM_WITH_VALUE(AuxVelocityLhs, 200)                                    \
    ENUM_ITEM_WITH_VALUE(VelocityLhs, 201)                                       \
      /*for pressureLhs see CBS */					         \
    ENUM_ITEM_WITH_VALUE(PressureGradientMatrix, 203)                            \
    ENUM_ITEM_WITH_VALUE(DivergenceMatrix, 204)					 \
    ENUM_ITEM_WITH_VALUE(VelocityLaplacianMatrix, 205)                           \
    ENUM_ITEM_WITH_VALUE(PressureLaplacianMatrix, 206)	                         \
    ENUM_ITEM_WITH_VALUE(StabilizationMassMatrix, 207)	                         \
    /* PFEM vectors */  \
    ENUM_ITEM_WITH_VALUE(PressureGradientVector, 208)		                 \
    ENUM_ITEM_WITH_VALUE(MassVelocityVector, 209)		                 \
    ENUM_ITEM_WITH_VALUE(MassAuxVelocityVector, 210)		                 \
    ENUM_ITEM_WITH_VALUE(LaplacePressureVector, 211)		                 \
    ENUM_ITEM_WITH_VALUE(LaplaceVelocityVector, 212)		                 \
    ENUM_ITEM_WITH_VALUE(DivergenceAuxVelocityVector, 213)                       \
    ENUM_ITEM_WITH_VALUE(DivergenceVelocityVector, 214) \
    /* MPM */ \
    ENUM_ITEM_WITH_VALUE(PermeabilityMatrix, 250) \
    /*ENUM_ITEM_WITH_VALUE(FluidMassBalancePressureContribution, 251)*/ \
    /*ENUM_ITEM_WITH_VALUE(BiotConstant, 252) */                        \
    /*ENUM_ITEM_WITH_VALUE(CompressibilityCoefficient, 253)*/           \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_StiffnessMatrix, 254) \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_PressureCouplingMatrix, 255) \
    ENUM_ITEM_WITH_VALUE(MassBalance_PermeabilityMatrix, 256) \
    ENUM_ITEM_WITH_VALUE(MassBalance_CompresibilityMatrix, 257) \
    ENUM_ITEM_WITH_VALUE(MassBalance_StressCouplingMatrix, 258) \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_StressResidual, 259) \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_PressureResidual, 260) \
    ENUM_ITEM_WITH_VALUE(MassBalance_StressRateResidual, 261) \
    ENUM_ITEM_WITH_VALUE(MassBalance_PressureResidual, 262) \
    ENUM_ITEM_WITH_VALUE(MassBalance_PressureRateResidual, 263) \
    /*ENUM_ITEM_WITH_VALUE(FluidViscosity, 264)*/               \
    /* MPM UPV */ \
    ENUM_ITEM_WITH_VALUE(MassBalance_UdotTerm, 265) \
    ENUM_ITEM_WITH_VALUE(MassBalance_VTerm, 266) \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_1, 267) \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_2, 268) \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_3, 269) \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_4, 270) \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_5, 271) \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_2, 272)       \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_3, 273)       \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_4, 274)           \
    ENUM_ITEM_WITH_VALUE(MassBalance_UdotTerm_Residual, 275)    \
    ENUM_ITEM_WITH_VALUE(MassBalance_VTerm_Residual, 276)       \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_1_Residual, 277) \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_2_Residual, 278) \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_3_Residual, 279) \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_4_Residual, 280) \
    ENUM_ITEM_WITH_VALUE(MomentumBalanceMatrix_5_Residual, 281) \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_2_Residual, 282)       \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_3_Residual, 283)       \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_4_Residual, 284)       \
    ENUM_ITEM_WITH_VALUE(InitialVolumeFraction, 285) \
    /* MPM ThermoMechanics*/ \
    ENUM_ITEM_WITH_VALUE(EnergyBalance_ConductivityMatrix, 300) \
    ENUM_ITEM_WITH_VALUE(EnergyBalance_CapacityMatrix, 301) \
    ENUM_ITEM_WITH_VALUE(EnergyBalance_StressCouplingMatrix, 302) \
    ENUM_ITEM_WITH_VALUE(EnergyBalance_Residual, 303) \
    ENUM_ITEM_WITH_VALUE(EnergyBalance_TemperatureResidual, 304) \
    ENUM_ITEM_WITH_VALUE(EnergyBalance_DSigmaDTMatrix, 305) \
    ENUM_ITEM_WITH_VALUE(InternalFluxVector, 306) \
    ENUM_ITEM_WITH_VALUE(MomentumBalance_ThermalCouplingMatrix, 307) \
    ENUM_ITEM_WITH_VALUE(EnergyBalance_ConvectionBCMatrix, 308) \
    ENUM_ITEM_WITH_VALUE(EnergyBalance_ConvectionBCResidual, 309) \
    ENUM_ITEM_WITH_VALUE(EnergyBalance_InternalSourceVector, 310) 


    


enum CharType {
    CharType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__CharTypeToString(CharType _value);
} // end namespace oofem
#endif // chartype_h
