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
    ENUM_ITEM_WITH_VALUE(UnknownCharType,0)                    \
    ENUM_ITEM_WITH_VALUE(StiffnessMatrix,1)                    \
    ENUM_ITEM_WITH_VALUE(TangentStiffnessMatrix,2)             \
    ENUM_ITEM_WITH_VALUE(SecantStiffnessMatrix,3)              \
    ENUM_ITEM_WITH_VALUE(ElasticStiffnessMatrix,4)             \
    ENUM_ITEM_WITH_VALUE(MassMatrix,5)                         \
    ENUM_ITEM_WITH_VALUE(LumpedMassMatrix,6)                   \
    ENUM_ITEM_WITH_VALUE(EffectiveStiffnessMatrix,7)           \
    ENUM_ITEM_WITH_VALUE(EffectiveMassMatrix,8)                \
    ENUM_ITEM_WITH_VALUE(ConductivityMatrix,9)                 \
    ENUM_ITEM_WITH_VALUE(CapacityMatrix,10)                    \
    ENUM_ITEM_WITH_VALUE(InitialStressMatrix,11)               \
    ENUM_ITEM_WITH_VALUE(HeatAndMoistureCharMatrix,12)         \
    /* CBS */                                                  \
    ENUM_ITEM_WITH_VALUE(IntermediateConvectionTerm,50)        \
    ENUM_ITEM_WITH_VALUE(IntermediateDiffusionTerm,51)         \
    ENUM_ITEM_WITH_VALUE(DensityRhsVelocityTerms,52)           \
    ENUM_ITEM_WITH_VALUE(DensityRhsPressureTerms,53)           \
    ENUM_ITEM_WITH_VALUE(DensityPrescribedTractionPressure,54) \
    ENUM_ITEM_WITH_VALUE(NumberOfNodalPrescribedTractionPressureContributions,55)   \
    ENUM_ITEM_WITH_VALUE(PressureLhs,56)                                            \
    ENUM_ITEM_WITH_VALUE(CorrectionRhs,57)                                          \
    ENUM_ITEM_WITH_VALUE(CriticalTimeStep,58)                                       \
    ENUM_ITEM_WITH_VALUE(PrescribedVelocityRhsVector,59)                            \
    ENUM_ITEM_WITH_VALUE(PrescribedDensityRhsVector,60)                             \
    /* SUPG/PSPG */                                                                      \
    ENUM_ITEM_WITH_VALUE(AccelerationTerm_MB,100)                                        \
    ENUM_ITEM_WITH_VALUE(AdvectionDerivativeTerm_MB,101)                                 \
    ENUM_ITEM_WITH_VALUE(DiffusionDerivativeTerm_MB,102)                                 \
    ENUM_ITEM_WITH_VALUE(SecantDiffusionDerivativeTerm_MB,103)                           \
    ENUM_ITEM_WITH_VALUE(TangentDiffusionDerivativeTerm_MB,104)                          \
    ENUM_ITEM_WITH_VALUE(InitialDiffusionDerivativeTerm_MB,105)                          \
    ENUM_ITEM_WITH_VALUE(PressureTerm_MB,106)                                            \
    ENUM_ITEM_WITH_VALUE(LSICStabilizationTerm_MB,107)                                   \
    ENUM_ITEM_WITH_VALUE(LinearAdvectionTerm_MC,108)                                     \
    ENUM_ITEM_WITH_VALUE(AdvectionTerm_MC,109)                                           \
    ENUM_ITEM_WITH_VALUE(AdvectionDerivativeTerm_MC,110)                                 \
    ENUM_ITEM_WITH_VALUE(AccelerationTerm_MC,111)                                        \
    ENUM_ITEM_WITH_VALUE(DiffusionDerivativeTerm_MC,112)                                 \
    ENUM_ITEM_WITH_VALUE(DiffusionTerm_MC,113)                                           \
    ENUM_ITEM_WITH_VALUE(PressureTerm_MC,114)                                            \
    ENUM_ITEM_WITH_VALUE(BCLhsTerm_MB,115)                                               \
    ENUM_ITEM_WITH_VALUE(BCLhsPressureTerm_MB,116)                                       \
    ENUM_ITEM_WITH_VALUE(BCLhsPressureTerm_MC,117)                                       \
    ENUM_ITEM_WITH_VALUE(BCRhsTerm_MB,118)                                               \
    ENUM_ITEM_WITH_VALUE(BCRhsTerm_MC,119)                                               \
    ENUM_ITEM_WITH_VALUE(AlgorithmicRhsTerm_MB,120)                                      \
    ENUM_ITEM_WITH_VALUE(AlgorithmicRhsTerm_MC,121)                                      \
    ENUM_ITEM_WITH_VALUE(AdvectionTerm_MB,122)                                           \
    ENUM_ITEM_WITH_VALUE(DiffusionTerm_MB,123)                                           \
    /* characteristic vectors */                                                         \
    ENUM_ITEM_WITH_VALUE(ExternalForcesVector,150)                                       \
    ENUM_ITEM_WITH_VALUE(InternalForcesVector,151)                                       \
    ENUM_ITEM_WITH_VALUE(LastEquilibratedInternalForcesVector,152)                       \
    ENUM_ITEM_WITH_VALUE(ElementBCTransportVector,153)                                   \
    ENUM_ITEM_WITH_VALUE(ElementInternalSourceVector,154)                                \
    ENUM_ITEM_WITH_VALUE(LHSBCMatrix,155) /* LHS due to Boundary Conditions (Transport problems) */                       \
    ENUM_ITEM_WITH_VALUE(NSTP_MidpointLhs,156) /* NonStationaryTransportProblem - LHS for midpoint discretization alg. */ \
    ENUM_ITEM_WITH_VALUE(NSTP_MidpointRhs,157) /* NonStationaryTransportProblem - RHS for midpoint discretization alg. */ \
    ENUM_ITEM_WITH_VALUE(IntSourceLHSMatrix,158) /* LHS due to material internal source (Transport problems) */           \
    ENUM_ITEM_WITH_VALUE(PrescribedRhsVector,159)

enum CharType {
    CharType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__CharTypeToString(CharType _value);
} // end namespace oofem
#endif // chartype_h
