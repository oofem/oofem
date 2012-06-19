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
#define CharType_DEF                            \
    ENUM_ITEM(UnknownCharType)                    \
    ENUM_ITEM(StiffnessMatrix)                    \
    ENUM_ITEM(TangentStiffnessMatrix)             \
    ENUM_ITEM(SecantStiffnessMatrix)              \
    ENUM_ITEM(ElasticStiffnessMatrix)             \
    ENUM_ITEM(MassMatrix)                         \
    ENUM_ITEM(LumpedMassMatrix)                   \
    ENUM_ITEM(EffectiveStiffnessMatrix)           \
    ENUM_ITEM(ConductivityMatrix)                 \
    ENUM_ITEM(CapacityMatrix)                     \
    ENUM_ITEM(InitialStressMatrix)                \
    ENUM_ITEM(HeatAndMoistureCharMatrix)          \
    /* CBS */ \
    ENUM_ITEM(IntermediateConvectionTerm)         \
    ENUM_ITEM(IntermediateDiffusionTerm)          \
    ENUM_ITEM(DensityRhsVelocityTerms)            \
    ENUM_ITEM(DensityRhsPressureTerms)            \
    ENUM_ITEM(DensityPrescribedTractionPressure)          \
    ENUM_ITEM(NumberOfNodalPrescribedTractionPressureContributions)       \
    ENUM_ITEM(PressureLhs)                                                \
    ENUM_ITEM(CorrectionRhs)                                              \
    ENUM_ITEM(CriticalTimeStep)                                           \
    ENUM_ITEM(PrescribedVelocityRhsVector)                                \
    ENUM_ITEM(PrescribedDensityRhsVector)                                 \
    /* SUPG/PSPG */ \
    ENUM_ITEM(AccelerationTerm_MB)                                        \
    ENUM_ITEM(AdvectionDerivativeTerm_MB)                                 \
    ENUM_ITEM(DiffusionDerivativeTerm_MB)                                 \
    ENUM_ITEM(SecantDiffusionDerivativeTerm_MB)                           \
    ENUM_ITEM(TangentDiffusionDerivativeTerm_MB)                          \
    ENUM_ITEM(InitialDiffusionDerivativeTerm_MB)                          \
    ENUM_ITEM(PressureTerm_MB)                                            \
    ENUM_ITEM(LSICStabilizationTerm_MB)                                   \
    ENUM_ITEM(LinearAdvectionTerm_MC)                                     \
    ENUM_ITEM(AdvectionTerm_MC)                                           \
    ENUM_ITEM(AdvectionDerivativeTerm_MC)                                 \
    ENUM_ITEM(AccelerationTerm_MC)                                        \
    ENUM_ITEM(DiffusionDerivativeTerm_MC)                                 \
    ENUM_ITEM(DiffusionTerm_MC)                                           \
    ENUM_ITEM(PressureTerm_MC)                                            \
    ENUM_ITEM(BCLhsTerm_MB)                                               \
    ENUM_ITEM(BCLhsPressureTerm_MB)                                       \
    ENUM_ITEM(BCRhsTerm_MB)                                               \
    ENUM_ITEM(BCRhsTerm_MC)                                               \
    ENUM_ITEM(AlgorithmicRhsTerm_MB)                                      \
    ENUM_ITEM(AlgorithmicRhsTerm_MC)                                      \
    ENUM_ITEM(AdvectionTerm_MB)                                           \
    ENUM_ITEM(DiffusionTerm_MB)                                           \
    /* characteristic vectors */ \
    ENUM_ITEM(ExternalForcesVector)                                       \
    ENUM_ITEM(InternalForcesVector)                                       \
    ENUM_ITEM(LastEquilibratedInternalForcesVector)                       \
    ENUM_ITEM(ElementBCTransportVector)                                   \
    ENUM_ITEM(ElementInternalSourceVector)                                \
    ENUM_ITEM(LHSBCMatrix) /* LHS due to Boundary Conditions (Transport problems) */  \
    ENUM_ITEM(NSTP_MidpointLhs) /* NonStationaryTransportProblem - LHS for midpoint discretization alg. */ \
    ENUM_ITEM(NSTP_MidpointRhs) /* NonStationaryTransportProblem - RHS for midpoint discretization alg. */ \
    ENUM_ITEM(IntSourceLHSMatrix) /* LHS due to material internal source (Transport problems) */  \
    ENUM_ITEM(PrescribedRhsVector)

enum CharType {
    CharType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__CharTypeToString(CharType _value);
} // end namespace oofem
#endif // chartype_h
