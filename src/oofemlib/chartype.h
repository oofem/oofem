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
    ENUM_ITEM_WITH_VALUE(DivergenceVelocityVector, 214)


enum CharType {
    CharType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__CharTypeToString(CharType _value);
} // end namespace oofem
#endif // chartype_h
