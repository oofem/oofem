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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef stateoperator_h
#define stateoperator_h

#include "enum.h"
#include "field.h"

#include <vector>

namespace oofem {

/**
 * Differential operator applied to a primary field when assembling the generalized state a
 * material consumes.
 *
 * Kept separate from the field identity so that the two compose: adding a field costs one
 * DofIDItem and no operator work, and adding an operator serves every field at once. Encoding the
 * pair as a single name instead would need one name per combination.
 */
#define ENUM_TYPE StateOperator
#define ENUM_DEF \
    ENUM_ITEM_WITH_VALUE(SO_Value, 0)              /* the field itself, N*u */ \
    ENUM_ITEM_WITH_VALUE(SO_Gradient, 1)           /* gradient of a scalar field, one row per spatial direction */ \
    ENUM_ITEM_WITH_VALUE(SO_SymmetricGradient, 2)  /* symmetric gradient of a vector field, i.e. strain */ \
    ENUM_ITEM_WITH_VALUE(SO_Divergence, 3)         /* divergence of a vector field */

#include "enum-impl.h"

/**
 * One entry of the generalized state layout a material advertises: which primary field, and what
 * is taken of it.
 *
 * The field is named by what it physically is, not by where its degrees of freedom happen to be
 * numbered, so a material never obliges a deck to renumber its dofs. FieldType is the core-level
 * vocabulary for this; the mpm module's Variable carries one, and a material in core could not
 * depend on an mpm-local type anyway.
 */
struct StateVariableSpec {
    FieldType field;
    StateOperator op;
};

/// Generalized state layout; see Material::giveStateVariableIDs.
typedef std::vector< StateVariableSpec > StateVariableLayout;

} // end namespace oofem
#endif // stateoperator_h
