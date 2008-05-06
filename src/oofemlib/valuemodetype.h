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
// FILE: valuemodetype.h
//

#ifndef valuemodetype_h
#define valuemodetype_h

#include "enumitem.h"

//
// following mode determines the mode of particular unknown
// which can be requested on DOF.
// particular DOF contain for example displacement type unknown,
// but we can request total value, increment of value or velocity or acceleration of
// this unknown. This has been done mainly in  order to  improve runtime checking
// ad Dof level.
//
// see also isUnknownModeIncrementalMode() function (cltypes.C)
// when adding new ValueModeType mode.
/**
 * Type representing the mode of UnknownType or CharType, or similar types.
 * Afore mentioned types usually describes the physical meaning of
 * value and ValueModeType provides the further necessary classification. For example "DisplacementVector"
 * value can be futher classified to be total displacement (TotalMode) or  velocity of
 * displacement (VelocityMode) an so on.
 */
#define ValueModeType_DEF \
    ENUM_ITEM_WITH_VALUE(VM_Unknown, 0) \
    ENUM_ITEM_WITH_VALUE(VM_Total, 1)              \
    ENUM_ITEM_WITH_VALUE(VM_Velocity, 2)           \
    ENUM_ITEM_WITH_VALUE(VM_Acceleration, 3)       \
    ENUM_ITEM_WITH_VALUE(VM_Incremental, 4)        \
    ENUM_ITEM_WITH_VALUE(VM_RhsTotal, 5)           \
    ENUM_ITEM_WITH_VALUE(VM_RhsIncremental, 6)     \
    ENUM_ITEM_WITH_VALUE(VM_RhsInitial, 7)

enum ValueModeType {
    ValueModeType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


char *__ValueModeTypeToString(ValueModeType _value);


#endif // valuemodetype_h

