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

#ifndef unknowntype_h
#define unknowntype_h
//#include "internalstatetype.h"
#include "enumitem.h"

namespace oofem {
#define UnknownType_DEF \
    ENUM_ITEM_WITH_VALUE(DisplacementVector, 1) \
    ENUM_ITEM_WITH_VALUE(GeneralizedDisplacementVector, 2) \
    ENUM_ITEM_WITH_VALUE(FluxVector, 3)                    \
    ENUM_ITEM_WITH_VALUE(VelocityVector, 4)                \
    ENUM_ITEM_WITH_VALUE(PressureVector, 5)                \
    ENUM_ITEM_WITH_VALUE(Temperature, 6)                   \
    ENUM_ITEM_WITH_VALUE(Humidity, 7)                      \
    ENUM_ITEM_WITH_VALUE(EigenVector, 8)                   \
    ENUM_ITEM_WITH_VALUE(DirectorField, 15) /* Vector field */ \
    ENUM_ITEM_WITH_VALUE(DeplanationFunction, 16)
/**
 * Type representing particular unknown (its physical meaning).
 */
enum UnknownType {
    UnknownType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h




const char *__UnknownTypeToString(UnknownType _value);
} // end namespace oofem
#endif // unknowntype_h
