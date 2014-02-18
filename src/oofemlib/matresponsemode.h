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

#ifndef matesponsemode_h
#define matesponsemode_h

#include "enumitem.h"

namespace oofem {
#define MatResponseMode_DEF \
    ENUM_ITEM(TangentStiffness) \
    ENUM_ITEM(SecantStiffness) \
    ENUM_ITEM(ElasticStiffness)                                           \
    ENUM_ITEM(Conductivity)  /* element level conductivity matrix */ \
    ENUM_ITEM(Conductivity_ww) /* material level conductivity submatrix */ \
    ENUM_ITEM(Conductivity_hh) /* material level conductivity submatrix */ \
    ENUM_ITEM(Conductivity_hw) /* material level conductivity submatrix */ \
    ENUM_ITEM(Conductivity_wh) /* material level conductivity submatrix */ \
    ENUM_ITEM(Capacity)                                                   \
    ENUM_ITEM(Capacity_ww) /* material level capacity submatrix */ \
    ENUM_ITEM(Capacity_hh) /* material level capacity submatrix */ \
    ENUM_ITEM(Capacity_hw) /* material level capacity submatrix */ \
    ENUM_ITEM(Capacity_wh) /* material level capacity submatrix */ \
    ENUM_ITEM(IntSource)                                                  \
    ENUM_ITEM(IntSource_ww) /* material level internal source submatrix - water source */ \
    ENUM_ITEM(IntSource_hh) /*  - heat source */ \
    ENUM_ITEM(IntSource_hw) /*  - heat source dependency on water content change */ \
    ENUM_ITEM(IntSource_wh) /*  - water source dependency on temperature change */

/**
 * Describes the character of characteristic material matrix.
 */
enum MatResponseMode {
    MatResponseMode_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__MatResponseModeToString(MatResponseMode _value);
} // end namespace oofem
#endif // matesponsemode_h
