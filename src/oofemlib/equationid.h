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

#ifndef equationid_h
#define equationid_h
#include "enumitem.h"

namespace oofem {
//EID_AuxMomentumBalance - Auxiliary equation for characteristic-split based methods
//EID_MomentumBalance_ConservationEquation - Coupled system

#define EquationID_DEF \
    ENUM_ITEM_WITH_VALUE(EID_Undefined, 0) \
    ENUM_ITEM_WITH_VALUE(EID_MomentumBalance, 1) \
    ENUM_ITEM_WITH_VALUE(EID_ConservationEquation, 3) \
    ENUM_ITEM_WITH_VALUE(EID_MomentumBalance_ConservationEquation, 4) \

/**
 * This type identifies the governing equation.
 */
enum EquationID {
    EquationID_DEF
};

// enum EquationID {
//     EID_MomentumBalance, ///< Momentum balance equation.
//     EID_ConservationEquation, ///< Conservation equation.
//     EID_MomentumBalance_ConservationEquation, ///< Coupled momentum balance and conservation equation.
// };
#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h

const char *__EquationIDToString(EquationID _value);
} // end namespace oofem
#endif // equationid_h
