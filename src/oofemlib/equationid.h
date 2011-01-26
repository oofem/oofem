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
// FILE: equationid.h
//

#ifndef equationid_h
#define equationid_h
#include "enumitem.h"

namespace oofem {
/**
 * This type identifies the governing equation
 */
//EID_AuxMomentumBalance - Auxiliary equation for characteristic-split based methods
//EID_MomentumBalance_ConservationEquation - Coupled system

#define EquationID_DEF \
    ENUM_ITEM_WITH_VALUE(EID_Undefined, 0) \
    ENUM_ITEM_WITH_VALUE(EID_MomentumBalance, 1) \
    ENUM_ITEM_WITH_VALUE(EID_AuxMomentumBalance, 2) \
    ENUM_ITEM_WITH_VALUE(EID_ConservationEquation, 3) \
    ENUM_ITEM_WITH_VALUE(EID_MomentumBalance_ConservationEquation, 4) \

enum EquationID {
    EquationID_DEF
};

// enum EquationID {
//     EID_MomentumBalance,
//     EID_AuxMomentumBalance,
//     EID_ConservationEquation,
//     EID_MomentumBalance_ConservationEquation,
// };
#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h

const char *__EquationIDToString(EquationID _value);
} // end namespace oofem
#endif // equationid_h

