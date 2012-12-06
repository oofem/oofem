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

#ifndef unknowntype_h
#define unknowntype_h

#include "enumitem.h"

namespace oofem {

#define UnknownType_DEF \
    ENUM_ITEM_WITH_VALUE(UnknownType_Unknown, 0) \
    ENUM_ITEM_WITH_VALUE(DisplacementVector, 1) \
    ENUM_ITEM_WITH_VALUE(GeneralizedDisplacementVector, 2) \
    ENUM_ITEM_WITH_VALUE(FluxVector, 3) \
    ENUM_ITEM_WITH_VALUE(VelocityVector, 4)                     \
    ENUM_ITEM_WITH_VALUE(PressureVector, 5)                     \
    ENUM_ITEM_WITH_VALUE(TemperatureVector, 6)                  \
    ENUM_ITEM_WITH_VALUE(EigenValue, 7)                         \
    ENUM_ITEM_WITH_VALUE(EigenVector, 8)                        \
    ENUM_ITEM_WITH_VALUE(TotalLoadLevel, 9)                     \
    ENUM_ITEM_WITH_VALUE(ReynoldsNumber, 10)                                             \
    ENUM_ITEM_WITH_VALUE(Theta_1, 11) /* CBS integration constant)*/ \
    ENUM_ITEM_WITH_VALUE(Theta_2, 12) /* CBS integration constant)*/ \
    ENUM_ITEM_WITH_VALUE(PrescribedTractionPressure, 13) /* CBS prescribed pressure due to applied traction)*/ \
    ENUM_ITEM_WITH_VALUE(InternalForcesEBENorm, 14)  /* Norm of nodal internal forces evaluated on element by element basis*/ \
    ENUM_ITEM_WITH_VALUE(DirectorField, 15) /* Vector field */
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
