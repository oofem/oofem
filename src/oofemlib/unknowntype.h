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
// FILE: unknowntype.h
//

#ifndef unknowntype_h
#define unknowntype_h

#include "enumitem.h"

/**
 * Type representing particular unknown (its physical meaning).
 */

#define UnknownType_DEF \
    ENUM_ITEM(UnknownType_Unknown) \
    ENUM_ITEM(DisplacementVector) \
    ENUM_ITEM(GeneralizedDisplacementVector) \
    ENUM_ITEM(FluxVector) \
    ENUM_ITEM(VelocityVector)                     \
    ENUM_ITEM(PressureVector)                     \
    ENUM_ITEM(TemperatureVector)                  \
    ENUM_ITEM(EigenValue)                         \
    ENUM_ITEM(EigenVector)                        \
    ENUM_ITEM(TotalLoadLevel)                     \
    ENUM_ITEM(ReynoldsNumber)                                             \
    ENUM_ITEM(Theta_1) /* CBS integration constan)*/ \
    ENUM_ITEM(Theta_2) /* CBS integration constan)*/ \
    ENUM_ITEM(PrescribedTractionPressure) /* CBS prescribed pressure due to applied tractio)*/

enum UnknownType {
    UnknownType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


char *__UnknownTypeToString(UnknownType _value);


#endif // unknowntype_h

