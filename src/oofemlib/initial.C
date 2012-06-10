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

#include "initial.h"
#include "cltypes.h"

namespace oofem {

double InitialCondition :: give(ValueModeType type)
// Returns the prescribed value of the kinematic unknown 'u'.
{
    char u;
    u =  cltypesGiveUnknownTypeModeKey(type);
    if ( this->hasConditionOn(u) ) {
        return initialValueDictionary.at(u);
    } else {
        return 0.;
    }
}


int InitialCondition :: hasConditionOn(int u)
// Returns True if the receiver submits the unknown 'u' to an initial
// condition, else returns False.
{
    return  ( initialValueDictionary.includes(u) );
}


int InitialCondition :: hasConditionOn(ValueModeType type)
// Returns True if the receiver submits the unknown 'u' to an initial
// condition, else returns False.
{
    char u = cltypesGiveUnknownTypeModeKey(type);
    return  ( initialValueDictionary.includes(u) );
}

void InitialCondition :: printYourself()
// Prints the receiver on screen.
{
    printf("Initial condition %d\ninitial values :\n", number);
    initialValueDictionary.printYourself();
}


IRResultType
InitialCondition :: initializeFrom(InputRecord *ir)
// Sets up the dictionary where the receiver stores the conditions it
// imposes.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, initialValueDictionary, IFT_InitialCondition_conditions, "conditions"); // Macro

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_GeneralBoundaryCondition_valType, "valtype"); // Macro
    valType = ( bcValType ) val;

    return IRRT_OK;
}


void
InitialCondition :: scale(ValueModeType type, double s)
{
    if ( this->hasConditionOn(type) ) {
        initialValueDictionary.at(type) *= s;
    }
}
} // end namespace oofem
