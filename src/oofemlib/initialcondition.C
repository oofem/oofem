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

#include "initialcondition.h"
#include "inputrecord.h"
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
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, initialValueDictionary, _IFT_InitialCondition_conditions);

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_InitialCondition_valType);
    valType = ( bcValType ) val;

    ///@todo Make these both not optional (and remove the old approach). Not done right now because it breaks backwards compatibility with input files.
    this->set = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->set, _IFT_InitialCondition_set);
    this->dofIDs.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, this->dofIDs, _IFT_InitialCondition_dofs);

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
