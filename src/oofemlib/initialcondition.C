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
#include "function.h"

namespace oofem {
double InitialCondition :: give(ValueModeType type, const FloatArray& coords)
// Returns the prescribed value of the kinematic unknown 'u'.
{
    if (this->compatibilityMode) {
        char u;
        u =  cltypesGiveUnknownTypeModeKey(type);
        if ( this->hasConditionOn(u) ) {
            return initialValueDictionary.at(u);
        } else {
            return 0.;
        }
    } else {
        int size = coords.giveSize();
        double x = (size > 0) ? coords.at(1) : 0.0;
        double y = (size > 1) ? coords.at(2) : 0.0;
        double z = (size > 2) ? coords.at(3) : 0.0;
        switch (type)
        {
        case VM_Total:
            return this->valueExpr.eval( { { "x", x }, {"y", y}, {"z",z} }, this->giveDomain() );
            break;
        case VM_Velocity:
            return this->velocityExpr.eval( { { "x", x }, {"y", y}, {"z",z} }, this->giveDomain() );
            break;
        case VM_Acceleration:
            return this->accelerationExpr.eval( { { "x", x }, {"y", y}, {"z",z} }, this->giveDomain() );
            break;
        default:
            return 0.0;
        }
    }
}


int InitialCondition :: hasConditionOn(int u)
// Returns True if the receiver submits the unknown 'u' to an initial
// condition, else returns False.
{
    if (this->compatibilityMode) {
        return  ( initialValueDictionary.includes(u) );
    } else {
        return true;
    }
}


int InitialCondition :: hasConditionOn(ValueModeType type)
// Returns True if the receiver submits the unknown 'u' to an initial
// condition, else returns False.
{
    if (this->compatibilityMode) {
        char u = cltypesGiveUnknownTypeModeKey(type);
        return  ( initialValueDictionary.includes(u) );
    } else {
        return true;
    }
}

void InitialCondition :: printYourself()
// Prints the receiver on screen.
{
    printf("Initial condition %d\ninitial values :\n", number);
    initialValueDictionary.printYourself();
}


void
InitialCondition :: initializeFrom(InputRecord &ir)
// Sets up the dictionary where the receiver stores the conditions it
// imposes.
{
    if ( ir.hasField(_IFT_InitialCondition_conditions) ) {
        // compatibility with old input files
        IR_GIVE_FIELD(ir, initialValueDictionary, _IFT_InitialCondition_conditions);
        this->compatibilityMode = true;
    } else {
        this->compatibilityMode = false;
        // new input file format
        valueExpr.setValue(0.0);
        if ( ir.hasField(_IFT_InitialCondition_f) ) {
            IR_GIVE_OPTIONAL_FIELD(ir, valueExpr, _IFT_InitialCondition_f);
        }
        velocityExpr.setValue(0.0);
        if ( ir.hasField(_IFT_InitialCondition_dfdt) ) {
            IR_GIVE_OPTIONAL_FIELD(ir, velocityExpr, _IFT_InitialCondition_dfdt);
        }
        accelerationExpr.setValue(0.0);
        if ( ir.hasField(_IFT_InitialCondition_d2fdt2) ) {
            IR_GIVE_OPTIONAL_FIELD(ir, accelerationExpr, _IFT_InitialCondition_d2fdt2);
        }
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_InitialCondition_valType);
    valType = ( bcValType ) val;

    ///@todo Make these both not optional (and remove the old approach). Not done right now because it breaks backwards compatibility with input files.
    this->set = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->set, _IFT_InitialCondition_set);
    this->dofIDs.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, this->dofIDs, _IFT_InitialCondition_dofs);
}


void
InitialCondition :: scale(ValueModeType type, double s)
{
    if (this->compatibilityMode) {
        if ( this->hasConditionOn(type) ) {
            initialValueDictionary.at(type) *= s;
        }
    } else {
        OOFEM_ERROR ("Not supoorted");
    }
}

} // end namespace oofem
