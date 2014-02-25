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

#include "boundaryload.h"
#include "function.h"
#include "floatarray.h"
#include "timestep.h"
#include "dynamicinputrecord.h"
#include "valuemodetype.h"

namespace oofem {
void
BoundaryLoad :: computeComponentArrayAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    // returns component array for elements which use direct formulae
    Load :: computeComponentArrayAt(answer, tStep, mode);
}


void
BoundaryLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode)
{
    // Evaluates the value at specific integration point
    int i, j, nSize;
    double value, factor;
    FloatArray N;

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        _error("computeValueAt: unknown mode");
    }

    answer.resize(this->nDofs);

    this->computeNArray(N, coords);
    nSize = N.giveSize();

    if ( ( this->componentArray.giveSize() / nSize ) != nDofs ) {
        _error("computeValueAt: componentArray size mismatch");
    }

    for ( i = 1; i <= nDofs; i++ ) {
        for ( value = 0., j = 1; j <= nSize; j++ ) {
            value += N.at(j) * this->componentArray.at(i + ( j - 1 ) * nDofs);
        }

        answer.at(i) = value;
    }

    // time distribution

    /*
     * factor = this -> giveTimeFunction() -> at(tStep->giveTime()) ;
     * if ((mode==VM_Incremental) && (!tStep->isTheFirstStep()))
     * //factor -= this->giveTimeFunction()->at(tStep->givePreviousStep()->giveTime()) ;
     * factor -= this->giveTimeFunction()->at(tStep->giveTime()-tStep->giveTimeIncrement());
     */
    factor = this->giveTimeFunction()->evaluate(tStep, mode);

    answer.times(factor);
}


IRResultType
BoundaryLoad :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = Load :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, nDofs, _IFT_BoundaryLoad_ndofs);

    int value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_BoundaryLoad_loadtype);
    lType = ( bcType ) value;

    value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_BoundaryLoad_cstype);
    coordSystemType = ( CoordSystType ) value;

    IR_GIVE_OPTIONAL_FIELD(ir, propertyDictionary, _IFT_BoundaryLoad_properties);
    IR_GIVE_OPTIONAL_FIELD(ir, propertyTimeFunctDictionary, _IFT_BoundaryLoad_propertyTimeFunctions);

    return result;
}


void
BoundaryLoad :: giveInputRecord(DynamicInputRecord &input)
{
    Load :: giveInputRecord(input);
    input.setField(this->nDofs, _IFT_BoundaryLoad_ndofs);
    input.setField(this->lType, _IFT_BoundaryLoad_loadtype);
    input.setField(this->coordSystemType, _IFT_BoundaryLoad_cstype);
    input.setField(this->propertyDictionary, _IFT_BoundaryLoad_properties);
    input.setField(this->propertyTimeFunctDictionary, _IFT_BoundaryLoad_propertyTimeFunctions);
}


double
BoundaryLoad :: giveProperty(int aProperty, TimeStep *tStep)
// Returns the value of the property aProperty (e.g. the area
// 'A') of the receiver.
{
    if ( propertyDictionary.includes(aProperty) ) {
        // check if time fuction registered under the same key
        if ( propertyTimeFunctDictionary.includes(aProperty) ) {
            return propertyDictionary.at(aProperty) * domain->giveFunction( propertyTimeFunctDictionary.at(aProperty) )->evaluate(tStep, VM_Total);
        } else {
            return propertyDictionary.at(aProperty);
        }
    } else {
        _error("give: property not defined");
    }

    return 0.0;
}
} // end namespace oofem
