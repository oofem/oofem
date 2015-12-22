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
#include "domain.h"

namespace oofem {

BoundaryLoad :: BoundaryLoad(int i, Domain * d) : Load(i, d), coordSystemType(CST_Global){

}

void
BoundaryLoad :: computeComponentArrayAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    // returns component array for elements which use direct formulae
    Load :: computeComponentArrayAt(answer, tStep, mode);
}


void
BoundaryLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    // Evaluates the value at specific integration point
    int nSize, nDofs;
    double factor;
    FloatArray N;

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("unknown mode");
    }

    this->computeNArray(N, coords);
    nSize = N.giveSize();

    nDofs = this->componentArray.giveSize() / nSize;

    answer.resize(nDofs);

    for ( int i = 1; i <= nDofs; i++ ) {
        double value = 0.;
        for ( int j = 1; j <= nSize; j++ ) {
            value += N.at(j) * this->componentArray.at(i + ( j - 1 ) * nDofs);
        }

        answer.at(i) = value;
    }

    // time distribution
    factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer.times(factor);
}


IRResultType
BoundaryLoad :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = Load :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    int dummy;
    IR_GIVE_OPTIONAL_FIELD(ir, dummy, "ndofs");

    int value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_BoundaryLoad_loadtype);
    lType = ( bcType ) value;

    value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_BoundaryLoad_cstype);
    coordSystemType = ( CoordSystType ) value;

    IR_GIVE_OPTIONAL_FIELD(ir, propertyDictionary, _IFT_BoundaryLoad_properties);
    IR_GIVE_OPTIONAL_FIELD(ir, propertyTimeFunctDictionary, _IFT_BoundaryLoad_propertyTimeFunctions);

    IR_GIVE_OPTIONAL_FIELD(ir, propertyMultExpr, _IFT_BoundaryLoad_propertyMultExpr);

    return IRRT_OK;
}


void
BoundaryLoad :: giveInputRecord(DynamicInputRecord &input)
{
    Load :: giveInputRecord(input);
    input.setField(this->lType, _IFT_BoundaryLoad_loadtype);
    input.setField(this->coordSystemType, _IFT_BoundaryLoad_cstype);
    input.setField(this->propertyDictionary, _IFT_BoundaryLoad_properties);
    input.setField(this->propertyTimeFunctDictionary, _IFT_BoundaryLoad_propertyTimeFunctions);
    input.setField(this->propertyMultExpr, _IFT_BoundaryLoad_propertyMultExpr);
}

void
BoundaryLoad :: setVariableState(int aVariable, double val)
{
    variableState.at(aVariable) = val;
}

double
BoundaryLoad :: giveVariableState(int aVariable)
{
    return variableState.at(aVariable);
}

double
BoundaryLoad :: giveProperty(int aProperty, TimeStep *tStep)
// Returns the value of the property aProperty (e.g. the area
// 'A') of the receiver.
{
    double answer;
    if ( propertyDictionary.includes(aProperty) ) {
        // check if time fuction registered under the same key
        if ( propertyTimeFunctDictionary.includes(aProperty) ) {
            answer = propertyDictionary.at(aProperty) * domain->giveFunction( (int)propertyTimeFunctDictionary.at(aProperty) )->evaluate(tStep, VM_Total);
        } else {
            answer = propertyDictionary.at(aProperty);
        }
    } else {
        OOFEM_ERROR("Property '%c' not defined", (char)aProperty);
    }

    if ( propertyMultExpr.isDefined() ){
        double x = giveVariableState('x');
        answer *= propertyMultExpr.eval( { { "x", x } }, this->giveDomain() );
    }
    return answer;
}
} // end namespace oofem
