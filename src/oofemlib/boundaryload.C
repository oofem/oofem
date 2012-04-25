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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

#include "boundaryload.h"
#include "loadtime.h"
#include "flotarry.h"
#include "timestep.h"

namespace oofem {
void
BoundaryLoad :: computeComponentArrayAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)  {
    // returns component array for elements which use direct formulae
    Load :: computeComponentArrayAt(answer, tStep, mode);
}


void
BoundaryLoad :: computeValueAt(FloatArray &answer, TimeStep *stepN, FloatArray &coords, ValueModeType mode)
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
     * factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
     * if ((mode==VM_Incremental) && (!stepN->isTheFirstStep()))
     * //factor -= this->giveLoadTimeFunction()->at(stepN->givePreviousStep()->giveTime()) ;
     * factor -= this->giveLoadTimeFunction()->at(stepN->giveTime()-stepN->giveTimeIncrement());
     */
    factor = this->giveLoadTimeFunction()->evaluate(stepN, mode);

    answer.times(factor);
}


IRResultType
BoundaryLoad :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = Load :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        IR_IOERR(giveClassName(), __proc, IFT_Unknown, "", ir, result);
    }

    IR_GIVE_FIELD(ir, nDofs, IFT_BoundaryLoad_ndofs, "ndofs"); // Macro

    int value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_BoundaryLoad_loadtype, "loadtype"); // Macro
    lType = ( bcType ) value;

    value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_BoundaryLoad_cstype, "cstype"); // Macro
    coordSystemType = ( BL_CoordSystType ) value;

    IR_GIVE_OPTIONAL_FIELD(ir, propertyDictionary, IFT_BoundaryLoad_properties, "properties"); // Macro

    return IRRT_OK;
}


int
BoundaryLoad :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];


    Load :: giveInputRecordString(str, keyword);
    sprintf(buff, " ndofs %d loadtype %d cstype %d", this->nDofs, ( int ) this->lType, ( int ) this->coordSystemType);
    str += buff;

    if ( this->propertyDictionary.giveSize() != 0 ) {
        std :: string helpStr;

        sprintf( buff, " properties %d", this->propertyDictionary.giveSize() );
        str += buff;
        this->propertyDictionary.formatAsString(helpStr);
        str += helpStr;
    }

    return 1;
}


double
BoundaryLoad :: giveProperty(int aProperty)
// Returns the value of the property aProperty (e.g. the area
// 'A') of the receiver.
{
    if ( propertyDictionary.includes(aProperty) ) {
        return propertyDictionary.at(aProperty);
    } else {
        _error("give: property not defined");
    }

    return 0.0;
}
} // end namespace oofem
