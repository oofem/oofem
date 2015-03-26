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

#include "generalboundarycondition.h"
#include "bcvaltype.h"
#include "function.h"
#include "timestep.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"
#include "domain.h"

namespace oofem {
GeneralBoundaryCondition :: GeneralBoundaryCondition(int n, Domain *d) : FEMComponent(n, d)
{
    timeFunction = 0;
    isImposedTimeFunction = 0;
    set = 0;
}


Function *GeneralBoundaryCondition :: giveTimeFunction()
// Returns the load-time function of the receiver. Reads its number in the
// data file if has not been done yet.
{
    if ( !timeFunction ) {
        OOFEM_ERROR("TimeFunction is not defined");
    }

    return domain->giveFunction(timeFunction);
}


IRResultType
GeneralBoundaryCondition :: initializeFrom(InputRecord *ir)
{
    IRResultType result;           // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, timeFunction, _IFT_GeneralBoundaryCondition_timeFunct);
    if ( timeFunction <= 0 ) {
        OOFEM_WARNING("bad TimeFunction id");
        return IRRT_BAD_FORMAT;
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_GeneralBoundaryCondition_valType);
    valType = ( bcValType ) val;

    dofs.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, dofs, _IFT_GeneralBoundaryCondition_dofs);

    isImposedTimeFunction = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, isImposedTimeFunction, _IFT_GeneralBoundaryCondition_isImposedTimeFunct);

    set = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, set, _IFT_GeneralBoundaryCondition_set);

    return IRRT_OK;
}

bool GeneralBoundaryCondition :: isImposed(TimeStep *tStep)
{
    // Returns a value of isImposedTimeFunction, indicating whether b.c. is imposed or not
    // in given time (nonzero indicates imposed b.c.).

    if ( isImposedTimeFunction ) {
        return ( domain->giveFunction(isImposedTimeFunction)->evaluateAtTime( tStep->giveIntrinsicTime() ) != 0. );
    } else {
        // zero value indicates default behavior -> b.c. is imposed
        // anytime
        return true;
    }
}


void
GeneralBoundaryCondition :: giveInputRecord(DynamicInputRecord &input)
{
    FEMComponent :: giveInputRecord(input);
    input.setField(this->timeFunction, _IFT_GeneralBoundaryCondition_timeFunct);

    if ( ( int ) this->giveBCValType() > 0 ) {
        input.setField(this->giveBCValType(), _IFT_GeneralBoundaryCondition_valType);
    }

    if ( this->giveDofIDs().giveSize() > 0 ) {
        input.setField(this->giveDofIDs(), _IFT_GeneralBoundaryCondition_dofs);
    }

    if ( this->isImposedTimeFunction > 0 ) {
        input.setField(this->isImposedTimeFunction, _IFT_GeneralBoundaryCondition_isImposedTimeFunct);
    }

    if ( this->giveSetNumber() > 0 ) {
        input.setField(this->giveSetNumber(), _IFT_GeneralBoundaryCondition_set);
    }
}

contextIOResultType
GeneralBoundaryCondition :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    if ( mode & CM_Definition ) {
        if ( !stream.write(timeFunction) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}


contextIOResultType
GeneralBoundaryCondition :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    if ( mode & CM_Definition ) {
        if ( !stream.read(timeFunction) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}
} // end namespace oofem
