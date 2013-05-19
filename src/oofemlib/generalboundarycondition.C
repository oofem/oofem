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

#include "generalboundarycondition.h"
#include "valuemodetype.h"
#include "bcvaltype.h"
#include "loadtimefunction.h"
#include "reinforcement.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"

namespace oofem {
GeneralBoundaryCondition :: GeneralBoundaryCondition(int n, Domain *d) : FEMComponent(n, d)
{
    loadTimeFunction = 0;
    isImposedTimeFunction = 0;
    set = 0;
}


LoadTimeFunction *GeneralBoundaryCondition :: giveLoadTimeFunction()
// Returns the load-time function of the receiver. Reads its number in the
// data file if has not been done yet.
{
    if ( !loadTimeFunction ) {
        _error("giveLoadTimeFunction: LoadTimeFunction is not defined");
    }

    return domain->giveLoadTimeFunction(loadTimeFunction);
}


IRResultType
GeneralBoundaryCondition :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;           // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, loadTimeFunction, _IFT_GeneralBoundaryCondition_LoadTimeFunct);
    if ( loadTimeFunction <= 0 ) {
        _error("initializeFrom: bad loadtimefunction id");
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_GeneralBoundaryCondition_valType);
    valType = ( bcValType ) val;

    defaultDofs.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, defaultDofs, _IFT_GeneralBoundaryCondition_defaultDofs);

    isImposedTimeFunction = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, isImposedTimeFunction, _IFT_GeneralBoundaryCondition_IsImposedTimeFunct);

    set = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, set, _IFT_GeneralBoundaryCondition_set);

    return IRRT_OK;
}

bool GeneralBoundaryCondition :: isImposed(TimeStep *tStep)
{
    // Returns a value of isImposedTimeFunction, indicating whether b.c. is imposed or not
    // in given time (nonzero indicates imposed b.c.).

    if ( isImposedTimeFunction ) {
        return ( domain->giveLoadTimeFunction(isImposedTimeFunction)->evaluate(tStep, VM_Total) != 0. );
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
    input.setField(this->loadTimeFunction, _IFT_GeneralBoundaryCondition_LoadTimeFunct);
    input.setField(this->isImposedTimeFunction, _IFT_GeneralBoundaryCondition_IsImposedTimeFunct);
}

contextIOResultType
GeneralBoundaryCondition :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    if ( mode & CM_Definition ) {
        if ( !stream->write(& loadTimeFunction, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}


contextIOResultType
GeneralBoundaryCondition :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    if ( mode & CM_Definition ) {
        if ( !stream->read(& loadTimeFunction, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}

} // end namespace oofem
