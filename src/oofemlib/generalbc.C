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

#include "generalbc.h"
#include "bcvaltype.h"

namespace oofem {
GeneralBoundaryCondition :: GeneralBoundaryCondition(int n, Domain *d) : FEMComponent(n, d)
{
    loadTimeFunction = 0;
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

    IR_GIVE_FIELD(ir, loadTimeFunction, IFT_GeneralBoundaryCondition_LoadTimeFunct, "loadtimefunction"); // Macro
    if ( loadTimeFunction <= 0 ) {
        _error("initializeFrom: bad loadtimefunction id");
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_GeneralBoundaryCondition_valType, "valtype"); // Macro
    valType = ( bcValType ) val;

    IR_GIVE_OPTIONAL_FIELD(ir, defaultDofs, IFT_GeneralBoundaryCondition_defaultDofs, "defaultdofs"); // Macro

    return IRRT_OK;
}


int
GeneralBoundaryCondition :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    FEMComponent :: giveInputRecordString(str, keyword);
    sprintf(buff, " loadtimefunction %d", this->loadTimeFunction);
    str += buff;

    return 1;
}
} // end namespace oofem
