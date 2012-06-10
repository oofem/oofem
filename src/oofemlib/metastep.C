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

#include "metastep.h"

namespace oofem {
MetaStep :: MetaStep(int n, EngngModel *e)
{
    this->number = n;
    this->eModel = e;
    this->numberOfSteps = 0;
    this->attributes = NULL;
}

MetaStep :: MetaStep(int n, EngngModel *e, int nsteps, InputRecord &attrib)
{
    this->number = n;
    this->eModel = e;
    this->numberOfSteps = nsteps;
    this->attributes = attrib.GiveCopy();
}


IRResultType
MetaStep :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, numberOfSteps, IFT_MetaStep_nsteps, "nsteps"); // Macro

    if ( attributes ) {
        delete attributes;
    }

    this->attributes = ir->GiveCopy();
    /*
     * this->readQuotedString (initString, "attributes", this->attributes, MetaStepAttrRecLenght);
     */
    return IRRT_OK;
}

int
MetaStep :: setStepBounds(int startStepNumber)
{
    sindex = startStepNumber;

    return sindex + numberOfSteps;
}

void
MetaStep :: setNumberOfSteps(int numberOfSteps)
{
    this->numberOfSteps = numberOfSteps;
}

int
MetaStep :: isStepValid(int solStepNumber)
{
    if ( ( solStepNumber >= sindex ) &&
        ( solStepNumber < ( sindex + numberOfSteps ) ) ) {
        return 1;
    }

    return 0;
}
} // end namespace oofem
