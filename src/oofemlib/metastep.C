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

MetaStep :: ~MetaStep()
{
    delete attributes;
}


IRResultType
MetaStep :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, numberOfSteps, _IFT_MetaStep_nsteps);

    delete attributes;
    this->attributes = ir->GiveCopy();

    return IRRT_OK;
}

int
MetaStep :: setStepBounds(int startStepNumber)
{
    sindex = startStepNumber;

    return sindex + numberOfSteps;
}

void
MetaStep :: setNumberOfSteps(int newNumberOfSteps)
{
    this->numberOfSteps = newNumberOfSteps;
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
