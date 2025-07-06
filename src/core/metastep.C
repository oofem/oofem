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
MetaStep :: MetaStep(int n, EngngModel *e) :
    eModel(e),
    numberOfSteps(0),
    number(n)
{}

MetaStep :: MetaStep(int n, EngngModel *e, int nsteps, InputRecord &attrib) :
    eModel(e),
    numberOfSteps(nsteps),
    attributes(attrib.clone()),
    number(n)
{}


void
MetaStep :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, numberOfSteps, _IFT_MetaStep_nsteps);

    this->attributes = ir.clone();
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
    return solStepNumber >= sindex && solStepNumber < ( sindex + numberOfSteps );
}
} // end namespace oofem
