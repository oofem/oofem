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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

// Milan ?????????????????
//#include "gpinitmodule.h"
// Milan ?????????????????

#include "dummyengngm.h"
#include "timestep.h"
#include "classfactory.h"


namespace oofem {
REGISTER_EngngModel(DummyEngngModel);

DummyEngngModel :: DummyEngngModel(int i, EngngModel *_master) : EngngModel (i, _master)
{
    ndomains = 1;
}

void
DummyEngngModel :: initializeFrom(InputRecord &ir)
{
    this->numberOfSteps = 1;
    this->nMetaSteps   = 0;
    this->suppressOutput = true;
    
}

TimeStep *DummyEngngModel :: giveNextStep()
{
    if ( !currentStep ) {
        // first step -> generate initial step
        //currentStep = std::make_unique<TimeStep>(*giveSolutionStepWhenIcApply());
        currentStep = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 1, 0., 1., 0);
    }
    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(*previousStep, 1.);

    return currentStep.get();
}

void DummyEngngModel :: solveYourselfAt(TimeStep *tStep)
{

}
} // end namespace oofem
