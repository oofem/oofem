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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#include "timestepcontroller.h"
#include "datastream.h"
#include "contextioerr.h"
#include "error.h"
#include "classfactory.h"
#include "engngm.h"
#include "function.h"

namespace oofem {


void
TimeStepController :: initializeFrom(InputRecord &ir)
{

  numberOfMetaSteps   = 0;
  IR_GIVE_OPTIONAL_FIELD(ir, numberOfMetaSteps, _IFT_TimeStepController_nmsteps);   
}


void
TimeStepController :: saveContext(DataStream &stream)
{

}


void
TimeStepController :: restoreContext(DataStream &stream)
{

}

 
TimeStep*
TimeStepController :: giveNextStep()
{
    int istep = eModel->giveNumberOfFirstStep();
    double totalTime = 1.;
    StateCounterType counter = 1;
    double dt = this->giveDeltaT();
    if ( !currentStep ) {
      currentStep = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), eModel, 1, 0., dt, 0);
	    //currentStep = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), eModel, 0, -dt, dt, 0);
    }

    auto mStepNum = currentStep->giveMetaStepNumber();
    previousStep = std :: move(currentStep);    
    dt = this->giveCurrentMetaStep()->giveDeltaT(istep, previousStep);
    totalTime = previousStep->giveTargetTime() + dt;
    istep =  previousStep->giveNumber() + 1;      
    counter = previousStep->giveSolutionStateCounter() + 1;

    //
    currentStep = std::make_unique<TimeStep>(istep, eModel, mStepNum, totalTime, dt, counter);
    //
    double intrinsicTime = currentStep->giveTargetTime();
    currentStep->setIntrinsicTime(intrinsicTime);
   
    return currentStep.get();
}


int
TimeStepController :: instanciateMetaSteps(DataReader &dr)
{
    double totalNumberOfSteps;
    // create meta steps
    metaStepList.clear();
    metaStepList.reserve(this->numberOfMetaSteps);
    for ( int i = 1; i <= this->numberOfMetaSteps; i++ ) {
        //MetaStep *mstep = new MetaStep(i, this);
        metaStepList.emplace_back(i, eModel);
    }

    // read problem domains
    for ( int i = 1; i <= this->numberOfMetaSteps; i++ ) {
        auto &ir = dr.giveInputRecord(DataReader :: IR_mstepRec, i);
        metaStepList[i-1].initializeFrom(ir);
	totalNumberOfSteps += metaStepList[i-1].giveNumberOfSteps();
    }

    //this->numberOfSteps = metaStepList.size();
    //this->recomputeFinalTime(totalNumberOfSteps);
    return 1;
}


MetaStep *
TimeStepController :: giveMetaStep(int i)
{
    if ( ( i > 0 ) && ( i <= this->numberOfMetaSteps ) ) {
        return &this->metaStepList[i-1];
    } else {
        OOFEM_ERROR("undefined metaStep (%d)", i);
    }

    return NULL;
}


int
TimeStepController :: instanciateDefaultMetaStep(InputRecord &ir)
{

  metaStepList.clear();
  this->numberOfMetaSteps = 1;
  metaStepList.reserve(1);
  //
  metaStepList.emplace_back(1, eModel);
  metaStepList[0].initializeFrom(ir);
  return 1;
}


void
TimeStepController :: initMetaStepAttributes(MetaStep *mStep)
{
    // update attributes
    eModel->updateAttributes(mStep); // virtual function
    // finish data acquiring
    mStep->giveAttributesRecord().finish();
}
  




void
TimeStepController :: postInitialize()
{
    int istep = eModel->giveNumberOfFirstStep(true);
    for ( auto &metaStep: metaStepList ) {
        istep = metaStep.setStepBounds(istep);
    }


}

void
TimeStepController :: reduceTimeStep()
{
  MetaStep *mS = this->giveCurrentMetaStep();
  mS->reduceTimeStep();
  this->currentStep->setTargetTime(this->currentStep->giveTargetTime()-this->currentStep->giveTimeIncrement()+mS->giveDeltaT());
  this->currentStep->setTimeIncrement(mS->giveDeltaT());
  
}

void
TimeStepController :: adaptTimeStep(int nIter, double targetTime)
{
  this->giveCurrentMetaStep()->adaptTimeStep(nIter, targetTime);
}




  
  
} // end namespace oofem
