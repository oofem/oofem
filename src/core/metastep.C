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

#include "metastep.h"
#include "function.h"

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
    timeStepReductionStrategyType = "NoReduction";
    // read time step reduction strategy type
    IR_GIVE_OPTIONAL_FIELD(ir, timeStepReductionStrategyType, _IFT_MetaStep_timeReductionStrategyType);
    // create and initialize time step reduction strategy
    timeStepReductionStrategy = classFactory.createTimeStepReductionStrategy(this->timeStepReductionStrategyType.c_str(),1);
    timeStepReductionStrategy->initializeFrom(ir);
    // default value of dT is 1;
    deltaT = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, _IFT_MetaStep_deltaT);
    prescribedTimes.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, finalTime, _IFT_MetaStep_finalT);
    // compute minDeltaT as min of deltaT and mindeltaT form timeStepReudctionStrategy
    minDeltaT = std::min(deltaT, timeStepReductionStrategy->giveMinDeltaT());
    
    if ( finalTime < 0 ) {
      OOFEM_ERROR("Final time of the analysis can't be negative");
    } else if(finalTime == 0) {
      IR_GIVE_OPTIONAL_FIELD(ir, numberOfSteps, _IFT_MetaStep_nsteps);
      if ( numberOfSteps < 0 ) {
	OOFEM_ERROR("Numer of steps has to be positive number");
      } else {
	dtFunction = 0;
	if ( ir.hasField(_IFT_MetaStep_dtFunction) ) {
	  IR_GIVE_FIELD(ir, this->dtFunction, _IFT_MetaStep_dtFunction);
	} else if ( ir.hasField(_IFT_MetaStep_prescribedTimes) ) {
	  IR_GIVE_OPTIONAL_FIELD(ir, prescribedTimes, _IFT_MetaStep_prescribedTimes);
	  if ( prescribedTimes.giveSize() > 0 ) {
	    numberOfSteps = prescribedTimes.giveSize();
	    finalTime = prescribedTimes.at(numberOfSteps);
	  } 
	} else {
	  //@todo: how to get dt
	  finalTime = numberOfSteps * deltaT + this->previousMetaStepFinalTime;
	}
      }
    } 
    
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


double
MetaStep :: giveDeltaT(int n, std::unique_ptr<TimeStep> &previousStep)
{
    //@todo: Handle the dTfunction and discrete times appropriately
    if ( giveDtFunction() ) {
      return giveDtFunction()->evaluateAtTime(n);
    }
  
    if ( discreteTimes.giveSize() > 0 ) {
      return this->giveDiscreteTime(n) - this->giveDiscreteTime(n - 1);
    }

    if ( this->prescribedTimes.giveSize() > 0 ) {
      return this->prescribedTimes.at(previousStep->giveNumber() + 1) - previousStep->giveTargetTime();
    }    

    return deltaT;
}


void
MetaStep :: reduceTimeStep()
{
    double dT = deltaT * this->timeStepReductionStrategy->giveNoConvergenceTimeIncrementReductionFactor();
    if(dT >= this->timeStepReductionStrategy->give_dTmin()) {
      this->deltaT =  dT;
    } else {
      OOFEM_ERROR("Minimum time increment reached");
    }
}
  
void
MetaStep :: adaptTimeStep(int niter, double targetTime)
{
  auto dT_iter = this->timeStepReductionStrategy->giveReqIterTimeIncrementAdaptationFactor(niter);
  auto dT_MatMax = this->timeStepReductionStrategy->giveMaterialTimeIncrementAdaptationFactorMax();
  auto dT_MatMin = this->timeStepReductionStrategy->giveMaterialTimeIncrementAdaptationFactorMin();
  auto dT = this->deltaT;

  if(dT_iter < 1) {
    dT *= std::min(dT_MatMin, dT_iter);
  } else {
    if(dT_MatMin < 1) {
      dT *= dT_MatMin;
    } else {
      dT *= std::max(dT_MatMax, dT_iter);
    }
  }

  double dtmin = this->timeStepReductionStrategy->give_dTmin();
  double dtmax = this->timeStepReductionStrategy->give_dTmax();

  if(dT > dtmax){
    OOFEM_WARNING("Setting time increment to the maximum value dTmax, %f", dtmax);
    dT = dtmax;    
  }else  if(dT < dtmin){
    OOFEM_WARNING("Setting time increment to the minimum value dTmin, %f", dtmin);
    dT = dtmin;    
  }

  if( targetTime + dT > this->finalTime ) {
    dT = finalTime - targetTime;
  }

  
  this->deltaT = dT;     

}


void
MetaStep :: setTimeStepReductionFactor(double tsrf)
{
  this->timeStepReductionStrategy->setTimeStepIncrementAdaptationFactor(tsrf);
}
  

  
} // end namespace oofem
