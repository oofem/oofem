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

#include "timestepreductionstrategy.h"
#include "classfactory.h"
namespace oofem {
REGISTER_TimeStepReductionStrategy(NoReductionStrategy);  
REGISTER_TimeStepReductionStrategy(SimpleReductionStrategy);

 NoReductionStrategy::NoReductionStrategy(int n):TimeStepReductionStrategy(n)
{
  this->reductionFlag = false;

}


  
  SimpleReductionStrategy::SimpleReductionStrategy(int n):TimeStepReductionStrategy(n)
{
  this->reductionFlag = true;

}



double
SimpleReductionStrategy :: giveMaterialTimeIncrementAdaptationFactorMax()
{
  return materialTimeStepAdaptationFactorMax;
}

double
SimpleReductionStrategy :: giveMaterialTimeIncrementAdaptationFactorMin()
{
  return materialTimeStepAdaptationFactorMin;
}

  
double
SimpleReductionStrategy :: giveReqIterTimeIncrementAdaptationFactor(int nIter)
{
  if(nIter == 0 || (nIter >= nMinRequiredIterations && nIter <= nMaxRequiredIterations)) {
    return 1.;
  } else {
    if(nIter < nMinRequiredIterations) {
      return 2.0;
	//return static_cast<double>(nIter)/static_cast<double>(nMinRequiredIterations);
      
    } else if(nIter > nMaxRequiredIterations) {
      return 0.5;
      //return static_cast<double>(nMinRequiredIterations)/static_cast<double>(nIter);
    }
  }
  return 1.;
}

double
SimpleReductionStrategy :: giveNoConvergenceTimeIncrementReductionFactor()
{
  return noConvergenceReductionFactor;
}

  

void
SimpleReductionStrategy :: initializeFrom(InputRecord &ir)
{
    noConvergenceReductionFactor = 0.5;
    IR_GIVE_OPTIONAL_FIELD(ir, noConvergenceReductionFactor, _IFT_SimpleReductionStrategy_noConvergenceReductionFactor);
    //default value 5 ?
    numberOfMaxTimeStepReductions = 5;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfMaxTimeStepReductions, _IFT_SimpleReductionStrategy_nMaxRestarts);
    IR_GIVE_FIELD(ir, deltaTmax, _IFT_SimpleReductionStrategy_deltaTmax);
    IR_GIVE_FIELD(ir, deltaTmin, _IFT_SimpleReductionStrategy_deltaTmin);
    IR_GIVE_OPTIONAL_FIELD(ir, nMinRequiredIterations, _IFT_SimpleReductionStrategy_minRequiredIter);
    IR_GIVE_OPTIONAL_FIELD(ir, nMaxRequiredIterations, _IFT_SimpleReductionStrategy_maxRequiredIter);
    //
    //IR_GIVE_OPTIONAL_FIELD(ir, tsaf, _IFT_SimpleReductionStrategy_tsaf);
}
  
  
}
