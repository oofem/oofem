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

#ifndef timestepreductionstrategy_h
#define timestepreductionstrategy_h

#include "timestepreductionstrategytype.h"
#include "timestep.h"
#include "classfactory.h"
#include <limits>

#define _IFT_NoReductionStrategy_Name "noreduction"

//
#define _IFT_SimpleReductionStrategy_Name "simplereduction"
#define _IFT_SimpleReductionStrategy_nMaxRestarts "nmaxrestarts"
#define _IFT_SimpleReductionStrategy_deltaTmax "dtmax"
#define _IFT_SimpleReductionStrategy_deltaTmin "dtmin"
#define _IFT_SimpleReductionStrategy_minRequiredIter "minrequirediter"
#define _IFT_SimpleReductionStrategy_maxRequiredIter "maxrequirediter"
#define _IFT_SimpleReductionStrategy_noConvergenceReductionFactor "ncrf"

namespace oofem {

/**
 * The class is responsible for time step reduction based on specific implementation. It can take into account number of global Newton's interations, convergence problems at the material or element level, etc.
 */
class OOFEM_EXPORT TimeStepReductionStrategy
{
protected:
int number;
double deltaTmax = 1;
double deltaTmin = 1;
int nMinRequiredIterations = 3;
int nMaxRequiredIterations = 8;
int nMaxRestarts = 5;
double noConvergenceReductionFactor;
double materialTimeStepAdaptationFactorMax = 1;
double materialTimeStepAdaptationFactorMin = 1;
bool reductionFlag;
public:
   TimeStepReductionStrategy(int n){this->number = n;}
   virtual ~TimeStepReductionStrategy() = default;
   virtual void initializeFrom(InputRecord &ir)  = 0;
   virtual TimeStepReductionStrategyType giveTimeStepReductionStrategyType() = 0;
   //virtual double giveMaxMaterialTimeIncrementAdaptationFactor(std::unique_ptr<TimeStep> &ts) = 0;
   virtual double giveMaterialTimeIncrementAdaptationFactorMax() = 0;
   virtual double giveMaterialTimeIncrementAdaptationFactorMin() = 0;
   virtual double giveReqIterTimeIncrementAdaptationFactor(int nIter) = 0;
   virtual double giveNoConvergenceTimeIncrementReductionFactor() = 0;
   virtual int giveNumberOfMaxTimeStepReductions() = 0;
   virtual void setTimeStepIncrementAdaptationFactor(double tsrf) =  0;
   virtual double giveMinDeltaT(){return deltaTmin;}
   
   int giveNumberOfMinRequiredIterations(){return nMinRequiredIterations;}
   int giveNumberOfMaxRequiredIterations(){return nMaxRequiredIterations;}
   virtual double give_dTmin(){return deltaTmin;}
   virtual double give_dTmax(){return deltaTmax;}
   
   virtual void updateYourself(TimeStep *tStep) {;}
   virtual bool giveReductionFlag() {return reductionFlag;};
};



class NoReductionStrategy : public TimeStepReductionStrategy
{
public:
  NoReductionStrategy(int n);
    void initializeFrom(InputRecord &ir) override {};
  TimeStepReductionStrategyType giveTimeStepReductionStrategyType() override{return TSRST_NoReduction;}
  double giveMaterialTimeIncrementAdaptationFactorMax() override {return 1.;}
  double giveMaterialTimeIncrementAdaptationFactorMin() override {return 1.;}
  double giveReqIterTimeIncrementAdaptationFactor(int nIter) override {return 1.;}
  double giveNoConvergenceTimeIncrementReductionFactor() override {return 1.;}
  //do nothing
  void setTimeStepIncrementAdaptationFactor(double tsrf) override {;}
 
  int giveNumberOfMaxTimeStepReductions() override {return 0;}
  double give_dTmin() override {return 0.;}
  double give_dTmax() override {return std::numeric_limits<double>::infinity();}
};


  
class SimpleReductionStrategy : public TimeStepReductionStrategy
{
private:
  int numberOfMaxTimeStepReductions;
public:
  SimpleReductionStrategy(int n);
    void initializeFrom(InputRecord &ir) override;
  TimeStepReductionStrategyType giveTimeStepReductionStrategyType() override{return TSRST_SimpleReduction;}
  double giveMaterialTimeIncrementAdaptationFactorMax() override;
  double giveMaterialTimeIncrementAdaptationFactorMin() override;
  double giveReqIterTimeIncrementAdaptationFactor(int nIter) override;
  double giveNoConvergenceTimeIncrementReductionFactor() override;
  int giveNumberOfMaxTimeStepReductions() override {return this->numberOfMaxTimeStepReductions;}
  //do nothing
  void setTimeStepIncrementAdaptationFactor(double tsrf) override {
    if(tsrf < 1){
      materialTimeStepAdaptationFactorMin = tsrf;
    } else {
      materialTimeStepAdaptationFactorMin = tsrf;
    }
  }
  virtual void updateYourself(TimeStep *tStep) override
  {
    materialTimeStepAdaptationFactorMin = materialTimeStepAdaptationFactorMax = 1.;
  }
  
};

  
} // end namespace oofem
#endif // timestepreductionstrategy_h

