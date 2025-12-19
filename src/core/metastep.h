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

#ifndef metastep_h
#define metastep_h

#include "oofemcfg.h"
#include "inputrecord.h"
#include "floatarray.h"
#include "timestepreductionstrategy.h"
#include <memory>

///@name Input fields for meta step
//@{
#define _IFT_MetaStep_Name "metastep"
#define _IFT_MetaStep_nsteps "nsteps"
#define _IFT_MetaStep_finalT "finalt"
#define _IFT_MetaStep_deltaT "deltat"
#define _IFT_MetaStep_prescribedTimes "prescribedtimes"
#define _IFT_MetaStep_dtFunction "dtfunction" ///< Function that determines size of time step.

#define _IFT_MetaStep_timeReductionStrategyType "treductiontype"

//@}

namespace oofem {
class EngngModel;
class TimeStepReductionStrategy;

/**
 * Class representing meta step. The meta step instance represent sequence of
 * solution steps (timeSteps). The meta step role is to describe the common
 * attributes related to solution steps it represent from the point of view of engineering model.
 * For example, meta step may represent series of solution steps, for which particular
 * solution control is used. The common attributes it represent depend on engineering model
 * representation. To store these dependent attributes, the metaStep record (currently string)
 * is read from input and is provided to engineering model upon request.
 *
 * The meta step maintains its number, the total number of steps it represent, time increment and
 * its e-model attributes.
 */
class OOFEM_EXPORT MetaStep
{
protected:
    /// Engineering model of receiver.
    EngngModel *eModel;
    /// Number of subsequent steps the receiver represent
    int numberOfSteps;
    /// Intrinsic time increment.
    double deltaT = 1;
    /// final time of the metastep
    double finalTime = 0;
    /// Engineering model attributes.
    std::shared_ptr<InputRecord> attributes;
    /// Start solution step number for which receiver is responsible.
    int sindex;
    /// Receiver number.
    int number;
    /// Specified times where the problem is solved    
    FloatArray prescribedTimes;
    /// Associated time function for time step increment.
    int dtFunction = 0;
    /// Specified times where the problem is solved merge with prescribedTimes
    FloatArray discreteTimes;
    /// Time step reduction strategy type - default no reduction of time step
    std :: string timeStepReductionStrategyType = "NoReductionStrategy";
    /// Time step reduction strategy type
    std :: unique_ptr<TimeStepReductionStrategy> timeStepReductionStrategy;
    /// Minimal value of the time step
    double minDeltaT;
    /// final time of the previous meta step
    double previousMetaStepFinalTime = 0;
    
public:
    /**
     * Constructor. Creates a new meta step.
     * @param n Meta step number.
     * @param e Reference to corresponding engineering model.
     */
    MetaStep(int n, EngngModel * e);
    MetaStep(int n, EngngModel * e, int nsteps, InputRecord & attrib);
    MetaStep(MetaStep &&ms) = default;

    /// Returns receiver's number.
    int giveNumber() { return number; }
    /// Returns number of Steps it represent.
    int giveNumberOfSteps() { return this->numberOfSteps; }
    /// Returns time increment.
    double giveTimeIncrement() { return this->deltaT; }
    /// Returns e-model attributes.
    InputRecord &giveAttributesRecord() { return *this->attributes; }
    /**
     * Instanciates the receiver from input record.
     */
    void initializeFrom(InputRecord &ir);
    /// Sets the receiver bounds according to given solution step number, returns end index.
    int setStepBounds(int startStepNumber);
    /// Sets the number of steps within the metastep.
    void setNumberOfSteps(int newNumberOfSteps);
    /// Tests if step number is maintained by receiver.
    int isStepValid(int solStepNumber);
    /// Returns the step relative number  to receiver.
    int giveStepRelativeNumber(int stepNumber) { return ( stepNumber - sindex + 1 ); }
    /// Returns first step number.
    int giveFirstStepNumber() { return sindex; }
    /// Returns last step number.
    int giveLastStepNumber() { return ( sindex + numberOfSteps - 1 ); }
    /// Returns class name of receiver.
    const char *giveClassName() const { return "MetaStep"; }
    /// Return maximal number of time step reductions
    int giveNumberOfMaxTimeStepReductions(){return timeStepReductionStrategy->giveNumberOfMaxTimeStepReductions();}
    /// Returns final time of the meta step
    double giveFinalTime(){return finalTime;}
    void setFinalTime(double ft){finalTime = ft;}
    void recomputeFinalTime(int nSteps) {finalTime =  nSteps * deltaT;}
    void setPreviousMetaStepFinalTime(double ft){previousMetaStepFinalTime = ft;}


    
    /** 
     * Returns time function for time step increment.
     * Used time function should provide step lengths as function of step number.
     * Initial step with number 0 is considered as [ -dt(0), 0 ], first step is [ 0, dt(1) ], ...
     */
    Function *giveDtFunction(){return NULL;}
    
    
    /// Returns the time step length.
    double giveDeltaT(){return this->deltaT;}
    double giveDeltaT(int n, std::unique_ptr<TimeStep> &previousStep);
   
    /// Set the time step length.
    void setDeltaT(double dT){this->deltaT = dT;}

    
    /// Return the minimal time step length set be the user
    double giveMinDeltaT(){return this->minDeltaT;}
    
    /**
     * Returns time for time step number n (array discreteTimes must be specified)
     */
    double giveDiscreteTime(int n){return discreteTimes.at(n);}
    void checkDiscreteTimes(){;}
    /// Reduce the time step length in case of no convergence.
    void reduceTimeStep();
    /// Adapt the time step length based on the number of iterations
    void adaptTimeStep(int nIter, double targetTime);
    /// Set the time step reduction factor used for reducing the time step length in case of no converfence
    void setTimeStepReductionFactor(double tStepRedFactor);
    //returns TimeStepReductionStrategy
    TimeStepReductionStrategy* giveTimeStepReductionStrategy() {return timeStepReductionStrategy.get();}
    

    
};
} // end namespace oofem
#endif // metastep_h
