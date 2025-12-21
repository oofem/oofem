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

#ifndef timestepcontroller_h
#define timestepcontroller_h


#include "timestep.h"
#include "metastep.h"

///@name Input fields
//@{
#define _IFT_TimeStepController_nmsteps "nmsteps"
#define _IFT_TimeStepController_alpha "alpha"

//@}

namespace oofem {
class DataReader;
  class MetaStep;

/**
 * This class is responsible for controlling the time step changes during the analysis based on the information from different levels
 * The class hold time step reduction strategy that defines how exactly is the step reduction/increa based on specific implementation. It can take into account number of global Newton's interations, convergence problems at the material or element level, etc.
 */
class OOFEM_EXPORT TimeStepController
{
protected:
    int numberOfMetaSteps = 0;
    int currentMetaStepNumber = 0;
    /// Current time step.
    std :: unique_ptr< TimeStep > currentStep;
    /// Previous time step.
    std :: unique_ptr< TimeStep > previousStep;
    ///
    EngngModel *eModel;
    /// Returns end of time interest (time corresponding to end of time integration).
    virtual double giveEndOfTimeOfInterest() { return 0.; }
    /// Returns the time step number, when initial conditions should apply.
    int giveNumberOfTimeStepWhenIcApply() { return 0; }
    /// Initial time from which the computation runs. Default is zero.
    double initT = 0.;
    /// List of problem metasteps.
    std :: vector< MetaStep > metaStepList;
    double alpha = 1.0; 
public:
    /**
     * Creates a new solution step.
     * @param n Solution step number.
     * @param e Reference to corresponding engng model.
     * @param mn Meta step number.
     * @param tt Intrinsic time.
     * @param dt Intrinsic time increment.
     * @param counter Solution state counter.
     * @param td Time discretization.
     */
     TimeStepController(EngngModel *e){eModel = e;}
    /**
     * Returns class name of receiver.
     * @return Pointer to s parameter filled with name.
     */
    const char *giveClassName() const { return "TimeStepController"; }


    /**
     * Returns the current time step.
     */
    virtual TimeStep *giveCurrentStep() {
      return currentStep.get();
      
    }
    /** Returns previous time step.
     *  @param force when set to true then previous step of receiver is returned instead of master (default)
     */
    virtual TimeStep *givePreviousStep() {
        return previousStep.get();
    }

    /// Returns next time step (next to current step) of receiver.
    //virtual TimeStep giveNextStep();
    virtual TimeStep* giveNextStep();


    TimeStep* generateNextStep(){return NULL;}

    

    
    void initializeFrom(InputRecord &ir);
    /**
     * Stores receiver state to output stream.
     * @param stream Output stream.
     * @exception ContextIOERR If error encountered.
     */
    void saveContext(DataStream &stream);
    /**
     * Restores the receiver state previously written in stream.
     * @see saveContext member function.
     */
    void restoreContext(DataStream &stream);
    //functions copied from engng
    /// Return number of meta steps.
    int giveNumberOfMetaSteps() { return numberOfMetaSteps; }
    /// Returns the i-th meta step.
    MetaStep *giveMetaStep(int i);

    int instanciateDefaultMetaStep(InputRecord &ir);
    int instanciateMetaSteps(DataReader &dr);
    void initMetaStepAttributes(MetaStep *mStep);
    int giveMetaStepNumber(){return 1;}
    void postInitialize();
    void reduceTimeStep();
    void adaptTimeStep(int nIter);
    //@todo improve: return number of steps of the current metastep
    int giveNumberOfSteps(){return this->giveCurrentMetaStep()->giveNumberOfSteps();}

    
    int giveFinalTime(){return 0;}
    double giveDeltaT(){return this->giveCurrentMetaStep()->giveDeltaT();}
    double giveMinDeltaT(){return this->giveCurrentMetaStep()->giveMinDeltaT();}
    void setDeltaT(double dT){return this->giveCurrentMetaStep()->setDeltaT(dT);}
    MetaStep *giveCurrentMetaStep(){return &metaStepList.at(currentMetaStepNumber);}
    void setCurrentMetaStepNumber(int smstep){currentMetaStepNumber = smstep;}

protected:
//  double finalTime = 0;


};
} // end namespace oofem
#endif // timestepcontroller_h

