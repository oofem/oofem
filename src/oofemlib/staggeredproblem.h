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

#ifndef staggeredproblem_h
#define staggeredproblem_h

#include "engngm.h"
#include "inputrecord.h"
#include "floatarray.h"

///@name Input fields for StaggeredProblem
//@{
#define _IFT_StaggeredProblem_Name "staggeredproblem"
#define _IFT_StaggeredProblem_deltat "deltat"
#define _IFT_StaggeredProblem_dtf "dtf"
#define _IFT_StaggeredProblem_timeDefinedByProb "timedefinedbyprob"
#define _IFT_StaggeredProblem_stepmultiplier "stepmultiplier"
//#define _IFT_StaggeredProblem_timeLag "timelag"
#define _IFT_StaggeredProblem_prescribedtimes "prescribedtimes"
#define _IFT_StaggeredProblem_prob1 "prob1"
#define _IFT_StaggeredProblem_prob2 "prob2"
#define _IFT_StaggeredProblem_coupling "coupling"
#define _IFT_StaggeredProblem_adaptiveStepLength "adaptivesteplength"
#define _IFT_StaggeredProblem_minsteplength "minsteplength"
#define _IFT_StaggeredProblem_maxsteplength "maxsteplength"
#define _IFT_StaggeredProblem_reqiterations "reqiterations"
#define _IFT_StaggeredProblem_endoftimeofinterest "endoftimeofinterest"
#define _IFT_StaggeredProblem_adaptivestepsince "adaptivestepsince"
//@}

namespace oofem {
class Function;

/**
 * Implementation of general sequence (staggered) problem. The problem consists in sequence of
 * low level problems (slaves) which are executed sequentially and where the results
 * of particular slave depends on the results of previous slaves in sequence.
 * Typical example is heat&mass transfer analysis followed by mechanical one, which
 * takes into account the temperature field from the first analysis.
 *
 * The sequence problem is represented by this class. It maintains list
 * of subsequent (slave) problems and it is executes the slave problems. It is responsible
 * for solution step generation and synchronization between slave problems.
 * The transfer of required state variables is done by mapping of corresponding variables
 * between problem domains. This allows to to transfer primary (nodal) values of one problem to
 * integration points of subsequent problem or to use completely different discretizations for
 * slave problems.
 *
 * Since the master problem is responsible for synchronization, it is responsible for
 * generation the solution steps. Therefore, the solution step specification, as well as
 * relevant meta step attributes are specified at master level.
 *
 * @note To avoid confusion,
 * the slaves are treated in so-called maintained mode. In this mode, the attributes and
 * meta step attributes are taken from the master. The local attributes, even if specified,
 * are ignored.
 */
class OOFEM_EXPORT StaggeredProblem : public EngngModel
{
protected:
    /// List of engineering models to solve sequentially.
    std :: vector< std :: unique_ptr< EngngModel > >emodelList;
    double deltaT;
    std :: vector< std :: string >inputStreamNames;
    /// Associated time function for time step increment
    int dtFunction;
    /**
     * Constant multiplier, optional input parameter. This parameter determines the ratio of
     * two consecutive time steps. Efficient for creep and relaxation analyses.
     */
    double stepMultiplier;

    /**
     * Time lag specifying how much is the second sub-problem delayed after the first one
     * during this period the second subproblem isn't solved at all.
     * Efficient tool for coupling structural problem with hydration.
     */
    //    double timeLag;

    /// Specified times where the problem is solved
    FloatArray discreteTimes;

    /// Optional parameter which specify problems to define load time functions
    int timeDefinedByProb;

    /// List of slave models to which this model is coupled
    IntArray coupledModels;
    bool adaptiveStepLength;
    /// adaptive time step length - minimum
    double minStepLength;
    /// adaptive time step length - maximum
    double maxStepLength;
    /// adaptive time step length - required (=optimum) number of iterations 
    double reqIterations;
    /// adaptive time step length applies after prescribed time
    double adaptiveStepSince;
    /**
     * alternative overriding the number of steps "nsteps" - necessary for time-driven analyses when
     * the appropriate number of steps is apriori unknow. If used, set "nsteps" to a high number e.g. 100000000
     */
    double endOfTimeOfInterest;

    double prevStepLength;
    double currentStepLength;
    

public:
    /**
     * Constructor. Creates an engineering model with number i belonging to domain d.
     */
    StaggeredProblem(int i, EngngModel * _master = NULL);
    /// Destructor.
    virtual ~StaggeredProblem();
    StaggeredProblem(const StaggeredProblem &) = delete;
    StaggeredProblem & operator=(const StaggeredProblem &) = delete;
    void setContextOutputMode(ContextOutputMode contextMode);
    void setUDContextOutputMode(int cStep);
    void setProblemMode(problemMode pmode);
    virtual void setRenumberFlag();

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);
    virtual int forceEquationNumbering();
    virtual void updateYourself(TimeStep *tStep);
    virtual void initializeYourself(TimeStep *tStep) { }
    virtual int initializeAdaptive(int tStepNumber) { return 0; }
    virtual void terminate(TimeStep *tStep);
    virtual void doStepOutput(TimeStep *tStep);

    virtual int instanciateYourself(DataReader *dr, InputRecord *ir, const char *outFileName, const char *desc);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void updateAttributes(MetaStep *mStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual void updateDomainLinks();

    void printYourself();
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual TimeStep *giveCurrentStep(bool force = false);
    virtual TimeStep *givePreviousStep(bool force = false);
    virtual TimeStep *giveSolutionStepWhenIcApply(bool force = false);
    virtual int giveNumberOfFirstStep(bool force = false);

    virtual TimeStep *giveNextStep();

    // identification
    virtual const char *giveClassName() const { return "StaggeredProblem"; }
    virtual const char *giveInputRecordName() const { return _IFT_StaggeredProblem_Name; }
    virtual int useNonlocalStiffnessOption() { return 0; }

    virtual fMode giveFormulation() { return UNKNOWN; }
    /**
     * Returns time function for time step increment.
     * Used time function should provide step lengths as function of step number.
     * Initial step with number 0 is considered as [ -dt(0), 0 ], first step is [ 0, dt(1) ], ...
     */
    Function *giveDtFunction();

    /**
     * Returns the timestep length for given step number n, initial step is number 0
     */
    double giveDeltaT(int n);

    /**
     * Returns time for time step number n (array discreteTimes must be specified)
     */
    double giveDiscreteTime(int n);

    /// Returns list of model number that this model is coupled with. Used for staggered approach.
    void giveCoupledModels(IntArray &answer) { answer = coupledModels; }

#ifdef __OOFEG
    virtual void drawYourself(oofegGraphicContext &gc);
    virtual void drawElements(oofegGraphicContext &gc);
    virtual void drawNodes(oofegGraphicContext &gc);
    virtual void showSparseMtrxStructure(int type, oofegGraphicContext &gc, TimeStep *tStep) { }
#endif

    virtual int checkProblemConsistency();

    virtual EngngModel *giveSlaveProblem(int i);
    virtual int giveNumberOfSlaveProblems() { return (int)inputStreamNames.size(); }
    virtual int instanciateDefaultMetaStep(InputRecord *ir);

protected:
    int instanciateSlaveProblems();
};
} // end namespace oofem
#endif // staggeredproblem_h
