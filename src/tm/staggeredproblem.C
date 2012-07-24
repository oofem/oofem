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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "staggeredproblem.h"
#include "engngm.h"
#include "timestep.h"
#include "loadtime.h"
#include "metastep.h"
#include "exportmodulemanager.h"
#include "mathfem.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "verbose.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
StaggeredProblem :: StaggeredProblem(int i, EngngModel *_master) : EngngModel(i, _master)
{
    ndomains = 1; // domain is needed to store the time step ltf
    nModels  = 2;
    emodelList = new AList< EngngModel >(nModels);
    inputStreamNames = new std :: string [ nModels ];

    dtTimeFunction = 0;
    stepMultiplier = 1.;
    timeDefinedByProb = 0;
}

StaggeredProblem ::  ~StaggeredProblem()
{
    delete emodelList;
    delete [] inputStreamNames;
}

///////////
int
StaggeredProblem :: instanciateYourself(DataReader *dr, InputRecord *ir, const char *dataOutputFileName, const char *desc)
{
    int result;
    result = EngngModel :: instanciateYourself(dr, ir, dataOutputFileName, desc);
    ir->finish();
    // instanciate slave problems
    result &= this->instanciateSlaveProblems();
    return result;
}

int
StaggeredProblem :: instanciateDefaultMetaStep(InputRecord *ir)
{
    if ( !timeDefinedByProb ) {
        EngngModel :: instanciateDefaultMetaStep(ir);
    }

    //there are no slave problems initiated so far, the overall metaStep will defined in a slave problem instantiation

    return 1;
}

int
StaggeredProblem :: instanciateSlaveProblems()
{
    int i;
    EngngModel *timeDefProb, *slaveProb;

    //first instantiate master problem if defined
    if ( timeDefinedByProb ) {
        OOFEMTXTDataReader dr( inputStreamNames [ timeDefinedByProb - 1 ].c_str() );
        timeDefProb = oofem :: InstanciateProblem(& dr, this->pMode, this->contextOutputMode, NULL);
        emodelList->put(timeDefinedByProb, timeDefProb);
    }

    for ( i = 1; i <= nModels; i++ ) {
        if ( emodelList->includes(i) ) {
            continue;
        }

        OOFEMTXTDataReader dr( inputStreamNames [ i - 1 ].c_str() );
        //the slave problem dictating time needs to have attribute master=NULL, other problems point to the dictating slave
        slaveProb = oofem :: InstanciateProblem(& dr, this->pMode, this->contextOutputMode, timeDefinedByProb != 0 ? timeDefProb : this);
        emodelList->put(i, slaveProb);
    }

    return 1;
}


IRResultType
StaggeredProblem :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    if ( ir->hasField(IFT_StaggeredProblem_deltat, "deltat") ) {
        EngngModel :: initializeFrom(ir);
        IR_GIVE_FIELD(ir, deltaT, IFT_StaggeredProblem_deltat, "deltat"); // Macro
        dtTimeFunction = 0;
    } else {
        IR_GIVE_FIELD(ir, timeDefinedByProb, IFT_StaggeredProblem_timeDefinedByProb, "timedefinedbyprob"); // Macro
    }

    if ( dtTimeFunction < 1 ) {
        ndomains = 0;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, dtTimeFunction, IFT_StaggeredProblem_dtf, "dtf"); // Macro

    IR_GIVE_OPTIONAL_FIELD(ir, stepMultiplier, IFT_StaggeredProblem_stepmultiplier, "stepmultiplier"); // Macro
    if ( stepMultiplier < 0 ) {
        _error("stepMultiplier must be > 0")
    }

    IR_GIVE_FIELD(ir, inputStreamNames [ 0 ], IFT_StaggeredProblem_prob1, "prob1"); // Macro
    IR_GIVE_FIELD(ir, inputStreamNames [ 1 ], IFT_StaggeredProblem_prob2, "prob2"); // Macro

    return IRRT_OK;
}


void
StaggeredProblem :: updateAttributes(MetaStep *mStep)
{
    const char *__proc = "updateAttributes"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                  // Required by IR_GIVE_FIELD macro

    InputRecord *ir = mStep->giveAttributesRecord();

    EngngModel :: updateAttributes(mStep);

    // update attributes of slaves
    for ( int i = 1; i <= nModels; i++ ) {
        this->giveSlaveProblem(i)->updateAttributes(mStep);
    }

    if ( !timeDefinedByProb ) {
        IR_GIVE_FIELD(ir, deltaT, IFT_StaggeredProblem_deltat, "deltat"); // Macro
    }
}

LoadTimeFunction *StaggeredProblem :: giveDtTimeFunction()
// Returns the load-time function of the receiver.
{
    if ( !dtTimeFunction || !ndomains ) {
        return NULL;
    }

    return giveDomain(1)->giveLoadTimeFunction(dtTimeFunction);
}

double
StaggeredProblem :: giveDeltaT(int n)
{
    if ( giveDtTimeFunction() ) {
        return deltaT * giveDtTimeFunction()->__at(n);
    }

    //in the first step the time increment is taken as the initial, user-specified value
    if ( stepMultiplier != 1 && currentStep != NULL ) {
        if ( currentStep->giveNumber() >= 2 ) {
            return ( currentStep->giveTargetTime() * ( stepMultiplier ) );
        }
    }

    return deltaT;
}

TimeStep *
StaggeredProblem :: giveSolutionStepWhenIcApply()
{
    if ( stepWhenIcApply == NULL ) {
        int inin = giveNumberOfTimeStepWhenIcApply();
        stepWhenIcApply = new TimeStep(inin, this, 0, -giveDeltaT(inin), giveDeltaT(inin), 0);
    }

    return stepWhenIcApply;
}

TimeStep *
StaggeredProblem :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    double totalTime = 0;
    StateCounterType counter = 1;

    if ( previousStep != NULL ) {
        delete previousStep;
        previousStep = NULL;
    }

    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        totalTime = currentStep->giveTargetTime() + this->giveDeltaT(istep);
        counter = currentStep->giveSolutionStateCounter() + 1;
    } else {
        TimeStep *newStep;
        // first step -> generate initial step
        newStep = giveSolutionStepWhenIcApply();
        currentStep = new TimeStep(*newStep);
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, totalTime, this->giveDeltaT(istep), counter);

    // time and dt variables are set eq to 0 for statics - has no meaning
    return currentStep;
}


void
StaggeredProblem :: solveYourself()
{
    if ( !timeDefinedByProb ) {
        EngngModel :: solveYourself();
    } else { //time dictated by slave problem
        int imstep, jstep, nTimeSteps;
        int smstep = 1, sjstep = 1;
        MetaStep *activeMStep;
        FILE *out = this->giveOutputStream();
        this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);

        EngngModel *sp = this->giveSlaveProblem(timeDefinedByProb);

        if ( sp->giveCurrentStep() ) {
            smstep = sp->giveCurrentStep()->giveMetaStepNumber();
            sjstep = sp->giveMetaStep(smstep)->giveStepRelativeNumber( sp->giveCurrentStep()->giveNumber() ) + 1;
        } else {
            nMetaSteps = sp->giveNumberOfMetaSteps();
        }

        for ( imstep = smstep; imstep <= nMetaSteps; imstep++, sjstep = 1 ) { //loop over meta steps
            activeMStep = sp->giveMetaStep(imstep);
            // update state according to new meta step in all slaves
            this->initMetaStepAttributes(activeMStep);

            nTimeSteps = activeMStep->giveNumberOfSteps();
            for ( jstep = sjstep; jstep <= nTimeSteps; jstep++ ) { //loop over time steps
                this->timer.startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
                this->timer.initTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
                sp->giveNextStep();

                // renumber equations if necessary. Ensure to call forceEquationNumbering() for staggered problems
                if ( this->requiresEquationRenumbering( sp->giveCurrentStep() ) || this->giveClassID() == classType(StaggeredProblemClass) ) {
                    this->forceEquationNumbering();
                }

                this->solveYourselfAt( sp->giveCurrentStep() );
                this->updateYourself( sp->giveCurrentStep() );
                this->terminate( sp->giveCurrentStep() );

                this->timer.stopTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
                double _steptime = this->timer.getUtime(EngngModelTimer :: EMTT_SolutionStepTimer);
                OOFEM_LOG_INFO("EngngModel info: user time consumed by solution step %d: %.2fs\n",
                               sp->giveCurrentStep()->giveNumber(), _steptime);

                fprintf(out, "\nUser time consumed by solution step %d: %.3f [s]\n\n",
                        sp->giveCurrentStep()->giveNumber(), _steptime);

#ifdef __PARALLEL_MODE
                if ( loadBalancingFlag ) {
                    this->balanceLoad( sp->giveCurrentStep() );
                }

#endif
            }
        }
    }
}

void
StaggeredProblem :: solveYourselfAt(TimeStep *stepN)
{
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %5d, time %e]\n", stepN->giveNumber(), stepN->giveTargetTime() );
#endif
    for ( int i = 1; i <= nModels; i++ ) {
        EngngModel *emodel = this->giveSlaveProblem(i);
        emodel->solveYourselfAt(stepN);
    }

    stepN->incrementStateCounter();
}

int
StaggeredProblem :: forceEquationNumbering()
{
    int neqs = 0;
    for ( int i = 1; i <= nModels; i++ ) {
        EngngModel *emodel = this->giveSlaveProblem(i);
        // renumber equations if necessary
        if ( emodel->requiresEquationRenumbering( emodel->giveCurrentStep() ) ) {
            neqs += emodel->forceEquationNumbering();
        }
    }

    return neqs;
}

void
StaggeredProblem :: updateYourself(TimeStep *stepN)
{
    for ( int i = 1; i <= nModels; i++ ) {
        this->giveSlaveProblem(i)->updateYourself(stepN);
    }

    EngngModel :: updateYourself(stepN);
}

void
StaggeredProblem :: terminate(TimeStep *tStep) {
    for ( int i = 1; i <= nModels; i++ ) {
        this->giveSlaveProblem(i)->terminate(tStep);
    }
}

void
StaggeredProblem :: doStepOutput(TimeStep *stepN)
{
    FILE *File = this->giveOutputStream();

    // print output
    this->printOutputAt(File, stepN);
    // export using export manager
    for ( int i = 1; i <= nModels; i++ ) {
        this->giveSlaveProblem(i)->giveExportModuleManager()->doOutput(stepN);
    }
}


void
StaggeredProblem :: printOutputAt(FILE *File, TimeStep *stepN)
{
    FILE *slaveFile;
    for ( int i = 1; i <= nModels; i++ ) {
        slaveFile = this->giveSlaveProblem(i)->giveOutputStream();
        this->giveSlaveProblem(i)->printOutputAt(slaveFile, stepN);
    }
}


contextIOResultType
StaggeredProblem :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    EngngModel :: saveContext(stream, mode, obj);
    for ( int i = 1; i <= nModels; i++ ) {
        this->giveSlaveProblem(i)->saveContext(stream, mode, obj);
    }

    return CIO_OK;
}


contextIOResultType
StaggeredProblem :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    EngngModel :: restoreContext(stream, mode, obj);
    for ( int i = 1; i <= nModels; i++ ) {
        this->giveSlaveProblem(i)->restoreContext(stream, mode, obj);
    }

    return CIO_OK;
}


EngngModel *
StaggeredProblem :: giveSlaveProblem(int i)
{
    if ( ( i > 0 ) && ( i <= this->nModels ) ) {
        return this->emodelList->at(i);
    } else {
        _error("giveSlaveProblem: Undefined problem");
    }

    return NULL;
}


int
StaggeredProblem :: checkConsistency()
{
    // check internal consistency
    // if success returns nonzero
    int result = 1;
    for ( int i = 1; i <= nModels; i++ ) {
        result &= this->giveSlaveProblem(i)->checkConsistency();
    }

    return result;
}

void
StaggeredProblem :: updateDomainLinks()
{
    for ( int i = 1; i <= nModels; i++ ) {
        this->giveSlaveProblem(i)->updateDomainLinks();
    }
}

void
StaggeredProblem :: setRenumberFlag()
{
    for ( int i = 1; i <= nModels; i++ ) {
        this->giveSlaveProblem(i)->setRenumberFlag();
    }
}

#ifdef __OOFEG
void StaggeredProblem :: drawYourself(oofegGraphicContext &context) {
    int ap = context.getActiveProblemIndx();
    if ( ( ap > 0 ) && ( ap <= nModels ) ) {
        this->giveSlaveProblem(ap)->drawYourself(context);
    }
}

void StaggeredProblem :: drawElements(oofegGraphicContext &context) {
    int ap = context.getActiveProblemIndx();
    if ( ( ap > 0 ) && ( ap <= nModels ) ) {
        this->giveSlaveProblem(ap)->drawElements(context);
    }
}

void StaggeredProblem :: drawNodes(oofegGraphicContext &context) {
    int ap = context.getActiveProblemIndx();
    if ( ( ap > 0 ) && ( ap <= nModels ) ) {
        this->giveSlaveProblem(ap)->drawNodes(context);
    }
}
#endif
} // end namespace oofem
