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

#include <iostream>
using namespace std;

#include "../sm/EngineeringModels/nlineardynamic.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/Elements/structuralelementevaluator.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "dof.h"
#include "error.h"
#include "verbose.h"
#include "sparsenonlinsystemnm.h"
#include "nrsolver.h"
#include "outputmanager.h"
#include "datastream.h"
#include "classfactory.h"
#include "timer.h"
#include "contextioerr.h"
#include "sparsemtrx.h"
#include "errorestimator.h"
#include "unknownnumberingscheme.h"
#include "assemblercallback.h"

#ifdef __PARALLEL_MODE
 #include "loadbalancer.h"
 #include "problemcomm.h"
 #include "processcomm.h"
#endif

namespace oofem {
REGISTER_EngngModel(NonLinearDynamic);

NonLinearDynamic :: NonLinearDynamic(int i, EngngModel *_master) : StructuralEngngModel(i, _master),
    totalDisplacement(), incrementOfDisplacement(), internalForces(), forcesVector()
{
    initialTimeDiscretization = TD_ThreePointBackward;

    ndomains                = 1;
    internalVarUpdateStamp  = 0;
    initFlag = commInitFlag = 1;
}


NonLinearDynamic :: ~NonLinearDynamic()
{
}


NumericalMethod *NonLinearDynamic :: giveNumericalMethod(MetaStep *mStep)
{
    if ( mStep == NULL ) {
        OOFEM_ERROR("undefined meta step");
    }

    if ( this->nMethod ) {
        this->nMethod->reinitialize();
    } else {
        this->nMethod.reset( new NRSolver(this->giveDomain(1), this) );
    }

    return this->nMethod.get();
}


void
NonLinearDynamic :: updateAttributes(MetaStep *mStep)
{
    IRResultType result;                     // Required by IR_GIVE_FIELD macro

    InputRecord *ir = mStep->giveAttributesRecord();

    StructuralEngngModel :: updateAttributes(mStep);

    deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, _IFT_NonLinearDynamic_deltat);
    if ( deltaT < 0. ) {
        OOFEM_ERROR("deltaT < 0");
    }

    eta = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, eta, _IFT_NonLinearDynamic_eta);

    delta = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, delta, _IFT_NonLinearDynamic_delta);
}


IRResultType
NonLinearDynamic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    result = StructuralEngngModel :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    nonlocalStiffnessFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nonlocalStiffnessFlag, _IFT_NonLinearDynamic_nonlocstiff);

    deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, _IFT_NonLinearDynamic_deltat);
    if ( deltaT < 0. ) {
        OOFEM_ERROR("deltaT < 0");
    }

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NonLinearDynamic_ddtScheme);
    initialTimeDiscretization = ( TimeDiscretizationType ) val;

    gamma = 0.5;
    beta = 0.25;              // Default Newmark parameters.
    if ( initialTimeDiscretization == TD_Newmark ) {
        OOFEM_LOG_INFO("Selecting Newmark-beta metod\n");
        IR_GIVE_OPTIONAL_FIELD(ir, gamma, _IFT_NonLinearDynamic_gamma);
        IR_GIVE_OPTIONAL_FIELD(ir, beta, _IFT_NonLinearDynamic_beta);
    } else if ( initialTimeDiscretization == TD_TwoPointBackward ) {
        OOFEM_LOG_INFO("Selecting Two-point Backward Euler method\n");
    } else if ( initialTimeDiscretization == TD_ThreePointBackward ) {
        OOFEM_LOG_INFO("Selecting Three-point Backward Euler metod\n");
    } else {
        OOFEM_WARNING("Time-stepping scheme not found!");
        return IRRT_BAD_FORMAT;
    }

    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, _IFT_NRSolver_manrmsteps);

#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        commBuff = new CommunicatorBuff(this->giveNumberOfProcesses(), CBT_static);
        communicator = new NodeCommunicator(this, commBuff, this->giveRank(),
                                            this->giveNumberOfProcesses());

        if ( ir->hasField(_IFT_NonLinearDynamic_nonlocalext) ) {
            nonlocalExt = 1;
            nonlocCommunicator = new ElementCommunicator(this, commBuff, this->giveRank(),
                                                         this->giveNumberOfProcesses());
        }
    }
#endif

    return IRRT_OK;
}


double NonLinearDynamic :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }
#endif

    if ( tStep != this->giveCurrentStep() ) {
        OOFEM_ERROR("unknown time step encountered");
        return 0.;
    }

    switch ( mode ) {
    case VM_Incremental:
        return incrementOfDisplacement.at(eq);

    case VM_Total:
        return totalDisplacement.at(eq);

    case VM_Velocity:
        return velocityVector.at(eq);

    case VM_Acceleration:
        return accelerationVector.at(eq);

    default:
        OOFEM_ERROR("Unknown is of undefined ValueModeType for this problem");
    }

    return 0.0;
}


TimeStep *NonLinearDynamic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    int mStepNum = 1;
    StateCounterType counter = 1;
    double deltaTtmp = deltaT;
    double totalTime = deltaT;
    TimeDiscretizationType td = initialTimeDiscretization;

    //Do not increase deltaT on microproblem
    if ( pScale == microScale ) {
        deltaTtmp = 0.;
    }

    if ( currentStep ) {
        totalTime = currentStep->giveTargetTime() + deltaTtmp;
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
        mStepNum = currentStep->giveMetaStepNumber();
        td = currentStep->giveTimeDiscretization();

        if ( !this->giveMetaStep(mStepNum)->isStepValid(istep) ) {
            mStepNum++;
            if ( mStepNum > nMetaSteps ) {
                OOFEM_ERROR("no next step available, mStepNum=%d > nMetaSteps=%d", mStepNum, nMetaSteps);
            }
        }
    }

    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(istep, this, mStepNum, totalTime, deltaTtmp, counter, td) );

    return currentStep.get();
}


void NonLinearDynamic :: solveYourself()
{
    if ( commInitFlag ) {
 #ifdef __VERBOSE_PARALLEL
        // Force equation numbering before setting up comm maps.
        int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());
        OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
 #endif
        if ( isParallel() ) {
            // Set up communication patterns.
            this->initializeCommMaps();
        }
        commInitFlag = 0;
    }

    StructuralEngngModel :: solveYourself();
}


void
NonLinearDynamic :: solveYourselfAt(TimeStep *tStep)
{
    if ( commInitFlag ) {
 #ifdef __VERBOSE_PARALLEL
        // Force equation numbering before setting up comm maps.
        int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());
        OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
 #endif
        if ( isParallel() ) {
            // Set up communication patterns.
            this->initializeCommMaps();
        }
        commInitFlag = 0;
    }

#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving       [Step number %8d, Time %15e]\n\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    proceedStep(1, tStep);
}

void
NonLinearDynamic :: terminate(TimeStep *tStep)
{
    this->doStepOutput(tStep);
    this->printReactionForces(tStep, 1);
    fflush( this->giveOutputStream() );
    this->saveStepContext(tStep);
}

void 
NonLinearDynamic :: initializeYourself(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());

    if ( tStep->isTheFirstStep() && initFlag ) {
        // Initialization
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        velocityVector.resize(neq);
        velocityVector.zero();
        accelerationVector.resize(neq);
        accelerationVector.zero();
        internalForces.resize(neq);
        internalForces.zero();
        previousIncrementOfDisplacement.resize(neq);
        previousIncrementOfDisplacement.zero();
        previousTotalDisplacement.resize(neq);
        previousTotalDisplacement.zero();
        previousVelocityVector.resize(neq);
        previousVelocityVector.zero();
        previousAccelerationVector.resize(neq);
        previousAccelerationVector.zero();
        previousInternalForces.resize(neq);
        previousInternalForces.zero();

        TimeStep *stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0,
                                                 -deltaT, deltaT, 0);

        // Considering initial conditions.
        for ( auto &node : domain->giveDofManagers() ) {
            for ( Dof *dof: *node ) {
                // Ask for initial values obtained from
                // bc (boundary conditions) and ic (initial conditions).
                if ( !dof->isPrimaryDof() ) {
                    continue;
                }

                int jj = dof->__giveEquationNumber();
                if ( jj ) {
                    totalDisplacement.at(jj)  = dof->giveUnknown(VM_Total, stepWhenIcApply);
                    velocityVector.at(jj)     = dof->giveUnknown(VM_Velocity, stepWhenIcApply);
                    accelerationVector.at(jj) = dof->giveUnknown(VM_Acceleration, stepWhenIcApply);
                }
            }
        }
        this->giveInternalForces(internalForces, true, 1, tStep);
    }
}

void
NonLinearDynamic :: proceedStep(int di, TimeStep *tStep)
{
    // creates system of governing eq's and solves them at given time step
    // first assemble problem at current time step

    int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());

    // Time-stepping constants
    this->determineConstants(tStep);



    if ( initFlag ) {
        // First assemble problem at current time step.
        // Option to take into account initial conditions.
        if ( !effectiveStiffnessMatrix ) {
            effectiveStiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
            massMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        }

        if ( !effectiveStiffnessMatrix || !massMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        if ( nonlocalStiffnessFlag ) {
            if ( !effectiveStiffnessMatrix->isAsymmetric() ) {
                OOFEM_ERROR("effectiveStiffnessMatrix does not support asymmetric storage");
            }
        }

        effectiveStiffnessMatrix->buildInternalStructure( this, di, EModelDefaultEquationNumbering() );
        massMatrix->buildInternalStructure( this, di, EModelDefaultEquationNumbering() );

        // Assemble mass matrix
        this->assemble( *massMatrix, tStep, MassMatrixAssembler(),
                       EModelDefaultEquationNumbering(), this->giveDomain(di) );

        // Initialize vectors
        help.resize(neq);
        help.zero();
        rhs.resize(neq);
        rhs.zero();
        rhs2.resize(neq);
        rhs2.zero();

        previousIncrementOfDisplacement = incrementOfDisplacement;
        previousTotalDisplacement = totalDisplacement;
        previousVelocityVector = velocityVector;
        previousAccelerationVector = accelerationVector;
        previousInternalForces = internalForces;

        forcesVector.resize(neq);
        forcesVector.zero();

        totIterations = 0;
        initFlag = 0;
    }

#ifdef VERBOSE
    OOFEM_LOG_DEBUG("Assembling load\n");
#endif

    // Assemble the external forces
    FloatArray loadVector;
    loadVector.resize( this->giveNumberOfDomainEquations( di, EModelDefaultEquationNumbering() ) );
    loadVector.zero();
    this->assembleVector( loadVector, tStep, ExternalForceAssembler(),
                         VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(di) );
    this->updateSharedDofManagers(loadVector, EModelDefaultEquationNumbering(), LoadExchangeTag);

    // Assembling the effective load vector
    for ( int i = 1; i <= neq; i++ ) {
        //help.beScaled(a2, previousVelocityVector);
        //help.add(a3, previousAccelerationVector);
        //help.add(a4 * eta, previousVelocityVector);
        //help.add(a5 * eta, previousAccelerationVector);
        //help.add(a6 * eta, previousIncrementOfDisplacement);
        help.at(i) = a2 * previousVelocityVector.at(i) + a3 *previousAccelerationVector.at(i)
        + eta * ( a4 * previousVelocityVector.at(i)
                 + a5 * previousAccelerationVector.at(i)
                 + a6 * previousIncrementOfDisplacement.at(i) );
    }

    massMatrix->times(help, rhs);

    if ( delta != 0 ) {
        //help.beScaled(a4 * delta, previousVelocityVector);
        //help.add(a5 * delta, previousAccelerationVector);
        //help.add(a6 * delta, previousIncrementOfDisplacement);
        for ( int i = 1; i <= neq; i++ ) {
            help.at(i) = delta * ( a4 * previousVelocityVector.at(i)
                                  + a5 * previousAccelerationVector.at(i)
                                  + a6 * previousIncrementOfDisplacement.at(i) );
        }
        this->timesMtrx(help, rhs2, TangentStiffnessMatrix, this->giveDomain(di), tStep);

        help.zero();
        rhs.add(rhs2);

        //this->assembleVector(rhs, tStep, MatrixProductAssembler(TangentAssembler(), help), VM_Total, 
        //                    EModelDefaultEquationNumbering(), this->giveDomain(1));
    }

    rhs.add(loadVector);
    rhs.subtract(previousInternalForces);

    //
    // Set-up numerical model.
    //
    this->giveNumericalMethod( this->giveCurrentMetaStep() );

    //
    // Call numerical model to solve problem.
    //
    double loadLevel = 1.0;

    if ( totIterations == 0 ) {
        incrementOfDisplacement.zero();
    }

    NM_Status numMetStatus = nMethod->solve(*effectiveStiffnessMatrix, rhs, NULL, NULL,
                                            totalDisplacement, incrementOfDisplacement, forcesVector,
                                            internalForcesEBENorm, loadLevel, SparseNonLinearSystemNM :: rlm_total, currentIterations, tStep);
    if ( !( numMetStatus & NM_Success ) ) {
        OOFEM_ERROR("NRSolver failed to solve problem");
    }

    rhs = previousVelocityVector;
    rhs2 = previousAccelerationVector;
    for ( int i = 1; i <= neq; i++ ) {
        accelerationVector.at(i) = a0 * incrementOfDisplacement.at(i) - a2 *rhs.at(i) - a3 *rhs2.at(i);
        velocityVector.at(i)     = a1 * incrementOfDisplacement.at(i) - a4 *rhs.at(i) - a5 *rhs2.at(i)
        - a6 *previousIncrementOfDisplacement.at(i);
    }
    totIterations += currentIterations;
}


void
NonLinearDynamic :: determineConstants(TimeStep *tStep)
{
    TimeDiscretizationType timeDiscretization = tStep->giveTimeDiscretization();

    if ( timeDiscretization == TD_Newmark ) {
        OOFEM_LOG_DEBUG("Solving using Newmark-beta method\n");
    } else if ( timeDiscretization == TD_TwoPointBackward ) {
        OOFEM_LOG_DEBUG("Solving using Backward Euler method\n");
    } else if ( timeDiscretization == TD_ThreePointBackward ) {
        OOFEM_LOG_DEBUG("Solving using Three-point Backward Euler method\n");
        if ( tStep->isTheFirstStep() ) {
            timeDiscretization = TD_TwoPointBackward;
        }
    } else {
        OOFEM_ERROR("Time-stepping scheme not found!\n")
    }

    deltaT = tStep->giveTimeIncrement();

    double dt2 = deltaT * deltaT;

    if ( timeDiscretization == TD_Newmark ) {
        a0 = 1. / ( beta * dt2 );
        a1 = gamma / ( beta * deltaT );
        a2 = 1. / ( beta * deltaT );
        a3 = 1. / ( 2. *  beta ) - 1.;
        a4 = ( gamma / beta ) - 1.;
        a5 = deltaT / 2. * ( gamma / beta - 2. );
        a6 = 0.;
    } else if ( timeDiscretization == TD_TwoPointBackward ) {
        a0 = 1. / dt2;
        a1 = 1. / deltaT;
        a2 = 1. / deltaT;
        a3 = 0.;
        a4 = 0.;
        a5 = 0.;
        a6 = 0.;
    } else if ( timeDiscretization == TD_ThreePointBackward ) {
        a0 = 2. / dt2;
        a1 = 3. / ( 2. * deltaT );
        a2 = 2. / deltaT;
        a3 = 0.;
        a4 = 0.;
        a5 = 0.;
        a6 = 1. / ( 2. * deltaT );
    }
}


void NonLinearDynamic :: updateYourself(TimeStep *tStep)
{
    totIterations = 0;

    previousVelocityVector = velocityVector;
    previousAccelerationVector = accelerationVector;
    previousIncrementOfDisplacement = incrementOfDisplacement;
    previousTotalDisplacement = totalDisplacement;
    previousInternalForces = internalForces;

    StructuralEngngModel :: updateYourself(tStep);
}

void NonLinearDynamic :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
//
// Updates some component, which is used by numerical method
// to newly reached state. used mainly by numerical method
// when new tanget stiffness is needed during finding
// of new equlibrium stage.
//
{
#ifdef TIME_REPORT
    Timer timer;
#endif

    switch ( cmpn ) {
    case NonLinearLhs:
        // Prevent assembly if already assembled ( totIterations > 0 )
        // Allow if MANRMSteps != 0
        if ( ( totIterations == 0 ) || MANRMSteps ) {
#ifdef VERBOSE
            OOFEM_LOG_DEBUG("Updating effective stiffness matrix\n");
#endif
#ifdef TIME_REPORT
            timer.startTimer();
#endif
#if 1
            effectiveStiffnessMatrix->zero();
            this->assemble(*effectiveStiffnessMatrix, tStep, EffectiveTangentAssembler(false, 1 + this->delta * a1,  this->a0 + this->eta * this->a1),
                           EModelDefaultEquationNumbering(), d);
#else
            this->assemble(effectiveStiffnessMatrix, tStep, TangentStiffnessMatrix,
                           EModelDefaultEquationNumbering(), d);
            effectiveStiffnessMatrix->times(1. + this->delta * a1);
            effectiveStiffnessMatrix->add(this->a0 + this->eta * this->a1, this->massMatrix);
#endif
#ifdef TIME_REPORT
            timer.stopTimer();
            OOFEM_LOG_DEBUG( "User time consumed by updating nonlinear LHS: %.2fs\n", timer.getUtime() );
#endif
        }
        break;
    case InternalRhs:
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Updating internal RHS\n");
#endif
        {
#ifdef TIME_REPORT
            timer.startTimer();
#endif
            if ( ( currentIterations != 0 ) || ( totIterations == 0 ) ) {
                this->giveInternalForces(internalForces, true, d->giveNumber(), tStep);

                // Updating the residual vector @ NR-solver
                help.beScaled(a0 + eta * a1, incrementOfDisplacement);

                massMatrix->times(help, rhs2);

                forcesVector = internalForces;
                forcesVector.add(rhs2);
                forcesVector.subtract(previousInternalForces);

                if ( delta != 0 ) {
                    help.beScaled(delta * a1, incrementOfDisplacement);
                    this->timesMtrx(help, rhs2, TangentStiffnessMatrix, this->giveDomain(1), tStep);
                    //this->assembleVector(rhs2, tStep, MatrixProductAssembler(TangentAssembler(), help), VM_Total, 
                    //                    EModelDefaultEquationNumbering(), this->giveDomain(1));

                    forcesVector.add(rhs2);
                }
            }
#ifdef TIME_REPORT
            timer.stopTimer();
            OOFEM_LOG_DEBUG( "User time consumed by updating internal RHS: %.2fs\n", timer.getUtime() );
#endif
        }
        break;
    default:
        OOFEM_ERROR("Unknown Type of component.");
    }
}

void
NonLinearDynamic :: printOutputAt(FILE *File, TimeStep *tStep)
{
    //FILE* File = this -> giveDomain() -> giveOutputStream() ;

    if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return; // Do not print even Solution step header
    }

    fprintf( File, "\n\nOutput for time %.3e, solution step number %d\n", tStep->giveTargetTime(), tStep->giveNumber() );
    fprintf(File, "Equilibrium reached in %d iterations\n\n", currentIterations);


    nMethod->printState(File);

    this->giveDomain(1)->giveOutputManager()->doDofManOutput(File, tStep);
    this->giveDomain(1)->giveOutputManager()->doElementOutput(File, tStep);
}

void
NonLinearDynamic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    static char dofchar[] = "dva";
    static ValueModeType dofmodes[] = {
        VM_Total, VM_Velocity, VM_Acceleration
    };

    iDof->printMultipleOutputAt(stream, tStep, dofchar, dofmodes, 3);
}

contextIOResultType NonLinearDynamic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file = NULL;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = incrementOfDisplacement.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = totalDisplacement.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    return CIO_OK;
}



contextIOResultType NonLinearDynamic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    int closeFlag = 0;
    int istep, iversion;
    contextIOResultType iores;
    FILE *file = NULL;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);
    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = incrementOfDisplacement.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = totalDisplacement.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    } // ensure consistent records

    return CIO_OK;
}


void
NonLinearDynamic :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();

    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
#ifdef __PARALLEL_MODE
    if ( this->giveLoadBalancer() ) {
        this->giveLoadBalancer()->setDomain( this->giveDomain(1) );
    }

#endif
}


void
NonLinearDynamic :: assemble(SparseMtrx &answer, TimeStep *tStep, const MatrixAssembler &ma,
                             const UnknownNumberingScheme &s, Domain *domain)
{
#ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
#endif

    EngngModel :: assemble(answer, tStep, ma, s, domain);

    if ( ( nonlocalStiffnessFlag ) && dynamic_cast< const TangentAssembler* >(&ma) ) {
        // Add nonlocal contribution.
        for ( auto &elem : domain->giveElements() ) {
            static_cast< NLStructuralElement * >( elem.get() )->addNonlocalStiffnessContributions(answer, s, tStep);
        }

        // Print storage statistics.
        answer.printStatistics();
    }

#ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_DEBUG( "NonLinearDynamic info: user time consumed by assembly: %.2fs\n", timer.getUtime() );
#endif
}

#ifdef __OOFEG
void
NonLinearDynamic :: showSparseMtrxStructure(int type, oofegGraphicContext &gc, TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    CharType ctype;

    if ( type != 1 ) {
        return;
    }

    ctype = TangentStiffnessMatrix;

    for ( auto &elem : domain->giveElements() ) {
        elem->showSparseMtrxStructure(ctype, gc, tStep);
    }

    for ( auto &elem : domain->giveElements() ) {
        elem->showExtendedSparseMtrxStructure(ctype, gc, tStep);
    }
}
#endif


void
NonLinearDynamic :: timesMtrx(FloatArray &vec, FloatArray &answer, CharType type, Domain *domain, TimeStep *tStep)
{
    int nelem = domain->giveNumberOfElements();
    //int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());
    int jj, kk, n;
    FloatMatrix charMtrx;
    IntArray loc;
    EModelDefaultEquationNumbering en;

    answer.resize( this->giveNumberOfDomainEquations(domain->giveNumber(), en) );
    answer.zero();
    for ( int i = 1; i <= nelem; i++ ) {
        Element *element = domain->giveElement(i);
        // Skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

        element->giveLocationArray(loc, en);
        element->giveCharacteristicMatrix(charMtrx, type, tStep);

        //
        // assemble it manually
        //
#ifdef DEBUG
        if ( loc.giveSize() != charMtrx.giveNumberOfRows() ) {
            OOFEM_ERROR("dimension mismatch");
        }

#endif

        n = loc.giveSize();

        for ( int j = 1; j <= n; j++ ) {
            jj = loc.at(j);
            if ( jj ) {
                for ( int k = 1; k <= n; k++ ) {
                    kk = loc.at(k);
                    if ( kk ) {
                        answer.at(jj) += charMtrx.at(j, k) * vec.at(kk);
                    }
                }
            }
        }
    }

    this->updateSharedDofManagers(answer, EModelDefaultEquationNumbering(), MassExchangeTag);
}


int
NonLinearDynamic :: estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType)
{
    int count = 0, pcount = 0;
    Domain *domain = this->giveDomain(1);

    if ( packUnpackType == 0 ) { ///@todo Fix this old ProblemCommMode__NODE_CUT value
        for ( int map: commMap ) {
            DofManager *dman = domain->giveDofManager( map );
            for ( Dof *dof: *dman ) {
                if ( dof->isPrimaryDof() && ( dof->__giveEquationNumber() ) ) {
                    count++;
                } else {
                    pcount++;
                }
            }
        }

        return ( buff.givePackSizeOfDouble(1) * max(count, pcount) );
    } else if ( packUnpackType == 1 ) {
        for ( int map: commMap ) {
            count += domain->giveElement( map )->estimatePackSize(buff);
        }

        return count;
    }

    return 0;
}


#ifdef __PARALLEL_MODE
LoadBalancer *
NonLinearDynamic :: giveLoadBalancer()
{
    if ( lb ) {
        return lb;
    }

    if ( loadBalancingFlag ) {
        lb = classFactory.createLoadBalancer( "parmetis", this->giveDomain(1) );
        return lb;
    } else {
        return NULL;
    }
}


LoadBalancerMonitor *
NonLinearDynamic :: giveLoadBalancerMonitor()
{
    if ( lbm ) {
        return lbm;
    }

    if ( loadBalancingFlag ) {
        lbm = classFactory.createLoadBalancerMonitor( "wallclock", this);
        return lbm;
    } else {
        return NULL;
    }
}
#endif

void
NonLinearDynamic :: packMigratingData(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    int ndofman = domain->giveNumberOfDofManagers(), _eq;

    for ( int idofman = 1; idofman <= ndofman; idofman++ ) {
        DofManager *_dm = domain->giveDofManager(idofman);
        for ( Dof *_dof: *_dm ) {
            if ( _dof->isPrimaryDof() ) {
                if ( ( _eq = _dof->__giveEquationNumber() ) ) {
                    // pack values in solution vectors
                    _dof->updateUnknownsDictionary( tStep, VM_Total, totalDisplacement.at(_eq) );
                }
            }
        }
    }
}

void
NonLinearDynamic :: unpackMigratingData(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    int ndofman = domain->giveNumberOfDofManagers(), _eq;
    //int myrank = this->giveRank();

    // resize target arrays
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
    totalDisplacement.resize(neq);
    incrementOfDisplacement.resize(neq);

    for ( int idofman = 1; idofman <= ndofman; idofman++ ) {
        DofManager *_dm = domain->giveDofManager(idofman);
        for ( Dof *_dof: *_dm ) {
            if ( _dof->isPrimaryDof() ) {
                if ( ( _eq = _dof->__giveEquationNumber() ) ) {
                    // pack values in solution vectors
                    totalDisplacement.at(_eq) = _dof->giveUnknownsDictionaryValue( tStep, VM_Total );
 #if 0
                    // debug print
                    if ( _dm->giveParallelMode() == DofManager_shared ) {
                        fprintf(stderr, "[%d] Shared: %d(%d) -> %d\n", myrank, idofman, idof, _eq);
                    } else {
                        fprintf(stderr, "[%d] Local : %d(%d) -> %d\n", myrank, idofman, idof, _eq);
                    }

 #endif 
                }
            }
        }
    }

    this->initializeCommMaps(true);
    nMethod->reinitialize();
    // reinitialize error estimator (if any)
    if ( this->giveDomainErrorEstimator(1) ) {
        this->giveDomainErrorEstimator(1)->reinitialize();
    }

    initFlag = true;
}

} // end namespace oofem
