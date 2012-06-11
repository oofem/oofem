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
//----------------------------------------------------------------------
#include <iostream>
using namespace std;
//----------------------------------------------------------------------

#include "nlineardynamic.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "error.h"
#include "verbose.h"
#include "sparsenonlinsystemnm.h"
#include "nrsolver.h"
#include "nlstructuralelement.h"
#include "structuralelementevaluator.h"
#include "outputmanager.h"
#include "datastream.h"
#include "usrdefsub.h"
#include "clock.h"
#include "contextioerr.h"
#include "sparsemtrx.h"
#include "errorestimator.h"

namespace oofem {
NonLinearDynamic :: NonLinearDynamic(int i, EngngModel *_master) : StructuralEngngModel(i, _master),
    totalDisplacement(), incrementOfDisplacement(), internalForces(), initialLoadVector(), incrementalLoadVector(),
    initialLoadVectorOfPrescribed(), incrementalLoadVectorOfPrescribed()
{
    ddtScheme    = euler3;
    numMetStatus = NM_None;

    ndomains                = 1;
    internalForcesEBENorm   = 0;
    internalVarUpdateStamp  = 0;
    initFlag = commInitFlag = 1;

    stiffnessMatrix  = NULL;
    nMethod          = NULL;

    refLoadInputMode = SparseNonLinearSystemNM :: rlm_total;
#ifdef __PARALLEL_MODE
    nonlocalExt = 0;

    communicator       = NULL;
    nonlocCommunicator = NULL;
    commBuff           = NULL;

    commMode = ProblemCommMode__NODE_CUT;
#endif
}


NonLinearDynamic :: ~NonLinearDynamic()
{
    if ( nMethod ) {
        delete nMethod;
    }
}


NumericalMethod *NonLinearDynamic :: giveNumericalMethod(MetaStep *mStep)
{
    //const char *__proc = "giveNumericalMethod"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                     // Required by IR_GIVE_FIELD macro

    if ( mStep == NULL ) {
        _error("giveNumericalMethod: undefined meta step");
    }

    SparseNonLinearSystemNM *nm = NULL;

    if ( nMethod ) {
        if ( nMethod->giveClassID() == NRSolverClass ) {
            nMethod->reinitialize();
            return nMethod;
        } else {
            delete nMethod;
        }
    }

    nm = ( SparseNonLinearSystemNM * ) new NRSolver(1, this->giveDomain(1), this, EID_MomentumBalance);
    nMethod = nm;

    return nm;
}


void
NonLinearDynamic :: updateAttributes(MetaStep *mStep)
{
    const char *__proc = "updateAttributes"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                  // Required by IR_GIVE_FIELD macro

    InputRecord *ir = mStep->giveAttributesRecord();

    StructuralEngngModel :: updateAttributes(mStep);

    deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, IFT_NonLinearDynamic_deltat, "deltat");
    if ( deltaT < 0. ) {
        _error("updateAttributes: deltaT < 0");
    }

    IR_GIVE_FIELD(ir, dumpingCoef, IFT_NonLinearDynamic_dumpcoef, "dumpcoef");

    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NonLinearDynamic_ddtScheme, "ddtscheme");

    ddtScheme = ( NonLinearDynamic_ddtScheme ) _val;

    if ( ddtScheme == newmark ) {
        OOFEM_LOG_INFO( "Selecting Newmark-beta metod\n" );
        IR_GIVE_FIELD(ir, alpha, IFT_NonLinearDynamic_alpha, "alpha");
        IR_GIVE_FIELD(ir, beta, IFT_NonLinearDynamic_beta, "beta");
    } else if ( ddtScheme == euler2 ) {
        OOFEM_LOG_INFO( "Selecting Two-point Backward Euler method\n" );
    } else if ( ddtScheme == euler3 ) {
        OOFEM_LOG_INFO( "Selecting Three-point Backward Euler metod\n" );
    } else {
        _error("NonLinearDynamic: Time-stepping scheme not found!\n");
    }

    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NonLinearDynamic_refloadmode, "refloadmode");
    refLoadInputMode = ( SparseNonLinearSystemNM :: referenceLoadInputModeType ) _val;
}


IRResultType
NonLinearDynamic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralEngngModel :: initializeFrom(ir);
    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_NonLinearDynamic_lstype, "lstype");
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_NonLinearDynamic_smtype, "smtype");
    sparseMtrxType = ( SparseMtrxType ) val;

    nonlocalStiffnessFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nonlocalStiffnessFlag, IFT_NonLinearDynamic_nonlocstiff, "nonlocstiff");

    deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, IFT_NonLinearDynamic_deltat, "deltat");
    if ( deltaT < 0. ) {
        _error("updateAttributes: deltaT < 0");
    }

#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        commBuff = new CommunicatorBuff(this->giveNumberOfProcesses(), CBT_static);
        communicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                               this->giveNumberOfProcesses(),
                                               this->commMode);

        if ( ir->hasField(IFT_NonLinearDynamic_nonlocalext, "nonlocalext") ) {
            nonlocalExt = 1;
            nonlocCommunicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                                         this->giveNumberOfProcesses(),
                                                         ProblemCommMode__REMOTE_ELEMENT_MODE);
        }
    }

#endif
    return IRRT_OK;
}


double NonLinearDynamic ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
                                                 TimeStep *tStep, Domain *d, Dof *dof)
{
    int eq = dof->__giveEquationNumber();
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }

    if ( tStep != this->giveCurrentStep() ) {
        _error("giveUnknownComponent: unknown time step encountered");
        return 0.;
    }

    if ( chc == EID_MomentumBalance ) {
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
            _error("giveUnknownComponent: Unknown is of undefined ValueModeType for this problem");
        }
    } else {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    return 0.0;
}


TimeStep *NonLinearDynamic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    int mstepNum = 1;
    StateCounterType counter = 1;
    double deltaTtmp = deltaT;
    double totalTime = deltaT;

    //Do not increase deltaT on microproblem
    if ( pScale == microScale ) {
        deltaTtmp = 0.;
    }

    delete previousStep;
    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTargetTime() + deltaTtmp;
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
        mstepNum = currentStep->giveMetaStepNumber();

        if ( !this->giveMetaStep(mstepNum)->isStepValid(istep) ) {
            mstepNum++;
            if ( mstepNum > nMetaSteps ) {
                OOFEM_ERROR3("giveNextStep: no next step available, mstepNum=%d > nMetaSteps=%d", mstepNum, nMetaSteps);
            }
        }
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, mstepNum, totalTime, deltaTtmp, counter);

    return currentStep;
}


void NonLinearDynamic :: solveYourself()
{
    if ( commInitFlag ) {
#ifdef __PARALLEL_MODE
 #ifdef __VERBOSE_PARALLEL
        // Force equation numbering before setting up comm maps.
        int neq = this->giveNumberOfEquations(EID_MomentumBalance);
        OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
 #endif
        if ( isParallel() ) {
            // Set up communication patterns.
            this->initializeCommMaps();
        }
#endif
        commInitFlag = 0;
    }

    StructuralEngngModel :: solveYourself();
}


void
NonLinearDynamic :: solveYourselfAt(TimeStep *tStep) {
    if ( commInitFlag ) {
#ifdef __PARALLEL_MODE
 #ifdef __VERBOSE_PARALLEL
        // Force equation numbering before setting up comm maps.
        int neq = this->giveNumberOfEquations(EID_MomentumBalance);
        OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
 #endif
        if ( isParallel() ) {
            // Set up communication patterns.
            this->initializeCommMaps();
        }
#endif
        commInitFlag = 0;
    }

#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving       [Step number %8d, Time %15e]\n\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    deltaT = tStep->giveTimeIncrement();

    proceedStep(1, tStep);
}

void
NonLinearDynamic :: terminate(TimeStep *tStep) {
    this->doStepOutput(tStep);
    this->printReactionForces(tStep, 1);
    this->saveStepContext(tStep);
}

void
NonLinearDynamic :: proceedStep(int di, TimeStep *tStep)
{
    // creates system of governing eq's and solves them at given time step
    // first assemble problem at current time step

    int neq = this->giveNumberOfEquations(EID_MomentumBalance);

    // Time-stepping constants
    double dt2 = deltaT * deltaT;

    if ( ddtScheme == newmark ) {
        OOFEM_LOG_DEBUG("Solving using Newmark-beta method\n");
        a0 = 1 / ( beta * dt2 );
        a1 = alpha / ( beta * deltaT );
        a2 = 1 / ( beta * deltaT );
        a3 = 1 / ( 2 *  beta ) - 1;
        a4 = ( alpha / beta ) - 1;
        a5 = deltaT / 2 * ( alpha / beta - 2 );
        a6 = 0;
    } else if ( ( ddtScheme == euler2 ) || ( tStep->giveNumber() == giveNumberOfFirstStep() ) ) {
        if ( ddtScheme != euler3 ) {
            OOFEM_LOG_DEBUG("Solving using Backward Euler method\n");
        } else {
            OOFEM_LOG_DEBUG("Solving initial step using Three-point Backward Euler method\n");
        }
        a0 = 1 / dt2;
        a1 = 1 / deltaT;
        a2 = 1 / deltaT;
        a3 = 0;
        a4 = 0;
        a5 = 0;
        a6 = 0;
    } else if ( ddtScheme == euler3 ) {
        OOFEM_LOG_DEBUG("Solving using Three-point Backward Euler method\n");
        a0 = 2 / dt2;
        a1 = 3 / ( 2 * deltaT );
        a2 = 2 / deltaT;
        a3 = 0;
        a4 = 0;
        a5 = 0;
        a6 = 1 / ( 2 * deltaT );
    } else {
        _error("NonLinearDynamic: Time-stepping scheme not found!\n")
    }

    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        // Initialization
        previousIncrementOfDisplacement.resize(neq);
        previousIncrementOfDisplacement.zero();
        previousTotalDisplacement.resize(neq);
        previousTotalDisplacement.zero();
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        previousInternalForces.resize(neq);
        previousInternalForces.zero();
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();
        velocityVector.resize(neq);
        velocityVector.zero();
        accelerationVector.resize(neq);
        accelerationVector.zero();

        TimeStep *stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0,
                                                 -deltaT, deltaT, 0);

        int nDofs, j, k, jj;
        int nman  = this->giveDomain(di)->giveNumberOfDofManagers();
        DofManager *node;
        Dof *iDof;

        // Considering initial conditions.
        for ( j = 1; j <= nman; j++ ) {
            node = this->giveDomain(di)->giveDofManager(j);
            nDofs = node->giveNumberOfDofs();

            for ( k = 1; k <= nDofs; k++ ) {
                // Ask for initial values obtained from
                // bc (boundary conditions) and ic (initial conditions).
                iDof  =  node->giveDof(k);
                if ( !iDof->isPrimaryDof() ) {
                    continue;
                }

                jj = iDof->__giveEquationNumber();
                if ( jj ) {
                    incrementOfDisplacement.at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Total, stepWhenIcApply);
                    velocityVector.at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Velocity, stepWhenIcApply);
                    accelerationVector.at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Acceleration, stepWhenIcApply);
                }
            }
        }
    } else {
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();
    }

    if ( initFlag ) {
        // First assemble problem at current time step.
        // Option to take into account initial conditions.
        if ( !stiffnessMatrix ) {
            stiffnessMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
        }

        if ( stiffnessMatrix == NULL ) {
            _error("proceedStep: sparse matrix creation failed");
        }

        if ( nonlocalStiffnessFlag ) {
            if ( !stiffnessMatrix->isAsymmetric() ) {
                _error("proceedStep: stiffnessMatrix does not support asymmetric storage");
            }
        }

        stiffnessMatrix->buildInternalStructure( this, di, EID_MomentumBalance, EModelDefaultEquationNumbering() );
        // Initialize vectors
        help.resize(neq);
        rhs.resize(neq);
        rhs2.resize(neq);
        internalForces.resize(neq);
        help.zero();
        rhs.zero();
        rhs2.zero();

        previousTotalDisplacement.resize(neq);
        for ( int i = 1; i <= neq; i++ ) {
            previousTotalDisplacement.at(i) = totalDisplacement.at(i);
        }
        initFlag = 0;
    }

#ifdef VERBOSE
    OOFEM_LOG_DEBUG("Assembling load\n");
#endif

    // Assemble the incremental reference load vector.
    this->assembleIncrementalReferenceLoadVectors(incrementalLoadVector, incrementalLoadVectorOfPrescribed,
                                                  refLoadInputMode, this->giveDomain(di), EID_MomentumBalance, tStep);

    // Assembling the effective load vector
    for ( int i = 1; i <= neq; i++ ) {
        help.at(i) = ( a2 * velocityVector.at(i) +
                       a3 * accelerationVector.at(i) ) +
                       dumpingCoef * (
                           a4 * velocityVector.at(i) +
                           a5 * accelerationVector.at(i) +
                           a6 * previousIncrementOfDisplacement.at(i)
                       );
    }

    this->timesMassMtrx(help, rhs2, this->giveDomain(di), tStep);

    for ( int i = 1; i <= neq; i++ ) {
        rhs2.at(i) += incrementalLoadVector.at(i) - previousInternalForces.at(i);
        totalDisplacement.at(i) = previousTotalDisplacement.at(i);
    }

    //
    // Set-up numerical model.
    //
    this->giveNumericalMethod( this->giveCurrentMetaStep() );

    //
    // Call numerical model to solve problem.
    //
    double loadLevel = 1.0;
    if ( initialLoadVector.isNotEmpty() ) {
        numMetStatus = nMethod->solve(stiffnessMatrix, & rhs2, & initialLoadVector,
                                      & totalDisplacement, & incrementOfDisplacement, & internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep);
    } else {
        numMetStatus = nMethod->solve(stiffnessMatrix, & rhs2, NULL,
                                      & totalDisplacement, & incrementOfDisplacement, & internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep);
    }

    OOFEM_LOG_INFO("Equilibrium reached in %d iterations\n", currentIterations);
}
//----------------------------------------------------------------------
void
NonLinearDynamic :: giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                                    CharType type, TimeStep *tStep, Domain *domain)
{
    // We don't directly call element ->GiveCharacteristicMatrix() function, because some
    // engngm classes may require special modification of base types supported on
    // element class level.

    if ( type == ModifiedStiffnessMatrix ) {
        Element *element;
        FloatMatrix charMtrx2, charMtrx3;

        element = ( NLStructuralElement * ) domain->giveElement(num);

        element->giveCharacteristicMatrix(answer, StiffnessMatrix, tStep);
        element->giveCharacteristicMatrix(charMtrx2, MassMatrix, tStep);
        element->giveCharacteristicMatrix(charMtrx3, MassMatrix, tStep);

        // Basically, performing Task.4
        answer.add(this->a0, charMtrx2);
        answer.add(this->a1*dumpingCoef, charMtrx3);
        return;
    } else {
        StructuralEngngModel :: giveElementCharacteristicMatrix(answer, num, type, tStep, domain);
    }
}
//----------------------------------------------------------------------

void NonLinearDynamic :: updateYourself(TimeStep *stepN)
{
    int neq =  this->giveNumberOfEquations(EID_MomentumBalance);

    for ( int i = 1; i <= neq; i++ ) {
        rhs.at(i)                             = velocityVector.at(i);
        rhs2.at(i)                            = accelerationVector.at(i);
        accelerationVector.at(i)              = a0 * incrementOfDisplacement.at(i) - a2 * rhs.at(i) - a3 * rhs2.at(i);
        velocityVector.at(i)                  = a1 * incrementOfDisplacement.at(i) - a4 * rhs.at(i) - a5 * rhs2.at(i) - a6 * previousIncrementOfDisplacement.at(i);
        previousIncrementOfDisplacement.at(i) = totalDisplacement.at(i) - previousTotalDisplacement.at(i);
        previousTotalDisplacement.at(i)       = totalDisplacement.at(i);
    }

    this->giveInternalForces(previousInternalForces, true, 1, stepN);

    // The following line is potentially serious performance leak.
    // The numerical method may compute their internal forces - thus causing
    // internal state to be updated, while checking equilibrium.
    // update internal state only if necessary
    this->updateInternalState(stepN);
    StructuralEngngModel :: updateYourself(stepN);
}

void NonLinearDynamic ::  updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
//
// Updates some component, which is used by numerical method
// to newly reached state. used mainly by numerical method
// when new tanget stiffness is needed during finding
// of new equlibrium stage.
//
{
    int neq =  this->giveNumberOfEquations(EID_MomentumBalance);

    switch ( cmpn ) {
    case NonLinearLhs:
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Updating nonlinear LHS\n");
#endif
        stiffnessMatrix->zero();
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling stiffness matrix\n");
#endif
        this->assemble(stiffnessMatrix, tStep, EID_MomentumBalance, ModifiedStiffnessMatrix,
                       EModelDefaultEquationNumbering(), d);
        break;
    case InternalRhs:
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Updating internal RHS\n");
#endif
        this->giveInternalForces(internalForces, true, 1, tStep);

        // Updating the residual vector @ NR-solver
        for ( int i = 1; i <= neq; i++ ) {
            help.at(i) = a0 * incrementOfDisplacement.at(i) + dumpingCoef *a1 *incrementOfDisplacement.at(i);
        }

        this->timesMassMtrx(help, rhs, this->giveDomain(1), tStep);

        for ( int i = 1; i <= neq; i++ ) {
            internalForces.at(i) += rhs.at(i) - previousInternalForces.at(i);
        }

        break;
    case NonLinearRhs_Incremental:
        this->assembleIncrementalReferenceLoadVectors(incrementalLoadVector, incrementalLoadVectorOfPrescribed,
                                                      refLoadInputMode, d, EID_MomentumBalance, tStep);
        break;
    default:
        _error("updateComponent: Unknown Type of component.");
    }
}

void
NonLinearDynamic :: printOutputAt(FILE *File, TimeStep *stepN)
{
    //FILE* File = this -> giveDomain() -> giveOutputStream() ;

    if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(stepN) ) {
        return; // Do not print even Solution step header
    }

    fprintf(File, "\n\nOutput for time % .3e, solution step number %d\n", stepN->giveTargetTime(), stepN->giveNumber() );
    fprintf(File, "Equilibrium reached in %d iterations\n\n", currentIterations);


    nMethod->printState(File);

    this->giveDomain(1)->giveOutputManager()->doDofManOutput(File, stepN);
    this->giveDomain(1)->giveOutputManager()->doElementOutput(File, stepN);
}

void
NonLinearDynamic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    static char dofchar[] = "dva";
    static ValueModeType dofmodes[] = {
        VM_Total, VM_Velocity, VM_Acceleration
    };

    iDof->printMultipleOutputAt(stream, atTime, dofchar, EID_MomentumBalance, dofmodes, 3);
}

contextIOResultType NonLinearDynamic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file;

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

    if ( ( iores = totalDisplacement.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = previousInternalForces.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = previousIncrementOfDisplacement.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = initialLoadVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = initialLoadVectorOfPrescribed.storeYourself(stream, mode) ) != CIO_OK ) {
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
    FILE *file;

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

    if ( ( iores = totalDisplacement.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = previousInternalForces.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = previousIncrementOfDisplacement.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = initialLoadVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = initialLoadVectorOfPrescribed.restoreYourself(stream, mode) ) != CIO_OK ) {
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


#ifdef __PETSC_MODULE
void
NonLinearDynamic :: initPetscContexts()
{
    petscContextList->growTo(ndomains);
    for ( int i = 1; i <= this->ndomains; i++ ) {
        petscContextList->put( i, new PetscContext(this, EID_MomentumBalance, false) ); // false == using local vectors.
    }
}
#endif

void
NonLinearDynamic :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID ut, CharType type,
                             const UnknownNumberingScheme &s, Domain *domain)
{
#ifdef TIME_REPORT
    oofem_timeval tstart;
    getUtime(tstart);
#endif

    EngngModel :: assemble(answer, tStep, ut, type, s, domain);

    if ( ( nonlocalStiffnessFlag ) && ( type == TangentStiffnessMatrix ) ) {
        // Add nonlocal contribution.
        int ielem, nelem = domain->giveNumberOfElements();
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            ( ( NLStructuralElement * ) ( domain->giveElement(ielem) ) )->addNonlocalStiffnessContributions(* answer, s, tStep);
        }

        // Print storage statistics.
        answer->printStatistics();
    }

#ifdef TIME_REPORT
    oofem_timeval tfin;
    getRelativeUtime(tfin, tstart);
    OOFEM_LOG_DEBUG( "NonLinearDynamic info: user time consumed by assembly: %.2fs\n",
                     ( double ) ( tfin.tv_sec + tfin.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif
}

#ifdef __OOFEG
void
NonLinearDynamic :: showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime)
{
    Domain *domain = this->giveDomain(1);
    CharType ctype;
    int i;

    if ( type != 1 ) {
        return;
    }

    ctype = TangentStiffnessMatrix;

    int nelems = domain->giveNumberOfElements();
    for ( i = 1; i <= nelems; i++ ) {
        domain->giveElement(i)->showSparseMtrxStructure(ctype, context, atTime);
    }

    for ( i = 1; i <= nelems; i++ ) {
        domain->giveElement(i)->showExtendedSparseMtrxStructure(ctype, context, atTime);
    }
}
#endif

void
NonLinearDynamic :: computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    if ( ( di == 1 ) && ( tStep == this->giveCurrentStep() ) ) {
        reactions = incrementalLoadVectorOfPrescribed;
        reactions.add(initialLoadVectorOfPrescribed);
    } else {
        _error("computeExternalLoadReactionContribution: unable to respond due to invalid solution step or domain");
    }
}

void
NonLinearDynamic :: assembleIncrementalReferenceLoadVectors(FloatArray &_incrementalLoadVector,
                                                            FloatArray &_incrementalLoadVectorOfPrescribed,
                                                            SparseNonLinearSystemNM :: referenceLoadInputModeType _refMode,
                                                            Domain *sourceDomain, EquationID ut, TimeStep *tStep)
{
    EModelDefaultEquationNumbering en;

    _incrementalLoadVector.resize( sourceDomain->giveEngngModel()->giveNumberOfEquations(EID_MomentumBalance) );
    _incrementalLoadVector.zero();
    _incrementalLoadVectorOfPrescribed.resize( sourceDomain->giveEngngModel()->giveNumberOfPrescribedEquations(EID_MomentumBalance) );
    _incrementalLoadVectorOfPrescribed.zero();

    if ( _refMode == SparseNonLinearSystemNM :: rlm_incremental ) {
        this->assembleVector(_incrementalLoadVector, tStep, ut, ExternalForcesVector,
                             VM_Incremental, EModelDefaultEquationNumbering(), sourceDomain);

        this->assembleVector(_incrementalLoadVectorOfPrescribed, tStep, ut, ExternalForcesVector,
                             VM_Incremental, EModelDefaultPrescribedEquationNumbering(), sourceDomain);
    } else {
        this->assembleVector(_incrementalLoadVector, tStep, ut, ExternalForcesVector,
                             VM_Total, EModelDefaultEquationNumbering(), sourceDomain);
        this->assembleVector(_incrementalLoadVectorOfPrescribed, tStep, ut, ExternalForcesVector,
                             VM_Total, EModelDefaultPrescribedEquationNumbering(), sourceDomain);
    }

#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers(_incrementalLoadVector, LoadExchangeTag);
#endif
}

void
NonLinearDynamic :: timesMassMtrx(FloatArray &vec, FloatArray &answer, Domain *domain, TimeStep *tStep)
{
    int nelem = domain->giveNumberOfElements();
    //int neq = this->giveNumberOfEquations(EID_MomentumBalance);
    int i, j, k, jj, kk, n;
    FloatMatrix charMtrx;
    IntArray loc;
    Element *element;
    EModelDefaultEquationNumbering en;

    //    answer.resize(neq);
    answer.zero();
    for ( i = 1; i <= nelem; i++ ) {
        element = domain->giveElement(i);
#ifdef __PARALLEL_MODE
        // Skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif

        element->giveLocationArray(loc, EID_MomentumBalance, en);
        element->giveCharacteristicMatrix(charMtrx, MassMatrix, tStep);

        //
        // assemble it manually
        //
#ifdef DEBUG
        if ( ( n = loc.giveSize() ) != charMtrx.giveNumberOfRows() ) {
            _error("solveYourselfAt : dimension mismatch");
        }

#endif

        n = loc.giveSize();

        for ( j = 1; j <= n; j++ ) {
            jj = loc.at(j);
            if ( jj ) {
                for ( k = 1; k <= n; k++ ) {
                    kk = loc.at(k);
                    if ( kk ) {
                        answer.at(jj) += charMtrx.at(j, k) * vec.at(kk);
                    }
                }
            }
        }
    }

#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers(answer, MassExchangeTag);
#endif
}

#ifdef __PARALLEL_MODE
int
NonLinearDynamic :: estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType)
{
    int mapSize = commMap.giveSize();
    int i, j, ndofs, count = 0, pcount = 0;
    IntArray locationArray;
    Domain *domain = this->giveDomain(1);
    DofManager *dman;
    Dof *jdof;

    if ( packUnpackType == ProblemCommMode__ELEMENT_CUT ) {
        for ( i = 1; i <= mapSize; i++ ) {
            count += domain->giveDofManager( commMap.at(i) )->giveNumberOfDofs();
        }

        return ( buff.givePackSize(MPI_DOUBLE, 1) * count );
    } else if ( packUnpackType == ProblemCommMode__NODE_CUT ) {
        for ( i = 1; i <= mapSize; i++ ) {
            ndofs = ( dman = domain->giveDofManager( commMap.at(i) ) )->giveNumberOfDofs();
            for ( j = 1; j <= ndofs; j++ ) {
                jdof = dman->giveDof(j);
                if ( jdof->isPrimaryDof() && ( jdof->__giveEquationNumber() ) ) {
                    count++;
                } else {
                    pcount++;
                }
            }
        }

        //printf ("\nestimated count is %d\n",count);
        return ( buff.givePackSize(MPI_DOUBLE, 1) * max(count, pcount) );
    } else if ( packUnpackType == ProblemCommMode__REMOTE_ELEMENT_MODE ) {
        for ( i = 1; i <= mapSize; i++ ) {
            count += domain->giveElement( commMap.at(i) )->estimatePackSize(buff);
        }

        return count;
    }

    return 0;
}


void
NonLinearDynamic :: initializeCommMaps(bool forceInit)
{
    // set up communication patterns
    communicator->setUpCommunicationMaps(this, true, forceInit);
    if ( nonlocalExt ) {
        nonlocCommunicator->setUpCommunicationMaps(this, true, forceInit);
    }
}

#endif

int
NonLinearDynamic :: checkConsistency()
{
    // check internal consistency
    // if success returns nonzero
    int i, nelem;
    Element *ePtr;
    NLStructuralElement *sePtr;
    StructuralElementEvaluator *see;
    Domain *domain = this->giveDomain(1);

    nelem = domain->giveNumberOfElements();
    // check for proper element type

    for ( i = 1; i <= nelem; i++ ) {
        ePtr = domain->giveElement(i);
        sePtr = dynamic_cast< NLStructuralElement * >( ePtr );
        see   = dynamic_cast< StructuralElementEvaluator * >( ePtr );

        if ( ( sePtr == NULL ) && ( see == NULL ) ) {
            _warning2("checkConsistency: element %d has no Structural support", i);
            return 0;
        }
    }

    EngngModel :: checkConsistency();

    return 1;
}

#ifdef __PARALLEL_MODE
LoadBalancer *
NonLinearDynamic :: giveLoadBalancer()
{
    if ( lb ) {
        return lb;
    }

    if ( loadBalancingFlag ) {
        lb = CreateUsrDefLoadBalancerOfType( ParmetisLoadBalancerClass, this->giveDomain(1) );
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
        lbm = CreateUsrDefLoadBalancerMonitorOfType(WallClockLoadBalancerMonitorClass, this);
        return lbm;
    } else {
        return NULL;
    }
}

void
NonLinearDynamic :: packMigratingData(TimeStep *atTime)
{
    Domain *domain = this->giveDomain(1);
    int ndofman =  domain->giveNumberOfDofManagers(), ndofs, idofman, idof, _eq;
    DofManager *_dm;
    Dof *_dof;
    bool initialLoadVectorEmpty = initialLoadVector.isEmpty();
    bool initialLoadVectorOfPrescribedEmpty = initialLoadVectorOfPrescribed.isEmpty();

    for ( idofman = 1; idofman <= ndofman; idofman++ ) {
        _dm = domain->giveDofManager(idofman);
        ndofs = _dm->giveNumberOfDofs();
        for ( idof = 1; idof <= ndofs; idof++ ) {
            _dof = _dm->giveDof(idof);
            if ( _dof->isPrimaryDof() ) {
                if ( ( _eq = _dof->__giveEquationNumber() ) ) {
                    // pack values in solution vectors
                    _dof->updateUnknownsDictionary( atTime, EID_MomentumBalance, VM_Total, totalDisplacement.at(_eq) );
                    if ( initialLoadVectorEmpty ) {
                        _dof->updateUnknownsDictionary(atTime, EID_MomentumBalance, VM_RhsInitial, 0.0);
                    } else {
                        _dof->updateUnknownsDictionary( atTime, EID_MomentumBalance, VM_RhsInitial, initialLoadVector.at(_eq) );
                    }

                    _dof->updateUnknownsDictionary( atTime, EID_MomentumBalance, VM_RhsIncremental, incrementalLoadVector.at(_eq) );
                } else if ( ( _eq = _dof->__givePrescribedEquationNumber() ) ) {
                    // pack values in prescribed solution vectors
                    if ( initialLoadVectorOfPrescribedEmpty ) {
                        _dof->updateUnknownsDictionary(atTime, EID_MomentumBalance, VM_RhsInitial, 0.0);
                    } else {
                        _dof->updateUnknownsDictionary( atTime, EID_MomentumBalance, VM_RhsInitial, initialLoadVectorOfPrescribed.at(_eq) );
                    }

                    _dof->updateUnknownsDictionary( atTime, EID_MomentumBalance, VM_RhsIncremental, incrementalLoadVectorOfPrescribed.at(_eq) );
                }
            } // end primary dof
        } // end dof loop
    } // end dofman loop
}

void
NonLinearDynamic :: unpackMigratingData(TimeStep *atTime)
{
    Domain *domain = this->giveDomain(1);
    int ndofman =  domain->giveNumberOfDofManagers(), ndofs, idofman, idof, _eq;
    //int myrank = this->giveRank();
    DofManager *_dm;
    Dof *_dof;

    // resize target arrays
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);
    totalDisplacement.resize(neq);
    incrementOfDisplacement.resize(neq);
    incrementalLoadVector.resize(neq);
    initialLoadVector.resize(neq);
    initialLoadVectorOfPrescribed.resize( giveNumberOfPrescribedEquations(EID_MomentumBalance) );
    incrementalLoadVectorOfPrescribed.resize( giveNumberOfPrescribedEquations(EID_MomentumBalance) );

    for ( idofman = 1; idofman <= ndofman; idofman++ ) {
        _dm = domain->giveDofManager(idofman);
        ndofs = _dm->giveNumberOfDofs();
        for ( idof = 1; idof <= ndofs; idof++ ) {
            _dof = _dm->giveDof(idof);
            if ( _dof->isPrimaryDof() ) {
                if ( ( _eq = _dof->__giveEquationNumber() ) ) {
                    // pack values in solution vectors
                    _dof->giveUnknownsDictionaryValue( atTime, EID_MomentumBalance, VM_Total, totalDisplacement.at(_eq) );
                    _dof->giveUnknownsDictionaryValue( atTime, EID_MomentumBalance, VM_RhsInitial, initialLoadVector.at(_eq) );
                    _dof->giveUnknownsDictionaryValue( atTime, EID_MomentumBalance, VM_RhsIncremental, incrementalLoadVector.at(_eq) );

 #if 0
                    // debug print
                    if ( _dm->giveParallelMode() == DofManager_shared ) {
                        fprintf(stderr, "[%d] Shared: %d(%d) -> %d\n", myrank, idofman, idof, _eq);
                    } else {
                        fprintf(stderr, "[%d] Local : %d(%d) -> %d\n", myrank, idofman, idof, _eq);
                    }

 #endif
                } else if ( ( _eq = _dof->__givePrescribedEquationNumber() ) ) {
                    // pack values in prescribed solution vectors
                    _dof->giveUnknownsDictionaryValue( atTime, EID_MomentumBalance, VM_RhsInitial, initialLoadVectorOfPrescribed.at(_eq) );
                    _dof->giveUnknownsDictionaryValue( atTime, EID_MomentumBalance, VM_RhsIncremental, incrementalLoadVectorOfPrescribed.at(_eq) );

 #if 0
                    // debug print
                    fprintf(stderr, "[%d] %d(%d) -> %d\n", myrank, idofman, idof, -_eq);
 #endif
                }
            } // end primary dof
        } // end dof loop
    } // end dofman loop

    this->initializeCommMaps(true);
    nMethod->reinitialize();
    // reinitialize error estimator (if any)
    if ( this->giveDomainErrorEstimator(1) ) {
        this->giveDomainErrorEstimator(1)->reinitialize();
    }

    initFlag = true;
}
#endif
} // end namespace oofem
