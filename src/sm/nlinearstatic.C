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

#include "nlinearstatic.h"
#include "structuralelement.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "error.h"
#include "verbose.h"
#include "sparsenonlinsystemnm.h"
#include "nrsolver.h"
#include "nrsolver2.h"
#include "calmls.h"
#include "outputmanager.h"
#include "datastream.h"
#include "usrdefsub.h"
#include "clock.h"
#include "contextioerr.h"
#include "sparsemtrx.h"
#include "errorestimator.h"
#include "mathfem.h"

namespace oofem {
NonLinearStatic :: NonLinearStatic(int i, EngngModel *_master) : LinearStatic(i, _master),
    totalDisplacement(), incrementOfDisplacement(), internalForces(), initialLoadVector(), incrementalLoadVector(),
    initialLoadVectorOfPrescribed(), incrementalLoadVectorOfPrescribed()
{
    //
    // constructor
    //
    prevStepLength = 0.;
    currentStepLength = 0.;
    internalForcesEBENorm = 0.0;
    loadLevel = cumulatedLoadLevel = 0.;
    mstepCumulateLoadLevelFlag = 0;
    numMetStatus = NM_None;
    stiffMode = nls_tangentStiffness; // default
    internalVarUpdateStamp = 0;
    initFlag = loadInitFlag = 1;
    controlMode = nls_indirectControl;
    refLoadInputMode = SparseNonLinearSystemNM :: rlm_total;
    nMethod = NULL;

#ifdef __PARALLEL_MODE
    commMode = ProblemCommMode__NODE_CUT;
    nonlocalExt = 0;
    communicator = nonlocCommunicator = NULL;
    commBuff = NULL;
#endif
}


NonLinearStatic :: ~NonLinearStatic()
{
    //
    // destructor
    //
    if ( nMethod ) {
        delete nMethod;
    }
}


NumericalMethod *NonLinearStatic :: giveNumericalMethod(MetaStep *mStep)
{
    const char *__proc = "giveNumericalMethod"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                     // Required by IR_GIVE_FIELD macro

    if ( mStep == NULL ) {
        _error("giveNumericalMethod: undefined meta step");
    }

    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD( ( mStep->giveAttributesRecord() ), _val, IFT_NonLinearStatic_controlmode, "controlmode" ); // Macro
    IR_GIVE_OPTIONAL_FIELD( ( mStep->giveAttributesRecord() ), _val, IFT_NonLinearStatic_controlmode, "controllmode" ); // for backward compatibility
    NonLinearStatic_controlType mode = ( NonLinearStatic_controlType ) _val;

    SparseNonLinearSystemNM *nm = NULL;
    if ( mode == nls_indirectControl ) {
        if ( nMethod ) {
            if ( nMethod->giveClassID() == CylindricalALMSolverClass ) {
                return nMethod;
            } else {
                delete nMethod;
            }
        }

        nm = ( SparseNonLinearSystemNM * ) new CylindricalALM(1, this->giveDomain(1), this, EID_MomentumBalance);
        nMethod = nm;
    } else if ( mode == nls_directControl ) {
        if ( nMethod ) {
            if ( nMethod->giveClassID() == NRSolverClass ) {
                return nMethod;
            } else {
                delete nMethod;
            }
        }

        nm = ( SparseNonLinearSystemNM * ) new NRSolver(1, this->giveDomain(1), this, EID_MomentumBalance);
        nMethod = nm;
    } else if ( mode == nls_directControl2 ) {
        if ( nMethod ) {
            if ( nMethod->giveClassID() == NRSolverClass ) {
                return nMethod;
            } else {
                delete nMethod;
            }
        }

        nm = ( SparseNonLinearSystemNM * ) new NRSolver2(1, this->giveDomain(1), this, EID_MomentumBalance);
        nMethod = nm;
    } else {
        _error("giveNumericalMethod: unsupported controlMode");
    }

    return nm;
}


void
NonLinearStatic :: updateAttributes(MetaStep *mStep)
{
    const char *__proc = "updateAttributes"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                  // Required by IR_GIVE_FIELD macro

    MetaStep *mStep1 = this->giveMetaStep( mStep->giveNumber() );//this line ensures correct input file in staggered problem
    InputRecord *ir = mStep1->giveAttributesRecord();

    LinearStatic :: updateAttributes(mStep1);

    /*
     * if ((mstep->giveFirstStepNumber() == atTime->giveNumber()) && hasString(initString, "fixload")) {
     * double factor;
     *
     * printf ("NonLinearStatic: fixed load level");
     * if (initialLoadVector.isEmpty()) initialLoadVector.resize(loadVector.giveSize());
     * if ((controlMode == nls_directControl) || (controlMode == nls_directControl2)) factor = 1.0;
     * else factor = loadLevel;
     * loadVector.times (factor);
     * initialLoadVector.add(loadVector);
     * loadVector.zero();
     * this->loadInitFlag = 1;
     * this->loadLevel = 0.0;
     * }
     */
    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NonLinearStatic_controlmode, "controlmode"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NonLinearStatic_controlmode, "controllmode"); // for backward compatibility
    controlMode = ( NonLinearStatic_controlType ) _val;

    deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, IFT_NonLinearStatic_deltat, "deltat"); // Macro
    if ( deltaT < 0. ) {
        _error("updateAttributes: deltaT < 0");
    }

    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NonLinearStatic_stiffmode, "stiffmode"); // Macro
    stiffMode = ( NonLinearStatic_stifnessMode ) _val;

    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NonLinearStatic_refloadmode, "refloadmode"); // Macro
    refLoadInputMode = ( SparseNonLinearSystemNM :: referenceLoadInputModeType ) _val;

    if ( ir->hasField(IFT_NonLinearStatic_keepll, "keepll") ) {
        mstepCumulateLoadLevelFlag = true;
    } else {
        mstepCumulateLoadLevelFlag = false;
    }

    // called just to mart filed as recognized, used later
    ir->hasField(IFT_NonLinearStatic_donotfixload, "donotfixload");
}


IRResultType
NonLinearStatic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LinearStatic :: initializeFrom(ir);
    nonlocalStiffnessFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nonlocalStiffnessFlag, IFT_NonLinearStatic_nonlocstiff, "nonlocstiff"); // Macro

#ifdef __PARALLEL_MODE
    //if (ir->hasField ("nodecutmode")) commMode = ProblemCommunicator::ProblemCommMode__NODE_CUT;
    //else if (ir->hasField ("elementcutmode")) commMode = ProblemCommunicator::ProblemCommMode__ELEMENT_CUT;
    //else _error ("instanciateFrom: ProblemCommunicator comm mode not specified");
    if ( isParallel() ) {
        //commBuff = new CommunicatorBuff (this->giveNumberOfProcesses(), CBT_dynamic);
        commBuff = new CommunicatorBuff(this->giveNumberOfProcesses(), CBT_static);
        communicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                               this->giveNumberOfProcesses(),
                                               this->commMode);

        if ( ir->hasField(IFT_NonLinearStatic_nonlocalext, "nonlocalext") ) {
            nonlocalExt = 1;
            nonlocCommunicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                                         this->giveNumberOfProcesses(),
                                                         ProblemCommMode__REMOTE_ELEMENT_MODE);
        }
    }
#endif

    return IRRT_OK;
}


double NonLinearStatic :: giveUnknownComponent(EquationID chc, ValueModeType mode,
                                               TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
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
            // return incrementOfDisplacement -> at(eq);
            // return nMethod-> giveUnknownComponent(IncrementOfSolution, eq);
            if ( incrementOfDisplacement.isNotEmpty() ) {
                return incrementOfDisplacement.at(eq);
            } else {
                return 0.;
            }

        case VM_Total:
            if ( totalDisplacement.isNotEmpty() ) {
                return totalDisplacement.at(eq);
            } else {
                return 0.;
            }

        default:
            _error("giveUnknownComponent: Unknown is of undefined ValueModeType for this problem");
        }
    } else {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    return 0.0;
}


double NonLinearStatic :: giveUnknownComponent(UnknownType chc, ValueModeType mode,
                                               TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }


    if ( tStep != this->giveCurrentStep() ) {
        _error("giveUnknownComponent: unknown time step encountered");
        return 0.;
    }

    if ( chc == TotalLoadLevel ) {
        return loadLevel;
    } else {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    return 0.0;
}


TimeStep *NonLinearStatic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    int mstepNum = 1;
    double totalTime = 0.0;
    StateCounterType counter = 1;
    double deltaTtmp = deltaT;

    //do not increase deltaT on microproblem
    if ( pScale == microScale ) {
        deltaTtmp = 0.;
    }

    if (previousStep != NULL){
        delete previousStep;
    }

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
    // dt variable are set eq to 0 for statics - has no meaning
    // *Wrong* It has meaning for viscoelastic materials.

    return currentStep;
}


void NonLinearStatic :: solveYourself()
{
#ifdef __PARALLEL_MODE
 #ifdef __VERBOSE_PARALLEL
    // force equation numbering before setting up comm maps
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);
    OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
 #endif

    // set up communication patterns
    this->initializeCommMaps();
    // init remote dofman list
    // this->initRemoteDofManList ();
#endif

    StructuralEngngModel :: solveYourself();
}


void
NonLinearStatic :: solveYourselfAt(TimeStep *tStep)
{
    proceedStep(1, tStep);
}


void
NonLinearStatic :: terminate(TimeStep *tStep)
{
    this->doStepOutput(tStep);
    this->printReactionForces(tStep, 1);
    // update load vectors before storing context
    this->updateLoadVectors(tStep);
    this->saveStepContext(tStep);
}


void
NonLinearStatic :: updateLoadVectors(TimeStep *stepN)
{
    MetaStep *mstep = this->giveMetaStep( stepN->giveMetaStepNumber() );
    bool isLastMetaStep = ( stepN->giveNumber() == mstep->giveLastStepNumber() );

    if ( controlMode == nls_indirectControl ) {
        //if ((stepN->giveNumber() == mstep->giveLastStepNumber()) && ir->hasField("fixload")) {
        if ( isLastMetaStep ) {
            if ( !mstep->giveAttributesRecord()->hasField(IFT_NonLinearStatic_donotfixload, "donotfixload") ) {
                OOFEM_LOG_INFO("Fixed load level\n");

                //update initialLoadVector
                if ( initialLoadVector.isEmpty() ) {
                    initialLoadVector.resize( incrementalLoadVector.giveSize() );
                }

                incrementalLoadVector.times(loadLevel);
                initialLoadVector.add(incrementalLoadVector);

                incrementalLoadVectorOfPrescribed.times(loadLevel);
                initialLoadVectorOfPrescribed.add(incrementalLoadVectorOfPrescribed);

                incrementalLoadVector.zero();
                incrementalLoadVectorOfPrescribed.zero();

                this->loadInitFlag = 1;
            }

            //if (!mstep->giveAttributesRecord()->hasField("keepll")) this->loadLevelInitFlag = 1;
        }
    } else { // direct control
        //update initialLoadVector after each step of direct control
        //(here the loading is not proportional)
        if ( initialLoadVector.isEmpty() ) {
            initialLoadVector.resize( incrementalLoadVector.giveSize() );
        }

        OOFEM_LOG_DEBUG("Fixed load level\n");

        incrementalLoadVector.times(loadLevel);
        initialLoadVector.add(incrementalLoadVector);

        incrementalLoadVectorOfPrescribed.times(loadLevel);
        initialLoadVectorOfPrescribed.add(incrementalLoadVectorOfPrescribed);

        incrementalLoadVector.zero();
        incrementalLoadVectorOfPrescribed.zero();

        this->loadInitFlag = 1;
    }


    // if (isLastMetaStep) {
    if ( isLastMetaStep && !mstep->giveAttributesRecord()->hasField(IFT_NonLinearStatic_donotfixload, "donotfixload") ) {
#ifdef VERBOSE
        OOFEM_LOG_INFO("Reseting load level\n");
#endif
        if ( mstepCumulateLoadLevelFlag ) {
            cumulatedLoadLevel += loadLevel;
        } else {
            cumulatedLoadLevel = 0.0;
        }

        this->loadLevel = 0.0;
    }
}


void
NonLinearStatic :: proceedStep(int di, TimeStep *tStep)
{
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step
    //

    if ( initFlag ) {
        //
        // first step  create space for stiffness Matrix
        //
        int neq = this->giveNumberOfEquations(EID_MomentumBalance);
        internalForces.resize(neq);
        internalForces.zero();

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
    }

#if 0
     if ((mstep->giveFirstStepNumber() == tStep->giveNumber())) {
#ifdef VERBOSE
        OOFEM_LOG_INFO("Resetting load level\n");
#endif
        if (mstepCumulateLoadLevelFlag) cumulatedLoadLevel += loadLevel;
        else cumulatedLoadLevel = 0.0;
        this->loadLevel = 0.0;
     }
#endif

    if ( loadInitFlag || ( controlMode == nls_directControl ) || ( controlMode == nls_directControl2 ) ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling reference load\n");
#endif
        //
        // assemble the incremental reference load vector
        //
        this->assembleIncrementalReferenceLoadVectors(incrementalLoadVector, incrementalLoadVectorOfPrescribed,
                                                      refLoadInputMode, this->giveDomain(di), EID_MomentumBalance, tStep);

        loadInitFlag = 0;
    }

    if ( tStep->giveNumber() == 1 ) {
        int neq = this->giveNumberOfEquations(EID_MomentumBalance);
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();
    }

    //
    //    ->   BEGINNING OF LOAD (OR DISPLACEMENT) STEP  <-
    //
    incrementOfDisplacement.zero();

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    //
    // call numerical model to solve arise problem
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving       [step number %5d.%d]\n\n", tStep->giveNumber(), tStep->giveVersion() );
#endif

    if ( initialLoadVector.isNotEmpty() ) {
        numMetStatus = nMethod->solve(stiffnessMatrix, & incrementalLoadVector, & initialLoadVector,
                                      & totalDisplacement, & incrementOfDisplacement, & internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep);
    } else {
        numMetStatus = nMethod->solve(stiffnessMatrix, & incrementalLoadVector, NULL,
                                      & totalDisplacement, & incrementOfDisplacement, & internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep);
    }

    ///@todo Use temporary variables. updateYourself() should set the final values, while proceedStep should be callable multiple times for each step (if necessary). / Mikael
    OOFEM_LOG_RELEVANT("Equilibrium reached at load level = %f in %d iterations\n", cumulatedLoadLevel + loadLevel, currentIterations);
    prevStepLength =  currentStepLength;
}


void
NonLinearStatic :: updateYourself(TimeStep *stepN)
{
    //
    // The following line is potentially serious performance leak.
    // The numerical method may compute their internal forces - thus causing
    // internal state to be updated, while checking equilibrium.
    // update internal state only if necessary
    this->updateInternalState(stepN);
    StructuralEngngModel :: updateYourself(stepN);
}


void
NonLinearStatic ::  updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
//
// updates some component, which is used by numerical method
// to newly reached state. used mainly by numerical method
// when new tangent stiffness is needed during finding
// of new equilibrium stage.
//
{
    switch ( cmpn ) {
    case NonLinearLhs:
        if ( stiffMode == nls_tangentStiffness ) {
            stiffnessMatrix->zero(); // zero stiffness matrix
#ifdef VERBOSE
            OOFEM_LOG_DEBUG("Assembling tangent stiffness matrix\n");
#endif
            this->assemble(stiffnessMatrix, tStep, EID_MomentumBalance, TangentStiffnessMatrix,
                           EModelDefaultEquationNumbering(), d);
        } else if ( ( stiffMode == nls_secantStiffness ) || ( stiffMode == nls_secantInitialStiffness && initFlag) ) {
#ifdef VERBOSE
            OOFEM_LOG_DEBUG("Assembling secant stiffness matrix\n");
#endif
            stiffnessMatrix->zero(); // zero stiffness matrix
            this->assemble(stiffnessMatrix, tStep, EID_MomentumBalance, SecantStiffnessMatrix,
                           EModelDefaultEquationNumbering(), d);
            initFlag = 0;
        } else if ( ( stiffMode == nls_elasticStiffness ) && ( initFlag ||
                ( this->giveMetaStep( tStep->giveMetaStepNumber() )->giveFirstStepNumber() == tStep->giveNumber() ) ) ) {
#ifdef VERBOSE
            OOFEM_LOG_DEBUG("Assembling elastic stiffness matrix\n");
#endif
            stiffnessMatrix->zero(); // zero stiffness matrix
            this->assemble( stiffnessMatrix, tStep, EID_MomentumBalance, ElasticStiffnessMatrix,
                            EModelDefaultEquationNumbering(), d);
            initFlag = 0;
        } else {
            // currently no action , this method is mainly intended to
            // assemble new tangent stiffness after each iteration
            // when secantStiffMode is on, we use the same stiffness
            // during iteration process
        }

        break;
    case InternalRhs:
        // update internalForces and internalForcesEBENorm concurrently
        this->giveInternalForces(internalForces, true, 1, tStep);
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
NonLinearStatic :: printOutputAt(FILE *File, TimeStep *stepN)
{
    if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(stepN) ) {
        return;                                                                      // do not print even Solution step header
    }

    fprintf( File, "\n\nOutput for time % .3e, solution step number %d\n", stepN->giveTargetTime(), stepN->giveNumber() );
    fprintf(File, "Reached load level : %20.6f in %d iterations\n\n",
            cumulatedLoadLevel + loadLevel, currentIterations);

    nMethod->printState(File);

    this->giveDomain(1)->giveOutputManager()->doDofManOutput(File, stepN);
    this->giveDomain(1)->giveOutputManager()->doElementOutput(File, stepN);
}


contextIOResultType
NonLinearStatic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //if ((iores = this->giveNumericalMethod(giveCurrentStep())->saveContext (stream)) != CIO_OK) THROW_CIOERR(iores);

    if ( ( iores = totalDisplacement.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = incrementOfDisplacement.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    int _cm = controlMode;
    if ( !stream->write(& _cm, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& loadLevel, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& cumulatedLoadLevel, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // store InitialLoadVector
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
    } // ensure consistent records

    return CIO_OK;
}


contextIOResultType
NonLinearStatic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    int closeFlag = 0;
    int istep, iversion;
    contextIOResultType iores;
    FILE *file;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);
    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    // save element context
    if ( ( iores = EngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //if ((iores = this->giveNumericalMethod(giveCurrentStep())->restoreContext (stream)) !=CIO_OK) THROW_CIOERR(iores);

    if ( ( iores = totalDisplacement.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = incrementOfDisplacement.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    int _cm;
    if ( !stream->read(& _cm, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    controlMode = ( NonLinearStatic_controlType ) _cm;
    if ( !stream->read(& loadLevel, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& cumulatedLoadLevel, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    // store InitialLoadVector
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
NonLinearStatic :: updateDomainLinks()
{
    LinearStatic :: updateDomainLinks();

    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
#ifdef __PARALLEL_MODE
    if ( this->giveLoadBalancer() ) {
        this->giveLoadBalancer()->setDomain( this->giveDomain(1) );
    }

#endif
}


#ifdef __PETSC_MODULE
void
NonLinearStatic :: initPetscContexts()
{
    petscContextList->growTo(ndomains);
    for ( int i = 1; i <= this->ndomains; i++ ) {
        petscContextList->put(i, new PetscContext(this, EID_MomentumBalance, false)); // false == using local vectors.
    }
}
#endif


void
NonLinearStatic :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID ut, CharType type,
                            const UnknownNumberingScheme &s, Domain *domain)
{
#ifdef TIME_REPORT
    oofem_timeval tstart;
    getUtime(tstart);
#endif

    LinearStatic :: assemble(answer, tStep, ut, type, s, domain);

    if ( ( nonlocalStiffnessFlag ) && ( type == TangentStiffnessMatrix ) ) {
        // add nonlocal contribution
        int ielem, nelem = domain->giveNumberOfElements();
        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            ( ( StructuralElement * ) ( domain->giveElement(ielem) ) )->addNonlocalStiffnessContributions(* answer, s, tStep);
        }

        // print storage statistics
        answer->printStatistics();
    }

#ifdef TIME_REPORT
    oofem_timeval tfin;
    getRelativeUtime(tfin, tstart);
    OOFEM_LOG_DEBUG( "NonLinearStatic: User time consumed by assembly: %.2fs\n",
                   ( double ) ( tfin.tv_sec + tfin.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif
}


#ifdef __OOFEG
void
NonLinearStatic :: showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime)
{
    Domain *domain = this->giveDomain(1);
    CharType ctype;
    int i;

    if ( type != 1 ) {
        return;
    }

    if ( stiffMode == nls_tangentStiffness ) {
        ctype = TangentStiffnessMatrix;
    } else if ( stiffMode == nls_secantStiffness ) {
        ctype = SecantStiffnessMatrix;
    } else {
        ctype = SecantStiffnessMatrix;
    }

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
NonLinearStatic :: computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    if ( ( di == 1 ) && ( tStep == this->giveCurrentStep() ) ) {
        reactions = incrementalLoadVectorOfPrescribed;
        reactions.times(loadLevel);
        reactions.add(initialLoadVectorOfPrescribed);
    } else {
        _error("computeExternalLoadReactionContribution: unable to respond due to invalid solution step or domain");
    }
}


void
NonLinearStatic :: assembleIncrementalReferenceLoadVectors(FloatArray &_incrementalLoadVector,
                                                           FloatArray &_incrementalLoadVectorOfPrescribed,
                                                           SparseNonLinearSystemNM :: referenceLoadInputModeType _refMode,
                                                           Domain *sourceDomain, EquationID ut, TimeStep *tStep)
{
    _incrementalLoadVector.resize( sourceDomain->giveEngngModel()->giveNumberOfEquations(EID_MomentumBalance) );
    _incrementalLoadVector.zero();
    _incrementalLoadVectorOfPrescribed.resize( sourceDomain->giveEngngModel()->giveNumberOfPrescribedEquations(EID_MomentumBalance) );
    _incrementalLoadVectorOfPrescribed.zero();

    if ( _refMode == SparseNonLinearSystemNM :: rlm_incremental ) {
        ///@todo This was almost definitely wrong before. It never seems to be used. Is this code even relevant?
        this->assembleVector( _incrementalLoadVector, tStep, ut, ExternalForcesVector,
                              VM_Incremental, EModelDefaultEquationNumbering(), sourceDomain);

        this->assembleVector( _incrementalLoadVectorOfPrescribed, tStep, ut, ExternalForcesVector,
                              VM_Incremental, EModelDefaultPrescribedEquationNumbering(), sourceDomain);
    } else {
        this->assembleVector( _incrementalLoadVector, tStep, ut, ExternalForcesVector,
                              VM_Total, EModelDefaultEquationNumbering(), sourceDomain);

        this->assembleVector( _incrementalLoadVectorOfPrescribed, tStep, ut, ExternalForcesVector,
                              VM_Total, EModelDefaultPrescribedEquationNumbering(), sourceDomain);
    }

#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers( _incrementalLoadVector, LoadExchangeTag  );
#endif
}


#ifdef __PARALLEL_MODE
int
NonLinearStatic :: estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType)
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


LoadBalancer *
NonLinearStatic :: giveLoadBalancer()
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
NonLinearStatic :: giveLoadBalancerMonitor()
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
NonLinearStatic :: packMigratingData(TimeStep *atTime)
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
NonLinearStatic :: unpackMigratingData(TimeStep *atTime)
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
