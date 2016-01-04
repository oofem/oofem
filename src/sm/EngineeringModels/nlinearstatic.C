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

#include "../sm/EngineeringModels/nlinearstatic.h"
#include "../sm/Elements/structuralelement.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "error.h"
#include "verbose.h"
#include "sparsenonlinsystemnm.h"
#include "nrsolver.h"
#include "calmls.h"
#include "outputmanager.h"
#include "datastream.h"
#include "classfactory.h"
#include "timer.h"
#include "contextioerr.h"
#include "sparsemtrx.h"
#include "errorestimator.h"
#include "mathfem.h"
#include "dofmanager.h"
#include "dof.h"
#include "unknownnumberingscheme.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "communicator.h"
 #include "loadbalancer.h"
#endif

namespace oofem {
REGISTER_EngngModel(NonLinearStatic);

NonLinearStatic :: NonLinearStatic(int i, EngngModel *_master) : LinearStatic(i, _master),
    totalDisplacement(), incrementOfDisplacement(), internalForces(), initialLoadVector(), incrementalLoadVector(),
    initialLoadVectorOfPrescribed(), incrementalLoadVectorOfPrescribed()
{
    //
    // constructor
    //
    prevStepLength = 0.;
    currentStepLength = 0.;
    loadLevel = cumulatedLoadLevel = 0.;
    mstepCumulateLoadLevelFlag = 0;
    numMetStatus = NM_None;
    stiffMode = nls_tangentStiffness; // default
    internalVarUpdateStamp = 0;
    initFlag = loadInitFlag = 1;
    controlMode = nls_indirectControl;
    refLoadInputMode = SparseNonLinearSystemNM :: rlm_total;
    nMethod = NULL;
    initialGuessType = IG_None;
}


NonLinearStatic :: ~NonLinearStatic()
{
    delete nMethod;
}


NumericalMethod *NonLinearStatic :: giveNumericalMethod(MetaStep *mStep)
{
    IRResultType result;                     // Required by IR_GIVE_FIELD macro

    if ( mStep == NULL ) {
        OOFEM_ERROR("undefined meta step");
    }

    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD( ( mStep->giveAttributesRecord() ), _val, _IFT_NonLinearStatic_controlmode );
    IR_GIVE_OPTIONAL_FIELD( ( mStep->giveAttributesRecord() ), _val, "controllmode" ); ///@todo If there is ever a major version change, remove this (for backward compatibility)
    NonLinearStatic_controlType mode = ( NonLinearStatic_controlType ) _val;

    if ( mode == nls_indirectControl ) {
        if ( nMethod ) {
            if ( dynamic_cast< CylindricalALM * >(nMethod) ) {
                return nMethod;
            } else {
                delete nMethod;
            }
        }

        this->nMethod = new CylindricalALM(this->giveDomain(1), this);
    } else if ( mode == nls_directControl ) {
        if ( nMethod ) {
            if ( dynamic_cast< NRSolver * >(nMethod) ) {
                return nMethod;
            } else {
                delete nMethod;
            }
        }

        this->nMethod = new NRSolver(this->giveDomain(1), this);
    } else {
        OOFEM_ERROR("unsupported controlMode");
    }

    return this->nMethod;
}


void
NonLinearStatic :: updateAttributes(MetaStep *mStep)
{
    IRResultType result;                  // Required by IR_GIVE_FIELD macro

    MetaStep *mStep1 = this->giveMetaStep( mStep->giveNumber() ); //this line ensures correct input file in staggered problem
    InputRecord *ir = mStep1->giveAttributesRecord();

    LinearStatic :: updateAttributes(mStep1);

    /*
     * if ((mstep->giveFirstStepNumber() == tStep->giveNumber()) && hasString(initString, "fixload")) {
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
    int _val = nls_indirectControl;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_NonLinearStatic_controlmode);
    IR_GIVE_OPTIONAL_FIELD(ir, _val, "controllmode"); /// @todo If there is ever a major version change, remove this (for backward compatibility)
    this->controlMode = ( NonLinearStatic_controlType ) _val;

    _val = IG_None;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_EngngModel_initialGuess);
    this->initialGuessType = ( InitialGuess ) _val;

    deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, _IFT_NonLinearStatic_deltat);
    if ( deltaT < 0. ) {
        OOFEM_ERROR("deltaT < 0");
    }

    _val = nls_tangentStiffness;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_NonLinearStatic_stiffmode);
    this->stiffMode = ( NonLinearStatic_stiffnessMode ) _val;

    _val = SparseNonLinearSystemNM :: rlm_total;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_NonLinearStatic_refloadmode);
    this->refLoadInputMode = ( SparseNonLinearSystemNM :: referenceLoadInputModeType ) _val;

    mstepCumulateLoadLevelFlag = ir->hasField(_IFT_NonLinearStatic_keepll);

    // called just to mark field as recognized, used later
    ir->hasField(_IFT_NonLinearStatic_donotfixload);
}


IRResultType
NonLinearStatic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = LinearStatic :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    nonlocalStiffnessFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nonlocalStiffnessFlag, _IFT_NonLinearStatic_nonlocstiff);

    updateElasticStiffnessFlag = false;
    if ( ir->hasField(_IFT_NonLinearStatic_updateElasticStiffnessFlag) ) {
      updateElasticStiffnessFlag = true;
    }
    
#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        //commBuff = new CommunicatorBuff (this->giveNumberOfProcesses(), CBT_dynamic);
        commBuff = new CommunicatorBuff(this->giveNumberOfProcesses(), CBT_static);
        communicator = new NodeCommunicator(this, commBuff, this->giveRank(),
                                            this->giveNumberOfProcesses());

        if ( ir->hasField(_IFT_NonLinearStatic_nonlocalext) ) {
            nonlocalExt = 1;
            nonlocCommunicator = new ElementCommunicator(this, commBuff, this->giveRank(),
                                                         this->giveNumberOfProcesses());
        }
    }
#endif

    return IRRT_OK;
}


double NonLinearStatic :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
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

    case VM_Velocity:
        if ( incrementOfDisplacement.isNotEmpty() ) {
            return incrementOfDisplacement.at(eq) / tStep->giveTimeIncrement();
        } else {
            return 0.;
        }

    default:
        OOFEM_ERROR("Unknown is of undefined ValueModeType for this problem");
    }

    return 0.0;
}

TimeStep *NonLinearStatic :: giveSolutionStepWhenIcApply(bool force)
{
  if ( master && (!force)) {
    return master->giveSolutionStepWhenIcApply();
  } else {
    if ( !stepWhenIcApply ) {
        int inin = giveNumberOfTimeStepWhenIcApply();
	//        int nFirst = giveNumberOfFirstStep();
        stepWhenIcApply.reset(new TimeStep(inin, this, 0, -deltaT, deltaT, 0));
    }

    return stepWhenIcApply.get();
  }
}


TimeStep *NonLinearStatic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    int mStepNum = 1;
    double totalTime = 0.0;
    StateCounterType counter = 1;
    double deltaTtmp = deltaT;

    //do not increase deltaT on microproblem
    if ( pScale == microScale ) {
        deltaTtmp = 0.;
    }

    if ( currentStep ) {
        totalTime = currentStep->giveTargetTime() + deltaTtmp;
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
        mStepNum = currentStep->giveMetaStepNumber();

        if ( !this->giveMetaStep(mStepNum)->isStepValid(istep) ) {
            mStepNum++;
            if ( mStepNum > nMetaSteps ) {
                OOFEM_ERROR("no next step available, mStepNum=%d > nMetaSteps=%d", mStepNum, nMetaSteps);
            }
        }
    } else {
        // first step -> generate initial step
        TimeStep *newStep = giveSolutionStepWhenIcApply();
        currentStep.reset(new TimeStep(*newStep));
    }

    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(istep, this, mStepNum, totalTime, deltaTtmp, counter) );
    // dt variable are set eq to 0 for statics - has no meaning
    // *Wrong* It has meaning for viscoelastic materials.

    return currentStep.get();
}


void NonLinearStatic :: solveYourself()
{
    if ( this->isParallel() ) {
 #ifdef __VERBOSE_PARALLEL
        // force equation numbering before setting up comm maps
        int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
        OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
 #endif

        // set up communication patterns
        this->initializeCommMaps();
        // init remote dofman list
        // this->initRemoteDofManList ();
    }
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
    fflush( this->giveOutputStream() );
    this->updateLoadVectors(tStep);
    this->saveStepContext(tStep);
}


void
NonLinearStatic :: updateLoadVectors(TimeStep *tStep)
{
    MetaStep *mstep = this->giveMetaStep( tStep->giveMetaStepNumber() );
    bool isLastMetaStep = ( tStep->giveNumber() == mstep->giveLastStepNumber() );

    if ( controlMode == nls_indirectControl ) {
        //if ((tStep->giveNumber() == mstep->giveLastStepNumber()) && ir->hasField("fixload")) {
        if ( isLastMetaStep ) {
            if ( !mstep->giveAttributesRecord()->hasField(_IFT_NonLinearStatic_donotfixload) ) {
                OOFEM_LOG_INFO("Fixed load level\n");

                //update initialLoadVector
                initialLoadVector.add(loadLevel, incrementalLoadVector);

                initialLoadVectorOfPrescribed.add(loadLevel, incrementalLoadVectorOfPrescribed);

                incrementalLoadVector.zero();
                incrementalLoadVectorOfPrescribed.zero();

                this->loadInitFlag = 1;
            }

            //if (!mstep->giveAttributesRecord()->hasField("keepll")) this->loadLevelInitFlag = 1;
        }
    } else { // direct control
        //update initialLoadVector after each step of direct control
        //(here the loading is not proportional)
        OOFEM_LOG_DEBUG("Fixed load level\n");

        initialLoadVector.add(loadLevel, incrementalLoadVector);

        initialLoadVectorOfPrescribed.add(loadLevel, incrementalLoadVectorOfPrescribed);

        incrementalLoadVector.zero();
        incrementalLoadVectorOfPrescribed.zero();

        this->loadInitFlag = 1;
    }


    // if (isLastMetaStep) {
    if ( isLastMetaStep && !mstep->giveAttributesRecord()->hasField(_IFT_NonLinearStatic_donotfixload) ) {
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
        if ( !stiffnessMatrix ) {
            stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        }

        if ( !stiffnessMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        if ( nonlocalStiffnessFlag ) {
            if ( !stiffnessMatrix->isAsymmetric() ) {
                OOFEM_ERROR("stiffnessMatrix does not support asymmetric storage");
            }
        }

        stiffnessMatrix->buildInternalStructure( this, di, EModelDefaultEquationNumbering() );
    }

#if 0
    if ( ( mstep->giveFirstStepNumber() == tStep->giveNumber() ) ) {
 #ifdef VERBOSE
        OOFEM_LOG_INFO("Resetting load level\n");
 #endif
        if ( mstepCumulateLoadLevelFlag ) {
            cumulatedLoadLevel += loadLevel;
        } else {
            cumulatedLoadLevel = 0.0;
        }
        this->loadLevel = 0.0;
    }
#endif

    if ( loadInitFlag || controlMode == nls_directControl ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling reference load\n");
#endif
        //
        // assemble the incremental reference load vector
        //
        this->assembleIncrementalReferenceLoadVectors(incrementalLoadVector, incrementalLoadVectorOfPrescribed,
                                                      refLoadInputMode, this->giveDomain(di), tStep);

        loadInitFlag = 0;
    }

    if ( tStep->giveNumber() == 1 ) {
        int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();
    }

    //
    //    ->   BEGINNING OF LOAD (OR DISPLACEMENT) STEP  <-
    //

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    //
    // call numerical model to solve arise problem
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving       [step number %5d.%d, time = %e]\n\n", tStep->giveNumber(), tStep->giveVersion(), tStep->giveIntrinsicTime() );
#endif

    FloatArray extrapolatedForces;
    FloatArray *extrapolatedForcesPtr = &extrapolatedForces;
    if ( this->initialGuessType == IG_Tangent ) {
#ifdef VERBOSE
        OOFEM_LOG_RELEVANT("Computing initial guess\n");
#endif
        this->assembleExtrapolatedForces( extrapolatedForces, tStep, TangentStiffnessMatrix, this->giveDomain(di) );
        extrapolatedForces.negated();

        this->updateComponent( tStep, NonLinearLhs, this->giveDomain(di) );
        SparseLinearSystemNM *linSolver = nMethod->giveLinearSolver();
        OOFEM_LOG_RELEVANT("solving for increment\n");
        linSolver->solve(*stiffnessMatrix, extrapolatedForces, incrementOfDisplacement);
        OOFEM_LOG_RELEVANT("initial guess found\n");
        totalDisplacement.add(incrementOfDisplacement);
    } else if ( this->initialGuessType == IG_Original ) {
        incrementOfDisplacement.zero();
        this->assembleExtrapolatedForces( extrapolatedForces, tStep, ElasticStiffnessMatrix, this->giveDomain(di) );
        extrapolatedForces.negated();
        
    } else if ( this->initialGuessType != IG_None ) {
        OOFEM_ERROR("Initial guess type: %d not supported", initialGuessType);
    } else {
        incrementOfDisplacement.zero();
        extrapolatedForcesPtr = NULL;
    }

    //totalDisplacement.printYourself();
    if ( initialLoadVector.isNotEmpty() ) {
      numMetStatus = nMethod->solve(*stiffnessMatrix, incrementalLoadVector, & initialLoadVector, extrapolatedForcesPtr,
                                      totalDisplacement, incrementOfDisplacement, internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep);
    } else {
      numMetStatus = nMethod->solve(*stiffnessMatrix, incrementalLoadVector, NULL, extrapolatedForcesPtr,
                                      totalDisplacement, incrementOfDisplacement, internalForces,
                                      internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, tStep);
    }
    ///@todo Martin: ta bort!!!
    //this->updateComponent(tStep, NonLinearLhs, this->giveDomain(di));
    ///@todo Use temporary variables. updateYourself() should set the final values, while proceedStep should be callable multiple times for each step (if necessary). / Mikael
    OOFEM_LOG_RELEVANT("Equilibrium reached at load level = %f in %d iterations\n", cumulatedLoadLevel + loadLevel, currentIterations);
    prevStepLength =  currentStepLength;
}

void
NonLinearStatic :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
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
            this->assemble(*stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                           EModelDefaultEquationNumbering(), d);
        } else if ( ( stiffMode == nls_secantStiffness ) || ( stiffMode == nls_secantInitialStiffness && initFlag ) ) {
#ifdef VERBOSE
            OOFEM_LOG_DEBUG("Assembling secant stiffness matrix\n");
#endif
            stiffnessMatrix->zero(); // zero stiffness matrix
            this->assemble(*stiffnessMatrix, tStep, TangentAssembler(SecantStiffness),
                           EModelDefaultEquationNumbering(), d);
            initFlag = 0;
        } else if ( ( stiffMode == nls_elasticStiffness ) && ( initFlag ||
                                                              ( this->giveMetaStep( tStep->giveMetaStepNumber() )->giveFirstStepNumber() == tStep->giveNumber() ) || (updateElasticStiffnessFlag) ) ) {
#ifdef VERBOSE
            OOFEM_LOG_DEBUG("Assembling elastic stiffness matrix\n");
#endif
            stiffnessMatrix->zero(); // zero stiffness matrix
            this->assemble(*stiffnessMatrix, tStep, TangentAssembler(ElasticStiffness),
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
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Updating internal forces\n");
#endif
        // update internalForces and internalForcesEBENorm concurrently
        this->giveInternalForces(internalForces, true, d->giveNumber(), tStep);
        break;

    default:
        OOFEM_ERROR("Unknown Type of component.");
    }
}


void
NonLinearStatic :: printOutputAt(FILE *File, TimeStep *tStep)
{
    if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return;                                                                      // do not print even Solution step header
    }

    fprintf( File, "\n\nOutput for time %.3e, solution step number %d\n", tStep->giveTargetTime(), tStep->giveNumber() );
    fprintf(File, "Reached load level : %20.6f in %d iterations\n\n",
            cumulatedLoadLevel + loadLevel, currentIterations);

    nMethod->printState(File);

    this->giveDomain(1)->giveOutputManager()->doDofManOutput(File, tStep);
    this->giveDomain(1)->giveOutputManager()->doElementOutput(File, tStep);
}


contextIOResultType
NonLinearStatic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file = NULL;

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

    if ( ( iores = totalDisplacement.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = incrementOfDisplacement.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    int _cm = controlMode;
    if ( !stream->write(_cm) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(loadLevel) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(cumulatedLoadLevel) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // store InitialLoadVector
    if ( ( iores = initialLoadVector.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = initialLoadVectorOfPrescribed.storeYourself(*stream) ) != CIO_OK ) {
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
    FILE *file = NULL;

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

    if ( ( iores = totalDisplacement.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = incrementOfDisplacement.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    int _cm;
    if ( !stream->read(_cm) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    controlMode = ( NonLinearStatic_controlType ) _cm;
    if ( !stream->read(loadLevel) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(cumulatedLoadLevel) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    // store InitialLoadVector
    if ( ( iores = initialLoadVector.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = initialLoadVectorOfPrescribed.restoreYourself(*stream) ) != CIO_OK ) {
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


void
NonLinearStatic :: assemble(SparseMtrx &answer, TimeStep *tStep, const MatrixAssembler &ma,
                            const UnknownNumberingScheme &s, Domain *domain)
{
#ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
#endif

    LinearStatic :: assemble(answer, tStep, ma, s, domain);

    if ( ( nonlocalStiffnessFlag ) && dynamic_cast< const TangentAssembler* >(&ma) ) {
        // add nonlocal contribution
        for ( auto &elem : domain->giveElements() ) {
            static_cast< StructuralElement * >( elem.get() )->addNonlocalStiffnessContributions(answer, s, tStep);
        }

        // print storage statistics
        answer.printStatistics();
    }

#ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_DEBUG( "NonLinearStatic: User time consumed by assembly: %.2fs\n", timer.getUtime() );
#endif
}


#ifdef __OOFEG
void
NonLinearStatic :: showSparseMtrxStructure(int type, oofegGraphicContext &gc, TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    CharType ctype;

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

    for ( auto &elem : domain->giveElements() ) {
        elem->showSparseMtrxStructure(ctype, gc, tStep);
    }

    for ( auto &elem : domain->giveElements() ) {
        elem->showExtendedSparseMtrxStructure(ctype, gc, tStep);
    }
}
#endif


void
NonLinearStatic :: computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    if ( ( di == 1 ) && ( tStep == this->giveCurrentStep() ) ) {
        reactions = initialLoadVectorOfPrescribed;
        reactions.add(loadLevel, incrementalLoadVectorOfPrescribed);
    } else {
        OOFEM_ERROR("unable to respond due to invalid solution step or domain");
    }
}


void
NonLinearStatic :: assembleIncrementalReferenceLoadVectors(FloatArray &_incrementalLoadVector,
                                                           FloatArray &_incrementalLoadVectorOfPrescribed,
                                                           SparseNonLinearSystemNM :: referenceLoadInputModeType _refMode,
                                                           Domain *sourceDomain, TimeStep *tStep)
{
    _incrementalLoadVector.resize( sourceDomain->giveEngngModel()->giveNumberOfDomainEquations( sourceDomain->giveNumber(), EModelDefaultEquationNumbering() ) );
    _incrementalLoadVector.zero();
    _incrementalLoadVectorOfPrescribed.resize( sourceDomain->giveEngngModel()->giveNumberOfDomainEquations( sourceDomain->giveNumber(), EModelDefaultPrescribedEquationNumbering() ) );
    _incrementalLoadVectorOfPrescribed.zero();

    if ( _refMode == SparseNonLinearSystemNM :: rlm_incremental ) {
        ///@todo This was almost definitely wrong before. It never seems to be used. Is this code even relevant?
        this->assembleVector(_incrementalLoadVector, tStep, ExternalForceAssembler(),
                             VM_Incremental, EModelDefaultEquationNumbering(), sourceDomain);

        this->assembleVector(_incrementalLoadVectorOfPrescribed, tStep, ExternalForceAssembler(),
                             VM_Incremental, EModelDefaultPrescribedEquationNumbering(), sourceDomain);
    } else {
        this->assembleVector(_incrementalLoadVector, tStep, ExternalForceAssembler(),
                             VM_Total, EModelDefaultEquationNumbering(), sourceDomain);

        this->assembleVector(_incrementalLoadVectorOfPrescribed, tStep, ExternalForceAssembler(),
                             VM_Total, EModelDefaultPrescribedEquationNumbering(), sourceDomain);
    }

    this->updateSharedDofManagers(_incrementalLoadVector, EModelDefaultEquationNumbering(), LoadExchangeTag);
}


int
NonLinearStatic :: estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType)
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

        //printf ("\nestimated count is %d\n",count);
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
NonLinearStatic :: giveLoadBalancer()
{
    if ( lb ) {
        return lb;
    }

    if ( loadBalancingFlag ) {
        ///@todo Make the name possibly optional (but currently, there is just one choice, "parmetis")
        lb = classFactory.createLoadBalancer( "parmetis", this->giveDomain(1) );
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
        lbm = classFactory.createLoadBalancerMonitor( "wallclock", this);
        return lbm;
    } else {
        return NULL;
    }
}
#endif

void
NonLinearStatic :: packMigratingData(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    int ndofman = domain->giveNumberOfDofManagers();
    bool initialLoadVectorEmpty = initialLoadVector.isEmpty();
    bool initialLoadVectorOfPrescribedEmpty = initialLoadVectorOfPrescribed.isEmpty();

    for ( int idofman = 1; idofman <= ndofman; idofman++ ) {
        DofManager *_dm = domain->giveDofManager(idofman);
        for ( Dof *_dof: *_dm ) {
            if ( _dof->isPrimaryDof() ) {
                int _eq;
                if ( ( _eq = _dof->__giveEquationNumber() ) ) {
                    // pack values in solution vectors
                    _dof->updateUnknownsDictionary( tStep, VM_Total, totalDisplacement.at(_eq) );
                    if ( initialLoadVectorEmpty ) {
                        _dof->updateUnknownsDictionary(tStep, VM_RhsInitial, 0.0);
                    } else {
                        _dof->updateUnknownsDictionary( tStep, VM_RhsInitial, initialLoadVector.at(_eq) );
                    }

                    _dof->updateUnknownsDictionary( tStep, VM_RhsIncremental, incrementalLoadVector.at(_eq) );
                } else if ( ( _eq = _dof->__givePrescribedEquationNumber() ) ) {
                    // pack values in prescribed solution vectors
                    if ( initialLoadVectorOfPrescribedEmpty ) {
                        _dof->updateUnknownsDictionary(tStep, VM_RhsInitial, 0.0);
                    } else {
                        _dof->updateUnknownsDictionary( tStep, VM_RhsInitial, initialLoadVectorOfPrescribed.at(_eq) );
                    }

                    _dof->updateUnknownsDictionary( tStep, VM_RhsIncremental, incrementalLoadVectorOfPrescribed.at(_eq) );
                }
            } // end primary dof
        } // end dof loop
    } // end dofman loop
}


void
NonLinearStatic :: unpackMigratingData(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    int ndofman = domain->giveNumberOfDofManagers();
    //int myrank = this->giveRank();

    // resize target arrays
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
    totalDisplacement.resize(neq);
    incrementOfDisplacement.resize(neq);
    incrementalLoadVector.resize(neq);
    initialLoadVector.resize(neq);
    initialLoadVectorOfPrescribed.resize( giveNumberOfDomainEquations( 1, EModelDefaultPrescribedEquationNumbering() ) );
    incrementalLoadVectorOfPrescribed.resize( giveNumberOfDomainEquations( 1, EModelDefaultPrescribedEquationNumbering() ) );

    for ( int idofman = 1; idofman <= ndofman; idofman++ ) {
        DofManager *_dm = domain->giveDofManager(idofman);
        for ( Dof *_dof: *_dm ) {
            if ( _dof->isPrimaryDof() ) {
                int _eq;
                if ( ( _eq = _dof->__giveEquationNumber() ) ) {
                    // pack values in solution vectors
                    totalDisplacement.at(_eq) = _dof->giveUnknownsDictionaryValue( tStep, VM_Total );
                    initialLoadVector.at(_eq) = _dof->giveUnknownsDictionaryValue( tStep, VM_RhsInitial );
                    incrementalLoadVector.at(_eq) = _dof->giveUnknownsDictionaryValue( tStep, VM_RhsIncremental );

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
                    initialLoadVectorOfPrescribed.at(_eq) = _dof->giveUnknownsDictionaryValue( tStep, VM_RhsInitial );
                    incrementalLoadVectorOfPrescribed.at(_eq) = _dof->giveUnknownsDictionaryValue( tStep, VM_RhsIncremental );

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

} // end namespace oofem
