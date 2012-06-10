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

#include "supg.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "initial.h"
#include "maskedprimaryfield.h"
#include "verbose.h"
#include "supgelement.h"
#include "usrdefsub.h"
#include "mathfem.h"
#include "dofdistributedprimaryfield.h"
#include "leplic.h"
#include "levelsetpcs.h"
#include "datastream.h"
#include "loadtime.h"
#include "contextioerr.h"
#ifdef TIME_REPORT
 #include "clock.h"
#endif

namespace oofem {
/* define if implicit interface update required */
//#define SUPG_IMPLICIT_INTERFACE


NumericalMethod *SUPG :: giveNumericalMethod(MetaStep *mStep)
{
    if ( nMethod ) {
        return nMethod;
    }

    nMethod = CreateUsrDefSparseLinSolver(solverType, 1, this->giveDomain(1), this);
    if ( nMethod == NULL ) {
        _error("giveNumericalMethod: linear solver creation failed");
    }

    return nMethod;
}

IRResultType
SUPG :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    EngngModel :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, rtolv, IFT_SUPG_rtolv, "rtolv"); // Macro
    atolv = 1.e-15;
    IR_GIVE_OPTIONAL_FIELD(ir, atolv, IFT_SUPG_atolv, "atolv"); // Macro


    int __val = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, __val, IFT_SUPG_stopmaxiter, "stopmaxiter"); // Macro
    if ( __val ) {
        stopmaxiter = true;
    } else {
        stopmaxiter = false;
    }

    maxiter = 200;
    IR_GIVE_OPTIONAL_FIELD(ir, maxiter, IFT_SUPG_maxiter, "maxiter"); // Macro

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_SUPG_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_SUPG_smtype, "smtype"); // Macro
    sparseMtrxType = ( SparseMtrxType ) val;

    IR_GIVE_FIELD(ir, deltaT, IFT_SUPG_deltat, "deltat"); // Macro
    deltaTLTF = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaTLTF, IFT_SUPG_deltatltf, "deltatltf");

    IR_GIVE_OPTIONAL_FIELD(ir, consistentMassFlag, IFT_SUPG_cmflag, "cmflag");

    alpha = 0.5;
    IR_GIVE_OPTIONAL_FIELD(ir, alpha, IFT_SUPG_alpha, "alpha");

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_SUPG_scaleflag, "scaleflag");
    equationScalingFlag = val;
    if ( equationScalingFlag ) {
        IR_GIVE_FIELD(ir, lscale, IFT_SUPG_lscale, "lscale"); // Macro
        IR_GIVE_FIELD(ir, uscale, IFT_SUPG_uscale, "uscale"); // Macro
        IR_GIVE_FIELD(ir, dscale, IFT_SUPG_dscale, "dscale"); // Macro
        double vref = 1.0; // reference viscosity
        Re = dscale * uscale * lscale / vref;
    } else {
        lscale = uscale = dscale = 1.0;
        Re = 1.0;
    }

    if ( requiresUnknownsDictionaryUpdate() ) {
        VelocityPressureField = new DofDistributedPrimaryField(this, 1, FT_VelocityPressure, EID_MomentumBalance_ConservationEquation, 1);
    } else {
        VelocityPressureField = new PrimaryField(this, 1, FT_VelocityPressure, EID_MomentumBalance_ConservationEquation, 1);
    }

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_SUPG_miflag, "miflag");
    if ( val == 1 ) {
        this->materialInterface = new LEPlic( 1, this->giveDomain(1) );
        this->materialInterface->initializeFrom(ir);
        // export velocity field
        FieldManager *fm = this->giveContext()->giveFieldManager();
        IntArray mask(3);
        mask.at(1) = V_u; mask.at(2) = V_v; mask.at(3) = V_w;
        MaskedPrimaryField* _velocityField = new MaskedPrimaryField (FT_Velocity, this->VelocityPressureField, mask);
        fm->registerField(_velocityField, FT_Velocity, true);

        //fsflag = 0;
        //IR_GIVE_OPTIONAL_FIELD (ir, fsflag, IFT_SUPG_fsflag, "fsflag");
    } else if ( val == 2 ) {
        // positive coefficient scheme level set alg
        this->materialInterface = new LevelSetPCS( 1, this->giveDomain(1) );
        this->materialInterface->initializeFrom(ir);
    }

    return IRRT_OK;
}


double
SUPG :: giveUnknownComponent(EquationID chc, ValueModeType mode,
                             TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displaacement, velocity of equation eq
// This function translates this request to numerical method language
{

    if ( this->requiresUnknownsDictionaryUpdate() ) {
        int hash = this->giveUnknownDictHashIndx(chc, mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
        } else {

        int eq = dof->__giveEquationNumber();
        if ( eq == 0 ) {
            _error("giveUnknownComponent: invalid equation number");
        }

        if ( chc == EID_ConservationEquation ) { // pressures
            return VelocityPressureField->giveUnknownValue(dof, mode, tStep);
        } else if ( chc == EID_MomentumBalance ) { // velocities & accelerations
            if ( mode == VM_Acceleration ) {
                return accelerationVector.at(eq);
                /*
                * if (tStep->isTheCurrentTimeStep()) {
                * return accelerationVector.at(eq);
                * } else if (tStep == givePreviousStep()) {
                * return previousAccelerationVector.at(eq);
                * } else _error ("giveUnknownComponent: VM_Acceleration history past previous step not supported");
                */
            } else if ( mode == VM_Velocity ) {
                return VelocityPressureField->giveUnknownValue(dof, VM_Total, tStep);
            } else {
                return VelocityPressureField->giveUnknownValue(dof, mode, tStep);
            }
        } else {
            _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
            return 0.;
        }
    }
    return 0;
}


double
SUPG :: giveUnknownComponent(UnknownType chc, ValueModeType mode,
                             TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displaacement, velocity of equation eq
// This function translates this request to numerical method language
{
    if ( chc == ReynoldsNumber ) {
        if ( equationScalingFlag ) {
            return this->Re;
        } else {
            return 1.0;
        }
    } else {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
    }

    return 0;
}


TimeStep *
SUPG :: giveSolutionStepWhenIcApply()
{
    if ( stepWhenIcApply == NULL ) {
        double dt = deltaT / this->giveVariableScale(VST_Time);

        stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0,
                                       0.0, dt, 0);
    }

    return stepWhenIcApply;
}

TimeStep *
SUPG :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    int i, nelem;
    double totalTime = 0;
    double dt = deltaT;
    StateCounterType counter = 1;
    delete previousStep;

    Domain *domain = this->giveDomain(1);
    nelem = domain->giveNumberOfElements();

    if ( currentStep == NULL ) {
        // first step -> generate initial step
        currentStep = new TimeStep( *giveSolutionStepWhenIcApply() );
    } else {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;

    // FORCE EQUATION NUMBERING
    this->giveNumberOfEquations(EID_MomentumBalance_ConservationEquation);

    if ( deltaTLTF ) {
        dt *= domain->giveLoadTimeFunction(deltaTLTF)->__at(istep);
    }

    // check for critical time step
    for ( i = 1; i <= nelem; i++ ) {
        dt = min( dt, ( ( SUPGElement * ) domain->giveElement(i) )->computeCriticalTimeStep(previousStep) );
    }

    if ( materialInterface ) {
        dt = min( dt, materialInterface->computeCriticalTimeStep(previousStep) );
    }

    // dt *= 0.6;
    dt /= this->giveVariableScale(VST_Time);

    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTargetTime() + dt;
    }

    currentStep = new TimeStep(istep, this, 1, totalTime, dt, counter);

    OOFEM_LOG_INFO( "SolutionStep %d : t = %e, dt = %e\n", istep, totalTime * this->giveVariableScale(VST_Time), dt * this->giveVariableScale(VST_Time) );
    // time and dt variables are set eq to 0 for staics - has no meaning
    return currentStep;
}


void
SUPG :: solveYourselfAt(TimeStep *tStep)
{
    int i, nite = 0;
    int neq =  this->giveNumberOfEquations(EID_MomentumBalance_ConservationEquation);
    double _absErrResid, err, avn = 0.0, aivn = 0.0, rnorm;
    FloatArray *solutionVector = NULL, *prevSolutionVector = NULL;
    FloatArray rhs(neq);
    int j,jj,ndofs, nnodes = this->giveDomain(1)->giveNumberOfDofManagers();
    double val,rnorm_mb, rnorm_mc;
    DofManager* inode;
    Dof* jDof;
    DofIDItem type;


    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        TimeStep *stepWhenIcApply = tStep->givePreviousStep();
        if ( materialInterface ) {
            materialInterface->initialize();
        }

        this->applyIC(stepWhenIcApply);
        //if (this->fsflag) this->updateDofManActivityMap(tStep);
    }

    if ( !requiresUnknownsDictionaryUpdate() ) {
        VelocityPressureField->advanceSolution(tStep);
        solutionVector = VelocityPressureField->giveSolutionVector(tStep);
        prevSolutionVector = VelocityPressureField->giveSolutionVector( tStep->givePreviousStep() );
    }

    if ( initFlag ) {
        if ( !requiresUnknownsDictionaryUpdate() ) {
            //previousAccelerationVector.resize(neq);
            accelerationVector.resize(neq);
            solutionVector->resize(neq);
            prevSolutionVector->resize(neq);
        }

        incrementalSolutionVector.resize(neq);

        lhs = CreateUsrDefSparseMtrx(sparseMtrxType);
        if ( lhs == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        lhs->buildInternalStructure( this, 1, EID_MomentumBalance_ConservationEquation, EModelDefaultEquationNumbering() );

        if ( materialInterface ) {
            this->updateElementsForNewInterfacePosition(tStep);
        }

        initFlag = 0;
    } else if ( requiresUnknownsDictionaryUpdate() ) {
        // rebuild lhs structure and resize solution vector
        incrementalSolutionVector.resize(neq);
        lhs->buildInternalStructure( this, 1, EID_MomentumBalance_ConservationEquation, EModelDefaultEquationNumbering() );
    }


    if ( !requiresUnknownsDictionaryUpdate() ) {
        for ( i = 1; i <= neq; i++ ) {
            solutionVector->at(i) = prevSolutionVector->at(i);
        }
    }

    //previousAccelerationVector=accelerationVector;

    // evaluate element supg and sppg stabilization coeffs
    this->evaluateElementStabilizationCoeffs( tStep);

    //
    // predictor
    //

    if ( requiresUnknownsDictionaryUpdate() ) {
        this->updateDofUnknownsDictionary_predictor(tStep);
    } else {
      this->updateSolutionVectors_predictor(*solutionVector, accelerationVector, tStep);
    }

    if ( tStep->giveNumber() != 1 ) {
      if ( materialInterface ) {
            //if (this->fsflag) updateDofManVals(tStep);
#ifdef SUPG_IMPLICIT_INTERFACE
 #ifdef TIME_REPORT
            oofem_timeval tstart;
            :: getUtime(tstart);
 #endif
            materialInterface->updatePosition( this->giveCurrentStep() );
            updateElementsForNewInterfacePosition(tStep);
 #ifdef TIME_REPORT
            oofem_timeval ut;
            :: getRelativeUtime(ut, tstart);
            OOFEM_LOG_INFO( "SUPG info: user time consumed by updating interfaces: %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
 #endif
#else
            //updateElementsForNewInterfacePosition (tStep);
#endif
            //if (this->fsflag) this->updateDofManActivityMap(tStep);
        }
    }

    // evaluate element supg and sppg stabilization coeffs
    // this -> evaluateElementStabilizationCoeffs (tStep);

    //
    // assemble rhs (residual)
    //
    rhs.zero();
    this->assembleVectorFromElements( rhs, tStep, EID_MomentumBalance, BCRhsTerm_MB, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, BCRhsTerm_MC, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    // algoritmic rhs part (assembled by e-model (in giveCharComponent service) from various element contribs)
    this->assembleVectorFromElements( rhs, tStep, EID_MomentumBalance, AlgorithmicRhsTerm_MB, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, AlgorithmicRhsTerm_MC, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );

    //
    // corrector
    //
    OOFEM_LOG_INFO("Iteration  IncrSolErr      RelResidErr     AbsResidErr\n_______________________________________________________\n");
    do {
        nite++;
        //
        // Assemble lhs
        //
        //if (nite == 1)
        {
            // momentum balance part
            lhs->zero();
            this->assemble( lhs, tStep, EID_MomentumBalance, AccelerationTerm_MB,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( lhs, tStep, EID_MomentumBalance, AdvectionDerivativeTerm_MB,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            if ( 1 ) { //if ((nite > 5)) // && (rnorm < 1.e4))
                this->assemble( lhs, tStep, EID_MomentumBalance, TangentDiffusionDerivativeTerm_MB,
                               EModelDefaultEquationNumbering(), this->giveDomain(1) );
            } else {
                this->assemble( lhs, tStep, EID_MomentumBalance, SecantDiffusionDerivativeTerm_MB,
                               EModelDefaultEquationNumbering(), this->giveDomain(1) );
            }

            this->assemble( lhs, tStep, EID_MomentumBalance, EID_ConservationEquation, PressureTerm_MB,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( lhs, tStep, EID_MomentumBalance, LSICStabilizationTerm_MB,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( lhs, tStep, EID_MomentumBalance, BCLhsTerm_MB,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( lhs, tStep, EID_MomentumBalance, EID_ConservationEquation, BCLhsPressureTerm_MB,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            // conservation eq part
            this->assemble( lhs, tStep, EID_ConservationEquation, EID_MomentumBalance, LinearAdvectionTerm_MC,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( lhs, tStep, EID_ConservationEquation, EID_MomentumBalance, AdvectionDerivativeTerm_MC,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( lhs, tStep, EID_ConservationEquation, EID_MomentumBalance, AccelerationTerm_MC,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( lhs, tStep, EID_ConservationEquation, EID_MomentumBalance, DiffusionDerivativeTerm_MC,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( lhs, tStep, EID_ConservationEquation, EID_ConservationEquation, PressureTerm_MC,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }
        //if (this->fsflag) this->imposeAmbientPressureInOuterNodes(lhs,&rhs,tStep);

#if 1
        nMethod->solve(lhs, & rhs, & incrementalSolutionVector);
#else

        SparseMtrx *__lhs = lhs->GiveCopy();

        nMethod->solve(lhs, & rhs, & incrementalSolutionVector);

        // check solver
        FloatArray __rhs(neq);
        double __absErr = 0., __relErr = 0.;
        __lhs->times(incrementalSolutionVector, __rhs);
        for ( i = 1; i <= neq; i++ ) {
            if ( fabs( __rhs.at(i) - rhs.at(i) ) > __absErr ) {
                __absErr = fabs( __rhs.at(i) - rhs.at(i) );
            }

            if ( fabs( ( __rhs.at(i) - rhs.at(i) ) / rhs.at(i) ) > __relErr ) {
                __relErr = fabs( ( __rhs.at(i) - rhs.at(i) ) / rhs.at(i) );
            }
        }

        OOFEM_LOG_INFO("SUPG: solver check: absoluteError %e, relativeError %e\n", __absErr, __relErr);
        delete(__lhs);
#endif



        if ( requiresUnknownsDictionaryUpdate() ) {
            this->updateDofUnknownsDictionary_corrector(tStep);
        } else {

        //update
        this->updateSolutionVectors(*solutionVector, accelerationVector, incrementalSolutionVector, tStep);
        avn = accelerationVector.computeSquaredNorm();
        aivn = incrementalSolutionVector.computeSquaredNorm();
    }   // end update

#if 0
 #ifdef SUPG_IMPLICIT_INTERFACE
        if ( materialInterface ) {
  #ifdef TIME_REPORT
            oofem_timeval tstart;
            :: getUtime(tstart);
  #endif
            //if (this->fsflag) updateDofManVals(tStep);
            materialInterface->updatePosition( this->giveCurrentStep() );
            updateElementsForNewInterfacePosition(tStep);
  #ifdef TIME_REPORT
            oofem_timeval ut;
            :: getRelativeUtime(ut, tstart);
            OOFEM_LOG_INFO( "SUPG info: user time consumed by updating interfaces: %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
  #endif
            //if (this->fsflag) this->updateDofManActivityMap(tStep);
        }

 #endif
#else
        if ( materialInterface ) {
            //if (this->fsflag) updateDofManVals(tStep);
            //if (this->fsflag) this->updateDofManActivityMap(tStep);
        }

#endif

        // check convergence and repeat iteration if desired
        if ( avn < 1.e-12 ) {
            err = aivn;
        } else {
            err = sqrt(aivn / avn);
        }

        //
        // assemble rhs (residual)
        //
        rhs.zero();
        this->assembleVectorFromElements( rhs, tStep, EID_MomentumBalance, BCRhsTerm_MB, VM_Total,
                                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
        this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, BCRhsTerm_MC, VM_Total,
                                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
        // algoritmic rhs part (assembled by e-model (in giveCharComponent service) from various element contribs)
        this->assembleVectorFromElements( rhs, tStep, EID_MomentumBalance, AlgorithmicRhsTerm_MB, VM_Total,
                                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
        this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, AlgorithmicRhsTerm_MC, VM_Total,
                                         EModelDefaultEquationNumbering(), this->giveDomain(1) );

        // check convergence and repeat iteration if desired
        rnorm_mb = rnorm_mc = 0.0;
        for ( i = 1; i <= nnodes; i++ ) {
            inode = this->giveDomain(1)->giveDofManager(i);
            ndofs = inode->giveNumberOfDofs();
            for ( j = 1; j <= ndofs; j++ ) {
                jDof  =  inode->giveDof(j);
                type  =  jDof->giveDofID();
                if ((jj = jDof->__giveEquationNumber())) {
                    val = rhs.at(jj);
                    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                        rnorm_mb+=val*val;
                    } else {
                        rnorm_mc+=val*val;
                    }
                }
            }
        }
        rnorm_mb=sqrt(rnorm_mb);
        rnorm_mc=sqrt(rnorm_mc);

        rnorm = rhs.computeNorm();

        _absErrResid = 0.0;
        for ( i = 1; i <= neq; i++ ) {
            _absErrResid = max( _absErrResid, fabs( rhs.at(i) ) );
        }

        if ( requiresUnknownsDictionaryUpdate() ) {
            OOFEM_LOG_INFO("%-10d       n/a       %-15e %-15e\n", nite, rnorm, _absErrResid);
        } else {
            OOFEM_LOG_INFO("%-10d %-15e %-15e (%-15e, %-15e) %-15e\n", nite, err, rnorm, rnorm_mb, rnorm_mc, _absErrResid);
        }

        if ( 0 ) {
            // evaluate element supg and sppg stabilization coeffs
            this->evaluateElementStabilizationCoeffs(tStep);
            //
            // assemble rhs (residual)
            //
            rhs.zero();
            this->assembleVectorFromElements( rhs, tStep, EID_MomentumBalance, BCRhsTerm_MB, VM_Total,
                                             EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, BCRhsTerm_MC, VM_Total,
                                             EModelDefaultEquationNumbering(), this->giveDomain(1) );
            // algoritmic rhs part (assembled by e-model (in giveCharComponent service) from various element contribs)
            this->assembleVectorFromElements( rhs, tStep, EID_MomentumBalance, AlgorithmicRhsTerm_MB, VM_Total,
                            EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, AlgorithmicRhsTerm_MC, VM_Total,
                            EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }

    } while ( ( rnorm > rtolv ) && ( _absErrResid > atolv ) && ( nite <= maxiter ) );

    if ( nite <= maxiter ) {
        OOFEM_LOG_INFO("SUPG info: number of iterations: %d\n", nite);
    } else {
        OOFEM_WARNING2("SUPG info: Convergence not reached, number of iterations: %d\n", nite);
        if ( stopmaxiter ) {
            exit(1);
        }
    }

#ifndef SUPG_IMPLICIT_INTERFACE
    if ( materialInterface ) {
 #ifdef TIME_REPORT
        oofem_timeval tstart;
        getUtime(tstart);
 #endif
        materialInterface->updatePosition( this->giveCurrentStep() );
        updateElementsForNewInterfacePosition(tStep);
 #ifdef TIME_REPORT
        oofem_timeval ut;
        getRelativeUtime(ut, tstart);
        OOFEM_LOG_INFO( "SUPG info: user time consumed by updating interfaces: %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
 #endif
        //if (this->fsflag) this->updateDofManActivityMap(tStep);
    }

#endif

    // update solution state counter
    tStep->incrementStateCounter();
}


void
SUPG :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    EngngModel :: updateYourself(stepN);
    if ( materialInterface ) {
        materialInterface->updateYourself(stepN);
    }

    //previousSolutionVector = solutionVector;
}



void
SUPG :: updateInternalState(TimeStep *stepN)
{
    int j, nnodes;
    Domain *domain;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);

        nnodes = domain->giveNumberOfDofManagers();
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( j = 1; j <= nnodes; j++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(j), stepN);
            }
        }

        int nelem = domain->giveNumberOfElements();
        for ( j = 1; j <= nelem; j++ ) {
            domain->giveElement(j)->updateInternalState(stepN);
        }
    }
}


int
SUPG :: forceEquationNumbering(int id)
{
    // Necessary to number DOFs in special order to guarantee that Skyline matrix factorization to work.

    int i, j, k, innodes, nnodes, nelem, nbc, ndofs;
    Element *elem;
    DofManager *dman;
    Dof *jDof;
    DofIDItem type;
    Domain *domain = this->giveDomain(id);
    TimeStep *currStep = this->giveCurrentStep();
    IntArray loc;

    this->domainNeqs.at(id) = 0;
    this->domainPrescribedNeqs.at(id) = 0;

    nnodes = domain->giveNumberOfDofManagers();
    nelem  = domain->giveNumberOfElements();
    nbc    = domain->giveNumberOfBoundaryConditions();

    // First velocity.
    for ( i = 1; i <= nnodes; i++ ) {
        dman = domain->giveDofManager(i);
        ndofs = dman->giveNumberOfDofs();
        for ( j = 1; j <= ndofs; j++ ) {
            jDof = dman->giveDof(j);
            type = jDof->giveDofID();
            if ( (type == V_u) || (type == V_v) || (type == V_w) ) {
                jDof->askNewEquationNumber(currStep);
            }
        }
    }
    for ( i = 1; i <= nelem; ++i ) {
        elem = domain->giveElement(i);
        innodes = elem->giveNumberOfInternalDofManagers();
        for ( k = 1; k <= innodes; k++ ) {
            dman = elem->giveInternalDofManager(k);
            ndofs = dman->giveNumberOfDofs();
            for ( j = 1; j <= ndofs; j++ ) {
                jDof = dman->giveDof(j);
                type = jDof->giveDofID();
                if ( (type == V_u) || (type == V_v) || (type == V_w) ) {
                    jDof->askNewEquationNumber(currStep);
                }
            }
        }
    }
    for ( i = 1; i <= nbc; ++i ) {
        GeneralBoundaryCondition *bc = domain->giveBc(i);
        innodes = bc->giveNumberOfInternalDofManagers();
        for ( k = 1; k <= innodes; k++ ) {
            dman = bc->giveInternalDofManager(k);
            ndofs = dman->giveNumberOfDofs();
            for ( j = 1; j <= ndofs; j++ ) {
                jDof = dman->giveDof(j);
                type = jDof->giveDofID();
                if ( (type == V_u) || (type == V_v) || (type == V_w) ) {
                    jDof->askNewEquationNumber(currStep);
                }
            }
        }
    }
    // Then the rest
    for ( i = 1; i <= nnodes; i++ ) {
        dman = domain->giveDofManager(i);
        ndofs = dman->giveNumberOfDofs();
        for ( j = 1; j <= ndofs; j++ ) {
            jDof = dman->giveDof(j);
            type = jDof->giveDofID();
            if ( !( (type == V_u) || (type == V_v) || (type == V_w) ) ) {
                jDof->askNewEquationNumber(currStep);
            }
        }
    }
    for ( i = 1; i <= nelem; ++i ) {
        elem = domain->giveElement(i);
        innodes = elem->giveNumberOfInternalDofManagers();
        for ( k = 1; k <= innodes; k++ ) {
            dman = elem->giveInternalDofManager(k);
            ndofs = dman->giveNumberOfDofs();
            for ( j = 1; j <= ndofs; j++ ) {
                jDof = dman->giveDof(j);
                type = jDof->giveDofID();
                if ( !( type == V_u || type == V_v || type == V_w ) ) {
                    jDof->askNewEquationNumber(currStep);
                }
            }
        }
    }
    for ( i = 1; i <= nbc; ++i ) {
        GeneralBoundaryCondition *bc = domain->giveBc(i);
        innodes = bc->giveNumberOfInternalDofManagers();
        for ( k = 1; k <= innodes; k++ ) {
            dman = bc->giveInternalDofManager(k);
            ndofs = dman->giveNumberOfDofs();
            for ( j = 1; j <= ndofs; j++ ) {
                jDof = dman->giveDof(j);
                type = jDof->giveDofID();
                if ( !( type == V_u || type == V_v || type == V_w ) ) {
                    jDof->askNewEquationNumber(currStep);
                }
            }
        }
    }

    return domainNeqs.at(id);
}


contextIOResultType
SUPG :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
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

    if ( ( iores = VelocityPressureField->saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( materialInterface ) {
        if ( ( iores = materialInterface->saveContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                        // ensure consistent records

    return CIO_OK;
}



contextIOResultType
SUPG :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    int istep, iversion;
    FILE *file;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = VelocityPressureField->restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( materialInterface ) {
        if ( ( iores = materialInterface->restoreContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                        // ensure consistent records

    return CIO_OK;
}


int
SUPG :: checkConsistency()
{
    // check internal consistency
    // if success returns nonzero
    int i, nelem;
    Element *ePtr;
    SUPGElement *sePtr;
    GeneralBoundaryCondition *bcPtr;
    InitialCondition *icPtr;
    Domain *domain = this->giveDomain(1);

    nelem = domain->giveNumberOfElements();
    // check for proper element type

    for ( i = 1; i <= nelem; i++ ) {
        ePtr = domain->giveElement(i);
        sePtr = dynamic_cast< SUPGElement * >(ePtr);
        if ( sePtr == NULL ) {
            _warning2("Element %d has no SUPG base", i);
            return 0;
        }
    }

    EngngModel :: checkConsistency();


    // scale boundary and initial conditions
    if ( equationScalingFlag ) {
        int nbc = domain->giveNumberOfBoundaryConditions();
        for ( i = 1; i <= nbc; i++ ) {
            bcPtr = domain->giveBc(i);
            if ( bcPtr->giveBCValType() == VelocityBVT ) {
                bcPtr->scale(1. / uscale);
            } else if ( bcPtr->giveBCValType() == PressureBVT ) {
                bcPtr->scale( 1. / this->giveVariableScale(VST_Pressure) );
            } else if ( bcPtr->giveBCValType() == ForceLoadBVT ) {
                bcPtr->scale( 1. / this->giveVariableScale(VST_Force) );
            } else {
                _error("checkConsistency: unknown bc/ic type\n");
            }
        }

        int nic = domain->giveNumberOfInitialConditions();
        for ( i = 1; i <= nic; i++ ) {
            icPtr = domain->giveIc(i);
            if ( icPtr->giveICValType() == VelocityBVT ) {
                icPtr->scale(VM_Total, 1. / uscale);
            } else if ( icPtr->giveICValType() == PressureBVT ) {
                icPtr->scale( VM_Total, 1. / this->giveVariableScale(VST_Pressure) );
            } else {
                _error("checkConsistency: unknown bc/ic type\n");
            }
        }
    }

    return 1;
}


void
SUPG :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}

void
SUPG :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    double pscale = ( dscale * uscale * uscale );

    DofIDItem type = iDof->giveDofID();
    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
        iDof->printSingleOutputAt(stream, atTime, 'v', EID_MomentumBalance, VM_Total, uscale);
    } else if ( type == P_f ) {
        iDof->printSingleOutputAt(stream, atTime, 'p', EID_ConservationEquation, VM_Total, pscale);
    } else {
        _error("printDofOutputAt: unsupported dof type");
    }
}

void
SUPG :: applyIC(TimeStep *stepWhenIcApply)
{
    Domain *domain = this->giveDomain(1);
    int neq =  this->giveNumberOfEquations(EID_MomentumBalance_ConservationEquation);
    FloatArray *vp_vector = NULL;

#ifdef VERBOSE
    OOFEM_LOG_INFO("Applying initial conditions\n");
#endif
    int nDofs, j, k, jj;
    int nman  = domain->giveNumberOfDofManagers();
    DofManager *node;
    Dof *iDof;
    DofIDItem type;

    if ( !requiresUnknownsDictionaryUpdate() ) {
        VelocityPressureField->advanceSolution(stepWhenIcApply);
        vp_vector = VelocityPressureField->giveSolutionVector(stepWhenIcApply);
        vp_vector->resize(neq);
        vp_vector->zero();
    }

    accelerationVector.resize(neq);
    accelerationVector.zero();

    if ( !requiresUnknownsDictionaryUpdate() ) {
        for ( j = 1; j <= nman; j++ ) {
            node = domain->giveDofManager(j);
            nDofs = node->giveNumberOfDofs();

            for ( k = 1; k <= nDofs; k++ ) {
                // ask for initial values obtained from
                // bc (boundary conditions) and ic (initial conditions)
                iDof  =  node->giveDof(k);
                if ( !iDof->isPrimaryDof() ) {
                    continue;
                }

                jj = iDof->__giveEquationNumber();
                type = iDof->giveDofID();
                if ( jj ) {
                    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                        vp_vector->at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Total, stepWhenIcApply);
                    } else {
                        vp_vector->at(jj) = iDof->giveUnknown(EID_ConservationEquation, VM_Total, stepWhenIcApply);
                    }
                }
            }
        }
    }

    //ICs and BCs are scaled already in CheckConsistency
    //vp_vector->times(1./this->giveVariableScale(VST_Velocity));

    // update element state according to given ic
    int nelem = domain->giveNumberOfElements();
    SUPGElement *element;

    //this->initElementsForNewStep (stepWhenIcApply);
    if ( materialInterface ) {
        this->updateElementsForNewInterfacePosition(stepWhenIcApply);
    }

    for ( j = 1; j <= nelem; j++ ) {
        element = ( SUPGElement * ) domain->giveElement(j);
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
    }
}

int
SUPG :: giveNewEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return ++domainNeqs.at(domain);
    } else if ( id == P_f ) {
        return ++domainNeqs.at(domain);
    } else {
        _error("giveNewEquationNumber:: Unknown DofIDItem");
    }

    return 0;
}

int
SUPG :: giveNewPrescribedEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return ++domainPrescribedNeqs.at(domain);
    } else if ( id == P_f ) {
        return ++domainPrescribedNeqs.at(domain);
    } else {
        _error("giveNewPrescribedEquationNumber:: Unknown DofIDItem");
    }

    return 0;
}

int
SUPG :: giveNumberOfEquations(EquationID id) {
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    if ( ( id == EID_MomentumBalance ) || ( id == EID_ConservationEquation ) || ( id == EID_MomentumBalance_ConservationEquation ) ) {
        return numberOfEquations;
    } else {
        _error("giveNumberOfEquations: unknown equation id");
    }

    return 0;
}

int
SUPG :: giveNumberOfPrescribedEquations(EquationID id) {
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    if ( ( id == EID_MomentumBalance ) || ( id == EID_ConservationEquation ) || ( id == EID_MomentumBalance_ConservationEquation ) ) {
        return numberOfPrescribedEquations;
    } else {
        _error("giveNumberOfPrescribedEquations: unknown equation id");
    }

    return 0;
}

int
SUPG :: giveNumberOfDomainEquations(int d, EquationID id) {
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    if ( ( id == EID_MomentumBalance ) || ( id == EID_ConservationEquation ) || ( id == EID_MomentumBalance_ConservationEquation ) ) {
        return numberOfEquations;
    } else {
        _error("giveNumberOfDomainEquations: unknown equation id");
    }

    return 0;
}

int
SUPG :: giveNumberOfPrescribedDomainEquations(int d, EquationID id) {
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    if ( ( id == EID_MomentumBalance ) || ( id == EID_ConservationEquation ) || ( id == EID_MomentumBalance_ConservationEquation ) ) {
        return numberOfPrescribedEquations;
    } else {
        _error("giveNumberOfPrescribedDomainEquations: unknown equation id");
    }

    return 0;
}

double
SUPG :: giveVariableScale(VarScaleType varID)
{
    if ( varID == VST_Length ) {
        return this->lscale;
    } else if ( varID == VST_Velocity ) {
        return this->uscale;
    } else if ( varID == VST_Density ) {
        return this->dscale;
    } else if ( varID == VST_Time ) {
        return ( lscale / uscale );
    } else if ( varID == VST_Pressure ) {
        return ( dscale * uscale * uscale );
    } else if ( varID == VST_Force ) {
        return ( uscale * uscale / lscale );
    } else if ( varID == VST_Viscosity ) {
        return 1.0;
    } else {
        _error("giveVariableScale: unknown variable type");
    }

    return 0.0;
}

void
SUPG :: evaluateElementStabilizationCoeffs(TimeStep *atTime)
{
    Domain *domain = this->giveDomain(1);
    int i, nelem = domain->giveNumberOfElements();
    SUPGElement *ePtr;

    //printf ("#");
    for ( i = 1; i <= nelem; i++ ) {
        ePtr = ( SUPGElement * ) domain->giveElement(i);
        ePtr->updateStabilizationCoeffs(atTime);
    }
}

void
SUPG :: updateElementsForNewInterfacePosition(TimeStep *atTime)
{
    Domain *domain = this->giveDomain(1);
    int i, nelem = domain->giveNumberOfElements();
    SUPGElement *ePtr;

    printf("updating elements for interface position\n");


    for ( i = 1; i <= nelem; i++ ) {
        ePtr = ( SUPGElement * ) domain->giveElement(i);
        ePtr->updateElementForNewInterfacePosition(atTime);
    }
}


void
SUPG :: giveElementCharacteristicMatrix(FloatMatrix &answer, int num, CharType type, TimeStep *tStep, Domain *domain)
{
    EngngModel :: giveElementCharacteristicMatrix(answer, num, type, tStep, domain);
    if ( ( type == AdvectionDerivativeTerm_MB ) ||
        ( type == DiffusionDerivativeTerm_MB ) ||
        ( type == TangentDiffusionDerivativeTerm_MB ) ||
        ( type == SecantDiffusionDerivativeTerm_MB ) ||
        ( type == InitialDiffusionDerivativeTerm_MB ) ||
        (type == BCLhsTerm_MB )) {
        answer.times( alpha * tStep->giveTimeIncrement() );
    } else if ( type == LSICStabilizationTerm_MB ) {
        double coeff = lscale / ( dscale * uscale * uscale );
        answer.times(alpha * tStep->giveTimeIncrement() * coeff);
    } else if ( ( type == AdvectionDerivativeTerm_MC ) ||
               ( type == DiffusionDerivativeTerm_MC ) ) {
        answer.times( alpha * tStep->giveTimeIncrement() );
    } else if ( type == LinearAdvectionTerm_MC ) {
        //double coeff = uscale/lscale;
        double coeff = 1 / ( dscale * uscale );
        answer.times(alpha * tStep->giveTimeIncrement() * coeff);
    }
}


void
SUPG :: giveElementCharacteristicVector(FloatArray &answer, int num, CharType type, ValueModeType mode, TimeStep *tStep, Domain *domain)
{
    if ( type == AlgorithmicRhsTerm_MB ) {
        Element *eptr = domain->giveElement(num);
        FloatMatrix m1;
        FloatArray h, v, vp;
        answer.resize(6);
        // add (M+M_delta)*a
        eptr->giveCharacteristicMatrix(m1, AccelerationTerm_MB, tStep);
        eptr->computeVectorOf(EID_MomentumBalance, VM_Acceleration, tStep, v);
        answer.beProductOf(m1, v);
        // add advection terms N+N_delta
        eptr->giveCharacteristicVector(v, AdvectionTerm_MB, VM_Total, tStep);
        answer.add(v);
        // add diffusion terms (K+K_delta)v
        eptr->giveCharacteristicVector(v, DiffusionTerm_MB, VM_Total, tStep);
        answer.add(v);
        // add lsic stabilization term
        eptr->giveCharacteristicMatrix(m1, LSICStabilizationTerm_MB, tStep);
        m1.times( lscale / ( dscale * uscale * uscale ) );
        eptr->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.add(h);
        eptr->giveCharacteristicMatrix(m1, BCLhsTerm_MB, tStep);
        h.beProductOf(m1, v);
        answer.add(h);
        // add pressure term
        eptr->giveCharacteristicMatrix(m1, PressureTerm_MB, tStep);
        eptr->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        //eptr->computeVectorOfPrescribed (EID_ConservationEquation, VM_Total, tStep, vp);
        //h.beProductOf(m1,vp); // term due to prescribed pressure
        answer.add(h);
        eptr->giveCharacteristicMatrix(m1, BCLhsPressureTerm_MB, tStep);
        h.beProductOf(m1, v);
        answer.add(h);
        answer.negated();
    } else if ( type == AlgorithmicRhsTerm_MC ) {
        Element *eptr = domain->giveElement(num);
        FloatMatrix m1, m2;
        FloatArray h, v, vp;
        answer.resize(3);
        // G^T term - linear advection term
        eptr->giveCharacteristicMatrix(m1, LinearAdvectionTerm_MC, tStep);
        //m1.times(uscale/lscale);
        m1.times( 1. / ( dscale * uscale ) );
        eptr->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, v);
        //eptr->computeVectorOfPrescribed (EID_MomentumBalance, VM_Velocity, tStep, vp);
        answer.beProductOf(m1, v);
        //h.beProductOf(m1,vp); // term due to prescribed velocity
        //answer.add(h);
        // Diffusion term
        eptr->giveCharacteristicVector(v, DiffusionTerm_MC, VM_Total, tStep);
        answer.add(v);
        //h.beProductOf(m1,vp);// term due to prescribed velocity
        //answer.add(h);
        // AccelerationTerm
        eptr->giveCharacteristicMatrix(m1, AccelerationTerm_MC, tStep);
        eptr->computeVectorOf(EID_MomentumBalance, VM_Acceleration, tStep, v);
        h.beProductOf(m1, v);
        answer.add(h);
        // advection N term (nonlinear)
        eptr->giveCharacteristicVector(v, AdvectionTerm_MC, VM_Total, tStep);
        answer.add(v);
        // pressure term
        eptr->giveCharacteristicMatrix(m1, PressureTerm_MC, tStep);
        eptr->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, v);
        //eptr->computeVectorOfPrescribed (EID_ConservationEquation, VM_Total, tStep, vp);
        h.beProductOf(m1, v);
        answer.add(h);
        //h.beProductOf(m1,vp); // term due to prescribed pressure
        //answer.add(h);
        answer.negated();
    } else {
        EngngModel :: giveElementCharacteristicVector(answer, num, type, mode, tStep, domain);
    }
}




void
SUPG :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    // update DOF unknowns dictionary, where
    // unknowns are hold instead of keeping them in global unknowns
    // vectors in engng instancies
    // this is necessary, because during solution equation numbers for
    // particular DOFs may changed, and it is necesary to keep them
    // in DOF level.

    // empty -> all done in updateDofUnknownsDictionary_corrector (TimeStep* tStep);
}


void
SUPG :: updateDofUnknownsDictionary_predictor(TimeStep *tStep)
{
    int i, j, ndofs;
    Dof *iDof;
    DofManager *inode;
    DofIDItem type;
    double val, prev_val, accel;
    double deltaT = tStep->giveTimeIncrement();
    Domain *domain = this->giveDomain(1);

    int nnodes = domain->giveNumberOfDofManagers();

    if ( requiresUnknownsDictionaryUpdate() ) {
        for ( j = 1; j <= nnodes; j++ ) {
            inode = domain->giveDofManager(j);
            ndofs = inode->giveNumberOfDofs();
            for ( i = 1; i <= ndofs; i++ ) {
                iDof = inode->giveDof(i);
                type = iDof->giveDofID();
                if ( !iDof->hasBc(tStep) ) {
                    prev_val = iDof->giveUnknown( EID_MomentumBalance_ConservationEquation, VM_Total, tStep->givePreviousStep() );
                    accel   = iDof->giveUnknown( EID_MomentumBalance_ConservationEquation, VM_Acceleration, tStep->givePreviousStep() );
                    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                        val = prev_val +  deltaT * accel;
                    } else {
                        val = prev_val;
                    }

                    iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance_ConservationEquation, VM_Total, val); // velocity
                    iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance_ConservationEquation, VM_Acceleration, accel); // acceleration
                } else {
                    val = iDof->giveBcValue(VM_Total, tStep);
                    iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance_ConservationEquation, VM_Total, val); // velocity
                    //val = iDof -> giveBcValue (VM_Velocity,tStep) ; //velocity of velocity is acceleration
                    iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance_ConservationEquation, VM_Acceleration, 0.0); // acceleration
                }
            }
        }
    }
}


void
SUPG :: updateDofUnknownsDictionary_corrector(TimeStep *tStep)
{
    int i, j, ndofs;
    Dof *iDof;
    DofManager *inode;
    DofIDItem type;
    double val, prev_val;
    double deltaT = tStep->giveTimeIncrement();
    Domain *domain = this->giveDomain(1);
    int nnodes = domain->giveNumberOfDofManagers();


    if ( requiresUnknownsDictionaryUpdate() ) {
        for ( j = 1; j <= nnodes; j++ ) {
            inode = domain->giveDofManager(j);
            ndofs = inode->giveNumberOfDofs();
            for ( i = 1; i <= ndofs; i++ ) {
                iDof = inode->giveDof(i);
                type = iDof->giveDofID();
                if ( !iDof->hasBc(tStep) ) {
                    prev_val = iDof->giveUnknown(EID_MomentumBalance_ConservationEquation, VM_Total, tStep);
                    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                        val = prev_val +  deltaT *alpha *incrementalSolutionVector.at( iDof->__giveEquationNumber() );
                    } else {
                        val = prev_val +  incrementalSolutionVector.at( iDof->__giveEquationNumber() );
                    }

                    iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance_ConservationEquation, VM_Total, val); // velocity

                    prev_val = iDof->giveUnknown(EID_MomentumBalance_ConservationEquation, VM_Acceleration, tStep);
                    val = prev_val +  incrementalSolutionVector.at( iDof->__giveEquationNumber() );
                    iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance_ConservationEquation, VM_Acceleration, val); // acceleration
                }
            }
        }
    }
}

int
SUPG :: giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN)
{
    if ( ( stepN == this->giveCurrentStep() ) || ( stepN == this->givePreviousStep() ) ) {
        return ( stepN->giveNumber() % 2 ) * 100 + mode;
    } else {
        _error("giveUnknownDictHashIndx: unsupported solution step");
    }

    return 0;
}

void
SUPG:: updateSolutionVectors_predictor(FloatArray& solutionVector, FloatArray& accelerationVector, TimeStep* tStep)
{
    int j, k, jj,ji, nDofs, ndofman;
    double deltaT = tStep->giveTimeIncrement();
    DofManager *node, *dofman;
    Dof *iDof;
    DofIDItem type;
    Domain *domain = this->giveDomain(1);
    int nman =  this->giveDomain(1)->giveNumberOfDofManagers();
    Element *elem;


    for ( j = 1; j <= nman; j++ ) {
        node = domain->giveDofManager(j);
        nDofs = node->giveNumberOfDofs();

        for ( k = 1; k <= nDofs; k++ ) {
            iDof  =  node->giveDof(k);
            if ( !iDof->isPrimaryDof() ) {
                continue;
            }

            jj = iDof->__giveEquationNumber();
            type = iDof->giveDofID();

            if ( jj ) {
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) { // v = v + (deltaT)*a
                    solutionVector.at(jj) += deltaT * accelerationVector.at(jj);
                }
            }
        }
    }

    for (j=1; j<= domain->giveNumberOfElements(); j++) {
        elem = domain->giveElement(j);
        ndofman = elem->giveNumberOfInternalDofManagers();
        for (ji=1; ji<=ndofman; ji++) {
            dofman = elem->giveInternalDofManager(ji);
            nDofs = dofman->giveNumberOfDofs();
            for ( k = 1; k <= nDofs; k++ ) {
                iDof  =  dofman->giveDof(k);
                if ( !iDof->isPrimaryDof() ) {
                    continue;
                }

                jj = iDof->__giveEquationNumber();
                type = iDof->giveDofID();

                if ( jj ) {
                    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) { // v = v + (deltaT*alpha)*da
                        solutionVector.at(jj) += deltaT * accelerationVector.at(jj);;
                    }
                }
            }
        } // end loop over elem internal dofmans
    } // end loop over elems
}


void
SUPG:: updateSolutionVectors(FloatArray& solutionVector, FloatArray& accelerationVector, FloatArray& incrementalSolutionVector, TimeStep* tStep)
{
    int j, k, jj,ji, nDofs, ndofman;
    double deltaT = tStep->giveTimeIncrement();
    DofManager *node, *dofman;
    Dof *iDof;
    DofIDItem type;
    Domain *domain = this->giveDomain(1);
    int nman =  this->giveDomain(1)->giveNumberOfDofManagers();
    Element *elem;

    accelerationVector.add(incrementalSolutionVector);

    for ( j = 1; j <= nman; j++ ) {
        node = domain->giveDofManager(j);
        nDofs = node->giveNumberOfDofs();

        for ( k = 1; k <= nDofs; k++ ) {
            iDof  =  node->giveDof(k);
            if ( !iDof->isPrimaryDof() ) {
                continue;
            }

            jj = iDof->__giveEquationNumber();
            type = iDof->giveDofID();

            if ( jj ) {
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) { // v = v + (deltaT*alpha)*da
                    solutionVector.at(jj) += deltaT * alpha * incrementalSolutionVector.at(jj);
                } else {
                    solutionVector.at(jj) += incrementalSolutionVector.at(jj); // p = p + dp
                }
            }
        }
    } // end loop over dnam


    for (j=1; j<= domain->giveNumberOfElements(); j++) {
        elem = domain->giveElement(j);
        ndofman = elem->giveNumberOfInternalDofManagers();
        for (ji=1; ji<=ndofman; ji++) {
            dofman = elem->giveInternalDofManager(ji);
            nDofs = dofman->giveNumberOfDofs();
            for ( k = 1; k <= nDofs; k++ ) {
                iDof  =  dofman->giveDof(k);
                if ( !iDof->isPrimaryDof() ) {
                    continue;
                }

                jj = iDof->__giveEquationNumber();
                type = iDof->giveDofID();

                if ( jj ) {
                    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) { // v = v + (deltaT*alpha)*da
                        solutionVector.at(jj) += deltaT * alpha * incrementalSolutionVector.at(jj);
                    } else {
                        solutionVector.at(jj) += incrementalSolutionVector.at(jj); // p = p + dp
                    }
                }
            }
        } // end loop over elem internal dofmans
    } // end loop over elems
}


#ifdef __PETSC_MODULE
void
SUPG :: initPetscContexts()
{
    PetscContext *petscContext;
    petscContextList->growTo(ndomains);
    for ( int i = 1; i <= this->ndomains; i++ ) {
        petscContext =  new PetscContext(this, EID_MomentumBalance_ConservationEquation);
        petscContextList->put(i, petscContext);
    }
}
#endif


#define __VOF_TRESHOLD 1.e-5

#define __NODE_OUT 0
#define __NODE_IN 1
#define __NODE_OUT_ACTIVE 2
#define __NODE_INTERPOL_CANDIDATE 3

/*
 * void
 * SUPG::updateDofManActivityMap (TimeStep* tStep)
 * {
 * // loop over inactive nodes and test if some of its shared element has become active (temp_vof>0)
 * // then initialize its velocity and pressure
 *
 * Domain* domain = this->giveDomain(1);
 * Element* ielem;
 * int nnodes = domain->giveNumberOfDofManagers();
 * int nelems = domain->giveNumberOfElements();
 * int i, inode, ie;
 * //int inodes, _code;
 * LEPlicElementInterface *interface;
 * IntArray dofMask(2);
 * FloatArray vv;
 *
 * dofMask.at(1) = V_u;
 * dofMask.at(2) = V_v;
 *
 * printf ("SUPG::updateDofManActivityMap called\n");
 *
 * IntArray __old_DofManActivityMask (__DofManActivityMask);
 * __DofManActivityMask.resize(nnodes);
 * for (i=1; i<=nnodes;i++) __DofManActivityMask.at(i) = __NODE_OUT;
 *
 *
 * for (ie=1; ie<= nelems; ie++) {
 *  ielem = domain->giveElement(ie);
 *  if ((interface = (LEPlicElementInterface*) (ielem->giveInterface(LEPlicElementInterfaceType)))) {
 *    if (interface->giveVolumeFraction() > 0.9999) {
 *      // mark all its nodes as active
 *      nnodes=ielem->giveNumberOfNodes();
 *      for (i=1; i<=nnodes;i++) this->__DofManActivityMask.at(ielem->giveNode(i)->giveNumber())=__NODE_IN;
 *    } else if (interface->giveVolumeFraction() > 0.0) {
 *      nnodes=ielem->giveNumberOfNodes();
 *      for (i=1; i<=nnodes;i++) {
 *        inode = ielem->giveNode(i)->giveNumber();
 *        //if (__DofManActivityMask.at(inode) == __NODE_OUT) {
 *          FloatArray normal; double p,x,y;
 *          interface->giveTempInterfaceNormal(normal);
 *          p = interface->giveTempLineConstant();
 *
 *          x=domain->giveNode(inode)->giveCoordinate(1);
 *          y=domain->giveNode(inode)->giveCoordinate(2);
 *
 *          if ((normal.at(1)*x+normal.at(2)*y+p) >= 0.) {
 *            __DofManActivityMask.at(inode) = __NODE_IN;
 *          } else {
 *            if (__DofManActivityMask.at(inode) == __NODE_OUT) {
 *              if (interface->giveVolumeFraction() < 0.5) __DofManActivityMask.at(inode) = __NODE_INTERPOL_CANDIDATE;
 *              else __DofManActivityMask.at(inode) = __NODE_OUT_ACTIVE;
 *            }
 *          }
 *          //}
 *      }
 *    }
 *  }
 * }
 *
 * }
 *
 *
 * void
 * SUPG:: updateDofManVals (TimeStep* tStep)
 * {
 * Domain* domain = this->giveDomain(1);
 * int nnodes = domain->giveNumberOfDofManagers();
 * int inode;
 * const IntArray* ct;
 * IntArray dofMask(2);
 * FloatArray vv;
 * Dof* iDof;
 *
 * dofMask.at(1) = V_u;
 * dofMask.at(2) = V_v;
 *
 * printf ("SUPG::updateDofManVals called\n");
 * for (inode=1; inode <= nnodes; inode++) {
 *  if (__DofManActivityMask.at(inode) == __NODE_INTERPOL_CANDIDATE) {
 *    // loop over shared elements
 *    ct = domain->giveConnectivityTable()->giveDofManConnectivityArray(inode);
 *    // now we have to find active dofmans in the neighborhood to interpolate velocity
 *    double vx=0.,vy=0.;
 *    int _j, _k, _d, _cnt=0;
 *    Element* _je;
 *    for (_j=1; _j<=ct->giveSize(); _j++) {
 *      _je=domain->giveElement(ct->at(_j));
 *      for (_k=1; _k<=_je->giveNumberOfNodes(); _k++) {
 *        int _code = __DofManActivityMask.at(_je->giveNode(_k)->giveNumber());
 *        if ((_code == __NODE_IN) || (_code == __NODE_OUT_ACTIVE)) { // active node (not candidate)
 *          _je->giveNode(_k)->giveUnknownVector(vv, dofMask, EID_MomentumBalance, VM_Total, tStep);
 *          vx += vv.at(1); vy += vv.at(2);
 *          _cnt++;
 *        }
 *      }
 *    }
 *    if (_cnt) printf ("%d ", inode);
 *    vx /= _cnt; vy /= _cnt; // average
 *    int ndofs = domain->giveDofManager(inode)->giveNumberOfDofs();
 *    for (_d=1;_d<=ndofs;_d++)  {
 *      iDof = domain->giveDofManager(inode)->giveDof(_d);
 *      DofIDItem type = iDof->giveDofID();
 *      if (!iDof->hasBc(tStep)) {
 *        if (type == V_u)
 *          iDof->updateUnknownsDictionary (tStep, EID_MomentumBalance_ConservationEquation, VM_Total, vx); // velocity
 *        else if (type == V_v)
 *          iDof->updateUnknownsDictionary (tStep, EID_MomentumBalance_ConservationEquation, VM_Total, vy); // velocity
 *          else if (type == V_w) {
 *            _error ("3d not yet supported");
 *          } else if (type == P_f) {
 *          //iDof->updateUnknownsDictionary (tStep, EID_MomentumBalance_ConservationEquation, VM_Total, 0.0);
 *          } else {_error2 ("unknown DOF type encountered (%s)", __DofIDItemToString (type));}
 *      }
 *    }
 *  }
 * }
 * }
 *
 * void
 * SUPG::imposeAmbientPressureInOuterNodes(SparseMtrx* lhs, FloatArray* rhs, TimeStep* stepN)
 * {
 * Domain* domain = this->giveDomain(1);
 * int nnodes = domain->giveNumberOfDofManagers();
 * int inode, pdof, _code;
 *
 * IntArray inMask(nnodes); inMask.zero();
 * for (inode=1;inode<=nnodes;inode++) {
 *  _code = __DofManActivityMask.at(inode);
 *  if ((_code == __NODE_OUT_ACTIVE) || (_code == __NODE_INTERPOL_CANDIDATE) || (_code == __NODE_OUT)) {
 *    if ((pdof=domain->giveDofManager(inode)->findDofWithDofId(P_f))) {
 *      int _eq = domain->giveDofManager(inode)->giveDof(pdof)->giveEquationNumber();
 *      if (_eq) {
 *        if (fabs(lhs->at(_eq,_eq)) > 0.0) {
 *          lhs->at(_eq,_eq) *= 1.e6;
 *          rhs->at(_eq) = 0.0;
 *        } else {
 *          lhs->at(_eq,_eq) = 1.0;
 *          rhs->at(_eq) = 0.0;
 *        }
 *      }
 *    }
 *  }
 * }
 * }
 *
 * void
 * SUPG::__debug (TimeStep* atTime)
 * {
 * int i, in, id, nincr = 1000;
 * int neq =  this -> giveNumberOfEquations (EID_MomentumBalance_ConservationEquation);
 * IntArray loc;
 * SUPGElement* element = (SUPGElement*) giveDomain(1)->giveElement(1);
 * FloatArray vincr(6), F(6), fi(6), fprev(6), fapprox(6), *solutionVector;
 * FloatMatrix d(6,6);
 *
 * vincr.at(1) = 1.e-2;
 * vincr.at(2) = 1.e-7;
 * vincr.at(3) = 2.e-2;
 * element -> giveLocationArray (loc, EID_MomentumBalance);
 * VelocityPressureField->advanceSolution(atTime);
 * solutionVector = VelocityPressureField->giveSolutionVector(atTime);
 * solutionVector->resize(neq);
 * // restrict vincr
 * for (i=1; i<=6; i++)
 *  if (loc.at(i) <= 0) vincr.at(i) = 0.0;
 *
 *
 * for (in=1; in<=nincr; in++) {
 *  for (i=1; i<=6; i++)
 *    if ((id = loc.at(i))) solutionVector->at(id) += vincr.at(i);
 *  // compute F from derivative
 *  EngngModel::giveElementCharacteristicMatrix(d, 1, TangentDiffusionDerivativeTerm_MB, atTime, giveDomain(1));
 *  fi.beProductOf (d, vincr);
 *  fapprox = fprev;
 *  fapprox.add(fi);
 *  // compute true F
 *  this->giveElementCharacteristicVector(F, 1, DiffusionTerm_MB, VM_Total, atTime, giveDomain(1));
 *  // print
 *  printf ("%d  %e %e %e %e %e %e   %e %e %e %e %e %e\n", in, F.at(1), F.at(2), F.at(3), F.at(4), F.at(5), F.at(6),
 *          fapprox.at(1), fapprox.at(2), fapprox.at(3), fapprox.at(4), fapprox.at(5), fapprox.at(6));
 *  //    if (in==1) fapprox = F;
 *  fprev = fapprox;
 * }
 *
 * }
 *
 */
} // end namespace oofem
