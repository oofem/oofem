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

#include "cbs.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "initial.h"
#include "maskedprimaryfield.h"

#include "verbose.h"
#include "cbselement.h"
#include "usrdefsub.h"
#include "mathfem.h"
#include "datastream.h"
//<RESTRICTED_SECTION>
#include "leplic.h"
//</RESTRICTED_SECTION>
#ifdef TIME_REPORT
 #include "clock.h"
#endif
#include "contextioerr.h"

namespace oofem {

NumericalMethod *CBS :: giveNumericalMethod(MetaStep *mStep)
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
CBS :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    EngngModel :: initializeFrom(ir);
    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_CBS_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_CBS_smtype, "smtype"); // Macro
    sparseMtrxType = ( SparseMtrxType ) val;

    IR_GIVE_FIELD(ir, deltaT, IFT_CBS_deltat, "deltat"); // Macro
    minDeltaT = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, minDeltaT, IFT_CBS_mindeltat, "mindeltat"); // Macro

    IR_GIVE_OPTIONAL_FIELD(ir, consistentMassFlag, IFT_CBS_cmflag, "cmflag");

    theta [ 0 ] = theta [ 1 ] = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, theta [ 0 ], IFT_CBS_theta1, "theta1");
    IR_GIVE_OPTIONAL_FIELD(ir, theta [ 1 ], IFT_CBS_theta2, "theta2");

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_CBS_scaleflag, "scaleflag");
    equationScalingFlag = val;
    if ( equationScalingFlag ) {
        IR_GIVE_FIELD(ir, lscale, IFT_CBS_lscale, "lscale"); // Macro
        IR_GIVE_FIELD(ir, uscale, IFT_CBS_uscale, "uscale"); // Macro
        IR_GIVE_FIELD(ir, dscale, IFT_CBS_dscale, "dscale"); // Macro
        double vref = 1.0; // reference viscosity
        Re = dscale * uscale * lscale / vref;
    } else {
        lscale = uscale = dscale = 1.0;
        Re = 1.0;
    }

    //<RESTRICTED_SECTION>
    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_CBS_miflag, "miflag");
    if ( val ) {
        this->materialInterface = new LEPlic( 1, this->giveDomain(1) );
        // export velocity field
        FieldManager *fm = this->giveContext()->giveFieldManager();
        IntArray mask(3);
        mask.at(1) = V_u; mask.at(1) = V_v; mask.at(1) = V_w;
        MaskedPrimaryField* _velocityField = new MaskedPrimaryField (FT_Velocity, &this->VelocityField, mask);
        fm->registerField( _velocityField, FT_Velocity, true);
    }

    //</RESTRICTED_SECTION>

    return IRRT_OK;
}

double
CBS :: giveUnknownComponent(EquationID chc, ValueModeType mode,
                            TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displaacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }

    if ( chc == EID_ConservationEquation ) { // pressures
        return PressureField.giveUnknownValue(dof, mode, tStep);
    } else if ( chc == EID_MomentumBalance ) { // velocities
        return VelocityField.giveUnknownValue(dof, mode, tStep);
    } else if ( chc == EID_AuxMomentumBalance ) { // aux velocities
        switch ( mode ) {
        case VM_Incremental:
            if ( deltaAuxVelocity.isNotEmpty() ) {
                return deltaAuxVelocity.at(eq);
            } else {
                return 0.;
            }

        case VM_Total:
            _error("giveUnknownComponent: Unknown is of undefined ValueModeType for this problem");
        default:
            _error("giveUnknownComponent: Unknown is of undefined ValueModeType for this problem");
        }
    } else {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    return 0;
}


double
CBS :: giveUnknownComponent(UnknownType chc, ValueModeType mode,
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
    } else if ( chc == Theta_1 ) {
        return this->theta [ 0 ];
    } else if ( chc == Theta_2 ) {
        return this->theta [ 1 ];
    } else if ( chc == PrescribedTractionPressure ) {
        if ( mode == VM_Total ) {
            int eq = dof->__givePrescribedEquationNumber();
            if ( eq ) {
                return prescribedTractionPressure.at(eq);
            } else {
                _error("giveUnknownComponent: prescribed traction pressure requested for dof with no BC");
            }
        } else {
            _error("giveUnknownComponent: only total values supported for PrescribedTractionPressure");
        }
    } else {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
    }

    return 0;
}


TimeStep *
CBS :: giveSolutionStepWhenIcApply()
{
    if ( stepWhenIcApply == NULL ) {
        /*
         * stepWhenIcApply = new TimeStep (giveNumberOfTimeStepWhenIcApply(),this,0,
         * -deltaT,deltaT,0);
         */
        stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0,
                                       0.0, deltaT, 0);
    }

    return stepWhenIcApply;
}

TimeStep *
CBS :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    int i, nelem;
    double totalTime = 0;
    double dt = deltaT;
    StateCounterType counter = 1;
    delete previousStep;

    if ( currentStep == NULL ) {
        // first step -> generate initial step
        currentStep = new TimeStep( *giveSolutionStepWhenIcApply() );
    } else {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;

    // FORCE EQUATION NUMBERING
    this->giveNumberOfEquations(EID_MomentumBalance);
    Domain *domain = this->giveDomain(1);
    nelem = domain->giveNumberOfElements();
    // check for critical time step
    for ( i = 1; i <= nelem; i++ ) {
        dt = min( dt, ( ( CBSElement * ) domain->giveElement(i) )->computeCriticalTimeStep(previousStep) );
    }

    dt *= 0.6;
    dt = max(dt, minDeltaT);
    dt /= this->giveVariableScale(VST_Time);

    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTargetTime() + dt;
    }

    currentStep = new TimeStep(istep, this, 1, totalTime, dt, counter);
    // time and dt variables are set eq to 0 for staics - has no meaning

    OOFEM_LOG_INFO( "SolutionStep %d : t = %e, dt = %e\n", istep, totalTime * this->giveVariableScale(VST_Time), dt * this->giveVariableScale(VST_Time) );

    return currentStep;
}

void
CBS :: solveYourselfAt(TimeStep *tStep)
{
    int i;
    int momneq =  this->giveNumberOfEquations(EID_MomentumBalance);
    int presneq =  this->giveNumberOfEquations(EID_ConservationEquation);
    int presneq_prescribed = this->giveNumberOfPrescribedEquations(EID_ConservationEquation);
    double deltaT = tStep->giveTimeIncrement();

    FloatArray rhs(momneq);


    if ( initFlag ) {
        deltaAuxVelocity.resize(momneq);

        prescribedTractionPressure.resize(presneq_prescribed);
        nodalPrescribedTractionPressureConnectivity.resize(presneq_prescribed);
        this->assembleVectorFromElements( nodalPrescribedTractionPressureConnectivity, tStep, EID_ConservationEquation,
                                          NumberOfNodalPrescribedTractionPressureContributions, VM_Total,
                                          EModelDefaultPrescribedEquationNumbering(), this->giveDomain(1) );


        lhs = CreateUsrDefSparseMtrx(sparseMtrxType);
        if ( lhs == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        lhs->buildInternalStructure( this, 1, EID_ConservationEquation, EModelDefaultEquationNumbering() );

        this->assemble( lhs, stepWhenIcApply, EID_ConservationEquation, PressureLhs,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        lhs->times(deltaT * theta [ 0 ] * theta [ 1 ]);

        if ( consistentMassFlag ) {
            mss = CreateUsrDefSparseMtrx(sparseMtrxType);
            if ( mss == NULL ) {
                _error("solveYourselfAt: sparse matrix creation failed");
            }

            mss->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );
            this->assemble( mss, stepWhenIcApply, EID_MomentumBalance, MassMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        } else {
            mm.resize(momneq);
            mm.zero();
            this->assembleVectorFromElements( mm, tStep, EID_MomentumBalance, LumpedMassMatrix, VM_Total,
                                             EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }

        //<RESTRICTED_SECTION>
        // init material interface
        if ( materialInterface ) {
            materialInterface->initialize();
        }

        //</RESTRICTED_SECTION>
        initFlag = 0;
    }
    //<RESTRICTED_SECTION>
    else if ( materialInterface ) {
        lhs->zero();
        this->assemble( lhs, stepWhenIcApply, EID_ConservationEquation, PressureLhs,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        lhs->times(deltaT * theta [ 0 ] * theta [ 1 ]);

        if ( consistentMassFlag ) {
            mss->zero();
            this->assemble( mss, stepWhenIcApply, EID_MomentumBalance, MassMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        } else {
            mm.zero();
            this->assembleVectorFromElements( mm, tStep, EID_MomentumBalance, LumpedMassMatrix, VM_Total,
                                             EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }
    }

    //</RESTRICTED_SECTION>

    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        TimeStep *stepWhenIcApply = tStep->givePreviousStep();
        this->applyIC(stepWhenIcApply);
    }

    VelocityField.advanceSolution(tStep);
    PressureField.advanceSolution(tStep);
    FloatArray *velocityVector = VelocityField.giveSolutionVector(tStep);
    FloatArray *prevVelocityVector = VelocityField.giveSolutionVector( tStep->givePreviousStep() );
    FloatArray *pressureVector = PressureField.giveSolutionVector(tStep);

    /* STEP 1 - calculates auxiliary velocities*/
    rhs.zero();
    this->assembleVectorFromElements( rhs, tStep, EID_AuxMomentumBalance, IntermediateConvectionTerm, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleVectorFromElements( rhs, tStep, EID_AuxMomentumBalance, IntermediateDiffusionTerm, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    //this->assembleVectorFromElements(mm, tStep, EID_AuxMomentumBalance, LumpedMassMatrix, VM_Total, this->giveDomain(1));

    if ( consistentMassFlag ) {
        rhs.times(deltaT);
        this->assembleVectorFromElements( rhs, tStep, EID_AuxMomentumBalance, PrescribedVelocityRhsVector, VM_Incremental,
                                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
        nMethod->solve(mss, & rhs, & deltaAuxVelocity);
    } else {
        for ( i = 1; i <= momneq; i++ ) {
            deltaAuxVelocity.at(i) = deltaT * rhs.at(i) / mm.at(i);
        }
    }

    /* STEP 2 - calculates pressure (implicit solver) */
    rhs.resize(presneq);
    rhs.zero();
    this->assembleVectorFromElements( prescribedTractionPressure, tStep, EID_ConservationEquation,
                                      DensityPrescribedTractionPressure, VM_Total,
                                      EModelDefaultPrescribedEquationNumbering(), this->giveDomain(1) );
    for ( i = 1; i <= presneq_prescribed; i++ ) {
        prescribedTractionPressure.at(i) /= nodalPrescribedTractionPressureConnectivity.at(i);
    }

    //prescribedTractionPressure.printYourself();
    this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, DensityRhsVelocityTerms, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, DensityRhsPressureTerms, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    pressureVector->resize(presneq);
    nMethod->solve(lhs, & rhs, pressureVector);
    pressureVector->times(this->theta [ 1 ]);
    pressureVector->add(*PressureField.giveSolutionVector( tStep->givePreviousStep() ) );

    /* STEP 3 - velocity correction step */
    rhs.resize(momneq);
    rhs.zero();
    velocityVector->resize(momneq);
    this->assembleVectorFromElements( rhs, tStep, EID_MomentumBalance, CorrectionRhs, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    if ( consistentMassFlag ) {
        rhs.times(deltaT);
        //this->assembleVectorFromElements(rhs, tStep, EID_MomentumBalance, PrescribedRhsVector, VM_Incremental, this->giveDomain(1));
        nMethod->solve(mss, & rhs, velocityVector);
        velocityVector->add(deltaAuxVelocity);
        velocityVector->add(*prevVelocityVector);
    } else {
        for ( i = 1; i <= momneq; i++ ) {
            velocityVector->at(i) = prevVelocityVector->at(i) + deltaAuxVelocity.at(i) + deltaT *rhs.at(i) / mm.at(i);
        }
    }

    // update solution state counter
    tStep->incrementStateCounter();

    //<RESTRICTED_SECTION>
    if ( materialInterface ) {
#ifdef TIME_REPORT
        oofem_timeval tstart;
        getUtime(tstart);
#endif
        materialInterface->updatePosition( this->giveCurrentStep() );
#ifdef TIME_REPORT
        oofem_timeval ut;
        getRelativeUtime(ut, tstart);
        OOFEM_LOG_INFO( "CBS info: user time consumed by updating interfaces: %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif
    }

    //</RESTRICTED_SECTION>

}


void
CBS :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    EngngModel :: updateYourself(stepN);
    //<RESTRICTED_SECTION>
    if ( materialInterface ) {
        materialInterface->updateYourself(stepN);
    }

    //</RESTRICTED_SECTION>
    //previousSolutionVector = solutionVector;
}


void
CBS :: updateInternalState(TimeStep *stepN)
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


contextIOResultType
CBS :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( ( iores = PressureField.saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = VelocityField.saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = prescribedTractionPressure.storeYourself(stream, mode) ) != CIO_OK ) {
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
CBS :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( ( iores = PressureField.restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = VelocityField.restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = prescribedTractionPressure.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    } // ensure consistent records

    return CIO_OK;
}


int
CBS :: checkConsistency()
{
    // check internal consistency
    // if success returns nonzero
    int i, nelem;
    Element *ePtr;
    CBSElement *sePtr;
    GeneralBoundaryCondition *bcPtr;
    InitialCondition *icPtr;
    Domain *domain = this->giveDomain(1);

    nelem = domain->giveNumberOfElements();
    // check for proper element type

    for ( i = 1; i <= nelem; i++ ) {
        ePtr = domain->giveElement(i);
        sePtr = dynamic_cast< CBSElement * >(ePtr);
        if ( sePtr == NULL ) {
            _warning2("Element %d has no CBS base", i);
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
CBS :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}


void
CBS :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
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
CBS :: applyIC(TimeStep *stepWhenIcApply)
{
    Domain *domain = this->giveDomain(1);
    int mbneq =  this->giveNumberOfEquations(EID_MomentumBalance);
    int pdneq =  this->giveNumberOfEquations(EID_ConservationEquation);
    FloatArray *velocityVector, *pressureVector;

#ifdef VERBOSE
    OOFEM_LOG_INFO("Applying initial conditions\n");
#endif
    int nDofs, j, k, jj;
    int nman  = domain->giveNumberOfDofManagers();
    DofManager *node;
    Dof *iDof;
    DofIDItem type;

    VelocityField.advanceSolution(stepWhenIcApply);
    velocityVector = VelocityField.giveSolutionVector(stepWhenIcApply);
    velocityVector->resize(mbneq);
    velocityVector->zero();

    PressureField.advanceSolution(stepWhenIcApply);
    pressureVector = PressureField.giveSolutionVector(stepWhenIcApply);
    pressureVector->resize(pdneq);
    pressureVector->zero();


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
                    velocityVector->at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Total, stepWhenIcApply);
                } else {
                    pressureVector->at(jj) = iDof->giveUnknown(EID_ConservationEquation, VM_Total, stepWhenIcApply);
                }
            }
        }
    }

    // update element state according to given ic
    int nelem = domain->giveNumberOfElements();
    CBSElement *element;

    for ( j = 1; j <= nelem; j++ ) {
        element = ( CBSElement * ) domain->giveElement(j);
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
    }
}


int
CBS :: giveNewEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return ++numberOfMomentumEqs;
    } else if ( id == P_f ) {
        return ++numberOfConservationEqs;
    } else {
        _error("giveNewEquationNumber:: Unknown DofIDItem");
    }

    return 0;
}


int
CBS :: giveNewPrescribedEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return ++numberOfPrescribedMomentumEqs;
    } else if ( id == P_f ) {
        return ++numberOfPrescribedConservationEqs;
    } else {
        _error("giveNewPrescribedEquationNumber:: Unknown DofIDItem");
    }

    return 0;
}


int
CBS :: giveNumberOfEquations(EquationID id)
{
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    if ( ( id == EID_MomentumBalance ) || ( id == EID_AuxMomentumBalance ) ) {
        return numberOfMomentumEqs;
    } else if ( id == EID_ConservationEquation ) {
        return numberOfConservationEqs;
    } else {
        _error("giveNumberOfEquations: unknown equation id");
    }

    return 0;
}


int CBS :: giveNumberOfPrescribedEquations(EquationID id)
{
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    if ( ( id == EID_MomentumBalance ) || ( id == EID_AuxMomentumBalance ) ) {
        return numberOfPrescribedMomentumEqs;
    } else if ( id == EID_ConservationEquation ) {
        return numberOfPrescribedConservationEqs;
    } else {
        _error("giveNumberOfPrescribedEquations: unknown equation id");
    }

    return 0;
}


int CBS :: giveNumberOfDomainEquations(int d, EquationID id)
{
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    if ( ( id == EID_MomentumBalance ) || ( id == EID_AuxMomentumBalance ) ) {
        return numberOfMomentumEqs;
    } else if ( id == EID_ConservationEquation ) {
        return numberOfConservationEqs;
    } else {
        _error("giveNumberOfDomainEquations: unknown equation id");
    }

    return 0;
}


int CBS :: giveNumberOfPrescribedDomainEquations(int d, EquationID id)
{
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    if ( ( id == EID_MomentumBalance ) || ( id == EID_AuxMomentumBalance ) ) {
        return numberOfPrescribedMomentumEqs;
    } else if ( id == EID_ConservationEquation ) {
        return numberOfPrescribedConservationEqs;
    } else {
        _error("giveNumberOfPrescribedDomainEquations: unknown equation id");
    }

    return 0;
}


double CBS :: giveVariableScale(VarScaleType varID)
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

/*
 * CBS :: printOutputAt (FILE * File,TimeStep* stepN)
 * {
 * //FILE* File = this -> giveDomain() -> giveOutputStream() ;
 * int domCount = 0, idomain;
 * Domain* domain;
 *
 * // fprintf (File,"\nOutput for time step number %d \n\n",stepN->giveNumber());
 * for (idomain = 1; idomain <= this->ndomains; idomain++) {
 * domain= this->giveDomain(idomain);
 * domCount += domain->giveOutputManager()->testTimeStepOutput (stepN);
 * }
 * if (domCount == 0) return;  // do not print even Solution step header
 *
 * fprintf (File,"\nOutput for time % .8e \n\n",stepN->giveTime()/this->giveVariableScale (VST_Time));
 * for (idomain = 1; idomain <= this->ndomains; idomain++) {
 *
 * domain= this->giveDomain(idomain);
 * fprintf (File,"\nOutput for domain %3d\n",domain->giveNumber());
 *
 * domain->giveOutputManager()->doDofManOutput  (File, stepN);
 * domain->giveOutputManager()->doElementOutput (File, stepN);
 * }
 * }
 *
 */
} // end namespace oofem
