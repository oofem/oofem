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

#include "cbs.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "initialcondition.h"
#include "maskedprimaryfield.h"
#include "verbose.h"
#include "cbselement.h"
#include "classfactory.h"
#include "mathfem.h"
#include "datastream.h"
//<RESTRICTED_SECTION>
#include "leplic.h"
//</RESTRICTED_SECTION>
#ifdef TIME_REPORT
 #include "timer.h"
#endif
#include "contextioerr.h"

namespace oofem {
REGISTER_EngngModel(CBS);


void NumberOfNodalPrescribedTractionPressureAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    static_cast< CBSElement & >( element ).computeNumberOfNodalPrescribedTractionPressureContributions(vec, tStep);
}

void IntermediateConvectionDiffusionAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    FloatArray vec_p;
    static_cast< CBSElement & >( element ).computeConvectionTermsI(vec, tStep);
    static_cast< CBSElement & >( element ).computeDiffusionTermsI(vec_p, tStep);
    vec.add(vec_p);
}

void PrescribedVelocityRhsAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    static_cast< CBSElement & >( element ).computePrescribedTermsI(vec, tStep);
    
}

void DensityPrescribedTractionPressureAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    static_cast< CBSElement & >( element ).computePrescribedTractionPressure(vec, tStep);
}

void DensityRhsAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    FloatArray vec_p;
    static_cast< CBSElement & >( element ).computeDensityRhsVelocityTerms(vec, tStep);
    static_cast< CBSElement & >( element ).computeDensityRhsPressureTerms(vec_p, tStep);
    vec.add(vec_p);
}

void CorrectionRhsAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    static_cast< CBSElement & >( element ).computeCorrectionRhs(vec, tStep);
    
}

void PressureLhsAssembler :: matrixFromElement(FloatMatrix &answer, Element &element, TimeStep *tStep) const
{
    static_cast< CBSElement & >( element ).computePressureLhs(answer, tStep);
}

void PressureLhsAssembler :: locationFromElement(IntArray& loc, Element& element, const UnknownNumberingScheme& s, IntArray* dofIds) const
{
    element.giveLocationArray(loc, {P_f}, s, dofIds);
}



CBS :: CBS(int i, EngngModel* _master) : FluidModel ( i, _master ),
    PressureField ( this, 1, FT_Pressure, 1 ),
    VelocityField ( this, 1, FT_Velocity, 1 ),
    vnum ( false ), vnumPrescribed ( true ), pnum ( false ), pnumPrescribed ( true )
{
    initFlag = 1;
    ndomains = 1;
    consistentMassFlag = 0;
    equationScalingFlag = false;
    lscale = uscale = dscale = 1.0;
}

CBS :: ~CBS()
{
}

NumericalMethod *CBS :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod.reset( classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this) );
        if ( !nMethod ) {
            OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
        }
    }
    return nMethod.get();
}

IRResultType
CBS :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    IR_GIVE_FIELD(ir, deltaT, _IFT_CBS_deltat);
    minDeltaT = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, minDeltaT, _IFT_CBS_mindeltat);

    IR_GIVE_OPTIONAL_FIELD(ir, consistentMassFlag, _IFT_CBS_cmflag);

    theta1 = theta2 = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, theta1, _IFT_CBS_theta1);
    IR_GIVE_OPTIONAL_FIELD(ir, theta2, _IFT_CBS_theta2);

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_CBS_scaleflag);
    equationScalingFlag = val > 0;
    if ( equationScalingFlag ) {
        IR_GIVE_FIELD(ir, lscale, _IFT_CBS_lscale);
        IR_GIVE_FIELD(ir, uscale, _IFT_CBS_uscale);
        IR_GIVE_FIELD(ir, dscale, _IFT_CBS_dscale);
        double vref = 1.0; // reference viscosity
        Re = dscale * uscale * lscale / vref;
    } else {
        lscale = uscale = dscale = 1.0;
        Re = 1.0;
    }

    //<RESTRICTED_SECTION>
    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_CBS_miflag);
    if ( val ) {
        this->materialInterface.reset( new LEPlic( 1, this->giveDomain(1) ) );
        // export velocity field
        FieldManager *fm = this->giveContext()->giveFieldManager();
        IntArray mask = {V_u, V_v, V_w};

        //std::shared_ptr<Field> _velocityField = make_shared<MaskedPrimaryField>(FT_Velocity, &this->VelocityField, mask);
        std :: shared_ptr< Field > _velocityField( new MaskedPrimaryField ( FT_Velocity, &this->VelocityField, mask ) );
        fm->registerField(_velocityField, FT_Velocity);
    }
    //</RESTRICTED_SECTION>

    return EngngModel :: initializeFrom(ir);
}


double
CBS :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of dof
{
#ifdef DEBUG
    if ( dof->__giveEquationNumber() == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }

#endif

    if ( dof->giveDofID() == P_f ) { // pressures
        return PressureField.giveUnknownValue(dof, mode, tStep);
    } else { // velocities
        return VelocityField.giveUnknownValue(dof, mode, tStep);
    }
}


double
CBS :: giveReynoldsNumber()
{
    return equationScalingFlag ? this->Re : 1.0;
}


double CBS :: giveTheta1() { return this->theta1; }

double CBS :: giveTheta2() { return this->theta2; }

double
CBS :: giveTractionPressure(Dof *dof)
{
    ///@todo This should just disappear completely.
    int eq = dof->__givePrescribedEquationNumber();
    if ( eq ) {
        return prescribedTractionPressure.at(eq);
    } else {
        OOFEM_ERROR("prescribed traction pressure requested for dof with no BC");
    }

    return 0;
}


TimeStep *
CBS :: giveSolutionStepWhenIcApply(bool force)
{
  if ( master && (!force)) {
    return master->giveSolutionStepWhenIcApply();
  } else {
    if ( !stepWhenIcApply ) {
        stepWhenIcApply.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0,
                                       0.0, deltaT, 0) );
    }

    return stepWhenIcApply.get();
  }
}

TimeStep *
CBS :: giveNextStep()
{
    double dt = deltaT;
    if ( !currentStep ) {
        // first step -> generate initial step
        currentStep.reset( new TimeStep( *giveSolutionStepWhenIcApply() ) );
    }

    previousStep = std :: move(currentStep);

    Domain *domain = this->giveDomain(1);
    // check for critical time step
    for ( auto &elem : domain->giveElements() ) {
        dt = min( dt, static_cast< CBSElement & >( *elem ).computeCriticalTimeStep(previousStep.get()) );
    }

    dt *= 0.6;
    dt = max(dt, minDeltaT);
    dt /= this->giveVariableScale(VST_Time);

    currentStep.reset( new TimeStep(*previousStep, dt) );

    OOFEM_LOG_INFO( "SolutionStep %d : t = %e, dt = %e\n", currentStep->giveNumber(),
                    currentStep->giveTargetTime() * this->giveVariableScale(VST_Time), dt * this->giveVariableScale(VST_Time) );

    return currentStep.get();
}

void
CBS :: solveYourselfAt(TimeStep *tStep)
{
    int momneq = this->giveNumberOfDomainEquations(1, vnum);
    int presneq = this->giveNumberOfDomainEquations(1, pnum);
    int presneq_prescribed = this->giveNumberOfDomainEquations(1, pnumPrescribed);
    double deltaT = tStep->giveTimeIncrement();

    FloatArray rhs(momneq);

    if ( initFlag ) {
        deltaAuxVelocity.resize(momneq);

        nodalPrescribedTractionPressureConnectivity.resize(presneq_prescribed);
        nodalPrescribedTractionPressureConnectivity.zero();
        this->assembleVectorFromElements( nodalPrescribedTractionPressureConnectivity, tStep,
                                         NumberOfNodalPrescribedTractionPressureAssembler(), VM_Total,
                                         pnumPrescribed, this->giveDomain(1) );


        lhs.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        if ( !lhs ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        lhs->buildInternalStructure(this, 1, pnum);

        this->assemble( *lhs, stepWhenIcApply.get(), PressureLhsAssembler(),
                       pnum, this->giveDomain(1) );
        lhs->times(deltaT * theta1 * theta2);

        if ( consistentMassFlag ) {
            mss.reset( classFactory.createSparseMtrx(sparseMtrxType) );
            if ( !mss ) {
                OOFEM_ERROR("sparse matrix creation failed");
            }

            mss->buildInternalStructure(this, 1, vnum);
            this->assemble( *mss, stepWhenIcApply.get(), MassMatrixAssembler(),
                           vnum, this->giveDomain(1) );
        } else {
            mm.resize(momneq);
            mm.zero();
            this->assembleVectorFromElements( mm, tStep, LumpedMassVectorAssembler(), VM_Total,
                                             vnum, this->giveDomain(1) );
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
        this->assemble( *lhs, stepWhenIcApply.get(), PressureLhsAssembler(),
                       pnum, this->giveDomain(1) );
        lhs->times(deltaT * theta1 * theta2);

        if ( consistentMassFlag ) {
            mss->zero();
            this->assemble( *mss, stepWhenIcApply.get(), MassMatrixAssembler(),
                           vnum, this->giveDomain(1) );
        } else {
            mm.zero();
            this->assembleVectorFromElements( mm, tStep, LumpedMassVectorAssembler(), VM_Total,
                                             vnum, this->giveDomain(1) );
        }
    }

    //</RESTRICTED_SECTION>

    if ( tStep->isTheFirstStep() ) {
        TimeStep *stepWhenIcApply = tStep->givePreviousStep();
        this->applyIC(stepWhenIcApply);
    }

    VelocityField.advanceSolution(tStep);
    PressureField.advanceSolution(tStep);
    FloatArray *velocityVector = VelocityField.giveSolutionVector(tStep);
    FloatArray *prevVelocityVector = VelocityField.giveSolutionVector( tStep->givePreviousStep() );
    FloatArray *pressureVector = PressureField.giveSolutionVector(tStep);
    FloatArray *prevPressureVector = PressureField.giveSolutionVector( tStep->givePreviousStep() );

    velocityVector->resize(momneq);
    pressureVector->resize(presneq);

    /* STEP 1 - calculates auxiliary velocities*/
    rhs.zero();
    // Depends on old v:
    this->assembleVectorFromElements( rhs, tStep, IntermediateConvectionDiffusionAssembler(), VM_Total, vnum, this->giveDomain(1) );
    //this->assembleVectorFromElements(mm, tStep, LumpedMassVectorAssembler(), VM_Total, this->giveDomain(1));

    if ( consistentMassFlag ) {
        rhs.times(deltaT);
        // Depends on prescribed v
        this->assembleVectorFromElements( rhs, tStep, PrescribedVelocityRhsAssembler(), VM_Total, vnum, this->giveDomain(1) );
        nMethod->solve(*mss, rhs, deltaAuxVelocity);
    } else {
        for ( int i = 1; i <= momneq; i++ ) {
            deltaAuxVelocity.at(i) = deltaT * rhs.at(i) / mm.at(i);
        }
    }

    /* STEP 2 - calculates pressure (implicit solver) */
    this->prescribedTractionPressure.resize(presneq_prescribed);
    this->prescribedTractionPressure.zero();
    this->assembleVectorFromElements( prescribedTractionPressure, tStep,
                                     DensityPrescribedTractionPressureAssembler(), VM_Total,
                                     pnumPrescribed, this->giveDomain(1) );
    for ( int i = 1; i <= presneq_prescribed; i++ ) {
        prescribedTractionPressure.at(i) /= nodalPrescribedTractionPressureConnectivity.at(i);
    }

    // DensityRhsVelocityTerms needs this: Current velocity without correction;
    * velocityVector = * prevVelocityVector;
    velocityVector->add(this->theta1, deltaAuxVelocity);

    // Depends on old V + deltaAuxV * theta1 and p:
    rhs.resize(presneq);
    rhs.zero();
    this->assembleVectorFromElements( rhs, tStep, DensityRhsAssembler(), VM_Total,
                                     pnum, this->giveDomain(1) );
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    nMethod->solve(*lhs, rhs, *pressureVector);
    pressureVector->times(this->theta2);
    pressureVector->add(* prevPressureVector);

    /* STEP 3 - velocity correction step */
    rhs.resize(momneq);
    rhs.zero();
    // Depends on p:
    this->assembleVectorFromElements( rhs, tStep, CorrectionRhsAssembler(), VM_Total,
                                     vnum, this->giveDomain(1) );
    if ( consistentMassFlag ) {
        rhs.times(deltaT);
        //this->assembleVectorFromElements(rhs, tStep, PrescribedRhsAssembler(), VM_Incremental, vnum, this->giveDomain(1));
        nMethod->solve(*mss, rhs, *velocityVector);
        velocityVector->add(deltaAuxVelocity);
        velocityVector->add(* prevVelocityVector);
    } else {
        for ( int i = 1; i <= momneq; i++ ) {
            velocityVector->at(i) = prevVelocityVector->at(i) + deltaAuxVelocity.at(i) + deltaT *rhs.at(i) / mm.at(i);
        }
    }

    // update solution state counter
    tStep->incrementStateCounter();

    //<RESTRICTED_SECTION>
    if ( materialInterface ) {
#ifdef TIME_REPORT
        Timer timer;
        timer.startTimer();
#endif
        materialInterface->updatePosition( this->giveCurrentStep() );
#ifdef TIME_REPORT
        timer.stopTimer();
        OOFEM_LOG_INFO( "CBS info: user time consumed by updating interfaces: %.2fs\n", timer.getUtime() );
#endif
    }

    //</RESTRICTED_SECTION>
}


void
CBS :: updateYourself(TimeStep *tStep)
{
    this->updateInternalState(tStep);
    EngngModel :: updateYourself(tStep);
    //<RESTRICTED_SECTION>
    if ( materialInterface ) {
        materialInterface->updateYourself(tStep);
    }

    //</RESTRICTED_SECTION>
    //previousSolutionVector = solutionVector;
}


void
CBS :: updateInternalState(TimeStep *tStep)
{
    for ( auto &domain: domainList ) {
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( auto &dman : domain->giveDofManagers() ) {
                this->updateDofUnknownsDictionary(dman.get(), tStep);
            }
        }

        for ( auto &elem : domain->giveElements() ) {
            elem->updateInternalState(tStep);
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

    if ( ( iores = PressureField.saveContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = VelocityField.saveContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = prescribedTractionPressure.storeYourself(*stream) ) != CIO_OK ) {
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
    FILE *file = NULL;

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

    if ( ( iores = PressureField.restoreContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = VelocityField.restoreContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = prescribedTractionPressure.restoreYourself(*stream) ) != CIO_OK ) {
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
    Domain *domain = this->giveDomain(1);

    // check for proper element type
    for ( auto &elem : domain->giveElements() ) {
        if ( !dynamic_cast< CBSElement * >( elem.get() ) ) {
            OOFEM_WARNING("Element %d has no CBS base", elem->giveLabel());
            return 0;
        }
    }

    EngngModel :: checkConsistency();


    // scale boundary and initial conditions
    if ( equationScalingFlag ) {
        for ( auto &bc: domain->giveBcs() ) {
            if ( bc->giveBCValType() == VelocityBVT ) {
                bc->scale(1. / uscale);
            } else if ( bc->giveBCValType() == PressureBVT ) {
                bc->scale( 1. / this->giveVariableScale(VST_Pressure) );
            } else if ( bc->giveBCValType() == ForceLoadBVT ) {
                bc->scale( 1. / this->giveVariableScale(VST_Force) );
            } else {
                OOFEM_WARNING("unknown bc/ic type");
                return 0;
            }
        }

        for ( auto &ic : domain->giveIcs() ) {
            if ( ic->giveICValType() == VelocityBVT ) {
                ic->scale(VM_Total, 1. / uscale);
            } else if ( ic->giveICValType() == PressureBVT ) {
                ic->scale( VM_Total, 1. / this->giveVariableScale(VST_Pressure) );
            } else {
                OOFEM_WARNING("unknown bc/ic type");
                return 0;
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
CBS :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    double pscale = ( dscale * uscale * uscale );

    DofIDItem type = iDof->giveDofID();
    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
        iDof->printSingleOutputAt(stream, tStep, 'd', VM_Total, uscale);
    } else if ( type == P_f ) {
        iDof->printSingleOutputAt(stream, tStep, 'd', VM_Total, pscale);
    } else {
        OOFEM_ERROR("unsupported dof type");
    }
}


void
CBS :: applyIC(TimeStep *stepWhenIcApply)
{
    Domain *domain = this->giveDomain(1);
    int mbneq = this->giveNumberOfDomainEquations(1, vnum);
    int pdneq = this->giveNumberOfDomainEquations(1, pnum);
    FloatArray *velocityVector, *pressureVector;

#ifdef VERBOSE
    OOFEM_LOG_INFO("Applying initial conditions\n");
#endif

    VelocityField.advanceSolution(stepWhenIcApply);
    velocityVector = VelocityField.giveSolutionVector(stepWhenIcApply);
    velocityVector->resize(mbneq);
    velocityVector->zero();

    PressureField.advanceSolution(stepWhenIcApply);
    pressureVector = PressureField.giveSolutionVector(stepWhenIcApply);
    pressureVector->resize(pdneq);
    pressureVector->zero();

    for ( auto &node : domain->giveDofManagers() ) {
        for ( Dof *iDof: *node ) {
            // ask for initial values obtained from
            // bc (boundary conditions) and ic (initial conditions)
            if ( !iDof->isPrimaryDof() ) {
                continue;
            }

            int jj = iDof->__giveEquationNumber();
            if ( jj ) {
                DofIDItem type = iDof->giveDofID();
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                    velocityVector->at(jj) = iDof->giveUnknown(VM_Total, stepWhenIcApply);
                } else {
                    pressureVector->at(jj) = iDof->giveUnknown(VM_Total, stepWhenIcApply);
                }
            }
        }
    }

    // update element state according to given ic
    for ( auto &elem : domain->giveElements() ) {
        CBSElement *element = static_cast< CBSElement * >( elem.get() );
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
    }
}


int
CBS :: giveNewEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return this->vnum.askNewEquationNumber();
    } else if ( id == P_f ) {
        return this->pnum.askNewEquationNumber();
    } else {
        OOFEM_ERROR("Unknown DofIDItem");
    }

    return 0;
}


int
CBS :: giveNewPrescribedEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return this->vnumPrescribed.askNewEquationNumber();
    } else if ( id == P_f ) {
        return this->pnumPrescribed.askNewEquationNumber();
    } else {
        OOFEM_ERROR("Unknown DofIDItem");
    }

    return 0;
}


int CBS :: giveNumberOfDomainEquations(int d, const UnknownNumberingScheme &num)
{
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    return num.giveRequiredNumberOfDomainEquation();
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
        OOFEM_ERROR("unknown variable type");
    }

    return 0.0;
}

#if 0
void CBS :: printOutputAt(FILE *File, TimeStep *tStep)
{
    //FILE* File = this->giveDomain()->giveOutputStream();
    int domCount = 0;
    // fprintf (File,"\nOutput for time step number %d \n\n",tStep->giveNumber());
    for ( auto &domain: this->domainList ) {
        domCount += domain->giveOutputManager()->testTimeStepOutput(tStep);
    }

    if ( domCount == 0 ) {
        return;                 // do not print even Solution step header
    }

    fprintf( File, "\nOutput for time % .8e \n\n", tStep->giveTime() / this->giveVariableScale(VST_Time) );
    for ( auto &domain: this->domainList ) {
        fprintf( File, "\nOutput for domain %3d\n", domain->giveNumber() );
        domain->giveOutputManager()->doDofManOutput(File, tStep);
        domain->giveOutputManager()->doElementOutput(File, tStep);
    }
}
#endif
} // end namespace oofem
