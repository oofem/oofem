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

#include "supg.h"
#include "sparsemtrx.h"
#include "nrsolver.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "initialcondition.h"
#include "maskedprimaryfield.h"
#include "verbose.h"
#include "supgelement.h"
#include "classfactory.h"
#include "mathfem.h"
#include "dofdistributedprimaryfield.h"
#include "leplic.h"
#include "levelsetpcs.h"
#include "datastream.h"
#include "function.h"
#include "contextioerr.h"
#include "timer.h"
#include "unknownnumberingscheme.h"
#include "boundarycondition.h"

namespace oofem {
/* define if implicit interface update required */
//#define SUPG_IMPLICIT_INTERFACE

REGISTER_EngngModel(SUPG);


SUPGInternalForceAssembler :: SUPGInternalForceAssembler(double l, double d, double u) : 
    VectorAssembler(), lscale(l), dscale(d), uscale(u)
{}

void SUPGInternalForceAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    SUPGElement *eptr = static_cast< SUPGElement * >( &element );
    FloatMatrix m1;
    FloatArray vacc, vtot, ptot, h;
    IntArray vloc, ploc;

    int size = eptr->computeNumberOfDofs();
    eptr->giveLocalVelocityDofMap(vloc);
    eptr->giveLocalPressureDofMap(ploc);
    vec.resize(size);
    vec.zero();

    ///@todo Computing v, p, a, is repeated in almost all the methods. We should pass these instead.
    /// I would do so right now, but I'm not sure which code to trust, as they all do different things.
    /// Some elements also use the prev solution for upwind stabilitzation, but, for example, quad10_supg does not.
    eptr->computeVectorOfVelocities(VM_Velocity, tStep, vacc);
    eptr->computeVectorOfVelocities(VM_Intermediate, tStep, vtot);
    eptr->computeVectorOfPressures(VM_Total, tStep, ptot);

    // MB contributions:
    // add (M+M_delta)*a
    eptr->computeAccelerationTerm_MB(m1, tStep);
    h.beProductOf(m1, vacc);
    vec.assemble(h, vloc);
    // add advection terms N+N_delta
    eptr->computeAdvectionTerm_MB(h, tStep);
    vec.assemble(h, vloc);
    // add diffusion terms (K+K_delta)h
    eptr->computeDiffusionTerm_MB(h, tStep);
    vec.assemble(h, vloc);
    // add lsic stabilization term
    eptr->computeLSICStabilizationTerm_MB(m1, tStep);
    m1.times( lscale / ( dscale * uscale * uscale ) );
    h.beProductOf(m1, vtot);
    vec.assemble(h, vloc);
    eptr->computeBCLhsTerm_MB(m1, tStep);
    if ( m1.isNotEmpty() ) {
        h.beProductOf(m1, vtot);
        vec.assemble(h, vloc);
    }

    // add pressure term
    eptr->computePressureTerm_MB(m1, tStep);
    h.beProductOf(m1, ptot); // term due to prescribed pressure
    vec.assemble(h, vloc);
    eptr->computeBCLhsPressureTerm_MB(m1, tStep);
    if ( m1.isNotEmpty() ) {
        h.beProductOf(m1, ptot);
        vec.assemble(h, vloc);
    }

    // MC contributions:
    // G^T term - linear advection term
    eptr->computeLinearAdvectionTerm_MC(m1, tStep);
    //m1.times(uscale/lscale);
    m1.times( 1. / ( dscale * uscale ) );
    eptr->computeVectorOfVelocities(VM_Intermediate, tStep, vtot);
    h.beProductOf(m1, vtot); // term due to prescribed velocity
    vec.assemble(h, ploc);
    // Diffusion term
    eptr->computeDiffusionTerm_MC(h, tStep);
    vec.assemble(h, ploc);
    // AccelerationTerm
    eptr->computeAccelerationTerm_MC(m1, tStep);
    h.beProductOf(m1, vacc);
    vec.assemble(h, ploc);
    eptr->computeBCLhsPressureTerm_MC(m1, tStep);
    if ( m1.isNotEmpty() ) {
        h.beProductOf(m1, vacc);
        vec.assemble(h, ploc);
    }

    // advection N term (nonlinear)
    eptr->computeAdvectionTerm_MC(h, tStep);
    vec.assemble(h, ploc);
    // pressure term
    eptr->computePressureTerm_MC(m1, tStep);
    h.beProductOf(m1, ptot);// term due to prescribed pressure
    vec.assemble(h, ploc);
};

SUPGTangentAssembler :: SUPGTangentAssembler(MatResponseMode m, double l, double d, double u, double a) : 
    MatrixAssembler(), rmode(m), lscale(l), dscale(d), uscale(u), alpha(a)
{}


void SUPGTangentAssembler :: matrixFromElement(FloatMatrix &answer, Element &el, TimeStep *tStep) const
{
    SUPGElement *element = static_cast< SUPGElement * >( &el );
    IntArray vloc, ploc;
    FloatMatrix h;
    int size = element->computeNumberOfDofs();
    element->giveLocalVelocityDofMap(vloc);
    element->giveLocalPressureDofMap(ploc);
    answer.resize(size, size);
    answer.zero();

    element->computeAccelerationTerm_MB(h, tStep);
    answer.assemble(h, vloc);
    element->computeAdvectionDerivativeTerm_MB(h, tStep);
    h.times( alpha * tStep->giveTimeIncrement() );
    answer.assemble(h, vloc);
    element->computeDiffusionDerivativeTerm_MB(h, rmode, tStep);

    h.times( alpha * tStep->giveTimeIncrement() );
    answer.assemble(h, vloc);
    element->computePressureTerm_MB(h, tStep);
    answer.assemble(h, vloc, ploc);
    element->computeLSICStabilizationTerm_MB(h, tStep);
    h.times( alpha * tStep->giveTimeIncrement() * lscale / ( dscale * uscale * uscale ) );
    answer.assemble(h, vloc);
    element->computeBCLhsTerm_MB(h, tStep);
    if ( h.isNotEmpty() ) {
        h.times( alpha * tStep->giveTimeIncrement() );
        answer.assemble(h, vloc);
    }

    element->computeBCLhsPressureTerm_MB(h, tStep);
    if ( h.isNotEmpty() ) {
        answer.assemble(h, vloc, ploc);
    }

    // conservation eq part
    element->computeLinearAdvectionTerm_MC(h, tStep);
    h.times( alpha * tStep->giveTimeIncrement() * 1.0 / ( dscale * uscale ) );
    answer.assemble(h, ploc, vloc);
    element->computeAdvectionDerivativeTerm_MC(h, tStep);
    h.times( alpha * tStep->giveTimeIncrement() );
    answer.assemble(h, ploc, vloc);
    element->computeAccelerationTerm_MC(h, tStep);
    answer.assemble(h, ploc, vloc);
    element->computeBCLhsPressureTerm_MC(h, tStep);
    if ( h.isNotEmpty() ) {
        answer.assemble(h, ploc, vloc);
    }

    element->computeDiffusionDerivativeTerm_MC(h, tStep);
    h.times( alpha * tStep->giveTimeIncrement() );
    answer.assemble(h, ploc, vloc);
    element->computePressureTerm_MC(h, tStep);
    answer.assemble(h, ploc);
}


SUPG :: SUPG(int i, EngngModel * _master) : FluidModel(i, _master)
{
    initFlag = 1;
    ndomains = 1;
    consistentMassFlag = 0;
    equationScalingFlag = false;
    lscale = uscale = dscale = 1.0;
}


SUPG :: ~SUPG()
{
}


NumericalMethod *SUPG :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !this->nMethod ) {
        this->nMethod.reset( classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this) ); 
        if ( !this->nMethod ) { 
            OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
        }
    }
    return this->nMethod.get();
}


IRResultType
SUPG :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = FluidModel :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    IR_GIVE_FIELD(ir, rtolv, _IFT_SUPG_rtolv);
    atolv = 1.e-15;
    IR_GIVE_OPTIONAL_FIELD(ir, atolv, _IFT_SUPG_atolv);


    stopmaxiter = ir->hasField(_IFT_SUPG_stopmaxiter);

    maxiter = 200;
    IR_GIVE_OPTIONAL_FIELD(ir, maxiter, _IFT_SUPG_maxiter);

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    IR_GIVE_FIELD(ir, deltaT, _IFT_SUPG_deltat);
    deltaTF = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaTF, _IFT_SUPG_deltatFunction);

    IR_GIVE_OPTIONAL_FIELD(ir, consistentMassFlag, _IFT_SUPG_cmflag);

    alpha = 0.5;
    IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_SUPG_alpha);

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_SUPG_scaleflag);
    equationScalingFlag = val > 0;
    if ( equationScalingFlag ) {
        IR_GIVE_FIELD(ir, lscale, _IFT_SUPG_lscale);
        IR_GIVE_FIELD(ir, uscale, _IFT_SUPG_uscale);
        IR_GIVE_FIELD(ir, dscale, _IFT_SUPG_dscale);
        double vref = 1.0; // reference viscosity
        Re = dscale * uscale * lscale / vref;
    } else {
        lscale = uscale = dscale = 1.0;
        Re = 1.0;
    }

    velocityPressureField = std::make_unique<DofDistributedPrimaryField>(this, 1, FT_VelocityPressure, 3, this->alpha);

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_SUPG_miflag);
    if ( val == 1 ) {
        this->materialInterface = std::make_unique<LEPlic>( 1, this->giveDomain(1) );
        this->materialInterface->initializeFrom(ir);
        // export velocity field
        FieldManager *fm = this->giveContext()->giveFieldManager();
        IntArray mask = {V_u, V_v, V_w};

        std :: shared_ptr< Field > _velocityField = std::make_shared<MaskedPrimaryField>( FT_Velocity, this->velocityPressureField.get(), mask );
        fm->registerField(_velocityField, FT_Velocity);

        //fsflag = 0;
        //IR_GIVE_OPTIONAL_FIELD (ir, fsflag, _IFT_SUPG_fsflag, "fsflag");
    } else if ( val == 2 ) {
        // positive coefficient scheme level set alg
        this->materialInterface = std::make_unique<LevelSetPCS>( 1, this->giveDomain(1) );
        this->materialInterface->initializeFrom(ir);
    }

    return IRRT_OK;
}


double
SUPG :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return this->velocityPressureField->giveUnknownValue(dof, mode, tStep);
}


void
SUPG :: computeInitialGuess(FloatArray &solution, TimeStep *tStep, Domain *d)
{
    // This is a bit narly, and I would welcome a nicer approach. The root problem is that the solution vector here is not stored
    // As part of the discretized solution. A good initial guess here is that the acceleration is constant:
    // a_{n} = a_{n-1} = (v_{n-1} - v_{n-2})/dt_{n-1}.
    // In addition, the solution vector includes other DOFs than velocities, which may have any other possible initial guess.
    // dp_{n} = dp_{n-1} = p_{n-1} - p_{n-2}
    // It might be worth simply assembling an old tangent and computing an initial guess. Such approach could handle changing BCs
    // which the current initial guess (constant acceleration) doesn't give good results fo. It would also let us ignore this whole mess.

    ///@todo This approach is broken because unkownsdicts can't handle more than 2 time-steps, because givePreviousStep only works 1 level.

    EModelDefaultEquationNumbering s;
    TimeStep *oldStep = tStep->givePreviousStep();

    auto a_to_v = [&solution, &s, &oldStep, this](DofManager &dman) {
        for ( Dof *dof: dman ) {
            DofIDItem id = dof->giveDofID();
            if ( !dof->isPrimaryDof() ) continue;
            int eqNum = dof->giveEquationNumber(s);
            if ( eqNum ) {
                ///@todo Couldn't/shouldn't SUPG elements just scale the elements with dt to obtain dv directly?
                // Then this whole thing could just be sol = old_sol + dsol
                if ( id == V_u || id == V_v || id == V_w ) {
                    solution.at(eqNum) = this->velocityPressureField->giveUnknownValue(dof, VM_Velocity, oldStep);
                } else {
                    solution.at(eqNum) = this->velocityPressureField->giveUnknownValue(dof, VM_Incremental, oldStep);
                }
            }
        }
    };

    for ( auto &node : d->giveDofManagers() ) {
        a_to_v(*node);
    }

    for ( auto &elem : d->giveElements() ) {
        int ndman = elem->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            a_to_v(*elem->giveInternalDofManager(i));
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        int ndman = bc->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            a_to_v(*bc->giveInternalDofManager(i));
        }
    }

    this->updateSolution(solution, tStep, d);
}


void SUPG :: updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d)
{
    // Two options here: Overload DofDistributedPrimaryField and implement the "mixed" time integration (w.r.t to pressure vs velocities).
    // It might also be nicer to just keep pressure and velocity separately.
    // Third alternative, compute the VM_Total value here, and use that to update the field.
    // This requires us to convert solutionVector (which contains [a, p] into [v, p]. Reminder:
    // Reminder, alpha doesn't come into play here: a = (v1 - v0)/dt --->   v1 = v0 + a*dt
    // alpha only affects VM_Intermediate, and is computed by velocityPressureField dynamically based on VM_Total.
    FloatArray totalSolution(solutionVector.giveSize());
    EModelDefaultEquationNumbering s;
    double dt = tStep->giveTimeIncrement();
    FloatArray prevSolutionVector;
    this->velocityPressureField->initialize(VM_Total, tStep->givePreviousStep(), prevSolutionVector, s);

    auto a_to_v = [&totalSolution, &prevSolutionVector, &solutionVector, &s, dt](DofManager &dman) {
        for ( Dof *dof: dman ) {
            DofIDItem id = dof->giveDofID();
            if ( !dof->isPrimaryDof() ) continue;
            int eqNum = dof->giveEquationNumber(s);
            if ( eqNum ) {
                ///@todo Couldn't/shouldn't SUPG elements just scale the elements with dt to obtain dv directly?
                // Then this whole thing could just be sol = old_sol + dsol
                if ( id == V_u || id == V_v || id == V_w ) {
                    totalSolution.at(eqNum) = prevSolutionVector.at(eqNum) + solutionVector.at(eqNum) * dt;
                } else {
                    totalSolution.at(eqNum) = prevSolutionVector.at(eqNum) + solutionVector.at(eqNum);
                }
            }
        }
    };

    for ( auto &node : d->giveDofManagers() ) {
        a_to_v(*node);
    }

    for ( auto &elem : d->giveElements() ) {
        int ndman = elem->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            a_to_v(*elem->giveInternalDofManager(i));
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        int ndman = bc->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            a_to_v(*bc->giveInternalDofManager(i));
        }
    }

    this->velocityPressureField->update(VM_Total, tStep, totalSolution, s);
    // update element stabilization
    for ( auto &elem : d->giveElements() ) {
        static_cast< FMElement * >( elem.get() )->updateStabilizationCoeffs(tStep);
    }
}


void SUPG :: updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm)
{
    answer.zero();
    this->assembleVector(answer, tStep, SUPGInternalForceAssembler(lscale, dscale, uscale), VM_Total,
                         EModelDefaultEquationNumbering(), d, eNorm);
    this->updateSharedDofManagers(answer, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
}


void SUPG :: updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d)
{
    mat.zero();
    //if ( 1 ) { //if ((nite > 5)) // && (rnorm < 1.e4))
    this->assemble( mat, tStep, SUPGTangentAssembler(TangentStiffness, lscale, dscale, uscale, alpha),
                    EModelDefaultEquationNumbering(), d );
    // } else {
    //     this->assemble( lhs, tStep, SUPGTangentAssembler(SecantStiffness),
    //                     EModelDefaultEquationNumbering(), d );
    // }
}


double
SUPG :: giveReynoldsNumber()
{
    if ( equationScalingFlag ) {
        return this->Re;
    } else {
        return 1.0;
    }
}


TimeStep *
SUPG :: giveSolutionStepWhenIcApply(bool force)
{
    if ( master && (!force)) {
        return master->giveSolutionStepWhenIcApply();
    } else {
        if ( !stepWhenIcApply ) {
            double dt = deltaT / this->giveVariableScale(VST_Time);

            stepWhenIcApply = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 0, 0.0, dt, 0);
        }

        return stepWhenIcApply.get();
    }
}


TimeStep *
SUPG :: giveNextStep()
{
    double dt = deltaT;

    Domain *domain = this->giveDomain(1);

    if ( !currentStep ) {
        // first step -> generate initial step
        currentStep = std::make_unique<TimeStep>( *giveSolutionStepWhenIcApply() );
    }

    previousStep = std :: move(currentStep);

    if ( deltaTF ) {
        dt *= domain->giveFunction(deltaTF)->evaluateAtTime(previousStep->giveNumber() + 1);
    }

    // check for critical time step
    for ( auto &elem : domain->giveElements() ) {
        dt = min( dt, static_cast< SUPGElement & >( *elem ).computeCriticalTimeStep(previousStep.get()) );
    }

    if ( materialInterface ) {
        dt = min( dt, materialInterface->computeCriticalTimeStep(previousStep.get()) );
    }

    // dt *= 0.6;
    dt /= this->giveVariableScale(VST_Time);

    currentStep = std::make_unique<TimeStep>(*previousStep, dt);

    OOFEM_LOG_INFO( "SolutionStep %d : t = %e, dt = %e\n", currentStep->giveNumber(), 
                    currentStep->giveTargetTime() * this->giveVariableScale(VST_Time), dt * this->giveVariableScale(VST_Time) );
    return currentStep.get();
}


void
SUPG :: solveYourselfAt(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
    FloatArray externalForces(neq);
    FloatArray internalForces(neq);
    FloatArray eNorm;

    if ( tStep->isTheFirstStep() ) {
        TimeStep *stepWhenIcApply = tStep->givePreviousStep();
        if ( materialInterface ) {
            materialInterface->initialize();
        }

        this->applyIC(stepWhenIcApply);
        //if (this->fsflag) this->updateDofManActivityMap(tStep);
    }

    velocityPressureField->advanceSolution(tStep);

    if ( !lhs ) {
        lhs.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        if ( !lhs ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        lhs->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    }

    if ( initFlag ) {

        if ( materialInterface ) {
            this->updateElementsForNewInterfacePosition(tStep);
        }

        initFlag = 0;
    }

    // evaluate element supg and sppg stabilization coeffs
    this->evaluateElementStabilizationCoeffs(tStep);

    if ( tStep->giveNumber() != 1 ) {
        if ( materialInterface ) {
            //if (this->fsflag) updateDofManVals(tStep);
#ifdef SUPG_IMPLICIT_INTERFACE
 #ifdef TIME_REPORT
            Timer timer;
            timer.startTimer();
 #endif
            materialInterface->updatePosition( this->giveCurrentStep() );
            updateElementsForNewInterfacePosition(tStep);
 #ifdef TIME_REPORT
            timer.stopTimer();
            OOFEM_LOG_INFO("SUPG info: user time consumed by updating interfaces: %.2fs\n", timer.getUtime();
                           );
 #endif
#else
            //updateElementsForNewInterfacePosition (tStep);
#endif
            //if (this->fsflag) this->updateDofManActivityMap(tStep);
        }
    }

    // Initial guess. Note, we should get this from the field, relying only on old values.
    // But there are shortcomings of fields that needs to be sorted out first. Until then, we keep the old solution vector cached.
    // Unfortunately, with a cached vector, we have no option to start from 0 if we renumber eqs.
    FloatArray incrementalSolutionVector(neq);
#if 0
    FloatArray solutionVector(neq);
    this->computeInitialGuess(solutionVector, tStep, domain);
#else
    if ( solutionVector.giveSize() != neq ) {
        solutionVector.resize(neq);
        solutionVector.zero();
    }
    this->updateSolution(solutionVector, tStep, domain);
#endif

    externalForces.zero();
    this->assembleVector(externalForces, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), domain);
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(),  LoadExchangeTag);

#if 0
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    this->initMetaStepAttributes( this->giveCurrentMetaStep() );
    double loadLevel;
    int currentIterations;
    this->updateInternalRHS( internalForces, tStep, domain, &eNorm );
    NM_Status status = this->nMethod->solve(*this->lhs,
                                            externalForces,
                                            NULL,
                                            solutionVector,
                                            incrementalSolutionVector,
                                            internalForces,
                                            eNorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total,
                                            currentIterations,
                                            tStep);

    if ( !( status & NM_Success ) ) {
        OOFEM_ERROR("No success in solving problem at time step", tStep->giveNumber());
    }
#else
    int nite = 0;
    double _absErrResid, err, avn = 0.0, aivn = 0.0, rnorm;
    FloatArray rhs;

    // algoritmic rhs part (assembled by e-model (in giveCharComponent service) from various element contribs)
    this->updateInternalRHS(internalForces, tStep, domain, &eNorm);

    rhs.beDifferenceOf(externalForces, internalForces);

    //OOFEM_LOG_INFO("Iteration  IncrSolErr      RelResidErr     AbsResidErr\n_______________________________________________________\n");
    OOFEM_LOG_INFO("Iteration  IncrSolErr      RelResidErr     (Rel_MB      ,Rel_MC      ) AbsResidErr\n__________________________________________________________________________________\n");
    do {
        nite++;
        //
        // Assemble lhs
        //
        //if (nite == 1)
        {
            // momentum balance part
            lhs->zero();
            if ( 1 ) { //if ((nite > 5)) // && (rnorm < 1.e4))
                this->assemble( *lhs, tStep, SUPGTangentAssembler(TangentStiffness, lscale, dscale, uscale, alpha),
                               EModelDefaultEquationNumbering(), this->giveDomain(1) );
            } else {
                this->assemble( *lhs, tStep, SUPGTangentAssembler(SecantStiffness, lscale, dscale, uscale, alpha),
                               EModelDefaultEquationNumbering(), this->giveDomain(1) );
            }
        }
        //if (this->fsflag) this->imposeAmbientPressureInOuterNodes(lhs,&rhs,tStep);

        nMethod->solve(*lhs, rhs, incrementalSolutionVector);
        solutionVector.add(incrementalSolutionVector);
        this->updateSolution(solutionVector, tStep, domain);

#if 0
 #ifdef SUPG_IMPLICIT_INTERFACE
        if ( materialInterface ) {
  #ifdef TIME_REPORT
            Timer timer;
            timer.startTimer();
  #endif
            //if (this->fsflag) updateDofManVals(tStep);
            materialInterface->updatePosition( this->giveCurrentStep() );
            updateElementsForNewInterfacePosition(tStep);
  #ifdef TIME_REPORT
            timer.stopTimer();
            OOFEM_LOG_INFO( "SUPG info: user time consumed by updating interfaces: %.2fs\n", timer.getUtime() );
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
        internalForces.zero();
        this->assembleVector(internalForces, tStep, SUPGInternalForceAssembler(lscale, dscale, uscale), VM_Total,
                             EModelDefaultEquationNumbering(), domain );
        this->updateSharedDofManagers(internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
        rhs.beDifferenceOf(externalForces, internalForces);

        // check convergence and repeat iteration if desired
        double rnorm_mb = 0.0, rnorm_mc = 0.0;
        for ( auto &inode : domain->giveDofManagers() ) {
            for ( Dof *dof: *inode ) {
                int jj;
                if ( ( jj = dof->__giveEquationNumber() ) ) {
                    DofIDItem type = dof->giveDofID();
                    double val = rhs.at(jj);
                    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                        rnorm_mb += val * val;
                    } else {
                        rnorm_mc += val * val;
                    }
                }
            }
        }

        rnorm_mb = sqrt(rnorm_mb);
        rnorm_mc = sqrt(rnorm_mc);

        rnorm = rhs.computeNorm();

        _absErrResid = 0.0;
        for ( int i = 1; i <= neq; i++ ) {
            _absErrResid = max( _absErrResid, fabs( rhs.at(i) ) );
        }

        //if ( requiresUnknownsDictionaryUpdate() ) {
        //    OOFEM_LOG_INFO("%-10d       n/a       %-15e %-15e\n", nite, rnorm, _absErrResid);
        //} else {
        OOFEM_LOG_INFO("%-10d %-15e %-15e (%-10e,%-10e) %-15e\n", nite, err, rnorm, rnorm_mb, rnorm_mc, _absErrResid);
        //}

        if ( 0 ) {
            // evaluate element supg and sppg stabilization coeffs
            this->evaluateElementStabilizationCoeffs(tStep);
            //
            // assemble rhs (residual)
            //
            internalForces.zero();
            this->assembleVector( internalForces, tStep, SUPGInternalForceAssembler(lscale, dscale, uscale), VM_Total,
                                 EModelDefaultEquationNumbering(), domain );
            this->updateSharedDofManagers(internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
            rhs.beDifferenceOf(externalForces, internalForces);
        }
    } while ( ( rnorm > rtolv ) && ( _absErrResid > atolv ) && ( nite <= maxiter ) );

    if ( nite <= maxiter ) {
        OOFEM_LOG_INFO("SUPG info: number of iterations: %d\n", nite);
    } else {
        OOFEM_WARNING("Convergence not reached, number of iterations: %d\n", nite);
        if ( stopmaxiter ) {
            exit(1);
        }
    }
#endif

#ifndef SUPG_IMPLICIT_INTERFACE
    if ( materialInterface ) {
 #ifdef TIME_REPORT
        Timer timer;
        timer.startTimer();
 #endif
        materialInterface->updatePosition( this->giveCurrentStep() );
        updateElementsForNewInterfacePosition(tStep);
 #ifdef TIME_REPORT
        timer.stopTimer();
        OOFEM_LOG_INFO( "SUPG info: user time consumed by updating interfaces: %.2fs\n", timer.getUtime() );
 #endif
        //if (this->fsflag) this->updateDofManActivityMap(tStep);
    }

#endif

    // update solution state counter
    tStep->incrementStateCounter();
}


int
SUPG :: forceEquationNumbering()
{
    lhs = nullptr;
    solutionVector.clear();
    return FluidModel::forceEquationNumbering();
}


bool
SUPG :: requiresEquationRenumbering(TimeStep *tStep)
{
    if ( tStep->isTheFirstStep() ) {
        return true;
    }
    // Check if Dirichlet b.c.s has changed.
    Domain *d = this->giveDomain(1);
    for ( auto &gbc : d->giveBcs() ) {
        BoundaryCondition *bc = dynamic_cast< BoundaryCondition * >(gbc.get());
        //ActiveBoundaryCondition *active_bc = dynamic_cast< ActiveBoundaryCondition * >(gbc.get());
        // We only need to consider Dirichlet b.c.s
        if ( bc ) {
            // Check of the dirichlet b.c. has changed in the last step (if so we need to renumber)
            if ( gbc->isImposed(tStep) != gbc->isImposed(tStep->givePreviousStep()) ) {
                return true;
            }
        }
    }
    return false;
}


void
SUPG :: updateYourself(TimeStep *tStep)
{
    this->updateInternalState(tStep);
    EngngModel :: updateYourself(tStep);
    if ( materialInterface ) {
        materialInterface->updateYourself(tStep);
    }

    //previousSolutionVector = solutionVector;
}


void
SUPG :: updateInternalState(TimeStep *tStep)
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
SUPG :: saveContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores;

    if ( ( iores = EngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    velocityPressureField->saveContext(stream);

    if ( materialInterface ) {
        if ( ( iores = materialInterface->saveContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType
SUPG :: restoreContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores;

    if ( ( iores = EngngModel :: restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    velocityPressureField->restoreContext(stream);

    if ( materialInterface ) {
        if ( ( iores = materialInterface->restoreContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


int
SUPG :: checkConsistency()
{
    // check internal consistency
    // if success returns nonzero
    Domain *domain = this->giveDomain(1);

    // check for proper element type
    for ( auto &elem : domain->giveElements() ) {
        if ( !dynamic_cast< SUPGElement * >( elem.get() ) ) {
            OOFEM_WARNING("Element %d has no SUPG base", elem->giveLabel());
            return 0;
        }
    }

    int ret = EngngModel :: checkConsistency();
    if ( ret == 0 ) {
        return 0;
    }

    // scale boundary and initial conditions
    if ( equationScalingFlag ) {
        for ( auto &bc : domain->giveBcs() ) {
            if ( bc->giveBCValType() == VelocityBVT ) {
                bc->scale(1. / uscale);
            } else if ( bc->giveBCValType() == PressureBVT ) {
                bc->scale( 1. / this->giveVariableScale(VST_Pressure) );
            } else if ( bc->giveBCValType() == ForceLoadBVT ) {
                bc->scale( 1. / this->giveVariableScale(VST_Force) );
            } else {
                OOFEM_ERROR("unknown bc/ic type");
            }
        }

        for ( auto &ic : domain->giveIcs() ) {
            if ( ic->giveICValType() == VelocityBVT ) {
                ic->scale(VM_Total, 1. / uscale);
            } else if ( ic->giveICValType() == PressureBVT ) {
                ic->scale( VM_Total, 1. / this->giveVariableScale(VST_Pressure) );
            } else {
                OOFEM_ERROR("unknown bc/ic type");
            }
        }
    }

    return 1;
}


void
SUPG :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->reinitialize();
    //this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}

void
SUPG :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
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
SUPG :: applyIC(TimeStep *stepWhenIcApply)
{
#ifdef VERBOSE
    OOFEM_LOG_INFO("Applying initial conditions\n");
#endif
    this->velocityPressureField->applyDefaultInitialCondition();

    // update element state according to given ic

    //this->initElementsForNewStep (stepWhenIcApply);
    if ( materialInterface ) {
        this->updateElementsForNewInterfacePosition(stepWhenIcApply);
    }
#if 0
    // Is this meaningful, or even needed?
    // We have to be *very* careful what the time integration and initial conditions actaully mean.
    // Internal states are computed in intermediate steps, which, doesn't exist until we proceed to the first step.
    // So, "stepWhenIcApply" makes no sense here.
    // And even then, just computing the internal forces will take care of this.
    Domain *domain = this->giveDomain(1);
    for ( auto &elem : domain->giveElements() ) {
        SUPGElement *element = static_cast< SUPGElement * >( elem.get() );
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
    }
#endif
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
        OOFEM_ERROR("unknown variable type");
    }

    return 0.0;
}


void
SUPG :: evaluateElementStabilizationCoeffs(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);

    for ( auto &elem : domain->giveElements() ) {
        FMElement *ePtr = static_cast< FMElement * >( elem.get() );
        ePtr->updateStabilizationCoeffs(tStep);
    }
}


void
SUPG :: updateElementsForNewInterfacePosition(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);

    OOFEM_LOG_DEBUG("SUPG :: updateElements - updating elements for interface position");

    for ( auto &elem : domain->giveElements() ) {
        SUPGElement *ePtr = static_cast< SUPGElement * >( elem.get() );
        ePtr->updateElementForNewInterfacePosition(tStep);
    }
}


int
SUPG :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % 2;
}

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
 *          _je->giveNode(_k)->giveUnknownVector(vv, dofMask, VM_Total, tStep);
 *          vx += vv.at(1); vy += vv.at(2);
 *          _cnt++;
 *        }
 *      }
 *    }
 *    if (_cnt) printf ("%d ", inode);
 *    vx /= _cnt; vy /= _cnt; // average
 *    for (Dof *iDof: domain->giveDofManager(inode))  {
 *      DofIDItem type = iDof->giveDofID();
 *      if (!iDof->hasBc(tStep)) {
 *        if (type == V_u)
 *          iDof->updateUnknownsDictionary (tStep, VM_Total, vx); // velocity
 *        else if (type == V_v)
 *          iDof->updateUnknownsDictionary (tStep, VM_Total, vy); // velocity
 *          else if (type == V_w) {
 *            OOFEM_ERROR("3d not yet supported");
 *          } else if (type == P_f) {
 *          //iDof->updateUnknownsDictionary (tStep, VM_Total, 0.0);
 *          } else {_error2 ("unknown DOF type encountered (%s)", __DofIDItemToString (type));}
 *      }
 *    }
 *  }
 * }
 * }
 *
 * void
 * SUPG::imposeAmbientPressureInOuterNodes(SparseMtrx* lhs, FloatArray* rhs, TimeStep* tStep)
 * {
 * Domain* domain = this->giveDomain(1);
 * int nnodes = domain->giveNumberOfDofManagers();
 * int inode, _code;
 * auto pdof;
 *
 * IntArray inMask(nnodes); inMask.zero();
 * for (inode=1;inode<=nnodes;inode++) {
 *  _code = __DofManActivityMask.at(inode);
 *  if ((_code == __NODE_OUT_ACTIVE) || (_code == __NODE_INTERPOL_CANDIDATE) || (_code == __NODE_OUT)) {
 *    if ((pdof=domain->giveDofManager(inode)->findDofWithDofId(P_f)) != domain->giveDofManager(inode)->end() ) {
 *      int _eq = (*pdof)->giveEquationNumber();
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
 * SUPG::__debug (TimeStep* tStep)
 * {
 * int i, in, id, nincr = 1000;
 * int neq =  this -> giveNumberOfDomainEquations (1, EModelDefaultEquationNumbering());
 * IntArray loc;
 * SUPGElement* element = (SUPGElement*) giveDomain(1)->giveElement(1);
 * FloatArray vincr(6), F(6), fi(6), fprev(6), fapprox(6), *solutionVector;
 * FloatMatrix d(6,6);
 *
 * vincr.at(1) = 1.e-2;
 * vincr.at(2) = 1.e-7;
 * vincr.at(3) = 2.e-2;
 * element -> giveLocationArray (loc);
 * velocityPressureField->advanceSolution(tStep);
 * solutionVector = velocityPressureField->giveSolutionVector(tStep);
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
 *  EngngModel::giveElementCharacteristicMatrix(d, 1, TangentDiffusionDerivativeTerm_MB, tStep, giveDomain(1));
 *  fi.beProductOf (d, vincr);
 *  fapprox = fprev;
 *  fapprox.add(fi);
 *  // compute true F
 *  this->giveElementCharacteristicVector(F, 1, DiffusionTerm_MB, VM_Total, tStep, giveDomain(1));
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
