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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "tm/EngineeringModels/transienttransportproblem.h"
#include "timestep.h"
#include "dofdistributedprimaryfield.h"
#include "maskedprimaryfield.h"
#include "intvarfield.h"
#include "tm/Elements/transportelement.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "nrsolver.h"
#include "unknownnumberingscheme.h"
#include "function.h"
#include "dofmanager.h"
#include "assemblercallback.h"
// Temporary:
#include "generalboundarycondition.h"
#include "boundarycondition.h"
#include "activebc.h"
#include "outputmanager.h"

namespace oofem {
REGISTER_EngngModel(TransientTransportProblem);

TransientTransportProblem :: TransientTransportProblem(int i, EngngModel *master) : EngngModel(i, master)
{
    ndomains = 1;
}


NumericalMethod *TransientTransportProblem :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod = std::make_unique<NRSolver>(this->giveDomain(1), this);
    }
    return nMethod.get();
}


void
TransientTransportProblem :: initializeFrom(InputRecord &ir)
{
    EngngModel :: initializeFrom(ir);

    int val = SMT_Skyline;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    this->sparseMtrxType = ( SparseMtrxType ) val;

    IR_GIVE_FIELD(ir, this->alpha, _IFT_TransientTransportProblem_alpha);

    if ( ir.hasField(_IFT_TransientTransportProblem_initt) ) {
        IR_GIVE_FIELD(ir, initT, _IFT_TransientTransportProblem_initt);
    }
    
    prescribedTimes.clear();
    dtFunction = 0;
    if ( ir.hasField(_IFT_TransientTransportProblem_dtFunction) ) {
        IR_GIVE_FIELD(ir, this->dtFunction, _IFT_TransientTransportProblem_dtFunction);
    } else if ( ir.hasField(_IFT_TransientTransportProblem_prescribedTimes) ) {
        IR_GIVE_FIELD(ir, this->prescribedTimes, _IFT_TransientTransportProblem_prescribedTimes);
    } else {
        IR_GIVE_FIELD(ir, this->deltaT, _IFT_TransientTransportProblem_deltaT);
    }

    this->keepTangent = ir.hasField(_IFT_TransientTransportProblem_keepTangent);

    this->lumped = ir.hasField(_IFT_TransientTransportProblem_lumped);

    field = std::make_unique<DofDistributedPrimaryField>(this, 1, FT_TransportProblemUnknowns, 2, this->alpha);

    // read field export flag
    exportFields.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, exportFields, _IFT_TransientTransportProblem_exportFields);
    if ( exportFields.giveSize() ) {
        FieldManager *fm = this->giveContext()->giveFieldManager();
        for ( int i = 1; i <= exportFields.giveSize(); i++ ) {
            if ( exportFields.at(i) == FT_Temperature ) {
                FieldPtr _temperatureField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->field.get(), {T_f} ) );
                fm->registerField( _temperatureField, ( FieldType ) exportFields.at(i) );
            } else if ( exportFields.at(i) == FT_HumidityConcentration ) {
                FieldPtr _concentrationField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->field.get(), {C_1} ) );
                fm->registerField( _concentrationField, ( FieldType ) exportFields.at(i) );
            }
        }
    }
    
//     InternalVariableField(IST_HydrationDegree, FT_Unknown, MMA_ClosestPoint, this->giveDomain(1));
}


double TransientTransportProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return this->field->giveUnknownValue(dof, mode, tStep);
}


double
TransientTransportProblem :: giveDeltaT(int n)
{
    if ( this->dtFunction ) {
        return this->giveDomain(1)->giveFunction(this->dtFunction)->evaluateAtTime(n);
    } else if ( this->prescribedTimes.giveSize() > 0 ) {
        return this->giveDiscreteTime(n) - this->giveDiscreteTime(n - 1);
    } else {
        return this->deltaT;
    }
}

double
TransientTransportProblem :: giveDiscreteTime(int iStep)
{
    if ( iStep > 0 && iStep <= this->prescribedTimes.giveSize() ) {
        return ( this->prescribedTimes.at(iStep) );
    } else if ( iStep == 0 ) {
        return initT;
    }

    OOFEM_ERROR("invalid iStep");
}


TimeStep *TransientTransportProblem :: giveNextStep()
{
    if ( !currentStep ) {
        // first step -> generate initial step
        currentStep = std::make_unique<TimeStep>( *giveSolutionStepWhenIcApply() );
    }

    double dt = this->giveDeltaT(currentStep->giveNumber()+1);
    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(*previousStep, dt);
    currentStep->setIntrinsicTime(previousStep->giveTargetTime() + alpha * dt);
    return currentStep.get();
}


TimeStep *TransientTransportProblem :: giveSolutionStepWhenIcApply(bool force)
{
    if ( master && (!force)) {
        return master->giveSolutionStepWhenIcApply();
    } else {
        if ( !stepWhenIcApply ) {
            double dt = this->giveDeltaT(1);
            stepWhenIcApply = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 0, this->initT, dt, 0);
            // The initial step goes from [-dt, 0], so the intrinsic time is at: -deltaT  + alpha*dt
            stepWhenIcApply->setIntrinsicTime(-dt + alpha * dt);
        }

        return stepWhenIcApply.get();
    }
}


void TransientTransportProblem :: solveYourselfAt(TimeStep *tStep)
{
    OOFEM_LOG_INFO( "\nSolving [step number %5d, time %e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
    
    Domain *d = this->giveDomain(1);
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    if ( tStep->isTheFirstStep() ) {
        this->applyIC();
    }

    field->advanceSolution(tStep);
    field->initialize(VM_Total, tStep, solution, EModelDefaultEquationNumbering());

    if ( !effectiveMatrix ) {
        effectiveMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        effectiveMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    }

    OOFEM_LOG_INFO("Assembling external forces\n");
    FloatArray externalForces(neq);
    externalForces.zero();
    this->assembleVector( externalForces, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), d );
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);

    // set-up numerical method
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    OOFEM_LOG_INFO("Solving for %d unknowns...\n", neq);

    internalForces.resize(neq);

    FloatArray incrementOfSolution;
    double loadLevel;
    int currentIterations;
    this->updateInternalRHS(this->internalForces, tStep, this->giveDomain(1), &this->eNorm); /// @todo Hack to ensure that internal RHS is evaluated before the tangent. This is not ideal, causing this to be evaluated twice for a linearproblem. We have to find a better way to handle this.
    this->nMethod->solve(*this->effectiveMatrix,
                         externalForces,
                         nullptr, // ignore
                         this->solution,
                         incrementOfSolution,
                         this->internalForces,
                         this->eNorm,
                         loadLevel, // ignore
                         SparseNonLinearSystemNM :: rlm_total, // ignore
                         currentIterations, // ignore
                         tStep);
}


void
TransientTransportProblem :: updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d)
{
    ///@todo NRSolver should report when the solution changes instead of doing it this way.
    this->field->update(VM_Total, tStep, solutionVector, EModelDefaultEquationNumbering());
    ///@todo Need to reset the boundary conditions properly since some "update" is doing strange
    /// things such as applying the (wrong) boundary conditions. This call will be removed when that code can be removed.
    this->field->applyBoundaryCondition(tStep);
}


void
TransientTransportProblem :: updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm)
{
    // F_eff = F(T^(k)) + C * dT/dt^(k)
    answer.zero();
    this->assembleVector(answer, tStep, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), d, eNorm);
    this->updateSharedDofManagers(answer, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
    if ( lumped ) {
        // Note, inertia contribution cannot be computed on element level when lumped mass matrices are used.
        FloatArray oldSolution, vel;
        this->field->initialize(VM_Total, tStep->givePreviousStep(), oldSolution, EModelDefaultEquationNumbering());
        vel.beDifferenceOf(solution, oldSolution);
        vel.times( 1./tStep->giveTimeIncrement() );
        FloatArray capacityDiag(vel.giveSize());
        this->assembleVector( capacityDiag, tStep, LumpedMassVectorAssembler(), VM_Total, EModelDefaultEquationNumbering(), d );
        for ( int i = 0; i < vel.giveSize(); ++i ) {
            answer[i] += capacityDiag[i] * vel[i];
        }
    } else {
        FloatArray tmp;
        this->assembleVector(answer, tStep, InertiaForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), d, & tmp);
        eNorm->add(tmp); ///@todo Fix this, assembleVector shouldn't zero eNorm inside the functions. / Mikael
    }
}


void
TransientTransportProblem :: updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d)
{
    // K_eff = (a*K + C/dt)
    if ( !this->keepTangent || !this->hasTangent ) {
        mat.zero();
        this->assemble(mat, tStep, EffectiveTangentAssembler(TangentStiffness, lumped, this->alpha, 1./tStep->giveTimeIncrement()),
                       EModelDefaultEquationNumbering(), d );
        this->hasTangent = true;
    }
}


void
TransientTransportProblem :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    // F(T) + C*dT/dt = Q, F(T)=(K_c+K_h)*T-R_q-R_h
    // Linearized:
    // F(T^(k)) + K*a*dT_1 = Q - C * dT/dt^(k) - C/dt * dT_1
    // Rearranged
    // (a*K + C/dt) * dT_1 = Q - (F(T^(k)) + C * dT/dt^(k))
    // K_eff        * dT_1 = Q - F_eff
    // Update:
    // T_1 += dT_1

    ///@todo NRSolver should report when the solution changes instead of doing it this way.
    this->field->update(VM_Total, tStep, solution, EModelDefaultEquationNumbering());
    ///@todo Need to reset the boundary conditions properly since some "update" is doing strange
    /// things such as applying the (wrong) boundary conditions. This call will be removed when that code can be removed.
    this->field->applyBoundaryCondition(tStep);

    if ( cmpn == InternalRhs ) {
        // F_eff = F(T^(k)) + C * dT/dt^(k)
        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, InternalForceAssembler(), VM_Total,
                             EModelDefaultEquationNumbering(), d, & this->eNorm);
        this->updateSharedDofManagers(this->internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
        if ( lumped ) {
            // Note, inertia contribution cannot be computed on element level when lumped mass matrices are used.
            FloatArray oldSolution, vel;
            this->field->initialize(VM_Total, tStep->givePreviousStep(), oldSolution, EModelDefaultEquationNumbering());
            vel.beDifferenceOf(solution, oldSolution);
            vel.times( 1./tStep->giveTimeIncrement() );
            FloatArray capacityDiag(vel.giveSize());
            this->assembleVector( capacityDiag, tStep, LumpedMassVectorAssembler(), VM_Total, EModelDefaultEquationNumbering(), d );
            for ( int i = 0; i < vel.giveSize(); ++i ) {
                this->internalForces[i] += capacityDiag[i] * vel[i];
            }
        } else {
            FloatArray tmp;
            this->assembleVector(this->internalForces, tStep, InertiaForceAssembler(), VM_Total,
                                EModelDefaultEquationNumbering(), d, & tmp);
            this->eNorm.add(tmp); ///@todo Fix this, assembleVector shouldn't zero eNorm inside the functions. / Mikael
        }

    } else if ( cmpn == NonLinearLhs ) {
        // K_eff = (a*K + C/dt)
        if ( !this->keepTangent || !this->hasTangent ) {
            this->effectiveMatrix->zero();
            this->assemble( *effectiveMatrix, tStep, EffectiveTangentAssembler(TangentStiffness, lumped, this->alpha, 1./tStep->giveTimeIncrement()),
                                                                               EModelDefaultEquationNumbering(), d );
            this->hasTangent = true;
        }
    } else {
        OOFEM_ERROR("Unknown component");
    }
}


void
TransientTransportProblem :: applyIC()
{
    Domain *domain = this->giveDomain(1);
    OOFEM_LOG_INFO("Applying initial conditions\n");

    this->field->applyDefaultInitialCondition();

    // set initial field IP values (needed by some nonlinear materials)
    TimeStep *s = this->giveSolutionStepWhenIcApply();
    for ( auto &elem : domain->giveElements() ) {
        TransportElement *element = static_cast< TransportElement * >( elem.get() );
        element->updateInternalState(s);
        element->updateYourself(s);
    }
}


bool
TransientTransportProblem :: requiresEquationRenumbering(TimeStep *tStep)
{
    ///@todo This method should be set as the default behavior instead of relying on a user specified flag. Then this function should be removed.
    if ( tStep->isTheFirstStep() ) {
        return true;
    }
    // Check if Dirichlet b.c.s has changed.
    Domain *d = this->giveDomain(1);
    for ( auto &gbc : d->giveBcs() ) {
        ActiveBoundaryCondition *active_bc = dynamic_cast< ActiveBoundaryCondition * >(gbc.get());
        BoundaryCondition *bc = dynamic_cast< BoundaryCondition * >(gbc.get());
        // We only need to consider Dirichlet b.c.s
        if ( bc || ( active_bc && ( active_bc->requiresActiveDofs() || active_bc->giveNumberOfInternalDofManagers() ) ) ) {
            // Check of the dirichlet b.c. has changed in the last step (if so we need to renumber)
            if ( gbc->isImposed(tStep) != gbc->isImposed(tStep->givePreviousStep()) ) {
                return true;
            }
        }
    }
    return false;
}

int
TransientTransportProblem :: forceEquationNumbering()
{
    this->effectiveMatrix = nullptr;
    return EngngModel :: forceEquationNumbering();
}


void
TransientTransportProblem :: printOutputAt(FILE *file, TimeStep *tStep)
{
    if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return; // Do not print even Solution step header
    }

    EngngModel :: printOutputAt(file, tStep);
    
    // computes and prints reaction forces in all supported or restrained dofs
    fprintf(file, "\n\n\tR E A C T I O N S  O U T P U T:\n\t_______________________________\n\n\n");
    
    IntArray restrDofMans, restrDofs, eqn;
    int di=1;
    //from StructuralEngngModel :: buildReactionTable(dofManMap, dofidMap, eqnMap, tStep, 1);
    // determine number of restrained dofs
    Domain *domain = this->giveDomain(di);
    int numRestrDofs = this->giveNumberOfDomainEquations( di, EModelDefaultPrescribedEquationNumbering() );
    int ndofMan = domain->giveNumberOfDofManagers();
    int rindex, count = 0;

    // initialize corresponding dofManagers and dofs for each restrained dof
    restrDofMans.resize(numRestrDofs);
    restrDofs.resize(numRestrDofs);
    eqn.resize(numRestrDofs);

    for ( int i = 1; i <= ndofMan; i++ ) {
        DofManager *inode = domain->giveDofManager(i);
        for ( Dof *jdof: *inode ) {
            if ( jdof->isPrimaryDof() && ( jdof->hasBc(tStep) ) ) { // skip slave dofs
                rindex = jdof->__givePrescribedEquationNumber();
                if ( rindex ) {
                    count++;
                    restrDofMans.at(count) = i;
                    restrDofs.at(count) = jdof->giveDofID();
                    eqn.at(count) = rindex;
                } else {
                    // NullDof has no equation number and no prescribed equation number
                    //_error("No prescribed equation number assigned to supported DOF");
                }
            }
        }
    }
    // Trim to size.
    restrDofMans.resizeWithValues(count);
    restrDofs.resizeWithValues(count);
    eqn.resizeWithValues(count);
    
    //restrDofMans.printYourself();
    //restrDofs.printYourself();
    //eqn.printYourself();
    
    //from StructuralEngngModel :: computeReaction()
    FloatArray internal, external;
    internal.resize( this->giveNumberOfDomainEquations( di, EModelDefaultPrescribedEquationNumbering() ) );
    internal.zero();
    
    // Add internal forces
    this->assembleVector( internal, tStep, InternalForceAssembler(), VM_Total,
                         EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di), &this->eNorm );
    // Subtract external loading
    external.resize(internal.giveSize());
    this->assembleVector( external, tStep, ExternalForceAssembler(), VM_Total,                   EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di), &this->eNorm );
    
    internal.subtract(external);
    //internal.printYourself();
    
    //
    // loop over reactions and print them
    //
    for ( int i = 1; i <= restrDofMans.giveSize(); i++ ) {
        if ( domain->giveOutputManager()->testDofManOutput(restrDofMans.at(i), tStep) ) {
            fprintf(file, "\tNode %8d iDof %2d reaction % .4e    [bc-id: %d]\n",
                    domain->giveDofManager( restrDofMans.at(i) )->giveLabel(),
                    restrDofs.at(i), internal.at( eqn.at(i) ),
                    domain->giveDofManager( restrDofMans.at(i) )->giveDofWithID( restrDofs.at(i) )->giveBcId() );
        }
    }
    
}


void
TransientTransportProblem :: updateYourself(TimeStep *tStep)
{
    EngngModel :: updateYourself(tStep);
}

void
TransientTransportProblem :: saveContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: saveContext(stream, mode);
    field->saveContext(stream);
}


void
TransientTransportProblem :: restoreContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: restoreContext(stream, mode);
    field->restoreContext(stream);
}


int
TransientTransportProblem :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % 2;
}


int
TransientTransportProblem :: requiresUnknownsDictionaryUpdate()
{
    return true;
}

int
TransientTransportProblem :: checkConsistency()
{
    // check for proper element type
    for ( auto &elem : this->giveDomain(1)->giveElements() ) {
        if ( !dynamic_cast< TransportElement * >( elem.get() ) ) {
            OOFEM_WARNING("Element %d has no TransportElement base", elem->giveLabel());
            return 0;
        }
    }

    return EngngModel :: checkConsistency();
}


void
TransientTransportProblem :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}


FieldPtr TransientTransportProblem::giveField(FieldType key, TimeStep *tStep)
{
    /* Note: the current implementation uses MaskedPrimaryField, that is automatically updated with the model progress, 
        so the returned field always refers to active solution step. 
    */

    if ( tStep != this->giveCurrentStep()) {
        OOFEM_ERROR("Unable to return field representation for non-current time step");
    }
    if ( key == FT_Temperature ) {
        return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{T_f} );
    } else if ( key == FT_HumidityConcentration ) {
        return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{C_1} );
    } else if ( key == FT_VOF ) {
        return this->giveContext()->giveFieldManager()->giveField(key);
    } else {
        return FieldPtr();
    }
}


} // end namespace oofem
