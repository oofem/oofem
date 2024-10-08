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

#include "stokesflow.h"
#include "fmelement.h"
#include "inputrecord.h"
#include "timestep.h"
#include "classfactory.h"
#include "domain.h"
#include "nrsolver.h"
#include "sparsenonlinsystemnm.h"
#include "meshqualityerrorestimator.h"
#include "topologydescription.h"
#include "parallelcontext.h"
#include "exportmodulemanager.h"
#include "dofdistributedprimaryfield.h"
#include "unknownnumberingscheme.h"

namespace oofem {
REGISTER_EngngModel(StokesFlow);

StokesFlow :: StokesFlow(int i, EngngModel *_master) : FluidModel(i, _master)
{
    this->ndomains = 1;
}


StokesFlow :: ~StokesFlow() { }


void StokesFlow :: initializeFrom(InputRecord &ir)
{
    FluidModel :: initializeFrom(ir);
    int val;

    val = ( int ) SMT_PetscMtrx;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    this->sparseMtrxType = ( SparseMtrxType ) val;

    val = ( int ) ST_Petsc;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    this->solverType = ( LinSystSolverType ) val;

    this->deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, _IFT_StokesFlow_deltat);

    this->velocityPressureField = std::make_unique<DofDistributedPrimaryField>(this, 1, FT_VelocityPressure, 1);
    this->stiffnessMatrix = nullptr;
    this->meshqualityee = nullptr;

    this->ts = TS_OK;

    this->maxdef = 25; ///@todo Deal with this parameter (set to some reasonable value by default now)
}


void StokesFlow :: solveYourselfAt(TimeStep *tStep)
{
    Domain *d = this->giveDomain(1);
    FloatArray externalForces;
    FloatArray incrementOfSolution;

    if ( d->giveNumberOfElements() == 0 && d->giveTopology() ) {
        d->giveTopology()->replaceFEMesh();
        this->meshqualityee = std::make_unique<MeshQualityErrorEstimator>( 1, d );
    }

    if ( d->giveTopology() && this->meshqualityee ) {
        // Check the quality of the deformed mesh.
        double meshdeformation = this->meshqualityee->giveValue(globalErrorEEV, tStep);
        if ( this->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO("StokesFlow :: solveYourselfAt - Mesh deformation at %e\n", meshdeformation);
        }

        if ( this->ts == TS_NeedsRemeshing || meshdeformation > this->maxdef ) {
            d->giveTopology()->replaceFEMesh();
            OOFEM_LOG_INFO( "StokesFlow :: solveYourselfAt - New mesh created (%d elements).\n", d->giveNumberOfElements() );
            /*meshdeformation =*/ this->meshqualityee->giveValue(globalErrorEEV, tStep);
            this->giveExportModuleManager()->initialize();
        }
    }

    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    velocityPressureField->advanceSolution(tStep);
    velocityPressureField->initialize(VM_Total, tStep, solutionVector, EModelDefaultEquationNumbering() );

    // Create "stiffness matrix"
    if ( !this->stiffnessMatrix ) {
        this->stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        if ( !this->stiffnessMatrix ) {
            OOFEM_ERROR("Couldn't create requested sparse matrix of type %d", sparseMtrxType);
        }

        this->stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    }

    incrementOfSolution.resize(neq);
    this->internalForces.resize(neq);

    // Build initial/external load (LoadVector)
    externalForces.resize(neq);
    externalForces.zero();
    this->assembleVector( externalForces, tStep, ExternalForceAssembler(), VM_Total,
                         EModelDefaultEquationNumbering(), d );
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);

    if ( this->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("StokesFlow :: solveYourselfAt - Solving step %d, metastep %d, (neq = %d)\n", tStep->giveNumber(), tStep->giveMetaStepNumber(), neq);
    }

    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    this->initMetaStepAttributes( this->giveCurrentMetaStep() );
#if 1
    double loadLevel;
    int currentIterations=0;
    this->updateInternalRHS( this->internalForces, tStep, d, &this->eNorm );
    ConvergedReason status = this->nMethod->solve(*this->stiffnessMatrix,
                                            externalForces,
                                            NULL,
                                            solutionVector,
                                            incrementOfSolution,
                                            this->internalForces,
                                            this->eNorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total,
                                            currentIterations,
                                            tStep);
#else
    this->updateInternalRHS( this->internalForces, tStep, d, nullptr );
    this->updateMatrix( this->stiffnessMatrix, tStep, d );
    this->internalForces.negated();
    this->internalForces.add(externalForces);
    ConvergedReason status = this->nMethod->giveLinearSolver()->solve(this->stiffnessMatrix.get(), & ( this->internalForces ), & solutionVector);
    this->updateSolution(solutionVector);
    this->updateInternalRHS( this->internalForces, tStep, d, nullptr );
#endif

    if (status != CR_CONVERGED ) {
        OOFEM_ERROR("No success in solving problem at time step", tStep->giveNumber());
    }
    tStep->numberOfIterations = currentIterations;
    tStep->convergedReason = status;
}

void StokesFlow :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    velocityPressureField->update(VM_Total, tStep, solutionVector, EModelDefaultEquationNumbering());

    // update element stabilization
    for ( auto &elem : d->giveElements() ) {
        static_cast< FMElement * >( elem.get() )->updateStabilizationCoeffs(tStep);
    }

    if ( cmpn == InternalRhs ) {
        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, InternalForceAssembler(), VM_Total,
                             EModelDefaultEquationNumbering(), d, & this->eNorm);
        this->updateSharedDofManagers(this->internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
        return;
    } else if ( cmpn == NonLinearLhs ) {
        this->stiffnessMatrix->zero();
        this->assemble(*stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                       EModelDefaultEquationNumbering(), d);
        return;
    } else {
        OOFEM_ERROR("Unknown component");
    }
}

void StokesFlow :: updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d)
{
    this->velocityPressureField->update(VM_Total, tStep, solutionVector, EModelDefaultEquationNumbering());
    // update element stabilization
    for ( auto &elem : d->giveElements() ) {
        static_cast< FMElement * >( elem.get() )->updateStabilizationCoeffs(tStep);
    }
}


void StokesFlow :: updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *enorm)
{
    answer.zero();
    this->assembleVector(answer, tStep, InternalForceAssembler(), VM_Total,
                         EModelDefaultEquationNumbering(), d, enorm);
    this->updateSharedDofManagers(answer, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
}


void StokesFlow :: updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d)
{
    mat.zero();
    this->assemble(mat, tStep, TangentAssembler(TangentStiffness),
                   EModelDefaultEquationNumbering(), d);
}


void StokesFlow :: updateYourself(TimeStep *tStep)
{
    this->updateInternalState(tStep);
    FluidModel :: updateYourself(tStep);
}


int StokesFlow :: forceEquationNumbering(int id)
{
    int neq = FluidModel :: forceEquationNumbering(id);
    this->equationNumberingCompleted = false;
    this->stiffnessMatrix = nullptr;
    return neq;
}


double StokesFlow :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return velocityPressureField->giveUnknownValue(dof, mode, tStep);
}


double StokesFlow :: giveReynoldsNumber()
{
    return 1.0;
}


int StokesFlow :: checkConsistency()
{
    Domain *domain = this->giveDomain(1);

    // check for proper element type
    for ( auto &elem : domain->giveElements() ) {
        if ( dynamic_cast< FMElement * >( elem.get() ) == NULL ) {
            OOFEM_WARNING("Element %d has no FMElement base", elem->giveLabel());
            return false;
        }
    }

    return EngngModel :: checkConsistency();
}

void StokesFlow :: updateInternalState(TimeStep *tStep)
{
    for ( auto &domain: domainList ) {
        if ( domain->giveTopology() ) {
            // Must be done before updating nodal positions
            this->ts = domain->giveTopology()->updateYourself(tStep);
        }

        for ( auto &elem : domain->giveElements() ) {
            elem->updateInternalState(tStep);
        }
    }
}

void StokesFlow :: doStepOutput(TimeStep *tStep)
{
    TopologyDescription *tp = this->giveDomain(1)->giveTopology();
    if ( tp ) {
        tp->doOutput(tStep);
    }

    FluidModel :: doStepOutput(tStep);
}

NumericalMethod *StokesFlow :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !this->nMethod ) {
        this->nMethod = std::make_unique<NRSolver>(this->giveDomain(1), this);
    }
    return this->nMethod.get();
}

TimeStep *StokesFlow :: giveNextStep()
{
    if ( !currentStep ) {
        // first step -> generate initial step
        //currentStep = std::make_unique<TimeStep>(*giveSolutionStepWhenIcApply());
        currentStep = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 1, 0., this->deltaT, 0);
    }
    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(*previousStep, this->deltaT);

    return currentStep.get();
}
} // end namespace oofem
