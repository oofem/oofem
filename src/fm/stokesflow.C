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
#include "primaryfield.h"

namespace oofem {
REGISTER_EngngModel(StokesFlow);

StokesFlow :: StokesFlow(int i, EngngModel *_master) : FluidModel(i, _master)
{
    this->nMethod = NULL;
    this->ndomains = 1;
    this->stiffnessMatrix = NULL;
    this->meshqualityee = NULL;
    this->velocityPressureField = NULL;
#ifdef __PARALLEL_MODE
    commMode = ProblemCommMode__NODE_CUT;
#endif
}

StokesFlow :: ~StokesFlow()
{
    delete this->velocityPressureField;
    delete this->nMethod;
    delete this->stiffnessMatrix;
    delete this->meshqualityee;
}

IRResultType StokesFlow :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    int val;

    val = ( int ) SMT_PetscMtrx;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    this->sparseMtrxType = ( SparseMtrxType ) val;

    val = ( int ) ST_Petsc;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    this->solverType = ( LinSystSolverType ) val;

    this->deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, _IFT_StokesFlow_deltat);

    delete this->velocityPressureField;
    this->velocityPressureField = new PrimaryField(this, 1, FT_VelocityPressure, EID_MomentumBalance_ConservationEquation, 1);
    delete this->stiffnessMatrix;
    this->stiffnessMatrix = NULL;
    delete this->meshqualityee;
    this->meshqualityee = NULL;

    this->ts = TS_OK;

    this->maxdef = 25; ///@todo Deal with this parameter (set to some reasonable value by default now)

    return EngngModel :: initializeFrom(ir);
}


void StokesFlow :: solveYourselfAt(TimeStep *tStep)
{
    FloatArray *solutionVector = NULL;

    if ( this->giveDomain(1)->giveNumberOfElements() == 0 && this->giveDomain(1)->giveTopology() ) {
        this->giveDomain(1)->giveTopology()->replaceFEMesh();
        this->meshqualityee = new MeshQualityErrorEstimator( 1, this->giveDomain(1) );
    }

    if ( this->giveDomain(1)->giveTopology() && this->meshqualityee ) {
        // Check the quality of the deformed mesh.
        double meshdeformation = this->meshqualityee->giveValue(globalErrorEEV, tStep);
        if ( this->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO("StokesFlow :: solveYourselfAt - Mesh deformation at %e\n", meshdeformation);
        }

        if ( this->ts == TS_NeedsRemeshing || meshdeformation > this->maxdef ) {
            this->giveDomain(1)->giveTopology()->replaceFEMesh();
            OOFEM_LOG_INFO( "StokesFlow :: solveYourselfAt - New mesh created (%d elements).\n", this->giveDomain(1)->giveNumberOfElements() );
            /*meshdeformation =*/ this->meshqualityee->giveValue(globalErrorEEV, tStep);
            this->giveExportModuleManager()->initialize();
        }
    }

    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    // Move solution space to current time step
    velocityPressureField->advanceSolution(tStep);

    // Point pointer SolutionVector to current solution in velocityPressureField
    solutionVector = velocityPressureField->giveSolutionVector(tStep);
    solutionVector->resize(neq);
    solutionVector->zero();

    // Create "stiffness matrix"
    if ( !this->stiffnessMatrix ) {
        this->stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        if ( !this->stiffnessMatrix ) {
            OOFEM_ERROR("Couldn't create requested sparse matrix of type %d", sparseMtrxType);
        }

        this->stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance_ConservationEquation, EModelDefaultEquationNumbering() );
    }

    this->incrementOfSolution.resize(neq);
    this->internalForces.resize(neq);

    // Build initial/external load (LoadVector)
    this->externalForces.resize(neq);
    this->externalForces.zero();
    this->assembleVector( this->externalForces, tStep, EID_MomentumBalance_ConservationEquation, ExternalForcesVector, VM_Total,
                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers(this->externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);
#endif

    if ( this->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("StokesFlow :: solveYourselfAt - Solving step %d, metastep %d, (neq = %d)\n", tStep->giveNumber(), tStep->giveMetaStepNumber(), neq);
    }

    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    this->initMetaStepAttributes( this->giveCurrentMetaStep() );
#if 1
    double loadLevel;
    int currentIterations;
    this->updateComponent( tStep, InternalRhs, this->giveDomain(1) );
    NM_Status status = this->nMethod->solve(this->stiffnessMatrix,
                                            & ( this->externalForces ),
                                            NULL,
                                            solutionVector,
                                            & ( this->incrementOfSolution ),
                                            & ( this->internalForces ),
                                            this->eNorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total,
                                            currentIterations,
                                            tStep);
#else
    this->updateComponent( tStep, InternalRhs, this->giveDomain(1) );
    this->updateComponent( tStep, NonLinearLhs, this->giveDomain(1) );
    this->internalForces.negated();
    this->internalForces.add(externalForces);
    NM_Status status = this->nMethod->giveLinearSolver()->solve(this->stiffnessMatrix, & ( this->internalForces ), solutionVector);
    this->updateComponent( tStep, NonLinearLhs, this->giveDomain(1) );
#endif

    if ( !( status & NM_Success ) ) {
        OOFEM_ERROR("No success in solving problem at time step", tStep->giveNumber());
    }


    // update element stabilization
    Domain *d = this->giveDomain(1);
    int nelem = d->giveNumberOfElements();
    for ( int i = 1; i <= nelem; ++i ) {
        static_cast< FMElement * >( d->giveElement(i) )->updateStabilizationCoeffs(tStep);
    }
}

void StokesFlow :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    // update element stabilization
    int nelem = d->giveNumberOfElements();
    for ( int i = 1; i <= nelem; ++i ) {
        static_cast< FMElement * >( d->giveElement(i) )->updateStabilizationCoeffs(tStep);
    }

    if ( cmpn == InternalRhs ) {
        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, EID_MomentumBalance_ConservationEquation, InternalForcesVector, VM_Total,
                             EModelDefaultEquationNumbering(), this->giveDomain(1), & this->eNorm);
#ifdef __PARALLEL_MODE
        this->updateSharedDofManagers(this->internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
#endif
        return;
    } else if ( cmpn == NonLinearLhs ) {
        this->stiffnessMatrix->zero();
        this->assemble(this->stiffnessMatrix, tStep, EID_MomentumBalance_ConservationEquation, StiffnessMatrix,
                       EModelDefaultEquationNumbering(), d);
        return;
    } else {
        OOFEM_ERROR("Unknown component");
    }
}

void StokesFlow :: updateYourself(TimeStep *tStep)
{
    this->updateInternalState(tStep);
    EngngModel :: updateYourself(tStep);
}

int StokesFlow :: forceEquationNumbering(int id)
{

  /*
    int neq = FluidModel :: forceEquationNumbering(id);
  */
    int neq = EngngModel::forceEquationNumbering(id);
    this->equationNumberingCompleted = false;
    if ( this->stiffnessMatrix ) {
        delete this->stiffnessMatrix;
        this->stiffnessMatrix = NULL;
    }

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


#ifdef __PARALLEL_MODE
void StokesFlow :: initParallelContexts()
{
    ParallelContext *parallelContext;
    parallelContextList->growTo(ndomains);
    for ( int i = 1; i <= this->ndomains; i++ ) {
        parallelContext =  new ParallelContext(this);
        parallelContextList->put(i, parallelContext);
    }
}
#endif

int StokesFlow :: checkConsistency()
{
    int nelem;
    FMElement *sePtr;
    Domain *domain = this->giveDomain(1);
    nelem = domain->giveNumberOfElements();

    // check for proper element type
    for ( int i = 1; i <= nelem; i++ ) {
        sePtr = dynamic_cast< FMElement * >( domain->giveElement(i) );
        if ( sePtr == NULL ) {
            OOFEM_WARNING("Element %d has no FMElement base", i);
            return false;
        }
    }

    return EngngModel :: checkConsistency();
}

void StokesFlow :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    DofIDItem type = iDof->giveDofID();
    ///@todo This won't work with slave dofs, xfem etc.
    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
        iDof->printSingleOutputAt(stream, tStep, 'v', VM_Total, 1);
    } else if ( type == P_f ) {
        iDof->printSingleOutputAt(stream, tStep, 'p', VM_Total, 1);
    } else {
        OOFEM_ERROR("unsupported dof type");
    }
}

void StokesFlow :: updateInternalState(TimeStep *tStep)
{
    Domain *domain;
    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        if ( domain->giveTopology() ) {
            // Must be done before updating nodal positions
            this->ts = domain->giveTopology()->updateYourself(tStep);
        }

        int nelem = domain->giveNumberOfElements();
        for ( int j = 1; j <= nelem; j++ ) {
            domain->giveElement(j)->updateInternalState(tStep);
        }
    }
}

void StokesFlow :: doStepOutput(TimeStep *tStep)
{
    TopologyDescription *tp = this->giveDomain(1)->giveTopology();
    if ( tp ) {
        tp->doOutput(tStep);
    }

    EngngModel :: doStepOutput(tStep);
}

NumericalMethod *StokesFlow :: giveNumericalMethod(MetaStep *mStep)
{
    if ( this->nMethod ) {
        return this->nMethod;
    }

    this->nMethod = new NRSolver(this->giveDomain(1), this);
    return this->nMethod;
}

TimeStep *StokesFlow :: giveNextStep()
{
    if ( previousStep ) {
        delete previousStep;
    }

    if ( currentStep == NULL ) {
        int istep = this->giveNumberOfFirstStep();
        // first step -> generate initial step
        previousStep = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, -this->deltaT, this->deltaT, 0);
        currentStep = new TimeStep(istep, this, 1, 0.0, this->deltaT, 1);
    } else {
        int istep =  currentStep->giveNumber() + 1;
        StateCounterType counter = currentStep->giveSolutionStateCounter() + 1;
        previousStep = currentStep;
        double dt = currentStep->giveTimeIncrement();
        double totalTime = currentStep->giveTargetTime() + dt;
        currentStep = new TimeStep(istep, this, 1, totalTime, dt, counter);
    }

    return currentStep;
}
} // end namespace oofem
