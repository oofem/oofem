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

#include "stokesflow.h"
#include "fmelement.h"
#include "inputrecord.h"
#include "timestep.h"
#include "usrdefsub.h"
#include "domain.h"
#include "nrsolver.h"
#include "sparsenonlinsystemnm.h"
#include "meshqualityerrorestimator.h"
#include "topologydescription.h"
#include "petsccontext.h"
#include "exportmodulemanager.h"
#include "primaryfield.h"

namespace oofem {
StokesFlow :: StokesFlow(int i, EngngModel *_master) : EngngModel(i, _master)
{
    this->nMethod = NULL;
    this->ndomains = 1;
    this->hasAdvanced = false;
    this->stiffnessMatrix = NULL;
    this->meshqualityee = NULL;
#ifdef __PARALLEL_MODE
    commMode = ProblemCommMode__NODE_CUT;
#endif
}

StokesFlow :: ~StokesFlow()
{
    delete this->velocityPressureField;
    delete this->nMethod;

    if ( this->stiffnessMatrix ) {
        delete this->stiffnessMatrix;
    }
    if ( this->meshqualityee ) {
        delete this->meshqualityee;
    }
}

IRResultType StokesFlow :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;
    int val;

    val = ( int ) SMT_PetscMtrx;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_StokesFlow_smtype, "smtype");
    this->sparseMtrxType = ( SparseMtrxType ) val;

    val = ( int ) ST_Petsc;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_SUPG_lstype, "lstype");
    this->solverType = ( LinSystSolverType ) val;

    this->deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, IFT_StokesFlow_deltat, "deltat");

    this->velocityPressureField = new PrimaryField(this, 1, FT_VelocityPressure, EID_MomentumBalance_ConservationEquation, 1);

    this->ts = TS_OK;

    this->maxdef = 25; ///@todo Deal with this parameter (set to some reasonable value by default now)

    return EngngModel :: initializeFrom(ir);
}

void StokesFlow :: solveYourselfAt(TimeStep *tStep)
{
    FloatArray *solutionVector = NULL;

    if ( this->giveDomain(1)->giveNumberOfElements() == 0 && this->giveDomain(1)->giveTopology() ) {
        this->giveDomain(1)->giveTopology()->replaceFEMesh();
        this->meshqualityee = new MeshQualityErrorEstimator(1, this->giveDomain(1));
    }

    if (this->giveDomain(1)->giveTopology() && this->meshqualityee) {
        // Check the quality of the deformed mesh.
        double meshdeformation = this->meshqualityee->giveValue(globalErrorEEV, tStep);
        OOFEM_LOG_INFO("StokesFlow :: solveYourselfAt - Mesh deformation at %e\n",meshdeformation);
        if (this->ts == TS_NeedsRemeshing || meshdeformation > this->maxdef) {
            this->giveDomain(1)->giveTopology()->replaceFEMesh();
            OOFEM_LOG_INFO("StokesFlow :: updateYourself - New mesh created (%d elements).\n",this->giveDomain(1)->giveNumberOfElements());
            meshdeformation = this->meshqualityee->giveValue(globalErrorEEV, tStep);
            this->giveExportModuleManager()->initialize();
        }
    }

    int neq = this->giveNumberOfEquations(EID_MomentumBalance_ConservationEquation);

    // Move solution space to current time step
    if ( !hasAdvanced ) {
        velocityPressureField->advanceSolution(tStep);
        hasAdvanced = true;
    }

    // Point pointer SolutionVector to current solution in velocityPressureField
    solutionVector = velocityPressureField->giveSolutionVector(tStep);
    solutionVector->resize(neq);
    solutionVector->zero();

    // Create "stiffness matrix"
    if ( !this->stiffnessMatrix ) {
        this->stiffnessMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
        if ( !this->stiffnessMatrix ) {
            OOFEM_ERROR2("StokesFlow :: solveYourselfAt - Couldn't create requested sparse matrix of type %d", sparseMtrxType);
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

    OOFEM_LOG_INFO("StokesFlow :: solveYourselfAt - Solving (neq = %d)\n", neq);

    this->giveNumericalMethod( this->giveCurrentMetaStep() );
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
                                            eNorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total, // Why this naming scheme? Should be RLM_Total, and ReferenceLoadInputModeType
                                            currentIterations,
                                            tStep);
#else
    this->updateComponent(tStep, InternalRhs, this->giveDomain(1));
    this->updateComponent(tStep, NonLinearLhs, this->giveDomain(1));
    this->internalForces.negated();
    this->internalForces.add(externalForces);
    NM_Status status = this->nMethod->giveLinearSolver()->solve(this->stiffnessMatrix, &(this->internalForces), solutionVector);
    this->updateComponent(tStep, NonLinearLhs, this->giveDomain(1));
#endif

    if ( !(status & NM_Success) ) {
        OOFEM_ERROR2( "No success in solving problem at time step", tStep->giveNumber() );
    }


    // update element stabilization
    Domain* d = this->giveDomain(1);
    int i, nelem = d->giveNumberOfElements();
    for ( i = 1; i <= nelem; ++i ) {
        ((FMElement*)d->giveElement(i))->updateStabilizationCoeffs(tStep);
    }
}

void StokesFlow :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    // update element stabilization
    int i, nelem = d->giveNumberOfElements();
    for ( i = 1; i <= nelem; ++i ) {
        ((FMElement*)d->giveElement(i))->updateStabilizationCoeffs(tStep);
    }

    if (cmpn == InternalRhs) {
        this->internalForces.zero();
        eNorm = this->assembleVector( this->internalForces, tStep, EID_MomentumBalance_ConservationEquation, InternalForcesVector, VM_Total,
                              EModelDefaultEquationNumbering(), this->giveDomain(1) );
        return;

    } else if (cmpn == NonLinearLhs) {
        this->stiffnessMatrix->zero();
        this->assemble(this->stiffnessMatrix, tStep, EID_MomentumBalance_ConservationEquation, StiffnessMatrix,
                EModelDefaultEquationNumbering(), d);
        return;

    } else {
        OOFEM_ERROR("StokesFlow::updateComponent - Unknown component");
    }
}

void StokesFlow :: updateYourself(TimeStep *tStep)
{
    hasAdvanced = false;
    this->updateInternalState(tStep);
    EngngModel :: updateYourself(tStep);
}

int StokesFlow :: forceEquationNumbering(int id)
{
    int neq = EngngModel :: forceEquationNumbering(id);

    this->equationNumberingCompleted = false;
    if ( this->stiffnessMatrix ) {
        delete this->stiffnessMatrix;
        this->stiffnessMatrix = NULL;
    }

    return neq;
}


double StokesFlow :: giveUnknownComponent(EquationID chc, ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    if ( ( chc == EID_ConservationEquation ) || ( chc == EID_MomentumBalance ) || ( chc == EID_MomentumBalance_ConservationEquation ) ) {
        return velocityPressureField->giveUnknownValue(dof, mode, tStep);
    } else {
        OOFEM_ERROR("giveUnknownComponent: Unknown is of undefined equation id for this problem");
    }

    return 0;
}

double StokesFlow::giveUnknownComponent(UnknownType ut, ValueModeType vmt, TimeStep *atTime, Domain *d, Dof *dof)
{
    if (ut==ReynoldsNumber)
        return 1.0;
    else
        return 0.0;
} // bp


#ifdef __PETSC_MODULE
void StokesFlow :: initPetscContexts()
{
    PetscContext *petscContext;
    petscContextList->growTo(ndomains);
    for ( int i = 1; i <= this->ndomains; i++ ) {
        petscContext =  new PetscContext(this, EID_MomentumBalance_ConservationEquation);
        petscContextList->put(i, petscContext);
    }
}
#endif

int StokesFlow :: checkConsistency()
{
    int i, nelem;
    FMElement *sePtr;
    Domain *domain = this->giveDomain(1);
    nelem = domain->giveNumberOfElements();

    // check for proper element type
    for ( i = 1; i <= nelem; i++ ) {
        sePtr = dynamic_cast< FMElement * >( domain->giveElement(i) );
        if ( sePtr == NULL ) {
            OOFEM_WARNING2("Element %d has no FMElement base", i);
            return false;
        }
    }

    return EngngModel :: checkConsistency();
}

void StokesFlow :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    DofIDItem type = iDof->giveDofID();
    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
        iDof->printSingleOutputAt(stream, atTime, 'v', EID_MomentumBalance, VM_Total, 1);
    } else if ( type == P_f ) {
        iDof->printSingleOutputAt(stream, atTime, 'p', EID_ConservationEquation, VM_Total, 1);
    } else {
        OOFEM_ERROR("printDofOutputAt: unsupported dof type");
    }
}

void StokesFlow :: updateInternalState(TimeStep *tStep)
{
    Domain *domain;
    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        if (domain->giveTopology()) {
            // Must be done before updating nodal positions
            this->ts = domain->giveTopology()->updateYourself(tStep);
        }
        int nelem = domain->giveNumberOfElements();
        for (int j = 1; j <= nelem; j++ ) {
            domain->giveElement(j)->updateInternalState(tStep);
        }
    }
}

void StokesFlow :: doStepOutput(TimeStep *tStep)
{
    TopologyDescription *tp = this->giveDomain(1)->giveTopology();
    if (tp) tp->doOutput(tStep);
    EngngModel :: doStepOutput(tStep);
}

NumericalMethod *StokesFlow :: giveNumericalMethod(MetaStep *mStep)
{
    if ( this->nMethod ) {
        return this->nMethod;
    }

    this->nMethod = new NRSolver(1, this->giveDomain(1), this, EID_MomentumBalance_ConservationEquation);
    if ( !nMethod ) {
        OOFEM_ERROR("giveNumericalMethod: numerical method creation failed");
    }
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
