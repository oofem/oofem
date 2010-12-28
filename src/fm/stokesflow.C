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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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
#include "masterdof.h"
#include "domain.h"
#include "nrsolver.h"

namespace oofem {
StokesFlow :: StokesFlow(int i, EngngModel *_master) : EngngModel(i, _master)
{
    this->nMethod = NULL;
    this->ndomains = 1;
    this->hasAdvanced = false;
    this->stiffnessMatrix = NULL;
}

StokesFlow :: ~StokesFlow()
{
    delete this->velocityPressureField;
    delete this->nMethod;

    if ( this->stiffnessMatrix ) {
        delete this->stiffnessMatrix;
    }
}

IRResultType StokesFlow :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;
    int val;

    EngngModel :: initializeFrom(ir);

    val = ( int ) SMT_PetscMtrx;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_StokesFlow_smtype, "smtype"); // Macro
    sparseMtrxType = ( SparseMtrxType ) val;

    this->deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, IFT_StokesFlow_deltat, "deltat"); // Macro

    velocityPressureField = new PrimaryField(this, 1, FT_VelocityPressure, EID_MomentumBalance_ConservationEquation, 1);

    return IRRT_OK;
}

void StokesFlow :: solveYourselfAt(TimeStep *tStep)
{
    FloatArray *solutionVector = NULL;

    int neq = this->giveNumberOfEquations(EID_MomentumBalance_ConservationEquation);

    // Move solution space to current timestep
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
        this->stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance_ConservationEquation, EModelDefaultEquationNumbering() );
    }

    // Build initial/external load (LoadVector)
    this->externalForces.resize(neq);
    this->externalForces.zero();
    this->assembleVectorFromElements( this->externalForces, tStep, EID_MomentumBalance_ConservationEquation, LoadVector, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );

    this->incrementOfSolution.resize(neq);
    this->internalForces.resize(neq);

#if 0
    this->giveNumericalMethod(tStep);
    double loadLevel, ebenorm;
    int currentIterations;
    this->updateComponent( tStep, InternalRhs, this->giveDomain(1) );
    this->updateComponent( tStep, NonLinearLhs, this->giveDomain(1) );
    FloatArray incrementalLoadVector(0); // Should be allowed to be null
    NM_Status status = this->nMethod->solve(this->stiffnessMatrix,
                                            & ( this->externalForces ),
                                            NULL,
                                            & incrementalLoadVector, //Can't be null?! I don't have it, yet a need to supply it.
                                            solutionVector,
                                            & ( this->incrementOfSolution ),
                                            & ( this->internalForces ),
                                            ebenorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total, // Why this naming scheme? Should be RLM_Total, and ReferenceLoadInputModeType
                                            currentIterations,
                                            tStep);
#else
    SparseLinearSystemNM *linMethod = CreateUsrDefSparseLinSolver(ST_Petsc, 1, this->giveDomain(1), this);
    this->updateComponent(tStep, InternalRhs, this->giveDomain(1));
    this->updateComponent(tStep, NonLinearLhs, this->giveDomain(1));
    this->internalForces.times(-1);
    this->internalForces.add(externalForces);
    NM_Status status = linMethod->solve(this->stiffnessMatrix, &(this->internalForces), solutionVector);
    delete(linMethod);
#endif

    if ( status == NM_NoSuccess ) {
        OOFEM_ERROR2( "No success in solving problem at time step", tStep->giveNumber() );
    }
}

void StokesFlow :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    switch ( cmpn ) {
    case InternalRhs:
        this->internalForces.zero();
        this->assembleVectorFromElements(this->internalForces, tStep, EID_MomentumBalance_ConservationEquation, NodalInternalForcesVector, VM_Total,
                                         EModelDefaultEquationNumbering(), d);
        return;

    case NonLinearLhs:
        this->stiffnessMatrix->zero();
        this->assemble(this->stiffnessMatrix, tStep, EID_MomentumBalance_ConservationEquation, StiffnessMatrix,
                       EModelDefaultEquationNumbering(), d);
        return;

    default:
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
    // Force numbering of equations. First all velocities on domain, then pressures.

    int i, j, ndofs, nnodes, nelem;
    DofManager *inode;
    Domain *domain = this->giveDomain(id);
    TimeStep *currStep = this->giveCurrentStep();
    IntArray loc;
    Dof *jDof;
    DofIDItem type;

    this->domainNeqs.at(id) = 0;
    this->domainPrescribedNeqs.at(id) = 0;

    nnodes = domain->giveNumberOfDofManagers();

    // number first velocities
    for ( i = 1; i <= nnodes; i++ ) {
        inode = domain->giveDofManager(i);
        ndofs = inode->giveNumberOfDofs();
        for ( j = 1; j <= ndofs; j++ ) {
            jDof  =  inode->giveDof(j);
            type  =  jDof->giveDofID();
            if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                jDof->askNewEquationNumber(currStep);
            }
        }
    }

    // and then pressures
    for ( i = 1; i <= nnodes; i++ ) {
        inode = domain->giveDofManager(i);
        ndofs = inode->giveNumberOfDofs();
        for ( j = 1; j <= ndofs; j++ ) {
            jDof  =  inode->giveDof(j);
            type  =  jDof->giveDofID();
            if ( type == P_f ) {
                jDof->askNewEquationNumber(currStep);
            }
        }
    }

    // invalidate element local copies of location arrays
    nelem = domain->giveNumberOfElements();
    for ( i = 1; i <= nelem; i++ ) {
        domain->giveElement(i)->invalidateLocationArray();
    }

    return domainNeqs.at(id);
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
            return 0;
        }
    }

    return EngngModel :: checkConsistency();
}

void StokesFlow :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    DofIDItem type = iDof->giveDofID();
    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
        iDof->printSingleOutputAt(stream, atTime, 'v', EID_MomentumBalance, VM_Total, 1);
    } else if ( ( type == P_f ) ) {
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
        int nelem = domain->giveNumberOfElements();
        for (int j = 1; j <= nelem; j++ ) {
            domain->giveElement(j)->updateInternalState(tStep);
        }
    }
}

NumericalMethod *StokesFlow :: giveNumericalMethod(TimeStep *tStep)
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
        currentStep = new TimeStep(istep, this, 1, 0.0, this->deltaT, 1);
    } else {
        int istep =  currentStep->giveNumber() + 1;
        StateCounterType counter = currentStep->giveSolutionStateCounter() + 1;
        previousStep = currentStep;
        double dt = currentStep->giveTimeIncrement();
        double totalTime = currentStep->giveTime() + dt;
        currentStep = new TimeStep(istep, this, 1, totalTime, dt, counter);
    }

    return currentStep;
}
} // end namespace oofem
