/*
 * DarcyFlow.C
 *
 *  Created on: Feb 25, 2010
 *      Author: carl
 */

#include "darcyflow.h"
#include "element.h"
#include "inputrecord.h"
#include "timestep.h"
#include "classfactory.h"
#include "sparselinsystemnm.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "nrsolver.h"
#include "primaryfield.h"
#include "unknownnumberingscheme.h"

#include <iostream>
#include <fstream>
#include <cstdio>

namespace oofem {
REGISTER_EngngModel(DarcyFlow);

DarcyFlow :: DarcyFlow(int i, EngngModel *_master) : EngngModel(i, _master)
{
    this->ndomains = 1;
    this->hasAdvanced = false;
}

DarcyFlow :: ~DarcyFlow()
{
}

IRResultType DarcyFlow :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    result = EngngModel :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    // Create solution space for pressure field
    PressureField.reset( new PrimaryField(this, 1, FT_Pressure, 1) );
    return IRRT_OK;
}

void DarcyFlow :: solveYourselfAt(TimeStep *tStep)
{
    /*
     * Assemble and solve system of equations as given timestep tStep.
     *
     */

    OOFEM_LOG_INFO("Parallelflag: %u\n", parallelFlag);

    FloatArray *solutionVector = NULL;
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    // Move solution space to current timestep
    if ( !hasAdvanced ) {
        PressureField->advanceSolution(tStep);
        hasAdvanced = true;
    }

    // Point pointer SolutionVector to current solution in VelocityPressureField
    solutionVector = PressureField->giveSolutionVector(tStep);
    solutionVector->resize(neq);
    solutionVector->zero();

    // Create "stiffness matrix"
    if ( !this->stiffnessMatrix ) {
        this->stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        this->stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    }


    // Build initial/external load (LoadVector)
    this->externalForces.resize(neq);
    this->externalForces.zero();
    this->assembleVectorFromElements( this->externalForces, tStep, ExternalForceAssembler(), VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->updateSharedDofManagers(this->externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);

    this->incrementOfSolution.resize(neq);
    this->internalForces.resize(neq);

    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    double loadLevel;
    int currentIterations;
    this->updateComponent( tStep, InternalRhs, this->giveDomain(1) );
    this->updateComponent( tStep, NonLinearLhs, this->giveDomain(1) );

    NM_Status status = this->nMethod->solve(*this->stiffnessMatrix,
                                            this->externalForces,
                                            NULL,
                                            NULL,
                                            * solutionVector,
                                            this->incrementOfSolution,
                                            this->internalForces,
                                            this->ebeNorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total, // Why this naming scheme? Should be RLM_Total, and ReferenceLoadInputModeType
                                            currentIterations,
                                            tStep);

    if ( status & NM_NoSuccess ) {
        OOFEM_ERROR("couldn't solve for time step %d\n", tStep->giveNumber() );
    }

#define DUMPMATRICES 0
#if DUMPMATRICES
    FloatMatrix LHS_backup;
    lhs->toFloatMatrix(LHS_backup);
    DumpMatricesToFile(& LHS_backup, & rhs, NULL);
#endif
}


void DarcyFlow :: DumpMatricesToFile(FloatMatrix *LHS, FloatArray *RHS, FloatArray *SolutionVector)
{
    FILE *rhsFile = fopen("RHS.txt", "w");
    // rhs.printYourself();

    for ( int i = 1; i <= RHS->giveSize(); i++ ) {
        fprintf( rhsFile, "%0.15e\n", RHS->at(i) );
    }
    fclose(rhsFile);

    FILE *lhsFile = fopen("LHS.txt", "w");

    for ( int i = 1; i <= this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ); i++ ) {
        for ( int j = 1; j <= this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ); j++ ) {
            fprintf( lhsFile, "%0.15e\t", LHS->at(i, j) );
        }
        fprintf(lhsFile, "\n");
    }

    fclose(lhsFile);

    if ( SolutionVector == NULL ) {
        return;
    }

    FILE *SolutionFile = fopen("SolutionVector.txt", "w");
    for ( int i = 1; i <= SolutionVector->giveSize(); i++ ) {
        fprintf( SolutionFile, "%0.15e\n", SolutionVector->at(i) );
    }
    fclose(SolutionFile);
}

void DarcyFlow :: updateYourself(TimeStep *tStep)
{
    EngngModel :: updateYourself(tStep);
}

double DarcyFlow :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return PressureField->giveUnknownValue(dof, mode, tStep);
}

void DarcyFlow :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    switch ( cmpn ) {
    case InternalRhs:
        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, InternalForceAssembler(), VM_Total,
                             EModelDefaultEquationNumbering(), d, & this->ebeNorm);
        this->updateSharedDofManagers(this->externalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
        break;

    case NonLinearLhs:

        this->stiffnessMatrix->zero();
        this->assemble( *this->stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        break;

    default:
        OOFEM_ERROR("Unknown component id (%d)", ( int ) cmpn);
    }
}

int DarcyFlow :: forceEquationNumbering(int id)  // Is this really needed???!?
{
    int neq = EngngModel :: forceEquationNumbering(id);

    this->equationNumberingCompleted = false;
    this->stiffnessMatrix.reset(NULL);

    return neq;
}

NumericalMethod *DarcyFlow :: giveNumericalMethod(MetaStep *mStep)
{
    /*
     * Returns pointer to NumericalMethod object. The solver used for StokesFlow is SparseLinearSystemNM.
     * If no solver has bee initialized, create one, otherwise, return the existing solver.
     */

    if ( !this->nMethod ) {
        this->nMethod.reset( new NRSolver(this->giveDomain(1), this) );
        if ( !nMethod ) {
            OOFEM_ERROR("numerical method creation failed");
        }
    }

    return this->nMethod.get();
}

TimeStep *DarcyFlow :: giveNextStep()
{
    if ( !currentStep ) {
        // first step -> generate initial step
        //currentStep.reset( new TimeStep(*giveSolutionStepWhenIcApply()) );
        currentStep.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 1, 0., 1.0, 0) );
    }
    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(*previousStep, 1.0) );

    return currentStep.get();
}

}
