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

#include "tm/EngineeringModels/stationarytransportproblem.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "dof.h"
#include "maskedprimaryfield.h"
#include "verbose.h"
#include "tm/Elements/transportelement.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "nrsolver.h"
#include "unknownnumberingscheme.h"
#include "dofdistributedprimaryfield.h"

namespace oofem {
REGISTER_EngngModel(StationaryTransportProblem);

StationaryTransportProblem :: StationaryTransportProblem(int i, EngngModel *_master = nullptr) : EngngModel(i, _master),
    nMethod(nullptr)
{
    ndomains = 1;
}


NumericalMethod *StationaryTransportProblem :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod = std::make_unique<NRSolver>(this->giveDomain(1), this);
    }
    return nMethod.get();
}


void
StationaryTransportProblem :: initializeFrom(InputRecord &ir)
{
    EngngModel :: initializeFrom(ir);

    int val = SMT_Skyline;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    this->sparseMtrxType = ( SparseMtrxType ) val;

    ///@todo Combine this option with structural problems, where it is possible to keep the secant tangent elastic tangent (or generally, the initial tangent) etc. One option should fit all common needs here.
    this->keepTangent = ir.hasField(_IFT_StationaryTransportProblem_keepTangent);

    // read field export flag
    IntArray exportFields;
    IR_GIVE_OPTIONAL_FIELD(ir, exportFields, _IFT_StationaryTransportProblem_exportfields);
    if ( exportFields.giveSize() ) {
        FieldManager *fm = this->giveContext()->giveFieldManager();
        for ( int i = 1; i <= exportFields.giveSize(); i++ ) {
            if ( exportFields.at(i) == FT_Temperature ) {
                FieldPtr _temperatureField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->UnknownsField.get(), {T_f} ) );
                fm->registerField( _temperatureField, ( FieldType ) exportFields.at(i) );
            } else if ( exportFields.at(i) == FT_HumidityConcentration ) {
                FieldPtr _concentrationField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->UnknownsField.get(), {C_1} ) );
                fm->registerField( _concentrationField, ( FieldType ) exportFields.at(i) );
            }
        }
    }

    if ( !UnknownsField ) { // can exist from nonstationary transport problem
        //UnknownsField = std::make_unique<DofDistributedPrimaryField>(this, 1, FT_TransportProblemUnknowns, 0);
        UnknownsField = std::make_unique<PrimaryField>(this, 1, FT_TransportProblemUnknowns, 0);
    }
}



double StationaryTransportProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
#ifdef DEBUG
    if ( dof->__giveEquationNumber() == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }
#endif
    if (mode == VM_TotalIntrinsic) mode = VM_Total;
    return UnknownsField->giveUnknownValue(dof, mode, tStep);
}


FieldPtr StationaryTransportProblem::giveField (FieldType key, TimeStep *tStep)
{
    /* Note: the current implementation uses MaskedPrimaryField, that is automatically updated with the model progress, 
        so the returned field always refers to active solution step. 
    */

    if ( tStep != this->giveCurrentStep() ) {
        OOFEM_ERROR("Unable to return field representation for non-current time step");
    }
    if ( key == FT_Temperature ) {
        FieldPtr _ptr ( new MaskedPrimaryField ( key, this->UnknownsField.get(), {T_f} ) );
        return _ptr;
    } else if ( key == FT_HumidityConcentration ) {
        FieldPtr _ptr ( new MaskedPrimaryField ( key, this->UnknownsField.get(), {C_1} ) );
        return _ptr;
    } else {
        return FieldPtr();
    }
}



TimeStep *StationaryTransportProblem :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(istep, this, 1, ( double ) istep, 0., counter);
    return currentStep.get();
}


void StationaryTransportProblem :: solveYourselfAt(TimeStep *tStep)
{


    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step
    UnknownsField->advanceSolution(tStep);
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    if ( tStep->giveNumber() == 1 ) {
        // allocate space for solution vector
        FloatArray *solutionVector = UnknownsField->giveSolutionVector(tStep);
        solutionVector->resize(neq);
        solutionVector->zero();

        conductivityMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        if ( !conductivityMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        conductivityMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
        if ( this->keepTangent ) {
            this->conductivityMatrix->zero();
            this->assemble( *conductivityMatrix, tStep, TangentAssembler(TangentStiffness),
                            EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }
    }

    Domain *domain = this->giveDomain(1);
    if ( tStep->isTheFirstStep() ) {
        // update element state according to given ic
        for ( auto &elem : domain->giveElements() ) {
            TransportElement *element = static_cast< TransportElement * >( elem.get() );
            element->updateInternalState(tStep);
            element->updateYourself(tStep);
        }
    }

    internalForces.resize(neq);

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling external forces\n");
#endif
    FloatArray externalForces(neq);
    externalForces.zero();
    this->assembleVector( externalForces, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);

    // set-up numerical method
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving for %d unknowns\n", neq);
#endif

    FloatArray incrementOfSolution;
    double loadLevel;
    int currentIterations;
    this->nMethod->solve(*this->conductivityMatrix,
                         externalForces,
                         NULL,
                         *UnknownsField->giveSolutionVector(tStep),
                         incrementOfSolution,
                         this->internalForces,
                         this->eNorm,
                         loadLevel, // Only relevant for incrementalBCLoadVector
                         SparseNonLinearSystemNM :: rlm_total,
                         currentIterations,
                         tStep);

    //nMethod->solve( *conductivityMatrix, rhsVector, *UnknownsField->giveSolutionVector(tStep) );
}


void
StationaryTransportProblem :: updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d)
{
    // No-op: currently uses "PrimaryField" and passes along the reference
    ///@todo this->field->update(VM_Total, tStep, solutionVector, EModelDefaultEquationNumbering());
}


void
StationaryTransportProblem :: updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *enorm)
{
    answer.zero();
    this->assembleVector(answer, tStep, InternalForceAssembler(), VM_Total,
                         EModelDefaultEquationNumbering(), this->giveDomain(1), enorm);
    this->updateSharedDofManagers(answer, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
}


void
StationaryTransportProblem :: updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d)
{
    mat.zero();
    this->assemble(mat, tStep, TangentAssembler(TangentStiffness), EModelDefaultEquationNumbering(), this->giveDomain(1) );
}


void
StationaryTransportProblem :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    if ( cmpn == InternalRhs ) {
        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, InternalForceAssembler(), VM_Total,
                             EModelDefaultEquationNumbering(), this->giveDomain(1), & this->eNorm);
        this->updateSharedDofManagers(this->internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
        return;
    } else if ( cmpn == NonLinearLhs ) {
        if ( !this->keepTangent ) {
            // Optimization for linear problems, we can keep the old matrix (which could save its factorization)
            this->conductivityMatrix->zero();
            this->assemble( *conductivityMatrix, tStep, TangentAssembler(TangentStiffness),
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }
        return;
    } else {
        OOFEM_ERROR("Unknown component");
    }
}

void
StationaryTransportProblem :: saveContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: saveContext(stream, mode);
    UnknownsField->saveContext(stream);
}


void
StationaryTransportProblem :: restoreContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: restoreContext(stream, mode);
    UnknownsField->restoreContext(stream);
}


int
StationaryTransportProblem :: checkConsistency()
{
    Domain *domain = this->giveDomain(1);

    // check for proper element type
    for ( auto &elem : domain->giveElements() ) {
        if ( !dynamic_cast< TransportElement * >( elem.get() ) ) {
            OOFEM_WARNING("Element %d has no TransportElement base", elem->giveLabel());
            return 0;
        }
    }

    return EngngModel :: checkConsistency();
}


void
StationaryTransportProblem :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}

} // end namespace oofem
