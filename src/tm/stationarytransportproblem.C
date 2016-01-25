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

#include "stationarytransportproblem.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "dof.h"
#include "maskedprimaryfield.h"
#include "intvarfield.h"
#include "verbose.h"
#include "transportelement.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "nrsolver.h"
#include "unknownnumberingscheme.h"

namespace oofem {
REGISTER_EngngModel(StationaryTransportProblem);

StationaryTransportProblem :: StationaryTransportProblem(int i, EngngModel *_master = NULL) : EngngModel(i, _master)
{
    ndomains = 1;
    nMethod = NULL;
}

StationaryTransportProblem :: ~StationaryTransportProblem()
{
    delete nMethod;
}

NumericalMethod *StationaryTransportProblem :: giveNumericalMethod(MetaStep *mStep)
{
    if ( nMethod ) {
        return nMethod;
    }

    nMethod = new NRSolver(this->giveDomain(1), this);
    return nMethod;
}

IRResultType
StationaryTransportProblem :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = EngngModel :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    int val = SMT_Skyline;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    this->sparseMtrxType = ( SparseMtrxType ) val;

    ///@todo Combine this option with structural problems, where it is possible to keep the secant tangent elastic tangent (or generally, the initial tangent) etc. One option should fit all common needs here.
    this->keepTangent = ir->hasField(_IFT_StationaryTransportProblem_keepTangent);

    // read field export flag
    IntArray exportFields;
    IR_GIVE_OPTIONAL_FIELD(ir, exportFields, _IFT_StationaryTransportProblem_exportfields);
    if ( exportFields.giveSize() ) {
        FieldManager *fm = this->giveContext()->giveFieldManager();
        for ( int i = 1; i <= exportFields.giveSize(); i++ ) {
            if ( exportFields.at(i) == FT_Temperature ) {
                FM_FieldPtr _temperatureField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->UnknownsField.get(), {T_f} ) );
                fm->registerField( _temperatureField, ( FieldType ) exportFields.at(i) );
            } else if ( exportFields.at(i) == FT_HumidityConcentration ) {
                FM_FieldPtr _concentrationField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->UnknownsField.get(), {C_1} ) );
                fm->registerField( _concentrationField, ( FieldType ) exportFields.at(i) );
            }
        }
    }

    if ( !UnknownsField ) { // can exist from nonstationary transport problem
        UnknownsField.reset( new PrimaryField(this, 1, FT_TransportProblemUnknowns, 0) );
    }

    return IRRT_OK;
}



double StationaryTransportProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
#ifdef DEBUG
    if ( dof->__giveEquationNumber() == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }
#endif
    return UnknownsField->giveUnknownValue(dof, mode, tStep);
}


EModelFieldPtr StationaryTransportProblem::giveField (FieldType key, TimeStep *tStep)
{
  /* Note: the current implementation uses MaskedPrimaryField, that is automatically updated with the model progress, 
     so the returned field always refers to active solution step. 
  */

  if ( tStep != this->giveCurrentStep()) {
    OOFEM_ERROR("Unable to return field representation for non-current time step");
  }
  if ( key == FT_Temperature ) {
    FM_FieldPtr _ptr ( new MaskedPrimaryField ( key, this->UnknownsField.get(), {T_f} ) );
    return _ptr;
  } else if ( key == FT_HumidityConcentration ) {
    FM_FieldPtr _ptr ( new MaskedPrimaryField ( key, this->UnknownsField.get(), {C_1} ) );
    return _ptr;
  } else {
    return FM_FieldPtr();
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
    currentStep.reset( new TimeStep(istep, this, 1, ( double ) istep, 0., counter) );
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

        conductivityMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        if ( conductivityMatrix == NULL ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        conductivityMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
        if ( this->keepTangent ) {
            this->conductivityMatrix->zero();
            this->assemble( *conductivityMatrix, tStep, TangentAssembler(TangentStiffness),
                            EModelDefaultEquationNumbering(), this->giveDomain(1) );
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
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    FloatArray incrementOfSolution;
    double loadLevel;
    int currentIterations;
    this->nMethod->solve(*this->conductivityMatrix,
                         externalForces,
                         NULL,
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
StationaryTransportProblem :: updateYourself(TimeStep *tStep)
{
    this->updateInternalState(tStep);
    EngngModel :: updateYourself(tStep);
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

contextIOResultType
StationaryTransportProblem :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( ( iores = UnknownsField->saveContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    // ensure consistent records

    return CIO_OK;
}



contextIOResultType
StationaryTransportProblem :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( ( iores = UnknownsField->restoreContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    // ensure consistent records

    return CIO_OK;
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


void
StationaryTransportProblem :: updateInternalState(TimeStep *tStep)
{
    ///@todo Remove this, unnecessary with solving as a nonlinear problem (left for now, since nonstationary problems might still need it)
    for ( auto &domain: domainList ) {
        for ( auto &elem : domain->giveElements() ) {
            elem->updateInternalState(tStep);
        }
    }
}

} // end namespace oofem
