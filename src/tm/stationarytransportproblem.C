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
#include "maskedprimaryfield.h"
#include "verbose.h"
#include "transportelement.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "nrsolver.h"

namespace oofem {
REGISTER_EngngModel(StationaryTransportProblem);

StationaryTransportProblem :: StationaryTransportProblem(int i, EngngModel *_master = NULL) : EngngModel(i, _master)
{
    UnknownsField = NULL;
    conductivityMatrix = NULL;
    ndomains = 1;
    nMethod = NULL;
}

StationaryTransportProblem :: ~StationaryTransportProblem()
{
    delete conductivityMatrix;
    delete nMethod;
    delete UnknownsField;
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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    EngngModel :: initializeFrom(ir);

    int val = SMT_Skyline;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    this->sparseMtrxType = ( SparseMtrxType ) val;

    ///@todo Combine this option with structural problems, where it is possible to keep the secant tangent elastic tangent (or generally, the initial tangent) etc. One option should fit all common needs here.
    this->keepTangent = ir->hasField(_IFT_StationaryTransportProblem_keepTangent);

    // read field export flag
    IntArray exportFields;
    exportFields.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, exportFields, _IFT_StationaryTransportProblem_exportfields);
    if ( exportFields.giveSize() ) {
        IntArray mask(1);
        FieldManager *fm = this->giveContext()->giveFieldManager();
        for ( int i = 1; i <= exportFields.giveSize(); i++ ) {
            if ( exportFields.at(i) == FT_Temperature ) {
                mask.at(1) = T_f;
#ifdef FIELDMANAGER_USE_SHARED_PTR
                //std::tr1::shared_ptr<Field> _temperatureField = make_shared<MaskedPrimaryField>(FT_Temperature, this->UnknownsField, mask);
                std :: tr1 :: shared_ptr< Field > _temperatureField( new MaskedPrimaryField ( FT_Temperature, this->UnknownsField, mask ) );
                fm->registerField( _temperatureField, ( FieldType ) exportFields.at(i) );
#else
                MaskedPrimaryField *_temperatureField = new MaskedPrimaryField(FT_Temperature, this->UnknownsField, mask);
                fm->registerField(_temperatureField, ( FieldType ) exportFields.at(i), true);
#endif
            } else if ( exportFields.at(i) == FT_HumidityConcentration ) {
                mask.at(1) = C_1;
#ifdef FIELDMANAGER_USE_SHARED_PTR
                //std::tr1::shared_ptr<Field> _temperatureField = make_shared<MaskedPrimaryField>(FT_Temperature, this->UnknownsField, mask);
                std :: tr1 :: shared_ptr< Field > _concentrationField( new MaskedPrimaryField ( FT_HumidityConcentration, this->UnknownsField, mask ) );
                fm->registerField( _concentrationField, ( FieldType ) exportFields.at(i) );
#else
                MaskedPrimaryField *_concentrationField = new MaskedPrimaryField(FT_HumidityConcentration, this->UnknownsField, mask);
                fm->registerField(_concentrationField, ( FieldType ) exportFields.at(i), true);
#endif
            }
        }
    }

    if ( UnknownsField == NULL ) { // can exist from nonstationary transport problem
        UnknownsField = new PrimaryField(this, 1, FT_TransportProblemUnknowns, EID_ConservationEquation, 0);
    }

    return IRRT_OK;
}



double StationaryTransportProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
#ifdef DEBUG
    if ( dof->__giveEquationNumber() == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }
#endif
    return UnknownsField->giveUnknownValue(dof, mode, tStep);
}


TimeStep *StationaryTransportProblem :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    //int mstep = 1;
    StateCounterType counter = 1;

    delete previousStep;

    if ( currentStep != NULL ) {
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, ( double ) istep, 0., counter);
    // time and dt variables are set eq to 0 for statics - has no meaning
    return currentStep;
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
        if ( conductivityMatrix == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        conductivityMatrix->buildInternalStructure( this, 1, EID_ConservationEquation, EModelDefaultEquationNumbering() );
        if ( this->keepTangent ) {
            this->conductivityMatrix->zero();
            this->assemble( conductivityMatrix, tStep, EID_ConservationEquation, ConductivityMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( conductivityMatrix, tStep, EID_ConservationEquation, LHSBCMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }
    }

    internalForces.resize(neq);

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling external forces\n");
#endif
    FloatArray externalForces(neq);
    externalForces.zero();
    this->assembleVector( externalForces, tStep, EID_ConservationEquation, ExternalForcesVector, VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );

    // set-up numerical method
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    FloatArray incrementOfSolution;
    double loadLevel;
    int currentIterations;
    this->nMethod->solve(this->conductivityMatrix,
                         & externalForces,
                         NULL,
                         UnknownsField->giveSolutionVector(tStep),
                         & incrementOfSolution,
                         & this->internalForces,
                         this->eNorm,
                         loadLevel, // Only relevant for incrementalBCLoadVector
                         SparseNonLinearSystemNM :: rlm_total,
                         currentIterations,
                         tStep);

    //nMethod->solve( conductivityMatrix, & rhsVector, UnknownsField->giveSolutionVector(tStep) );
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
        this->assembleVector(this->internalForces, tStep, EID_ConservationEquation, InternalForcesVector, VM_Total,
                             EModelDefaultEquationNumbering(), this->giveDomain(1), & this->eNorm);
        return;
    } else if ( cmpn == NonLinearLhs ) {
        if ( !this->keepTangent ) {
            // Optimization for linear problems, we can keep the old matrix (which could save its factorization)
            this->conductivityMatrix->zero();
            ///@todo We should use some problem-neutral names instead of "ConductivityMatrix" (and something nicer for LHSBCMatrix)
            this->assemble( conductivityMatrix, tStep, EID_ConservationEquation, ConductivityMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( conductivityMatrix, tStep, EID_ConservationEquation, LHSBCMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }
        return;
    } else {
        OOFEM_ERROR("StationaryTransportProblem::updateComponent - Unknown component");
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

    if ( ( iores = UnknownsField->saveContext(stream, mode) ) != CIO_OK ) {
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
    FILE *file;

    this->resolveCorrespondingtStepumber(istep, iversion, obj);

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

    if ( ( iores = UnknownsField->restoreContext(stream, mode) ) != CIO_OK ) {
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
    int nelem = domain->giveNumberOfElements();
    for ( int i = 1; i <= nelem; i++ ) {
        Element *ePtr = domain->giveElement(i);
        if ( !dynamic_cast< TransportElement * >(ePtr) ) {
            _warning2("Element %d has no TransportElement base", i);
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
StationaryTransportProblem :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    iDof->printSingleOutputAt(stream, tStep, 'f', VM_Total);
}


void
StationaryTransportProblem :: updateInternalState(TimeStep *tStep)
{
    ///@todo Remove this, unnecessary with solving as a nonlinear problem (left for now, since nonstationary problems might still need it)
    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        Domain *domain = this->giveDomain(idomain);
        int nelem = domain->giveNumberOfElements();
        for ( int j = 1; j <= nelem; j++ ) {
            domain->giveElement(j)->updateInternalState(tStep);
        }
    }
}


#ifdef __PETSC_MODULE
void
StationaryTransportProblem :: initPetscContexts()
{
    petscContextList->growTo(ndomains);
    for ( int i = 1; i <= this->ndomains; i++ ) {
        PetscContext *petscContext =  new PetscContext(this);
        petscContextList->put(i, petscContext);
    }
}
#endif
} // end namespace oofem
