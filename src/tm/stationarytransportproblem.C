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

#include "stationarytransportproblem.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "maskedprimaryfield.h"
#include "verbose.h"
#include "transportelement.h"
#include "usrdefsub.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {

NumericalMethod *StationaryTransportProblem :: giveNumericalMethod(MetaStep *mStep)
// only one has reason for LinearStatic
//     - SolutionOfLinearEquations

{
    if ( nMethod ) {
        return nMethod;
    }

    nMethod = CreateUsrDefSparseLinSolver(solverType, 1, this->giveDomain(1), this);
    if ( nMethod == NULL ) {
        _error("giveNumericalMethod: linear solver creation failed");
    }

    return nMethod;
}

IRResultType
StationaryTransportProblem :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    EngngModel :: initializeFrom(ir);
    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_StationaryTransportProblem_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_StationaryTransportProblem_smtype, "smtype"); // Macro
    sparseMtrxType = ( SparseMtrxType ) val;

    /* The following done in updateAttributes
     * if (this->giveNumericalMethod (giveCurrentStep())) nMethod -> instanciateFrom (ir);
     */

    // read field export flag
    IntArray exportFields;
    exportFields.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, exportFields, IFT_StationaryTransportProblem_exportfields, "exportfields");  // Macro
    if ( exportFields.giveSize() ) {
        IntArray mask(1);
        FieldManager *fm = this->giveContext()->giveFieldManager();
        for ( int i = 1; i <= exportFields.giveSize(); i++ ) {
            if ( exportFields.at(i) == FT_Temperature ) {
                mask.at(1) = T_f;
                MaskedPrimaryField *_temperatureField = new MaskedPrimaryField(FT_Temperature, &this->FluxField, mask);

                fm->registerField(_temperatureField, ( FieldType ) exportFields.at(i), true);
            } else if ( exportFields.at(i) == FT_HumidityConcentration ) {
                mask.at(1) = C_1;
                MaskedPrimaryField *_concentrationField = new MaskedPrimaryField(FT_HumidityConcentration, &this->FluxField, mask);

                fm->registerField(_concentrationField, ( FieldType ) exportFields.at(i), true);
            }
        }
    }


    return IRRT_OK;
}



double StationaryTransportProblem ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
                                                           TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displaacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }

    if ( chc == EID_ConservationEquation ) { // heat and mass concetration vector
        return FluxField.giveUnknownValue(dof, mode, tStep);
    } else {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
    }

    return 0.;
}


TimeStep *StationaryTransportProblem :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    //int mstep = 1;
    StateCounterType counter = 1;

    if (previousStep != NULL){
        delete previousStep;
    }

    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, ( double ) istep, 0., counter);
    // time and dt variables are set eq to 0 for staics - has no meaning
    return currentStep;
}


void StationaryTransportProblem :: solveYourselfAt(TimeStep *tStep) {
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step
    FluxField.advanceSolution(tStep);

    if ( tStep->giveNumber() == 1 ) {
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling conductivity matrix\n");
#endif

        // alocate space for solution vector
        FloatArray *solutionVector;
        solutionVector = FluxField.giveSolutionVector(tStep);
        solutionVector->resize( this->giveNumberOfEquations(EID_ConservationEquation) );
        solutionVector->zero();

        conductivityMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
        if ( conductivityMatrix == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        conductivityMatrix->buildInternalStructure( this, 1, EID_ConservationEquation, EModelDefaultEquationNumbering() );

        this->assemble( conductivityMatrix, tStep, EID_ConservationEquation, ConductivityMatrix,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        this->assemble( conductivityMatrix, tStep, EID_ConservationEquation, LHSBCMatrix,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
    }

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling rhs\n");
#endif

    //
    // assembling the element part of load vector
    //
    rhsVector.resize( this->giveNumberOfEquations(EID_ConservationEquation) );
    rhsVector.zero();

    this->assembleVectorFromElements( rhsVector, tStep, EID_ConservationEquation, ElementInternalSourceVector, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleVectorFromElements( rhsVector, tStep, EID_ConservationEquation, ElementBCTransportVector, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleDirichletBcRhsVector( rhsVector, tStep, EID_ConservationEquation, VM_Total, ConductivityMatrix,
                                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
    //
    // assembling the nodal part of load vector
    //
    this->assembleVectorFromDofManagers( rhsVector, tStep, EID_ConservationEquation, ExternalForcesVector, VM_Total,
                                        EModelDefaultEquationNumbering(), this->giveDomain(1) );

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveCurrentMetaStep() );

    //
    // call numerical model to solve arised problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    //nMethod -> solveYourselfAt(tStep);
    nMethod->solve( conductivityMatrix, & rhsVector, FluxField.giveSolutionVector(tStep) );
    // update solution state counter
    tStep->incrementStateCounter();
}

void
StationaryTransportProblem :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    EngngModel :: updateYourself(stepN);
}


contextIOResultType
StationaryTransportProblem :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    FILE *file;

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

    if ( ( iores = FluxField.saveContext(stream, mode) ) != CIO_OK ) {
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

    if ( ( iores = FluxField.restoreContext(stream, mode) ) != CIO_OK ) {
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
    // check internal consistency
    // if success returns nonzero
    int i, nelem;
    Element *ePtr;
    TransportElement *sePtr;
    Domain *domain = this->giveDomain(1);

    nelem = domain->giveNumberOfElements();
    // check for proper element type

    for ( i = 1; i <= nelem; i++ ) {
        ePtr = domain->giveElement(i);
        sePtr = dynamic_cast< TransportElement * >(ePtr);
        if ( sePtr == NULL ) {
            _warning2("Element %d has no TransportElement base", i);
            return 0;
        }
    }

    EngngModel :: checkConsistency();

    return 1;
}


void
StationaryTransportProblem :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}


void
StationaryTransportProblem :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    iDof->printSingleOutputAt(stream, atTime, 'f', EID_ConservationEquation, VM_Total);
}

void
StationaryTransportProblem :: assembleDirichletBcRhsVector(FloatArray &answer, TimeStep *tStep, EquationID ut,
                                                           ValueModeType mode, CharType lhsType,
                                                           const UnknownNumberingScheme &ns, Domain *d)
{
    int ielem;
    IntArray loc;
    Element *element;
    FloatArray rp, charVec;
    FloatMatrix s;

    int nelem = d->giveNumberOfElements();

    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = d->giveElement(ielem);

        element->computeVectorOfPrescribed(EID_ConservationEquation, mode, tStep, rp);
        if ( rp.containsOnlyZeroes() ) {
            continue;
        } else {
            this->giveElementCharacteristicMatrix(s, ielem, lhsType, tStep, d);
            charVec.beProductOf(s, rp);
            charVec.negated();

            element->giveLocationArray(loc, ut, ns);
            answer.assemble(charVec, loc);
        }
    }

    // end element loop
}


void
StationaryTransportProblem :: updateInternalState(TimeStep *stepN)
{
    int j, idomain, nelem;
    Domain *domain;

    for ( idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        nelem = domain->giveNumberOfElements();
        for ( j = 1; j <= nelem; j++ ) {
            domain->giveElement(j)->updateInternalState(stepN);
        }
    }
}



#ifdef __PETSC_MODULE
void
StationaryTransportProblem :: initPetscContexts()
{
    int i;
    PetscContext *petscContext;


    petscContextList->growTo(ndomains);
    for ( i = 0; i < this->ndomains; i++ ) {
        petscContext =  new PetscContext(this, EID_ConservationEquation);
        petscContextList->put(i + 1, petscContext);
    }
}
#endif
} // end namespace oofem
