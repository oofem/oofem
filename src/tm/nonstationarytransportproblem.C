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

#include "nonstationarytransportproblem.h"
#include "stationarytransportproblem.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "maskedprimaryfield.h"
#include "verbose.h"
#include "transportelement.h"
#include "usrdefsub.h"
#include "datastream.h"
#include "contextioerr.h"

#ifdef __CEMHYD_MODULE
 #include "cemhydmat.h"
#endif //__CEMHYD_MODULE

namespace oofem {

NonStationaryTransportProblem :: NonStationaryTransportProblem(int i, EngngModel *_master = NULL) : StationaryTransportProblem(i, _master)
{
    UnknownsField = NULL;
    lhs = NULL;
    nMethod = NULL;
    ndomains = 1;
    lumpedCapacityStab = 0;
    initT = 0.;
    deltaT = 0.;
    //exportFieldFlag = 0;
    dtTimeFunction = 0;
    internalVarUpdateStamp = 0;
    changingProblemSize = false;
}

NonStationaryTransportProblem :: ~NonStationaryTransportProblem()
{
    if ( lhs ) {
        delete lhs;
    }

    if ( nMethod ) {
        delete nMethod;
    }

    if ( UnknownsField ) {
       // delete UnknownsField;
    }
}


NumericalMethod *NonStationaryTransportProblem :: giveNumericalMethod(MetaStep *mStep)
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
NonStationaryTransportProblem :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    EngngModel :: initializeFrom(ir);
    
    if ( ir->hasField(IFT_NonStationaryTransportProblem_initt, "initt") ) {
        IR_GIVE_FIELD(ir, initT, IFT_NonStationaryTransportProblem_initt, "initt"); // Macro
    }

    if ( ir->hasField(IFT_NonStationaryTransportProblem_deltat, "deltat") ) {
        IR_GIVE_FIELD(ir, deltaT, IFT_NonStationaryTransportProblem_deltat, "deltat"); // Macro
    } else if ( ir->hasField(IFT_NonStationaryTransportProblem_deltat, "deltatfunction") ){
        IR_GIVE_FIELD(ir, dtTimeFunction, IFT_NonStationaryTransportProblem_dtf, "deltatfunction"); // Macro
    } else {
        OOFEM_ERROR("Time step not defined");
    }

    IR_GIVE_FIELD(ir, alpha, IFT_NonStationaryTransportProblem_alpha, "alpha"); // Macro
    /* The following done in updateAttributes
     * if (this->giveNumericalMethod (giveCurrentStep())) nMethod -> instanciateFrom (ir);
     */
    // read lumped capacity stabilization flag
    if ( ir->hasField(IFT_NonStationaryTransportProblem_lumpedcapa, "lumpedcapa") ) {
        lumpedCapacityStab = 1;
    }

    //secure equation renumbering, otherwise keep efficient algorithms
    if ( ir->hasField(IFT_NonStationaryTransportProblem_changingproblemsize, "changingproblemsize") ) {
        changingProblemSize = true;
        UnknownsField = new DofDistributedPrimaryField(this, 1, FT_TransportProblemUnknowns, EID_ConservationEquation, 1);
    } else {
        UnknownsField = new PrimaryField(this, 1, FT_TransportProblemUnknowns, EID_ConservationEquation, 1);
    }

    //read other input data from StationaryTransportProblem
    StationaryTransportProblem :: initializeFrom(ir);

    return IRRT_OK;
}


double NonStationaryTransportProblem :: giveUnknownComponent(EquationID type, ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
    if ( type != EID_ConservationEquation ) { // heat and mass concetration vector
        OOFEM_ERROR2( "giveUnknownComponent: EquationID %s is undefined for this problem", __EquationIDToString(type) );
        return 0.;
    }

    if ( this->requiresUnknownsDictionaryUpdate() ) {
        int hash = this->giveUnknownDictHashIndx(type, mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    }

    if ( dof->__giveEquationNumber() == 0 ) {
        OOFEM_ERROR2( "giveUnknownComponent: invalid equation number on DoF %d", dof->giveNumber() );
    }

    return UnknownsField->giveUnknownValue(dof, mode, tStep);
}


TimeStep *
NonStationaryTransportProblem :: giveSolutionStepWhenIcApply()
{
    if ( stepWhenIcApply == NULL ) {
        stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, this->initT-giveDeltaT(giveNumberOfFirstStep()), giveDeltaT(giveNumberOfFirstStep()), 0);
        //stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, -deltaT, deltaT, 0);
    }
    return stepWhenIcApply;
}


LoadTimeFunction *
NonStationaryTransportProblem :: giveDtTimeFunction()
// Returns the load-time function of the receiver.
{
    if ( !dtTimeFunction || !ndomains ) {
        return NULL;
    }

    return giveDomain(1)->giveLoadTimeFunction(dtTimeFunction);
}

double
NonStationaryTransportProblem :: giveDeltaT(int n)
{
    if ( giveDtTimeFunction() ) {
        return giveDtTimeFunction()->__at(n);
    }

    return deltaT;
}



TimeStep *
NonStationaryTransportProblem :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    double totalTime = this->initT;
    double intrinsicTime;
    StateCounterType counter = 1;
    delete previousStep;

    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        totalTime = currentStep->giveTargetTime() + giveDeltaT(istep);
        counter = currentStep->giveSolutionStateCounter() + 1;
    } else {
        // first step -> generate initial step
        currentStep = new TimeStep( *giveSolutionStepWhenIcApply() );
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, totalTime, this->giveDeltaT(istep), counter);
    //set intrinsic time to time of integration
    intrinsicTime = previousStep->giveTargetTime() + this->alpha *this->giveDeltaT(istep);
    currentStep->setIntrinsicTime(intrinsicTime);
    // time and dt variables are set eq to 0 for statics - has no meaning
    return currentStep;
}


void NonStationaryTransportProblem :: solveYourselfAt(TimeStep *tStep)
{
    // Creates system of governing eq's and solves them at given tStep
    // The solution is stored in UnknownsField. If the problem is growing/decreasing, the UnknownsField is projected on DoFs when needed.
    // If equations are not renumbered, the algorithm is efficient without projecting unknowns to DoFs (nodes).

    //Right hand side
    FloatArray rhs;

    int neq = this->giveNumberOfEquations(EID_ConservationEquation);
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    //Solution at the first time step needs history. Therefore, return back one time increment and create it.
    if ( tStep->isTheFirstStep() ) {
        this->giveSolutionStepWhenIcApply();

        bcRhs.resize(neq); //rhs vector from solution step i-1
        bcRhs.zero();

        this->applyIC(stepWhenIcApply);

        //edge or surface load on elements
        this->assembleVectorFromElements( bcRhs, stepWhenIcApply, EID_ConservationEquation, ElementBCTransportVector,
                                         VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add prescribed value, such as temperature, on nodes
        this->assembleDirichletBcRhsVector( bcRhs, stepWhenIcApply, EID_ConservationEquation, VM_Total,
                                           NSTP_MidpointLhs, EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add internal source vector on elements
        this->assembleVectorFromElements( bcRhs, stepWhenIcApply, EID_ConservationEquation, ElementInternalSourceVector,
                                         VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add nodal load
        this->assembleVectorFromDofManagers( bcRhs, stepWhenIcApply, EID_ConservationEquation, ExternalForcesVector,
                                            VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    }

    //Create a new lhs matrix if necessary
    if ( tStep->isTheFirstStep() || this->changingProblemSize ) {
        if ( lhs ) {
            delete lhs;
        }

        lhs = CreateUsrDefSparseMtrx(sparseMtrxType);
        if ( lhs == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        lhs->buildInternalStructure( this, 1, EID_ConservationEquation, EModelDefaultEquationNumbering() );

#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling conductivity and capacity matrices\n");
#endif

        this->assemble( lhs, stepWhenIcApply, EID_ConservationEquation, LHSBCMatrix,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        lhs->times(alpha);
        this->assemble( lhs, stepWhenIcApply, EID_ConservationEquation, NSTP_MidpointLhs,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
    }

    //obtain the last Rhs vector from DoFs directly
    if ( !tStep->isTheFirstStep() && this->changingProblemSize ) {
        UnknownsField->initialize(VM_RhsTotal, tStep, bcRhs);
    }

    //prepare position in UnknownsField to store the results
    UnknownsField->advanceSolution(tStep);

    FloatArray *solutionVector = UnknownsField->giveSolutionVector(tStep);
    solutionVector->resize(neq);
    solutionVector->zero();

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling rhs\n");
#endif
    // assembling load from elements
    rhs = bcRhs;
    rhs.times(1. - alpha);
    bcRhs.zero();
    this->assembleVectorFromElements( bcRhs, tStep, EID_ConservationEquation, ElementBCTransportVector,
                                     VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleDirichletBcRhsVector( bcRhs, tStep, EID_ConservationEquation, VM_Total, NSTP_MidpointLhs,
                                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleVectorFromElements( bcRhs, tStep, EID_ConservationEquation, ElementInternalSourceVector,
                                     VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );

    // assembling load from nodes
    this->assembleVectorFromDofManagers( bcRhs, tStep, EID_ConservationEquation, InternalForcesVector, VM_Total,
                                        EModelDefaultEquationNumbering(), this->giveDomain(1) );
    for ( int i = 1; i <= neq; i++ ) {
        rhs.at(i) += bcRhs.at(i) * alpha;
    }

    // add the rhs part depending on previous solution
    assembleAlgorithmicPartOfRhs( rhs, EID_ConservationEquation,
                                 EModelDefaultEquationNumbering(), tStep->givePreviousStep() );
    // set-up numerical model
    this->giveNumericalMethod( this->giveCurrentMetaStep() );

    //
    // call numerical model to solve arised problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif
    UnknownsField->giveSolutionVector(tStep)->resize(neq);
    nMethod->solve( lhs, & rhs, UnknownsField->giveSolutionVector(tStep) );
    // update solution state counter
    tStep->incrementStateCounter();
}

void
NonStationaryTransportProblem :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    EngngModel :: updateYourself(stepN);
}


void
NonStationaryTransportProblem :: updateInternalState(TimeStep *stepN)
{
    int j, idomain, nelem;
    Domain *domain;

    for ( idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);

        if ( requiresUnknownsDictionaryUpdate() ) {
            //update temperature vector
            UnknownsField->update( VM_Total, stepN, * ( this->UnknownsField->giveSolutionVector(stepN) ) );
            //update Rhs vector
            UnknownsField->update(VM_RhsTotal, stepN, bcRhs);
        }

        if ( internalVarUpdateStamp != stepN->giveSolutionStateCounter() ) {
            nelem = domain->giveNumberOfElements();
            for ( j = 1; j <= nelem; j++ ) {
                domain->giveElement(j)->updateInternalState(stepN);
            }

            internalVarUpdateStamp = stepN->giveSolutionStateCounter();
        }
    }
}

contextIOResultType
NonStationaryTransportProblem :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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
NonStationaryTransportProblem :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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
NonStationaryTransportProblem :: checkConsistency()
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
NonStationaryTransportProblem :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}

int
NonStationaryTransportProblem :: giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN) {
    if ( mode == VM_Total ) { //Nodal temperature
        return 0;
    } else if ( mode == VM_RhsTotal ) { //Nodal Rhs
        return 1;
    } else {
        _error2( "ValueModeType %s undefined", __ValueModeTypeToString(mode) );
    }

    return 0;
}


void
NonStationaryTransportProblem :: giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                                                 CharType type, TimeStep *tStep, Domain *domain)
{
    // we don't directly call element->GiveCharacteristicMatrix() function, because some
    // engngm classes may require special modification of base types supported on
    // element class level

    if ( ( type == NSTP_MidpointLhs ) || ( type == NSTP_MidpointRhs ) ) {
        Element *element;
        FloatMatrix charMtrx1, charMtrx2;

        element = domain->giveElement(num);
        element->giveCharacteristicMatrix(answer, ConductivityMatrix, tStep);
        element->giveCharacteristicMatrix(charMtrx2, CapacityMatrix, tStep);

        if ( lumpedCapacityStab ) {
            int i, j, size = charMtrx2.giveNumberOfRows();
            double s;
            for ( i = 1; i <= size; i++ ) {
                s = 0.0;
                for ( j = 1; j <= size; j++ ) {
                    s += charMtrx2.at(i, j);
                    charMtrx2.at(i, j) = 0.0;
                }

                charMtrx2.at(i, i) = s;
            }
        }

        if ( type == NSTP_MidpointLhs ) {
            answer.times(this->alpha);
            charMtrx2.times( 1. / tStep->giveTimeIncrement() );
        } else {
            answer.times(this->alpha - 1.0);
            charMtrx2.times( 1. / tStep->giveTimeIncrement() );
        }

        answer.add(charMtrx2);
        return;
    } else {
        EngngModel :: giveElementCharacteristicMatrix(answer, num, type, tStep, domain);
    }
}


void
NonStationaryTransportProblem :: assembleAlgorithmicPartOfRhs(FloatArray &answer, EquationID ut,
                                                              const UnknownNumberingScheme &s, TimeStep *tStep)
{
    int i;
    IntArray loc;
    FloatMatrix charMtrx, bcMtrx;
    FloatArray unknownVec, contrib;
    Element *element;

    Domain *domain = this->giveDomain(1);
    int nelem = domain->giveNumberOfElements();

    for ( i = 1; i <= nelem; i++ ) {
        element = domain->giveElement(i);
#ifdef __PARALLEL_MODE
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif
        element->giveLocationArray(loc, ut, s);
        this->giveElementCharacteristicMatrix(charMtrx, i, NSTP_MidpointRhs, tStep, domain);
        element->giveCharacteristicMatrix(bcMtrx, LHSBCMatrix, tStep);
        bcMtrx.times(this->alpha - 1.0);
        if ( bcMtrx.isNotEmpty() ) {
            charMtrx.add(bcMtrx);
        }

        if ( charMtrx.isNotEmpty() ) {
            element->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, unknownVec);
            contrib.beProductOf(charMtrx, unknownVec);
            answer.assemble(contrib, loc);
        }
    }
}


void
NonStationaryTransportProblem :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    iDof->printSingleOutputAt(stream, atTime, 'f', EID_ConservationEquation, VM_Total);
}

void
NonStationaryTransportProblem :: applyIC(TimeStep *stepWhenIcApply)
{
    Domain *domain = this->giveDomain(1);
    int neq =  this->giveNumberOfEquations(EID_ConservationEquation);
    FloatArray *solutionVector;
    double val;

#ifdef VERBOSE
    OOFEM_LOG_INFO("Applying initial conditions\n");
#endif
    int nDofs, j, k, jj;
    int nman  = domain->giveNumberOfDofManagers();
    DofManager *node;
    Dof *iDof;

    UnknownsField->advanceSolution(stepWhenIcApply);
    solutionVector = UnknownsField->giveSolutionVector(stepWhenIcApply);
    solutionVector->resize(neq);
    solutionVector->zero();

    for ( j = 1; j <= nman; j++ ) {
        node = domain->giveDofManager(j);
        nDofs = node->giveNumberOfDofs();

        for ( k = 1; k <= nDofs; k++ ) {
            // ask for initial values obtained from
            // bc (boundary conditions) and ic (initial conditions)
            iDof  =  node->giveDof(k);
            if ( !iDof->isPrimaryDof() ) {
                continue;
            }

            jj = iDof->__giveEquationNumber();
            if ( jj ) {
                val = iDof->giveUnknown(EID_ConservationEquation, VM_Total, stepWhenIcApply);
                solutionVector->at(jj) = val;
                //update in dictionary, if the problem is growing/decreasing
                if ( this->changingProblemSize ) {
                    iDof->updateUnknownsDictionary(stepWhenIcApply, EID_MomentumBalance, VM_Total, val);
                }
            }
        }
    }

#ifdef __CEMHYD_MODULE
    // Not relevant in linear case, but needed for CemhydMat for temperature averaging before solving balance equations
    // Update element state according to given ic
    int nelem = domain->giveNumberOfElements();
    TransportElement *element;
    CemhydMat *cem;
    for ( j = 1; j <= nelem; j++ ) {
        element = ( TransportElement * ) domain->giveElement(j);
        //assign status to each integration point on each element
        if ( element->giveMaterial()->giveClassID() == CemhydMatClass ) {
            element->giveMaterial()->initMaterial(element); //create microstructures and statuses on specific GPs
            element->updateInternalState(stepWhenIcApply);   //store temporary unequilibrated temperature
            element->updateYourself(stepWhenIcApply);   //store equilibrated temperature
            cem = ( CemhydMat * ) element->giveMaterial();
            cem->clearWeightTemperatureProductVolume(element);
            cem->storeWeightTemperatureProductVolume(element, stepWhenIcApply);
        }
    }

    //perform averaging on each material instance of CemhydMatClass
    int nmat = domain->giveNumberOfMaterialModels();
    for ( j = 1; j <= nmat; j++ ) {
        if ( domain->giveMaterial(j)->giveClassID() == CemhydMatClass ) {
            cem = ( CemhydMat * ) domain->giveMaterial(j);
            cem->averageTemperature();
        }
    }
#endif //__CEMHYD_MODULE
}


void
NonStationaryTransportProblem :: assembleDirichletBcRhsVector(FloatArray &answer, TimeStep *tStep, EquationID ut,
                                                              ValueModeType mode, CharType lhsType,
                                                              const UnknownNumberingScheme &ns, Domain *d)
{
    int ielem;
    IntArray loc;
    Element *element;
    FloatArray rp, charVec;
    FloatMatrix s, bcMtrx;

    int nelem = d->giveNumberOfElements();

    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = d->giveElement(ielem);

        element->computeVectorOfPrescribed(EID_ConservationEquation, mode, tStep, rp);
        if ( rp.containsOnlyZeroes() ) {
            continue;
        } else {
            this->giveElementCharacteristicMatrix(s, ielem, lhsType, tStep, d);
            element->giveCharacteristicMatrix(bcMtrx, LHSBCMatrix, tStep);
            s.add(bcMtrx);
            charVec.beProductOf(s, rp);
            charVec.negated();

            element->giveLocationArray(loc, ut, ns);
            answer.assemble(charVec, loc);
        }
    } // end element loop

}

// needed for CemhydMat
void
NonStationaryTransportProblem :: averageOverElements(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    int ielem, i;
    int nelem = domain->giveNumberOfElements();
    double dV;
    TransportElement *element;
    IntegrationRule *iRule;
    GaussPoint *gp;
    FloatArray vecTemperature;
    TransportMaterial *mat;




    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = ( TransportElement * ) domain->giveElement(ielem);
        mat = ( TransportMaterial * ) element->giveMaterial();
        if ( mat->giveClassID() == CemhydMatClass ) {
            iRule = element->giveDefaultIntegrationRulePtr();
            for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
                gp  = iRule->getIntegrationPoint(i);
                dV  = element->computeVolumeAround(gp);
                element->giveIPValue(vecTemperature, gp, IST_Temperature, tStep);
                //mat->IP_volume += dV;
                //mat->average_temp += vecState.at(1) * dV;
            }
        }
    }

    for ( i = 1; i <= domain->giveNumberOfMaterialModels(); i++ ) {
        mat = ( TransportMaterial * ) domain->giveMaterial(i);
        if ( mat->giveClassID() == CemhydMatClass ) {
            //mat->average_temp /= mat->IP_volume;
        }
    }
}


#ifdef __PETSC_MODULE
void
NonStationaryTransportProblem :: initPetscContexts()
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
