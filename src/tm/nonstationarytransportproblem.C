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
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "function.h"
#include "sparsenonlinsystemnm.h"
#include "unknownnumberingscheme.h"

#ifdef __CEMHYD_MODULE
 #include "cemhyd/cemhydmat.h"
#endif

namespace oofem {
REGISTER_EngngModel(NonStationaryTransportProblem);

void TransportExternalForceAssembler :: vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const
{
    TransportElement *telem = static_cast< TransportElement* >(&element);
    telem->computeBCVectorAt(vec, tStep, mode);
    FloatArray tmp;
    telem->computeInternalSourceRhsVectorAt(tmp, tStep, mode);
    vec.add(tmp);
}


MidpointLhsAssembler :: MidpointLhsAssembler(bool lumped, double alpha) : 
    MatrixAssembler(), lumped(lumped), alpha(alpha)
{}


void MidpointLhsAssembler :: matrixFromElement(FloatMatrix &answer, Element &el, TimeStep *tStep) const
{
    FloatMatrix capacity;
    el.giveCharacteristicMatrix(answer, TangentStiffnessMatrix, tStep);
    el.giveCharacteristicMatrix(capacity, this->lumped ? LumpedMassMatrix : MassMatrix, tStep);
    answer.times(this->alpha);
    answer.add(1. / tStep->giveTimeIncrement(), capacity);
}


void IntSourceLHSAssembler :: matrixFromElement(FloatMatrix &answer, Element &el, TimeStep *tStep) const
{
    static_cast< TransportElement * >( &el )->computeIntSourceLHSMatrix(answer, tStep);
}


NonStationaryTransportProblem :: NonStationaryTransportProblem(int i, EngngModel *_master = NULL) : StationaryTransportProblem(i, _master)
{
    ndomains = 1;
    lumpedCapacityStab = 0;
    initT = 0.;
    deltaT = 0.;
    dtFunction = 0;
    internalVarUpdateStamp = 0;
    changingProblemSize = false;
    solverType = ST_Direct;
}

NonStationaryTransportProblem :: ~NonStationaryTransportProblem()
{
}


NumericalMethod *NonStationaryTransportProblem :: giveNumericalMethod(MetaStep *mStep)
// only one has reason for LinearStatic
//     - SolutionOfLinearEquations

{
    if (!linSolver) 
        linSolver.reset( classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this) );
        if ( !linSolver ) {
            OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
        }
    return linSolver.get();
}

IRResultType
NonStationaryTransportProblem :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    result = EngngModel :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    if ( ir->hasField(_IFT_NonStationaryTransportProblem_initt) ) {
        IR_GIVE_FIELD(ir, initT, _IFT_NonStationaryTransportProblem_initt);
    }

    if ( ir->hasField(_IFT_NonStationaryTransportProblem_deltat) ) {
        IR_GIVE_FIELD(ir, deltaT, _IFT_NonStationaryTransportProblem_deltat);
    } else if ( ir->hasField(_IFT_NonStationaryTransportProblem_deltatfunction) ) {
        IR_GIVE_FIELD(ir, dtFunction, _IFT_NonStationaryTransportProblem_deltatfunction);
    } else if ( ir->hasField(_IFT_NonStationaryTransportProblem_prescribedtimes) ) {
        IR_GIVE_FIELD(ir, discreteTimes, _IFT_NonStationaryTransportProblem_prescribedtimes);
    } else {
        OOFEM_WARNING("Time step not defined");
        return IRRT_BAD_FORMAT;
    }

    IR_GIVE_FIELD(ir, alpha, _IFT_NonStationaryTransportProblem_alpha);
    /* The following done in updateAttributes
     * if (this->giveNumericalMethod (giveCurrentStep())) nMethod -> instanciateFrom (ir);
     */
    // read lumped capacity stabilization flag
    if ( ir->hasField(_IFT_NonStationaryTransportProblem_lumpedcapa) ) {
        lumpedCapacityStab = 1;
    }

    //secure equation renumbering, otherwise keep efficient algorithms
    if ( ir->hasField(_IFT_NonStationaryTransportProblem_changingproblemsize) ) {
        changingProblemSize = true;
        UnknownsField.reset( new DofDistributedPrimaryField(this, 1, FT_TransportProblemUnknowns, 1) );
    } else {
        UnknownsField.reset( new PrimaryField(this, 1, FT_TransportProblemUnknowns, 1) );
    }

    //read other input data from StationaryTransportProblem
    StationaryTransportProblem :: initializeFrom(ir);

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;


    return IRRT_OK;
}


double NonStationaryTransportProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        int hash = this->giveUnknownDictHashIndx(mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            OOFEM_ERROR("Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode));
        }
    }

    if ( dof->__giveEquationNumber() == 0 ) {
        OOFEM_ERROR("invalid equation number on DoF %d", dof->giveDofID());
    }

    return UnknownsField->giveUnknownValue(dof, mode, tStep);
}


TimeStep *
NonStationaryTransportProblem :: giveSolutionStepWhenIcApply(bool force)
{
    if ( master && (!force)) {
        return master->giveSolutionStepWhenIcApply();
    } else {
        if ( !stepWhenIcApply ) {
            stepWhenIcApply.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, this->initT - giveDeltaT ( giveNumberOfFirstStep() ), giveDeltaT ( giveNumberOfFirstStep() ), 0) );
            //stepWhenIcApply.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, -deltaT, deltaT, 0) );
        }

        return stepWhenIcApply.get();
    }
}


Function *
NonStationaryTransportProblem :: giveDtFunction()
// Returns the load-time function of the receiver.
{
    if ( !dtFunction ) {
        return NULL;
    }

    return giveDomain(1)->giveFunction(dtFunction);
}

double
NonStationaryTransportProblem :: giveDeltaT(int n)
{
    if ( giveDtFunction() ) {
        return giveDtFunction()->evaluateAtTime(n);
    }

    if ( discreteTimes.giveSize() > 0 ) {
        return this->giveDiscreteTime(n) - this->giveDiscreteTime(n - 1);
    }

    return deltaT;
}

double
NonStationaryTransportProblem :: giveDiscreteTime(int iStep)
{
    if ( ( iStep > 0 ) && ( iStep <= discreteTimes.giveSize() ) ) {
        return ( discreteTimes.at(iStep) );
    }

    if ( ( iStep == 0 ) && ( iStep <= discreteTimes.giveSize() ) ) {
        return ( initT );
    }

    OOFEM_ERROR("invalid iStep");
    return 0.0;
}


TimeStep *
NonStationaryTransportProblem :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    double totalTime = this->initT;
    double intrinsicTime;
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep =  currentStep->giveNumber() + 1;
        totalTime = currentStep->giveTargetTime() + giveDeltaT(istep);
        counter = currentStep->giveSolutionStateCounter() + 1;
    } else {
        // first step -> generate initial step
        currentStep.reset( new TimeStep( *giveSolutionStepWhenIcApply() ) );
    }

    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(istep, this, 1, totalTime, this->giveDeltaT ( istep ), counter) );
    //set intrinsic time to time of integration
    intrinsicTime = currentStep->giveTargetTime();
//     intrinsicTime = previousStep->giveTargetTime() + this->alpha *this->giveDeltaT(istep);
    currentStep->setIntrinsicTime(intrinsicTime);
    return currentStep.get();
}


void NonStationaryTransportProblem :: solveYourselfAt(TimeStep *tStep)
{
    // Creates system of governing eq's and solves them at given tStep
    // The solution is stored in UnknownsField. If the problem is growing/decreasing, the UnknownsField is projected on DoFs when needed.
    // If equations are not renumbered, the algorithm is efficient without projecting unknowns to DoFs (nodes).

    //Right hand side
    FloatArray rhs;
    TimeStep *icStep = this->giveSolutionStepWhenIcApply();

    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    //Solution at the first time step needs history. Therefore, return back one time increment and create it.
    if ( tStep->isTheFirstStep() ) {

        bcRhs.resize(neq); //rhs vector from solution step i-1
        bcRhs.zero();

        this->applyIC(icStep);

        //project initial conditions to have temporary temperature in integration points

        //edge or surface load on elements
        //add internal source vector on elements
        this->assembleVectorFromElements( bcRhs, icStep, TransportExternalForceAssembler(),
                                         VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add prescribed value, such as temperature, on nodes
        this->assembleDirichletBcRhsVector( bcRhs, icStep, VM_Total,
                                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add nodal load
        this->assembleVectorFromDofManagers( bcRhs, icStep, ExternalForceAssembler(),
                                            VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    }

    //Create a new lhs matrix if necessary
    if ( tStep->isTheFirstStep() || this->changingProblemSize ) {

        conductivityMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        if ( !conductivityMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        conductivityMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling conductivity and capacity matrices\n");
#endif

        //Add contribution of alpha*K+C/dt (where K has contributions from conductivity and neumann b.c.s)
        this->assemble( *conductivityMatrix, icStep, MidpointLhsAssembler(lumpedCapacityStab, alpha),
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
    }

    //get the previous Rhs vector
    if ( !tStep->isTheFirstStep() && this->changingProblemSize ) {
        UnknownsField->initialize( VM_RhsTotal, tStep, bcRhs, EModelDefaultEquationNumbering() );
    }

    //prepare position in UnknownsField to store the results
    FloatArray *solutionVector;
    UnknownsField->advanceSolution(tStep);
    solutionVector = UnknownsField->giveSolutionVector(tStep);
//     solutionVector->resize(neq);
//     solutionVector->zero();

    //Initialize and give solutionVector from previous solution
    //copy previous solution vector so we can use solution-dependent boundary conditions
    if ( changingProblemSize ) {
        if ( !tStep->isTheFirstStep() ) {
            //copy recent solution to previous position, copy from hash=0 to hash=1(previous)
            copyUnknownsInDictionary( VM_Total, tStep, tStep->givePreviousStep() );
        }
        UnknownsField->initialize( VM_Total, tStep->givePreviousStep(), *solutionVector, EModelDefaultEquationNumbering() );
    } else {
        //copy previous solution vector to actual
        *solutionVector = *UnknownsField->giveSolutionVector( tStep->givePreviousStep() );
    }

    ///@todo missing this->updateInternalState(& TauStep);

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling rhs\n");
#endif
    // assembling load from elements
    rhs = bcRhs;
    rhs.times(1. - alpha);
    bcRhs.zero();
    //boundary conditions evaluated at targetTime
    this->assembleVectorFromElements( bcRhs, tStep, TransportExternalForceAssembler(),
                                     VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleDirichletBcRhsVector( bcRhs, tStep, VM_Total,
                                       EModelDefaultEquationNumbering(), this->giveDomain(1) );

    // assembling load from nodes
    this->assembleVectorFromDofManagers( bcRhs, tStep, InternalForceAssembler(), VM_Total,
                                        EModelDefaultEquationNumbering(), this->giveDomain(1) );
    for ( int i = 1; i <= neq; i++ ) {
        rhs.at(i) += bcRhs.at(i) * alpha;
    }

    // add the rhs part depending on previous solution
    assembleAlgorithmicPartOfRhs( rhs, EModelDefaultEquationNumbering(), tStep->givePreviousStep() );
    // set-up numerical model
    this->giveNumericalMethod( this->giveCurrentMetaStep() );

    // call numerical model to solve arised problem
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif
//     UnknownsField->giveSolutionVector(tStep)->resize(neq);
    linSolver->solve(*conductivityMatrix, rhs, *UnknownsField->giveSolutionVector(tStep) );
    // update solution state counter
    tStep->incrementStateCounter();
}

void
NonStationaryTransportProblem :: updateYourself(TimeStep *tStep)
{
    this->updateInternalState(tStep);
    EngngModel :: updateYourself(tStep);

    ///@todo Find a cleaner way to do these cemhyd hacks
#ifdef __CEMHYD_MODULE
    for ( auto &domain: this->domainList ) {
        for ( int i = 1; i <= domain->giveNumberOfElements(); ++i ) {
            TransportElement *elem = static_cast< TransportElement * >( domain->giveElement(i) );
            //store temperature and associated volume on each GP before performing averaging
            CemhydMat *cem = dynamic_cast< CemhydMat * >( elem->giveMaterial() );
            if ( cem ) {
                cem->clearWeightTemperatureProductVolume(elem);
                cem->storeWeightTemperatureProductVolume(elem, tStep);
            }
        }
        //perform averaging on each material instance
        for ( int i = 1; i <= domain->giveNumberOfMaterialModels(); i++ ) {
            CemhydMat *cem = dynamic_cast< CemhydMat * >( domain->giveMaterial(i) );
            if ( cem ) {
                cem->averageTemperature();
            }
        }
    }
 #ifdef VERBOSE
    VERBOSE_PRINT0("Updated Materials ", 0)
 #endif
#endif
}

void
NonStationaryTransportProblem :: copyUnknownsInDictionary(ValueModeType mode, TimeStep *fromTime, TimeStep *toTime)
{
    Domain *domain = this->giveDomain(1);

    for ( auto &node : domain->giveDofManagers() ) {
        for ( Dof *dof: *node ) {
            double val = dof->giveUnknown(mode, fromTime);
            dof->updateUnknownsDictionary(toTime, mode, val);
        }
    }
}


void
NonStationaryTransportProblem :: updateInternalState(TimeStep *tStep)
{
    for ( auto &domain: domainList ) {
        if ( requiresUnknownsDictionaryUpdate() ) {
            //update temperature vector
            UnknownsField->update( VM_Total, tStep, * ( this->UnknownsField->giveSolutionVector(tStep) ), EModelDefaultEquationNumbering() );
            //update Rhs vector
            UnknownsField->update(VM_RhsTotal, tStep, bcRhs, EModelDefaultEquationNumbering());
        }

        if ( internalVarUpdateStamp != tStep->giveSolutionStateCounter() ) {
            for ( auto &elem : domain->giveElements() ) {
                elem->updateInternalState(tStep);
            }

            internalVarUpdateStamp = tStep->giveSolutionStateCounter();
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
NonStationaryTransportProblem :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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
NonStationaryTransportProblem :: checkConsistency()
{
    // check internal consistency
    // if success returns nonzero
    Domain *domain = this->giveDomain(1);

    // check for proper element type
    for ( auto &elem : domain->giveElements() ) {
        if ( !dynamic_cast< TransportElement * >( elem.get() ) ) {
            OOFEM_WARNING("Element %d has no TransportElement base", elem->giveLabel());
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
NonStationaryTransportProblem :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    if ( mode == VM_Total ) { //Nodal temperature
        return 0;
    } else if ( mode == VM_RhsTotal ) { //Nodal Rhs
        return 1;
    } else {
        OOFEM_ERROR("ValueModeType %s undefined", __ValueModeTypeToString(mode));
    }

    return 0;
}


void
NonStationaryTransportProblem :: assembleAlgorithmicPartOfRhs(FloatArray &answer,
                                                              const UnknownNumberingScheme &s, TimeStep *tStep)
{
    IntArray loc;
    FloatMatrix charMtrx, charMtrx2;
    FloatArray unknownVec, contrib, intSource;
    Element *element;

    Domain *domain = this->giveDomain(1);
    int nelem = domain->giveNumberOfElements();

    for ( int i = 1; i <= nelem; i++ ) {
        element = domain->giveElement(i);
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

        element->giveLocationArray(loc, s);
        //(alpha-1)*K+C/dt
        element->giveCharacteristicMatrix(charMtrx, TangentStiffnessMatrix, tStep);
        element->giveCharacteristicMatrix(charMtrx2, lumpedCapacityStab ? LumpedMassMatrix : MassMatrix, tStep);

        charMtrx.times(this->alpha - 1.0);
        charMtrx.add(1. / tStep->giveTimeIncrement(), charMtrx2);

        if ( charMtrx.isNotEmpty() ) {
            element->computeVectorOf(VM_Total, tStep, unknownVec);
            contrib.beProductOf(charMtrx, unknownVec);
            answer.assemble(contrib, loc);
        }
    }
}


void
NonStationaryTransportProblem :: applyIC(TimeStep *stepWhenIcApply)
{
    Domain *domain = this->giveDomain(1);
    int neq =  this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
    FloatArray *solutionVector;
    double val;

#ifdef VERBOSE
    OOFEM_LOG_INFO("Applying initial conditions\n");
#endif

    UnknownsField->advanceSolution(stepWhenIcApply);
    solutionVector = UnknownsField->giveSolutionVector(stepWhenIcApply);
    solutionVector->resize(neq);
    solutionVector->zero();

    for ( auto &node : domain->giveDofManagers() ) {

        for ( Dof *dof: *node ) {
            // ask for initial values obtained from
            // bc (boundary conditions) and ic (initial conditions)
            if ( !dof->isPrimaryDof() ) {
                continue;
            }

            int jj = dof->__giveEquationNumber();
            if ( jj ) {
                val = dof->giveUnknown(VM_Total, stepWhenIcApply);
                solutionVector->at(jj) = val;
                //update in dictionary, if the problem is growing/decreasing
                if ( this->changingProblemSize ) {
                    dof->updateUnknownsDictionary(stepWhenIcApply, VM_Total, val);
                }
            }
        }
    }


    //project initial temperature to integration points

    //     for ( int j = 1; j <= nelem; j++ ) {
    //         domain->giveElement(j)->updateInternalState(stepWhenIcApply);
    //     }

#ifdef __CEMHYD_MODULE
    // Not relevant in linear case, but needed for CemhydMat for temperature averaging before solving balance equations
    // Update element state according to given ic
    for ( auto &elem : domain->giveElements() ) {
        TransportElement *element = static_cast< TransportElement * >( elem.get() );
        CemhydMat *cem = dynamic_cast< CemhydMat * >( element->giveMaterial() );
        //assign status to each integration point on each element
        if ( cem ) {
            cem->initMaterial(element); //create microstructures and statuses on specific GPs
            element->updateInternalState(stepWhenIcApply);   //store temporary unequilibrated temperature
            element->updateYourself(stepWhenIcApply);   //store equilibrated temperature
            cem->clearWeightTemperatureProductVolume(element);
            cem->storeWeightTemperatureProductVolume(element, stepWhenIcApply);
        }
    }

    //perform averaging on each material instance of CemhydMatClass
    for ( auto &mat : domain->giveMaterials() ) {
        CemhydMat *cem = dynamic_cast< CemhydMat * >( mat.get() );
        if ( cem ) {
            cem->averageTemperature();
        }
    }

#endif //__CEMHYD_MODULE
}


void
NonStationaryTransportProblem :: assembleDirichletBcRhsVector(FloatArray &answer, TimeStep *tStep,
                                                              ValueModeType mode,
                                                              const UnknownNumberingScheme &ns, Domain *d)
{
    IntArray loc, dofids;
    FloatArray rp, charVec;
    FloatMatrix s;
    FloatMatrix capacity;

    int nelem = d->giveNumberOfElements();

    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *element = d->giveElement(ielem);

        element->giveElementDofIDMask(dofids);
        element->computeVectorOfPrescribed(dofids, mode, tStep, rp);
        if ( rp.containsOnlyZeroes() ) {
            continue;
        } else {
            element->giveCharacteristicMatrix(s, TangentStiffnessMatrix, tStep);
            element->giveCharacteristicMatrix(capacity, lumpedCapacityStab ? LumpedMassMatrix : MassMatrix, tStep);
            s.times(this->alpha);
            s.add(1. / tStep->giveTimeIncrement(), capacity);

            charVec.beProductOf(s, rp);
            charVec.negated();

            element->giveLocationArray(loc, ns);
            answer.assemble(charVec, loc);
        }
    } // end element loop
}

#ifdef __CEMHYD_MODULE
// needed for CemhydMat
void
NonStationaryTransportProblem :: averageOverElements(TimeStep *tStep)
{
    ///@todo Verify this, the function is completely unused.
    Domain *domain = this->giveDomain(1);
    FloatArray vecTemperature;

    for ( auto &elem : domain->giveElements() ) {
        TransportMaterial *mat = dynamic_cast< CemhydMat * >( elem->giveMaterial() );
        if ( mat ) {
            for ( GaussPoint *gp: *elem->giveDefaultIntegrationRulePtr() ) {
                elem->giveIPValue(vecTemperature, gp, IST_Temperature, tStep);
                //mat->IP_volume += dV;
                //mat->average_temp += vecState.at(1) * dV;
            }
        }
    }

    for ( auto &mat : domain->giveMaterials() ) {
        CemhydMat *cem = dynamic_cast< CemhydMat * >( mat.get() );
        if ( cem ) {
            //cem->average_temp /= mat->IP_volume;
        }
    }
}
#endif
} // end namespace oofem
