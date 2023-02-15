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

#include "mpmproblem.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "dictionary.h"
#include "verbose.h"
#include "classfactory.h"
#include "mathfem.h"
#include "assemblercallback.h"
#include "unknownnumberingscheme.h"

namespace oofem {
REGISTER_EngngModel(MPMProblem);

MPMLhsAssembler :: MPMLhsAssembler(CharType t, double c) : 
    MatrixAssembler(), type(t), factor(c)
{}


void MPMLhsAssembler :: matrixFromElement(FloatMatrix &answer, Element &el, TimeStep *tStep) const
{
    el.giveCharacteristicMatrix(answer, this->type, tStep);
    answer.times(this->factor);
}




MPMProblem :: MPMProblem(int i, EngngModel *_master = nullptr) : EngngModel(i, _master)
{
}

void
MPMProblem :: initializeFrom(InputRecord &ir)
{
    EngngModel :: initializeFrom(ir);

    if ( ir.hasField(_IFT_MPMProblem_initt) ) {
        IR_GIVE_FIELD(ir, initT, _IFT_MPMProblem_initt);
    }

    if ( ir.hasField(_IFT_MPMProblem_deltat) ) {
        IR_GIVE_FIELD(ir, deltaT, _IFT_MPMProblem_deltat);
    } else if ( ir.hasField(_IFT_MPMProblem_deltatfunction) ) {
        IR_GIVE_FIELD(ir, dtFunction, _IFT_MPMProblem_deltatfunction);
    } else if ( ir.hasField(_IFT_MPMProblem_prescribedtimes) ) {
        IR_GIVE_FIELD(ir, discreteTimes, _IFT_MPMProblem_prescribedtimes);
    } else {
        throw ValueInputException(ir, "none", "Time step not defined");
    }

    IR_GIVE_FIELD(ir, alpha, _IFT_MPMProblem_alpha);


    int val = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_MPMProblem_nsmax);
    nsmax = val;

    IR_GIVE_FIELD(ir, rtol, _IFT_MPMProblem_rtol);

    NR_Mode = nrsolverModifiedNRM;
    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, _IFT_MPMProblem_manrmsteps);
    if ( MANRMSteps > 0 ) {
        NR_Mode = nrsolverAccelNRM;
    } else {
        NR_Mode = nrsolverModifiedNRM;
    }
}

TimeStep *
MPMProblem :: giveNextStep()
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
        currentStep = std::make_unique<TimeStep>( *giveSolutionStepWhenIcApply() );
    }

    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(istep, this, 1, totalTime, this->giveDeltaT ( istep ), counter);
    //set intrinsic time to time of integration
    intrinsicTime = previousStep->giveTargetTime() + this->alpha*(currentStep->giveTargetTime()-previousStep->giveTargetTime());
    currentStep->setIntrinsicTime(intrinsicTime);
    return currentStep.get();
}


void MPMProblem :: solveYourselfAt(TimeStep *tStep)
{
    // creates system of governing eq's and solves them at given time step
    // first assemble problem at current time step

    // Right hand side
    FloatArray rhs;
    double solutionErr, incrementErr;
    int neq =  this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif
    //Delete lhs matrix and create a new one. This is necessary due to growing/decreasing number of equations.
    if ( tStep->isTheFirstStep() || this->changingProblemSize ) {

        jacobianMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        if ( !jacobianMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        jacobianMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling LHS matrix\n");
#endif
    }
    //create previous solution from IC or from previous tStep
    if ( tStep->isTheFirstStep() ) {
        if ( !stepWhenIcApply ) {
            stepWhenIcApply = std::make_unique<TimeStep>( *tStep->givePreviousStep() );
        }
        this->applyIC(stepWhenIcApply.get()); //insert solution to hash=1(previous), if changes in equation numbering
    }

    double dTTau = tStep->giveTimeIncrement();
    double Tau = tStep->giveTargetTime() - ( 1. - alpha ) * tStep->giveTimeIncrement();
    //Time step in which material laws are taken into account
    TimeStep TauStep(tStep->giveNumber(), this, tStep->giveMetaStepNumber(), Tau, dTTau, tStep->giveSolutionStateCounter() + 1);

    //Predictor
    FloatArray *solutionVector;
    UnknownsField->advanceSolution(tStep);
    solutionVector = UnknownsField->giveSolutionVector(tStep);

    //Initialize and give solutionVector from previous solution
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

    this->updateInternalState(& TauStep); //insert to hash=0(current), if changes in equation numbering

    FloatArray solutionVectorIncrement(neq);
    int nite = 0;

    OOFEM_LOG_INFO("Time            Iter       ResidNorm       IncrNorm\n__________________________________________________________\n");


    do {
        nite++;

        // Corrector
#ifdef VERBOSE
        // printf("\nAssembling conductivity and capacity matrices");
#endif

        if ( ( nite == 1 ) || ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
            jacobianMatrix->zero();
            MPMLhsAssembler conductivityAssembler(ConductivityMatrix, this->alpha);
            //Assembling left hand side - start with conductivity matrix
            this->assemble( *jacobianMatrix, & TauStep, conductivityAssembler,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            //Add capacity matrix
            MPMLhsAssembler capacityAssembler(CapacityMatrix, 1./tStep->giveTimeIncrement());
            this->assemble( *jacobianMatrix, & TauStep, capacityAssembler,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }

        rhs.resize(neq);
        rhs.zero();
        //edge or surface load on element
        //add internal source vector on elements
        //this->assembleVectorFromElements( rhs, tStep, TransportExternalForceAssembler(), VM_Total,
        //                                 EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add nodal load
        this->assembleVectorFromDofManagers( rhs, tStep, ExternalForceAssembler(), VM_Total,
                                            EModelDefaultEquationNumbering(), this->giveDomain(1) );

        // subtract the rhs part depending on previous solution
        assembleAlgorithmicPartOfRhs(rhs, EModelDefaultEquationNumbering(), tStep);
        // set-up numerical model
        this->giveNumericalMethod( this->giveCurrentMetaStep() );

        // call numerical model to solve arised problem
#ifdef VERBOSE
        //OOFEM_LOG_INFO("Solving ...\n");
#endif

        // compute norm of residuals from balance equations
        solutionErr = rhs.computeNorm();

        linSolver->solve(*jacobianMatrix, rhs, solutionVectorIncrement);
        solutionVector->add(solutionVectorIncrement);
        this->updateInternalState(tStep); //insert to hash=0(current), if changes in equation numbering
        // compute error in the solutionvector increment
        incrementErr = solutionVectorIncrement.computeNorm();

        // update solution state counter
        TauStep.incrementStateCounter();
        tStep->incrementStateCounter();

        OOFEM_LOG_INFO("%-15e %-10d %-15e %-15e\n", tStep->giveTargetTime(), nite, solutionErr, incrementErr);

        currentIterations = nite;

        if ( nite >= nsmax ) {
            OOFEM_ERROR("convergence not reached after %d iterations", nsmax);
        }
    } while ( ( fabs(solutionErr) > rtol ) || ( fabs(incrementErr) > rtol ) );
}


double MPMProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation
// This function translates this request to numerical method language
{
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        if ( mode == VM_Incremental ) { //get difference between current and previous time variable
            return dof->giveUnknowns()->at(0) - dof->giveUnknowns()->at(1);
        } else if ( mode == VM_TotalIntrinsic ) { // intrinsic value only for current step
            return this->alpha * dof->giveUnknowns()->at(0) + (1.-this->alpha) * dof->giveUnknowns()->at(1);
        }
        int hash = this->giveUnknownDictHashIndx(mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            OOFEM_ERROR("Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    }

    double t = tStep->giveTargetTime();
    TimeStep *previousStep = this->givePreviousStep(), *currentStep = this->giveCurrentStep();

    if ( dof->__giveEquationNumber() == 0 ) {
        OOFEM_ERROR("invalid equation number on DoF %d", dof->giveDofID() );
    }

    if ( ( t >= previousStep->giveTargetTime() ) && ( t <= currentStep->giveTargetTime() ) ) {
        ///@todo Shouldn't it be enough to just run this?
        //UnknownsField->giveUnknownValue(dof, mode, currentStep);
        double rtdt = UnknownsField->giveUnknownValue(dof, VM_Total, currentStep);
        double rt   = UnknownsField->giveUnknownValue(dof, VM_Total, previousStep);
        double psi = ( t - previousStep->giveTargetTime() ) / currentStep->giveTimeIncrement();
        if ( mode == VM_Velocity ) {
            return ( rtdt - rt ) / currentStep->giveTimeIncrement();
        } else if ( mode == VM_TotalIntrinsic ) {
            // only supported for current step
            return this->alpha * rtdt + ( 1. - this->alpha ) * rt; 
        } else if ( mode == VM_Total ) {
            return psi * rtdt + ( 1. - psi ) * rt;
        } else if ( mode == VM_Incremental ) {
            if ( previousStep->isIcApply() ) {
                return 0;
            } else {
                return ( rtdt - rt );
            }
        } else {
            OOFEM_ERROR("Unknown mode %s is undefined for this problem", __ValueModeTypeToString(mode) );
        }
    } else {
        OOFEM_ERROR("time value %f not within bounds %f and %f", t, previousStep->giveTargetTime(), currentStep->giveTargetTime() );
    }

    return 0.; // to make compiler happy;
}


void
MPMProblem :: applyIC(TimeStep *stepWhenIcApply)
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

    // update element state according to given ic
    for ( auto &elem : domain->giveElements() ) {
        Element *element = elem.get() ;
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
    }   
}

void
MPMProblem :: createPreviousSolutionInDofUnknownsDictionary(TimeStep *tStep)
{
    //Copy the last known temperature to be a previous solution
    for ( auto &domain: domainList ) {
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( auto &node : domain->giveDofManagers() ) {
                for ( Dof *dof: *node ) {
                    double val = dof->giveUnknown(VM_Total, tStep); //get number on hash=0(current)
                    dof->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, val);
                }
            }
        }
    }
}


void
MPMProblem :: updateYourself(TimeStep *tStep)
{
    //this->updateInternalState(tStep);
    //Set intrinsic time for a staggered problem here. This is important for materials such as hydratingconcretemat, who keep history of intrinsic times.
    for ( auto &domain: domainList ) {
        if ( requiresUnknownsDictionaryUpdate() ) {
            //update temperature vector
            UnknownsField->update( VM_Total, tStep, * ( this->UnknownsField->giveSolutionVector(tStep) ), EModelDefaultEquationNumbering() );
            //update Rhs vector
            //UnknownsField->update(VM_RhsTotal, tStep, bcRhs, EModelDefaultEquationNumbering());
        }

        if ( internalVarUpdateStamp != tStep->giveSolutionStateCounter() ) {
            for ( auto &elem : domain->giveElements() ) {
                elem->updateInternalState(tStep);
            }

            internalVarUpdateStamp = tStep->giveSolutionStateCounter();
        }
    }
    EngngModel :: updateYourself(tStep);
}


int
MPMProblem :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    if ( mode == VM_Total ) { //Nodal temperature
        if ( tStep->giveNumber() == this->giveCurrentStep()->giveNumber() ) { //current time
            return 0;
        } else if ( tStep->giveNumber() == this->giveCurrentStep()->giveNumber() - 1 ) { //previous time
            return 1;
        } else {
            OOFEM_ERROR("No history available at TimeStep %d = %f, called from TimeStep %d = %f", tStep->giveNumber(), tStep->giveTargetTime(), this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveTargetTime() );
        }
    } else {
        OOFEM_ERROR("ValueModeType %s undefined", __ValueModeTypeToString(mode));
    }

    return 0;
}

void
MPMProblem :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    // update DoF unknowns dictionary. Store the last and previous temperature only, see giveUnknownDictHashIndx
    for ( Dof *dof: *inode ) {
        int eqNum = dof->__giveEquationNumber();
        double val;
        if ( dof->hasBc(tStep) ) { // boundary condition
            val = dof->giveBcValue(VM_Total, tStep);
        } else {
            FloatArray *vect = this->UnknownsField->giveSolutionVector(tStep);
            val = vect->at(eqNum);
        }

        //update temperature, which is present in every node
        dof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}


void
MPMProblem :: updateInternalState(TimeStep *tStep)
{
    for ( auto &domain: domainList ) {
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( auto &dman : domain->giveDofManagers() ) {
                //update dictionary entry or add a new pair if the position is missing
                this->updateDofUnknownsDictionary(dman.get(), tStep);
            }
        }

        for ( auto &elem : domain->giveElements() ) {
            elem->updateInternalState(tStep);
        }
    }
}

void
MPMProblem :: assembleAlgorithmicPartOfRhs( FloatArray &answer,
                                            const UnknownNumberingScheme &ns, TimeStep *tStep)
{
    //
    // Computes right hand side on all nodes
    //
    IntArray loc;
    FloatArray contrib;
    Element *element;
    //TimeStep *previousStep = this->givePreviousStep(); //r_t
    //TimeStep *currentStep = this->giveCurrentStep(); //r_{t+\Delta t}. Note that *tStep is a Tau step between r_t and r_{t+\Delta t}

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

        if ( !element->isActivated(tStep) ) {
            continue;
        }

        element->giveLocationArray(loc, ns);

        element->giveCharacteristicVector(contrib, InternalForcesVector, VM_TotalIntrinsic, tStep);
        contrib.negated();
        answer.assemble(contrib, loc);
    }
}

double
MPMProblem :: giveDeltaT(int n)
{
    if ( giveDtFunction() ) {
        return giveDtFunction()->evaluateAtTime(n);
    }

    if ( discreteTimes.giveSize() > 0 ) {
        return this->giveDiscreteTime(n) - this->giveDiscreteTime(n - 1);
    }

    return deltaT;
}


void
MPMProblem :: copyUnknownsInDictionary(ValueModeType mode, TimeStep *fromTime, TimeStep *toTime)
{
    Domain *domain = this->giveDomain(1);

    for ( auto &node : domain->giveDofManagers() ) {
        for ( Dof *dof: *node ) {
            double val = dof->giveUnknown(mode, fromTime);
            dof->updateUnknownsDictionary(toTime, mode, val);
        }
    }
}

Function *
MPMProblem :: giveDtFunction()
// Returns the load-time function of the receiver.
{
    if ( !dtFunction ) {
        return NULL;
    }

    return giveDomain(1)->giveFunction(dtFunction);
}

double
MPMProblem :: giveDiscreteTime(int iStep)
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


} // end namespace oofem
