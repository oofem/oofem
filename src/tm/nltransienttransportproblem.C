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

#include "nltransienttransportproblem.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "verbose.h"
#include "transportelement.h"
#include "classfactory.h"
#include "mathfem.h"
#include "assemblercallback.h"
#include "unknownnumberingscheme.h"

namespace oofem {
REGISTER_EngngModel(NLTransientTransportProblem);

NLTransientTransportProblem :: NLTransientTransportProblem(int i, EngngModel *_master = NULL) : NonStationaryTransportProblem(i, _master)
{
    //constructor
}

NLTransientTransportProblem :: ~NLTransientTransportProblem()
{
    //destructor
}

IRResultType
NLTransientTransportProblem :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    result = NonStationaryTransportProblem :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    int val = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NLTransientTransportProblem_nsmax);
    nsmax = val;

    IR_GIVE_FIELD(ir, rtol, _IFT_NLTransientTransportProblem_rtol);

    NR_Mode = nrsolverModifiedNRM;
    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, _IFT_NLTransientTransportProblem_manrmsteps);
    if ( MANRMSteps > 0 ) {
        NR_Mode = nrsolverAccelNRM;
    } else {
        NR_Mode = nrsolverModifiedNRM;
    }

    return IRRT_OK;
}

TimeStep *
NLTransientTransportProblem :: giveNextStep()
{
    double intrinsicTime;
    NonStationaryTransportProblem :: giveNextStep();
    intrinsicTime = previousStep->giveTargetTime() + this->alpha*(currentStep->giveTargetTime()-previousStep->giveTargetTime());
    currentStep->setIntrinsicTime(intrinsicTime);
    return currentStep.get();
}


void NLTransientTransportProblem :: solveYourselfAt(TimeStep *tStep)
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

        conductivityMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        if ( !conductivityMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        conductivityMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling conductivity and capacity matrices\n");
#endif
    }

    //create previous solution from IC or from previous tStep
    if ( tStep->isTheFirstStep() ) {
        if ( !stepWhenIcApply ) {
            stepWhenIcApply.reset( new TimeStep( *tStep->givePreviousStep() ) );
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
            conductivityMatrix->zero();
            //Assembling left hand side - start with conductivity matrix
            this->assemble( *conductivityMatrix, & TauStep, IntSourceLHSAssembler(),
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            conductivityMatrix->times(alpha);
            //Add capacity matrix
            this->assemble( *conductivityMatrix, & TauStep, MidpointLhsAssembler(lumpedCapacityStab, alpha),
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }

        rhs.resize(neq);
        rhs.zero();
        //edge or surface load on element
        //add internal source vector on elements
        this->assembleVectorFromElements( rhs, & TauStep, TransportExternalForceAssembler(), VM_Total,
                                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add nodal load
        this->assembleVectorFromDofManagers( rhs, & TauStep, ExternalForceAssembler(), VM_Total,
                                            EModelDefaultEquationNumbering(), this->giveDomain(1) );

        // subtract the rhs part depending on previous solution
        assembleAlgorithmicPartOfRhs(rhs, EModelDefaultEquationNumbering(), & TauStep);
        // set-up numerical model
        this->giveNumericalMethod( this->giveCurrentMetaStep() );

        // call numerical model to solve arised problem
#ifdef VERBOSE
        //OOFEM_LOG_INFO("Solving ...\n");
#endif

        // compute norm of residuals from balance equations
        solutionErr = rhs.computeNorm();

        linSolver->solve(*conductivityMatrix, rhs, solutionVectorIncrement);
        solutionVector->add(solutionVectorIncrement);
        this->updateInternalState(& TauStep); //insert to hash=0(current), if changes in equation numbering
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


double NLTransientTransportProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation
// This function translates this request to numerical method language
{
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        if ( mode == VM_Incremental ) { //get difference between current and previous time variable
            return dof->giveUnknowns()->at(0) - dof->giveUnknowns()->at(1);
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
NLTransientTransportProblem :: applyIC(TimeStep *stepWhenIcApply)
{
    Domain *domain = this->giveDomain(1);

    NonStationaryTransportProblem :: applyIC(stepWhenIcApply);

    // update element state according to given ic
    for ( auto &elem : domain->giveElements() ) {
        TransportElement *element = static_cast< TransportElement * >( elem.get() );
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
    }
}

void
NLTransientTransportProblem :: createPreviousSolutionInDofUnknownsDictionary(TimeStep *tStep)
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
NLTransientTransportProblem :: updateYourself(TimeStep *tStep)
{
    //this->updateInternalState(tStep);
    //Set intrinsic time for a staggered problem here. This is important for materials such as hydratingconcretemat, who keep history of intrinsic times.
    double intrinsicTime = tStep->giveIntrinsicTime();
    double Tau = tStep->giveTargetTime() - ( 1. - alpha ) * tStep->giveTimeIncrement();
    tStep->setIntrinsicTime(Tau);
    NonStationaryTransportProblem :: updateYourself(tStep);
    //return back to original intrinsic time
    tStep->setIntrinsicTime(intrinsicTime);
}

int
NLTransientTransportProblem :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
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
NLTransientTransportProblem :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
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
NLTransientTransportProblem :: updateInternalState(TimeStep *tStep)
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
NLTransientTransportProblem :: assembleAlgorithmicPartOfRhs(FloatArray &answer,
                                                            const UnknownNumberingScheme &ns, TimeStep *tStep)
{
    //
    // Computes right hand side on all nodes
    //
    double t = tStep->giveTargetTime();
    IntArray loc;
    FloatMatrix charMtrxCond, charMtrxCap;
    FloatArray r, drdt, contrib, help;
    Element *element;
    TimeStep *previousStep = this->givePreviousStep(); //r_t
    TimeStep *currentStep = this->giveCurrentStep(); //r_{t+\Delta t}. Note that *tStep is a Tau step between r_t and r_{t+\Delta t}

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

        element->giveCharacteristicMatrix(charMtrxCond, TangentStiffnessMatrix, tStep);
        element->giveCharacteristicMatrix(charMtrxCap, CapacityMatrix, tStep);


        /*
         *  element -> computeVectorOf (VM_Total, tStep, r);
         *  element -> computeVectorOf (VM_Velocity, tStep, drdt);
         */

        if ( ( t >= previousStep->giveTargetTime() ) && ( t <= currentStep->giveTargetTime() ) ) {
            FloatArray rp, rc;
            element->computeVectorOf(VM_Total, currentStep, rc);
            element->computeVectorOf(VM_Total, previousStep, rp);

            //approximate derivative with a difference
            drdt.beDifferenceOf(rc, rp);
            drdt.times( 1. / currentStep->giveTimeIncrement() );
            //approximate current solution from linear interpolation
            rp.times(1 - alpha);
            rc.times(alpha);
            r = rc;
            r.add(rp);
        } else {
            OOFEM_ERROR("unsupported time value");
        }


        if ( lumpedCapacityStab ) {
            int size = charMtrxCap.giveNumberOfRows();
            double s;
            for ( int j = 1; j <= size; j++ ) {
                s = 0.0;
                for ( int k = 1; k <= size; k++ ) {
                    s += charMtrxCap.at(j, k);
                    charMtrxCap.at(j, k) = 0.0;
                }

                charMtrxCap.at(j, j) = s;
            }
        }

        help.beProductOf(charMtrxCap, drdt);

        contrib.beProductOf(charMtrxCond, r);
        contrib.add(help);
        contrib.negated();

        answer.assemble(contrib, loc);
    }
}
} // end namespace oofem
