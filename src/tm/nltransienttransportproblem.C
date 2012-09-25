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

#include "nltransienttransportproblem.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "verbose.h"
#include "transportelement.h"
#include "usrdefsub.h"
#include "mathfem.h"

namespace oofem {

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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    NonStationaryTransportProblem :: initializeFrom(ir);
    int val = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_NLTransientTransportProblem_nsmax, "nsmax"); // Macro
    nsmax = val;

    IR_GIVE_FIELD(ir, rtol, IFT_NLTransientTransportProblem_rtol, "rtol"); // Macro

    NR_Mode = nrsolverModifiedNRM;
    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, IFT_NLTransientTransportProblem_manrmsteps, "manrmsteps"); // Macro
    if ( MANRMSteps > 0 ) {
        NR_Mode = nrsolverAccelNRM;
    } else {
        NR_Mode = nrsolverModifiedNRM;
    }

    return IRRT_OK;
}



void NLTransientTransportProblem :: solveYourselfAt(TimeStep *tStep) {
    // creates system of governing eq's and solves them at given time step
    // first assemble problem at current time step

    // Right hand side
    FloatArray rhs;
    double solutionErr, incrementErr;
    int neq =  this->giveNumberOfEquations(EID_ConservationEquation);

#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif
    //Delete lhs matrix and create a new one. This is necessary due to growing/decreasing number of equations.
    if ( tStep->isTheFirstStep() || this->changingProblemSize ) {
        if ( conductivityMatrix ) {
            delete conductivityMatrix;
        }

        conductivityMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
        if ( conductivityMatrix == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        conductivityMatrix->buildInternalStructure( this, 1, EID_ConservationEquation, EModelDefaultEquationNumbering() );
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling conductivity and capacity matrices\n");
#endif
    }

    //create previous solution from IC or from previous tStep
    if ( tStep->isTheFirstStep() ) {
        if( stepWhenIcApply == NULL ){
            stepWhenIcApply = new TimeStep(*tStep->givePreviousStep());
        }
        this->applyIC(stepWhenIcApply); //insert solution to hash=1(previous), if changes in equation numbering
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
            copyUnknownsInDictionary( EID_ConservationEquation, VM_Total, tStep, tStep->givePreviousStep() );
        }

        UnknownsField->initialize(VM_Total, tStep->givePreviousStep(), * solutionVector);
    } else {
        //copy previous solution vector to actual
        * solutionVector = * UnknownsField->giveSolutionVector( tStep->givePreviousStep() );
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
            this->assemble( conductivityMatrix, & TauStep, EID_ConservationEquation, LHSBCMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( conductivityMatrix, & TauStep, EID_ConservationEquation, IntSourceLHSMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            conductivityMatrix->times(alpha);
            //Add capacity matrix
            this->assemble( conductivityMatrix, & TauStep, EID_ConservationEquation, NSTP_MidpointLhs,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }

        rhs.resize(neq);
        rhs.zero();
        //edge or surface load on element
        this->assembleVectorFromElements( rhs, & TauStep, EID_ConservationEquation, ElementBCTransportVector, VM_Total,
                                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add internal source vector on elements
        this->assembleVectorFromElements( rhs, & TauStep, EID_ConservationEquation, ElementInternalSourceVector, VM_Total,
                                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //add nodal load
        this->assembleVectorFromDofManagers( rhs, & TauStep, EID_ConservationEquation, ExternalForcesVector, VM_Total,
                                            EModelDefaultEquationNumbering(), this->giveDomain(1) );

        // subtract the rhs part depending on previous solution
        assembleAlgorithmicPartOfRhs(rhs, EID_ConservationEquation, EModelDefaultEquationNumbering(), & TauStep);
        // set-up numerical model
        this->giveNumericalMethod( this->giveCurrentMetaStep() );

        // call numerical model to solve arised problem
#ifdef VERBOSE
        //OOFEM_LOG_INFO("Solving ...\n");
#endif

        // compute norm of residuals from balance equations
        solutionErr = rhs.computeNorm();

        nMethod->solve(conductivityMatrix, & rhs, & solutionVectorIncrement);
        solutionVector->add(solutionVectorIncrement);
        this->updateInternalState(& TauStep); //insert to hash=0(current), if changes in equation numbering
        // compute error in the solutionvector increment
        incrementErr = solutionVectorIncrement.computeNorm();

        // update solution state counter
        TauStep.incrementStateCounter();
        tStep->incrementStateCounter();

        OOFEM_LOG_INFO("%-15e %-10d %-15e %-15e\n", tStep->giveTargetTime(), nite, solutionErr, incrementErr);

        if ( nite >= nsmax ) {
            _error2("convergence not reached after %d iterations", nsmax);
        }
    } while ( ( fabs(solutionErr) > rtol ) || ( fabs(incrementErr) > rtol ) );
}


double NLTransientTransportProblem ::  giveUnknownComponent(EquationID type, ValueModeType mode,
                                                            TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation
// This function translates this request to numerical method language
{
    if ( type != EID_ConservationEquation ) { // heat and mass concetration vector
        OOFEM_ERROR2( "giveUnknownComponent: EquationID %s is undefined for this problem", __EquationIDToString(type) );
        return 0.;
    }

    if ( this->requiresUnknownsDictionaryUpdate() ) {
        if (mode == VM_Incremental) { //get difference between current and previous time variable
            return dof->giveUnknowns()->at(0) - dof->giveUnknowns()->at(1);
        }
        int hash = this->giveUnknownDictHashIndx(type, mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    }

    double t = tStep->giveTargetTime();
    TimeStep *previousStep = this->givePreviousStep(), *currentStep = this->giveCurrentStep();

    if ( dof->__giveEquationNumber() == 0 ) {
        OOFEM_ERROR2( "giveUnknownComponent: invalid equation number on DoF %d", dof->giveNumber() );
    }

    if ( ( t >= previousStep->giveTargetTime() ) && ( t <= currentStep->giveTargetTime() ) ) {
        double rtdt = UnknownsField->giveUnknownValue(dof, VM_Total, currentStep);
        double rt   = UnknownsField->giveUnknownValue(dof, VM_Total, previousStep);
        double psi = ( t - previousStep->giveTargetTime() ) / currentStep->giveTimeIncrement();
        if ( mode == VM_Velocity ) {
            return ( rtdt - rt ) / currentStep->giveTimeIncrement();
        } else if ( mode == VM_Total ) {
            return psi * rtdt + ( 1. - psi ) * rt;
        } else {
            OOFEM_ERROR2( "giveUnknownComponent: Unknown mode %s is undefined for this problem", __ValueModeTypeToString(mode) );
        }
    } else {
        OOFEM_ERROR4( "giveUnknownComponent: time value %f not within bounds %f and %f", t, previousStep->giveTargetTime(), currentStep->giveTargetTime() );
    }

    return 0.; // to make compiler happy;
}


void
NLTransientTransportProblem :: applyIC(TimeStep *stepWhenIcApply)
{
    int j;
    Domain *domain = this->giveDomain(1);

    NonStationaryTransportProblem :: applyIC(stepWhenIcApply);

    // update element state according to given ic
    int nelem = domain->giveNumberOfElements();
    TransportElement *element;

    for ( j = 1; j <= nelem; j++ ) {
        element = ( TransportElement * ) domain->giveElement(j);
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
    }
}

void
NLTransientTransportProblem :: createPreviousSolutionInDofUnknownsDictionary(TimeStep *tStep) {
    //Copy the last known temperature to be a previous solution
    int i, nnodes, inode, nDofs;
    double val;
    Domain *domain;
    Dof *iDof;
    DofManager *node;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        nnodes = domain->giveNumberOfDofManagers();
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( inode = 1; inode <= nnodes; inode++ ) {
                node = domain->giveDofManager(inode);
                nDofs = node->giveNumberOfDofs();
                for ( i = 1; i <= nDofs; i++ ) {
                    iDof = node->giveDof(i);
                    val = iDof->giveUnknown(EID_ConservationEquation, VM_Total, tStep); //get number on hash=0(current)
                    iDof->updateUnknownsDictionary(tStep->givePreviousStep(), EID_MomentumBalance, VM_Total, val);
                }
            }
        }
    }
}


void
NLTransientTransportProblem :: updateYourself(TimeStep *stepN)
{
    //this->updateInternalState(stepN);
    NonStationaryTransportProblem :: updateYourself(stepN);
}

int
NLTransientTransportProblem :: giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN) {
    if ( mode == VM_Total ) { //Nodal temperature
        if ( stepN->giveNumber() == this->giveCurrentStep()->giveNumber() ) { //current time
            return 0;
        } else if ( stepN->giveNumber() == this->giveCurrentStep()->giveNumber() - 1 ) { //previous time
            return 1;
        } else {
            _error5( "No history available at TimeStep %d = %f, called from TimeStep %d = %f", stepN->giveNumber(), stepN->giveTargetTime(), this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveTargetTime() );
        }
    } else {
        _error2( "ValueModeType %s undefined", __ValueModeTypeToString(mode) );
    }

    return 0;
}

void
NLTransientTransportProblem :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    // update DoF unknowns dictionary. Store the last and previous temperature only, see giveUnknownDictHashIndx
    int i, ndofs = inode->giveNumberOfDofs();
    int eqNum;
    Dof *iDof;
    double val;
    FloatArray *vect;

    for ( i = 1; i <= ndofs; i++ ) {
        iDof = inode->giveDof(i);
        eqNum = iDof->__giveEquationNumber();
        if ( iDof->hasBc(tStep) ) { // boundary condition
            val = iDof->giveBcValue(VM_Total, tStep);
        } else {
            vect = this->UnknownsField->giveSolutionVector(tStep);
            val = vect->at(eqNum);
        }

        //update temperature, which is present in every node
        iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total, val);
    }
}


void
NLTransientTransportProblem :: copyUnknownsInDictionary(EquationID type, ValueModeType mode, TimeStep *fromTime, TimeStep *toTime) {
    int i, j, ndofs;
    double val;
    Domain *domain = this->giveDomain(1);
    int nnodes = domain->giveNumberOfDofManagers();
    DofManager *inode;
    Dof *iDof;

    for ( j = 1; j <= nnodes; j++ ) {
        inode = domain->giveDofManager(j);
        ndofs = inode->giveNumberOfDofs();
        for ( i = 1; i <= ndofs; i++ ) {
            iDof = inode->giveDof(i);
            val = iDof->giveUnknown(type, mode, fromTime);
            iDof->updateUnknownsDictionary(toTime, type, mode, val);
        }
    }
}


void
NLTransientTransportProblem :: updateInternalState(TimeStep *stepN)
{
    int j, nnodes;
    Domain *domain;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        nnodes = domain->giveNumberOfDofManagers();
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( j = 1; j <= nnodes; j++ ) {
                //update dictionary entry or add a new pair if the position is missing
                this->updateDofUnknownsDictionary(domain->giveDofManager(j), stepN);
            }
        }

        int nelem = domain->giveNumberOfElements();
        for ( j = 1; j <= nelem; j++ ) {
            domain->giveElement(j)->updateInternalState(stepN);
        }
    }
}

void
NLTransientTransportProblem :: assembleAlgorithmicPartOfRhs(FloatArray &answer, EquationID ut,
                                                            const UnknownNumberingScheme &ns, TimeStep *tStep)
{
    //
    // Computes right hand side on all nodes
    //
    int i;
    double t = tStep->giveTargetTime();
    IntArray loc;
    FloatMatrix charMtrxCond, charMtrxCap, bcMtrx;
    FloatArray r, drdt, contrib, help;
    Element *element;
    TimeStep *previousStep = this->givePreviousStep(); //r_t
    TimeStep *currentStep = this->giveCurrentStep(); //r_{t+\Delta t}. Note that *tStep is a Tau step between r_t and r_{t+\Delta t}

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
        if ( !element->isActivated(tStep) ) {
            continue;
        }
        
        element->giveLocationArray(loc, ut, ns);

        element->giveCharacteristicMatrix(charMtrxCond, ConductivityMatrix, tStep);
        element->giveCharacteristicMatrix(bcMtrx, LHSBCMatrix, tStep);
        element->giveCharacteristicMatrix(charMtrxCap, CapacityMatrix, tStep);


        /*
         *  element -> computeVectorOf (EID_ConservationEquation, VM_Total, tStep, r);
         *  element -> computeVectorOf (EID_ConservationEquation, VM_Velocity, tStep, drdt);
         */

        if ( ( t >= previousStep->giveTargetTime() ) && ( t <= currentStep->giveTargetTime() ) ) {
            FloatArray rp, rc;
            element->computeVectorOf(EID_ConservationEquation, VM_Total, currentStep, rc);
            element->computeVectorOf(EID_ConservationEquation, VM_Total, previousStep, rp);

            //approximate derivative with a difference
            drdt.beDifferenceOf(rc, rp);
            drdt.times( 1. / currentStep->giveTimeIncrement() );
            //approximate current solution from linear interpolation
            rp.times(1 - alpha);
            rc.times(alpha);
            r = rc;
            r.add(rp);
        } else {
            _error("assembleAlgorithmicPartOfRhs: unsupported time value");
        }


        if ( lumpedCapacityStab ) {
            int i, j, size = charMtrxCap.giveNumberOfRows();
            double s;
            for ( i = 1; i <= size; i++ ) {
                s = 0.0;
                for ( j = 1; j <= size; j++ ) {
                    s += charMtrxCap.at(i, j);
                    charMtrxCap.at(i, j) = 0.0;
                }

                charMtrxCap.at(i, i) = s;
            }
        }

        help.beProductOf(charMtrxCap, drdt);
        if ( bcMtrx.isNotEmpty() ) {
            charMtrxCond.add(bcMtrx);
        }

        contrib.beProductOf(charMtrxCond, r);
        contrib.add(help);
        contrib.negated();

        answer.assemble(contrib, loc);
    }
}
} // end namespace oofem
