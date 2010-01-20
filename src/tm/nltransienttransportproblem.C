/* $Header: /home/cvs/bp/oofem/tm/src/nltransienttransportproblem.C,v 1.2.4.1 2004/04/05 15:19:53 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#include "nummet.h"
#include "ldltfact.h"
#include "imlsolver.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "dofmanager.h"
#include "elementside.h"
#include "dof.h"

#include "verbose.h"
#include "conTable.h"
#include "transportelement.h"
#include "usrdefsub.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

namespace oofem {

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
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step
    double fluxErr, dispErr;
    int neq =  this->giveNumberOfEquations(EID_ConservationEquation);

    if ( initFlag ) {
        lhs = CreateUsrDefSparseMtrx(sparseMtrxType);
        if ( lhs == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        lhs->buildInternalStructure(this, 1, EID_ConservationEquation, EModelDefaultEquationNumbering());

        rhs.resize(neq);
        initFlag = 0;
    }

    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        TimeStep *stepWhenIcApply = tStep->givePreviousStep();

        this->applyIC(stepWhenIcApply);
    }

    double dTTau = tStep->giveTimeIncrement();
    double Tau = tStep->giveTime() - ( 1. - alpha ) * tStep->giveTimeIncrement();
    TimeStep TauStep(tStep->giveNumber(), this, tStep->giveMetaStepNumber(), Tau, dTTau, tStep->giveSolutionStateCounter() + 1);

    // predictor
    FloatArray *solutionVector, *prevSolutionVector;
    FluxField.advanceSolution(tStep);
    solutionVector = FluxField.giveSolutionVector(tStep);
    prevSolutionVector = FluxField.giveSolutionVector( tStep->givePreviousStep() );
    * solutionVector = * prevSolutionVector;
    /* smarter can be solutionVector = prevSolution + dt*D(prevSolution)/Dt*/
    this->updateInternalState(& TauStep);


    FloatArray solutionVectorIncrement(neq);
    int nite = 0;

    OOFEM_LOG_INFO("Time            Iter       ResidNorm       IncrNorm\n__________________________________________________________\n");


    do {
        nite++;

        // corrector
#ifdef VERBOSE
        // printf("\nAssembling conductivity and capacity matrices");
#endif

        if ( ( nite == 1 ) || ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
            lhs->zero();
            this->assemble( lhs, & TauStep, EID_ConservationEquation, LHSBCMatrix,
			    EModelDefaultEquationNumbering(), this->giveDomain(1) );
            this->assemble( lhs, & TauStep, EID_ConservationEquation, IntSourceLHSMatrix,
			    EModelDefaultEquationNumbering(), this->giveDomain(1) );
            lhs->times(alpha);
            this->assemble( lhs, & TauStep, EID_ConservationEquation, NSTP_MidpointLhs,
			    EModelDefaultEquationNumbering(), this->giveDomain(1) );
        }

#ifdef VERBOSE
        // printf("\nAssembling rhs");
#endif

        //
        // assembling the element part of load vector
        //
        rhs.zero();
        this->assembleVectorFromElements( rhs, & TauStep, EID_ConservationEquation, ElementBCTransportVector, VM_Total,
					  EModelDefaultEquationNumbering(), this->giveDomain(1) );
        // this->assembleDirichletBcRhsVector (rhs, &TauStep, VM_Total, NSTP_MidpointLhs, this->giveDomain(1));
        this->assembleVectorFromElements( rhs, & TauStep, EID_ConservationEquation, ElementInternalSourceVector, VM_Total,
					  EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //
        // assembling the nodal part of load vector
        //
        this->assembleVectorFromDofManagers( rhs, & TauStep, EID_ConservationEquation, NodalLoadVector, VM_Total,
					     EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //
        // add the rhs part depending on previous solution
        //
        assembleAlgorithmicPartOfRhs(rhs, EID_ConservationEquation, EModelDefaultEquationNumbering(), &TauStep, nite);
        //
        // set-up numerical model
        //
        this->giveNumericalMethod(tStep);

        //
        // call numerical model to solve arised problem
        //
#ifdef VERBOSE
        // printf("\nSolving ...");
#endif

        // compute norm of residual fluxes
        fluxErr = sqrt( dotProduct(rhs, rhs, neq) );

        //nMethod -> solveYourselfAt(tStep);
        nMethod->solve(lhs, & rhs, & solutionVectorIncrement);
        solutionVector->add(solutionVectorIncrement);
        this->updateInternalState(& TauStep);
        // compute error norms
        dispErr = sqrt( dotProduct(solutionVectorIncrement, solutionVectorIncrement, neq) );

        // update solution state counter
        TauStep.incrementStateCounter();
        tStep->incrementStateCounter();

        OOFEM_LOG_INFO("%-15e %-10d %-15e %-15e\n", tStep->giveTime(), nite, fluxErr, dispErr);

        if ( nite >= nsmax ) {
            _error2("convergence not reached after %d iterations", nsmax);
        }
    } while ( ( fabs(fluxErr) > rtol ) || ( fabs(dispErr) > rtol ) );

    //printf ("\n");
    // update nodes, elements, etc.
    this->updateYourself( this->giveCurrentStep() );
}


double NLTransientTransportProblem ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
                                                            TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displaacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
    double t = tStep->giveTime();

    TimeStep *previousStep = this->givePreviousStep(), *currentStep = this->giveCurrentStep();

    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }

    if ( chc != EID_ConservationEquation ) { // heat and mass concetration vector
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    if ( ( t >= previousStep->giveTime() ) && ( t <= currentStep->giveTime() ) ) {
        double rtdt = FluxField.giveUnknownValue(dof, VM_Total, currentStep);
        double rt   = FluxField.giveUnknownValue(dof, VM_Total, previousStep);
        double psi = ( t - previousStep->giveTime() ) / currentStep->giveTimeIncrement();
        if ( mode == VM_Velocity ) {
            return ( rtdt - rt ) / currentStep->giveTimeIncrement();
        } else if ( mode == VM_Total ) {
            return psi * rtdt + ( 1. - psi ) * rt;
        } else {
            _error("giveUnknownComponent: Unknown is of undefined mode for this problem");
        }
    } else {
        _error("giveUnknownComponent: unsupported value requested");
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
NLTransientTransportProblem :: updateYourself(TimeStep *stepN)
{
    //this->updateInternalState(stepN);
    NonStationaryTransportProblem :: updateYourself(stepN);
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
							    const UnknownNumberingScheme& ns, TimeStep *tStep, int nite)
{
    //
    // computes the real nodal fluxes on elements
    //
    int i;
    double t = tStep->giveTime();
    IntArray loc;
    FloatMatrix charMtrx, charMtrx2, bcMtrx;
    FloatArray r, drdt, contrib, help;
    Element *element;
    TimeStep *previousStep = this->givePreviousStep(), *currentStep = this->giveCurrentStep();

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
        element->giveLocationArray(loc, ut, ns);

        element->giveCharacteristicMatrix(charMtrx, ConductivityMatrix, tStep);
        element->giveCharacteristicMatrix(bcMtrx, LHSBCMatrix, tStep);
        element->giveCharacteristicMatrix(charMtrx2, CapacityMatrix, tStep);


        /*
         *  element -> computeVectorOf (EID_ConservationEquation, VM_Total, tStep, r);
         *  element -> computeVectorOf (EID_ConservationEquation, VM_Velocity, tStep, drdt);
         */

        if ( ( t >= previousStep->giveTime() ) && ( t <= currentStep->giveTime() ) ) {
            FloatArray rp, rc;
            element->computeVectorOf(EID_ConservationEquation, VM_Total, currentStep, rc);
            element->computeVectorOf(EID_ConservationEquation, VM_Total, previousStep, rp);

            drdt = rc;
            drdt.substract(rp);
            drdt.times( 1. / currentStep->giveTimeIncrement() );
            rp.times(1 - alpha);
            rc.times(alpha);
            r = rc;
            r.add(rp);
        } else {
            _error("assembleAlgorithmicPartOfRhs: unsupported time value");
        }


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

        help.beProductOf(charMtrx2, drdt);
        if ( bcMtrx.isNotEmpty() ) {
            charMtrx.plus(bcMtrx);
        }

        contrib.beProductOf(charMtrx, r);
        contrib.add(help);
        contrib.times(-1.0);

        answer.assemble(contrib, loc);
    }
}

} // end namespace oofem
