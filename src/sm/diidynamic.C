/* $Header: /home/cvs/bp/oofem/sm/src/diidynamic.C,v 1.5.4.1 2004/04/05 15:19:46 bp Exp $ */
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

//
// file DIIDynamic.cc
//

#include "diidynamic.h"
#include "nummet.h"
#include "ldltfact.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "elementside.h"
#include "dof.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#include "verbose.h"
#include "skyline.h"

namespace oofem {
NumericalMethod *DIIDynamic :: giveNumericalMethod(TimeStep *)
// only one has reason for DIIDynamic
//     - SolutionOfLinearEquations

{
    if ( nMethod ) {
        return nMethod;
    }

    SparseLinearSystemNM *nm;
    nm = ( SparseLinearSystemNM * ) new LDLTFactorization(1, this->giveDomain(1), this);
    nMethod = nm;
    return nm;
}

IRResultType
DIIDynamic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralEngngModel :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, deltaT, IFT_DIIDynamic_deltat, "deltat"); // Macro
    IR_GIVE_FIELD(ir, alpha, IFT_DIIDynamic_alpha, "alpha"); // Macro
    IR_GIVE_FIELD(ir, beta, IFT_DIIDynamic_beta, "beta"); // Macro
    IR_GIVE_FIELD(ir, Psi, IFT_DIIDynamic_psi, "psi"); // Macro
    if ( Psi < 1.1 ) {
        Psi = 1.0;
    }

    if ( ( Psi > 1.1 ) && ( Psi < 1.37 ) ) {
        Psi = 1.37;
    }

    return IRRT_OK;
}





double DIIDynamic ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
                                           TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displaacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }

    if ( tStep != this->giveCurrentStep() ) {
        _error("giveUnknownComponent: unknown time step encountered");
        return 0.;
    }

    if ( chc != EID_MomentumBalance ) {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    switch ( mode ) {
    case VM_Total:
        return displacementVector.at(eq);

    case VM_Velocity:
        return velocityVector.at(eq);

    case VM_Acceleration:
        return accelerationVector.at(eq);

    default:
        _error("giveUnknownComponent: Unknown is of undefined ValueModeType for this problem");
    }

    return 0.0;
}


TimeStep *DIIDynamic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    double totalTime = 0;
    StateCounterType counter = 1;

    delete previousStep;
    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTargetTime() + deltaT;
        istep     = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, totalTime, deltaT, counter);

    return currentStep;
}


void DIIDynamic :: solveYourselfAt(TimeStep *tStep) {
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step
    Domain *domain = this->giveDomain(1);
    FloatArray help;
    int neq =  this->giveNumberOfEquations(EID_MomentumBalance);
    int i;

    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        TimeStep *stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0,
                                                 -deltaT, deltaT, 0);
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling stiffness and mass matrices\n");
#endif

        //
        // first step  assemble stiffness and mass Matrices
        //
        /*
         * IntArray* mht = this -> GiveBanWidthVector ();
         * stiffnessMatrix = new Skyline ();
         * massMatrix      = new Skyline ();
         * stiffnessMatrix ->  checkSizeTowardsBanWidth (mht) ;
         * massMatrix      ->  checkSizeTowardsBanWidth (mht) ;
         * delete mht;
         */
        stiffnessMatrix = new Skyline();
        stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );
        massMatrix = stiffnessMatrix->GiveCopy();

        this->assemble(massMatrix, tStep, EID_MomentumBalance, MassMatrix, EModelDefaultEquationNumbering(), domain);

        //
        // determining starting displacemnts, velocities, accelerations
        //
        //previousLoadVector = new FloatArray (neq);
        //displacementVector = new FloatArray (neq);
        //velocityVector     = new FloatArray (neq);
        //accelerationVector = new FloatArray (neq);
        previousLoadVector.resize(neq);
        previousLoadVector.zero();
        displacementVector.resize(neq);
        displacementVector.zero();
        velocityVector.resize(neq);
        velocityVector.zero();
        accelerationVector.resize(neq);
        accelerationVector.zero();

        int nDofs, j, k, jj;
        int nman  = domain->giveNumberOfDofManagers();
        DofManager *node;
        Dof *iDof;

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
                    displacementVector.at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Total, stepWhenIcApply);
                    velocityVector.at(jj)     = iDof->giveUnknown(EID_MomentumBalance, VM_Velocity, stepWhenIcApply);
                    accelerationVector.at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Acceleration, stepWhenIcApply);
                }
            }
        }

        //
        // determining constants a0 .. a10
        //
        if ( Psi < 1.1 ) { // Newmark Method
            Psi = 1.;
            double dt2 = deltaT * deltaT;
            a0  = ( 4. + 2. * alpha * deltaT ) / ( dt2 + 2. * beta * deltaT );
            a1  = 4. / dt2 + 2 * ( alpha - a0 * beta ) / deltaT;
            a2  = 4. / deltaT + alpha - a0 * beta;
            a3  = 1.;
            a4  = 4. / ( dt2 + 2. * beta * deltaT );
            a5  = -a4;
            a6  = a5 * ( deltaT + beta );
            a7  = -1.;
            a8  = deltaT / 2.;
            a9  = dt2 / 4.;
            a10 = a9;
        } else {     // Wilson Method
            if ( Psi < 1.37 ) {
                Psi = 1.37;
            }

            double tau = Psi * deltaT;
            a0  = ( 6. + 3. * alpha * tau ) / ( tau * tau + 3. * beta * tau );
            a1  = 6. / ( tau * tau ) + 3. * ( alpha - beta * a0 ) / tau;
            a2  = 6. / tau + 2. * ( alpha - beta * a0 );
            a3  = 2. + tau * ( alpha - beta * a0 ) / 2.;
            a4  = 3. / ( ( 3. * beta + tau ) * tau );
            a5  = 3. * beta * a4 / tau - 6. / ( Psi * tau * tau );
            a6  = 2. * beta * a4 - 6. / ( Psi * tau );
            a7  = 1. - 3. / Psi + beta * tau * a4 / 2.;
            a8  = deltaT / 2.;
            a9  = deltaT * deltaT / 3.;
            a10 = deltaT * deltaT / 6.;
        }

        //
        // assemble LHS of problem ( K* = K + a0*M)
        //

        this->assemble(stiffnessMatrix, tStep, EID_MomentumBalance, DIIModifiedStiffnessMatrix,
                       EModelDefaultEquationNumbering(), domain);
        /*
         *  Element* element;
         *  IntArray* loc ;
         *
         *  int nelem = domain -> giveNumberOfElements ();
         *  FloatMatrix *charMtrx1,*charMtrx2;
         *  for ( i = 1; i <= nelem ; i++ ) {
         *    element = domain -> giveElement(i);
         *    loc = element -> giveLocationArray ();
         *    charMtrx1 = element -> GiveCharacteristicMatrix (StiffnessMatrix, tStep );
         *    charMtrx2 = element -> GiveCharacteristicMatrix (MassMatrix, tStep);
         *    charMtrx2 -> times(a0);
         *    charMtrx1 -> plus (charMtrx2);
         *    stiffnessMatrix ->  assemble (charMtrx1, loc) ;
         * delete charMtrx1;
         * delete charMtrx2;
         * }
         */
        delete stepWhenIcApply;
    }   // end of initializaton for first time step

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling load\n");
#endif

    //
    // assembling the element part of load vector
    //
    //loadVector = new FloatArray (this->giveNumberOfEquations());
    loadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
    loadVector.zero();

    this->assembleVectorFromElements(loadVector, tStep, EID_MomentumBalance, ElementForceLoadVector,
                                     VM_Total, EModelDefaultEquationNumbering(), domain);

    //
    // assembling the nodal part of load vector
    //
    this->assembleVectorFromDofManagers(loadVector, tStep, EID_MomentumBalance, NodalLoadVector,
                                        VM_Total, EModelDefaultEquationNumbering(), domain);

    //
    // assembling modified load vector
    //
    help.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
    for ( i = 1; i <= neq; i++ ) {
        help.at(i) = a1 * displacementVector.at(i) +
                     a2 *velocityVector.at(i)     +
                     a3 *accelerationVector.at(i);
    }

    massMatrix->times(help, rhs);
    //delete help;
    for ( i = 1; i <= neq; i++ ) {
        rhs.at(i) += previousLoadVector.at(i) +
                     Psi * ( loadVector.at(i) - previousLoadVector.at(i) );
        // store load vector
        previousLoadVector.at(i) = loadVector.at(i);
    }

    //
    // set-up numerical model
    //
    help.zero();
    /*
     * nMethod -> setSparseMtrxAsComponent ( LinearEquationLhs , stiffnessMatrix) ;
     * nMethod -> setFloatArrayAsComponent ( LinearEquationRhs , &rhs) ;
     * nMethod -> setFloatArrayAsComponent ( LinearEquationSolution, &help);
     */
    //
    // call numerical model to solve arised problem
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    //nMethod -> solveYourselfAt(tStep);
    nMethod->solve(stiffnessMatrix, & rhs, & help);

    //
    // compute displacements, velocities and accelerations at t+dt
    //
    for ( i = 1; i <= neq; i++ ) {
        rhs.at(i) = a4 * help.at(i) + a5 *displacementVector.at(i) +
                    a6 *velocityVector.at(i) +
                    a7 *accelerationVector.at(i);
        displacementVector.at(i) += deltaT * velocityVector.at(i) +
                                    a9 *accelerationVector.at(i) +
                                    a10 *rhs.at(i);
        velocityVector.at(i) += a8 * ( accelerationVector.at(i) + rhs.at(i) );
        accelerationVector.at(i) = rhs.at(i);
    }
}

void
DIIDynamic :: giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                              CharType type, TimeStep *tStep, Domain *domain)
{
    // we don't directlt call element ->GiveCharacteristicMatrix() function, because some
    // engngm classes may require special modification of base types supported on
    // element class level

    if ( type == DIIModifiedStiffnessMatrix ) {
        Element *element;
        // IntArray loc ;
        FloatMatrix charMtrx1, charMtrx2;

        element = domain->giveElement(num);
        // element -> giveLocationArray (loc);
        element->giveCharacteristicMatrix(answer, StiffnessMatrix, tStep);
        element->giveCharacteristicMatrix(charMtrx2, MassMatrix, tStep);
        charMtrx2.times(this->a0);
        answer.add(charMtrx2);

        return;
    } else {
        StructuralEngngModel :: giveElementCharacteristicMatrix(answer, num, type, tStep, domain);
    }
}


void DIIDynamic :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    StructuralEngngModel :: updateYourself(stepN);
}


void
DIIDynamic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    static char dofchar[] = "dva";
    static ValueModeType dofmodes[] = {
        VM_Total, VM_Velocity, VM_Acceleration
    };

    iDof->printMultipleOutputAt(stream, atTime, dofchar, EID_MomentumBalance, dofmodes, 3);
}
} // end namespace oofem
