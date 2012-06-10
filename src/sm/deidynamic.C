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

#include "deidynamic.h"
#include "timestep.h"
#include "dofmanager.h"
#include "element.h"
#include "dof.h"
#include "verbose.h"
#include "mathfem.h"

namespace oofem {
#define ZERO_MASS  1.E-10   // unit dependent !!!!

DEIDynamic :: ~DEIDynamic() { }

NumericalMethod *DEIDynamic :: giveNumericalMethod(MetaStep *mStep)
// only one has reason for DEIDynamic
//     - SolutionOfLinearEquations

{
    return NULL;  // not necessary here - diagonal matrix is used-simple inversion
}


IRResultType
DEIDynamic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralEngngModel :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, dumpingCoef, IFT_DEIDynamic_dumpcoef, "dumpcoef"); // C = dumpingCoef * M
    IR_GIVE_FIELD(ir, deltaT, IFT_DEIDynamic_deltat, "deltat"); // Macro

    return IRRT_OK;
}


double DEIDynamic ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
                                           TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// in time t
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


TimeStep *DEIDynamic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    double totalTime = 0.;
    StateCounterType counter = 1;

    delete previousStep;
    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTargetTime() + deltaT;
        istep     = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, totalTime, deltaT, counter);
    // time and dt variables are set eq to 0 for staics - has no meaning

    return currentStep;
}


void DEIDynamic :: solveYourselfAt(TimeStep *tStep)
{
    //
    // creates system of governing eq's and solves them at given time step
    //
    // this is an explicit problem: we assemble governing equating at time t
    // and solution is obtained for time t+dt
    //
    // first assemble problem at current time step to obtain results in following
    // time step.
    // and then print results for this step also.
    // for first time step we need special start code
    Domain *domain = this->giveDomain(1);
    int nelem = domain->giveNumberOfElements();
    int nman = domain->giveNumberOfDofManagers();
    IntArray loc;
    Element *element;
    DofManager *node;
    Dof *iDof;
    int nDofs, neq;
    int i, k, n, j, jj, kk, init = 0;
    double coeff, maxDt, maxOmi, maxOm = 0., maxOmEl, c1, c2, c3;
    FloatMatrix charMtrx, charMtrx2;
    FloatArray previousDisplacementVector;


    neq = this->giveNumberOfEquations(EID_MomentumBalance);
    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        init = 1;
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling mass matrix\n");
#endif

        //
        // first step  assemble mass Matrix
        //

        massMatrix.resize(neq);
        massMatrix.zero();
        EModelDefaultEquationNumbering dn;
        for ( i = 1; i <= nelem; i++ ) {
            element = domain->giveElement(i);
            element->giveLocationArray(loc, EID_MomentumBalance, dn);
            element->giveCharacteristicMatrix(charMtrx,  LumpedMassMatrix, tStep);
            // charMtrx.beLumpedOf(fullCharMtrx);
            element->giveCharacteristicMatrix(charMtrx2, StiffnessMatrix, tStep);

            //
            // assemble it manually
            //
#ifdef DEBUG
            if ( ( n = loc.giveSize() ) != charMtrx.giveNumberOfRows() ) {
                _error("solveYourselfAt : dimension mismatch");
            }

#endif

            n = loc.giveSize();

            maxOmEl = 0.;
            for ( j = 1; j <= n; j++ ) {
                if ( charMtrx.at(j, j) > ZERO_MASS ) {
                    maxOmi =  charMtrx2.at(j, j) / charMtrx.at(j, j);
                    if ( init ) {
                        maxOmEl = ( maxOmEl > maxOmi ) ? ( maxOmEl ) : ( maxOmi );
                    }
                }
            }

            maxOm = ( maxOm > maxOmEl ) ? ( maxOm ) : ( maxOmEl );

            for ( j = 1; j <= n; j++ ) {
                jj = loc.at(j);
                if ( ( jj ) && ( charMtrx.at(j, j) <= ZERO_MASS ) ) {
                    charMtrx.at(j, j) = charMtrx2.at(j, j) / maxOmEl;
                }
            }

            for ( j = 1; j <= n; j++ ) {
                jj = loc.at(j);
                if ( jj ) {
                    massMatrix.at(jj) += charMtrx.at(j, j);
                }
            }
        }

        // if init - try to determine the best deltaT
        if ( init ) {
            maxDt = 2 / sqrt(maxOm);
            if ( deltaT > maxDt ) {
                OOFEM_LOG_RELEVANT("DEIDynamic: deltaT reduced to %e\n", maxDt);
                deltaT = maxDt;
                tStep->setTimeIncrement(deltaT);
            }
        }


        //
        // special init step - compute displacements at tstep 0
        //
        displacementVector.resize(neq);
        displacementVector.zero();
        nextDisplacementVector.resize(neq);
        nextDisplacementVector.zero();
        velocityVector.resize(neq);
        velocityVector.zero();
        accelerationVector.resize(neq);
        accelerationVector.zero();


        for ( j = 1; j <= nman; j++ ) {
            node = domain->giveDofManager(j);
            nDofs = node->giveNumberOfDofs();

            for ( k = 1; k <= nDofs; k++ ) {
                // ask for initial values obtained from
                // bc (boundary conditions) and ic (initial conditions)
                // now we are setting initial cond. for step -1.
                iDof  =  node->giveDof(k);
                if ( !iDof->isPrimaryDof() ) {
                    continue;
                }

                jj = iDof->__giveEquationNumber();
                if ( jj ) {
                    nextDisplacementVector.at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                    // become displacementVector after init
                    velocityVector.at(jj)     = iDof->giveUnknown(EID_MomentumBalance, VM_Velocity, tStep);
                    // accelerationVector = iDof->giveUnknown(AccelerartionVector,tStep) ;
                }
            }
        }

        for ( j = 1; j <= neq; j++ ) {
            nextDisplacementVector.at(j) -= velocityVector.at(j) * ( deltaT );
        }

        return;
    } // end of init step

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling right hand side\n");
#endif


    c1 = ( 1. / ( deltaT * deltaT ) );
    c2 = ( 1. / ( 2. * deltaT ) );
    c3 = ( 2. / ( deltaT * deltaT ) );

    previousDisplacementVector = displacementVector;
    displacementVector         = nextDisplacementVector;

    //
    // assembling the element part of load vector
    //
    loadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
    loadVector.zero();
    this->assembleVector(loadVector, tStep, EID_MomentumBalance, ExternalForcesVector,
                         VM_Total, EModelDefaultEquationNumbering(), domain);

    //
    // assembling additional parts of right hand side
    //
    EModelDefaultEquationNumbering dn;
    for ( i = 1; i <= nelem; i++ ) {
        element = domain->giveElement(i);
        element->giveLocationArray(loc, EID_MomentumBalance, dn);
        element->giveCharacteristicMatrix(charMtrx, StiffnessMatrix, tStep);
        n = loc.giveSize();
        for ( j = 1; j <= n; j++ ) {
            jj = loc.at(j);
            if ( jj ) {
                for ( k = 1; k <= n; k++ ) {
                    kk = loc.at(k);
                    if ( kk ) {
                        loadVector.at(jj) -= charMtrx.at(j, k) * displacementVector.at(kk);
                    }
                }

                //
                // if init step - find minimum period of vibration in order to
                // determine maximal admissible time step
                //
                //maxOmi =  charMtrx.at(j,j)/massMatrix.at(jj) ;
                //if (init) maxOm = (maxOm > maxOmi) ? (maxOm) : (maxOmi) ;
            }
        }
    }



    for ( j = 1; j <= neq; j++ ) {
        coeff =  massMatrix.at(j);
        loadVector.at(j) += coeff * c3 * displacementVector.at(j) -
                            coeff * ( c1 - dumpingCoef * c2 ) *
                            previousDisplacementVector.at(j);
    }

    //
    // set-up numerical model
    //
    /* it is not necessary to call numerical method
     * approach used here is not good, but effective enough
     * inverse of diagonal mass matrix is done here
     */
    //
    // call numerical model to solve arose problem - done locally here
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif
    double prevD;

    for ( i = 1; i <= neq; i++ ) {
        prevD = previousDisplacementVector.at(i);
        nextDisplacementVector.at(i) = loadVector.at(i) /
                                       ( massMatrix.at(i) * ( c1 + dumpingCoef * c2 ) );
        velocityVector.at(i) = nextDisplacementVector.at(i) - prevD;
        accelerationVector.at(i) =
            nextDisplacementVector.at(i) -
            2. * displacementVector.at(i) + prevD;
    }

    accelerationVector.times(c1);
    velocityVector.times(c2);
}


void DEIDynamic :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    StructuralEngngModel :: updateYourself(stepN);
}


void DEIDynamic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    static char dofchar[] = "dva";
    static ValueModeType dofmodes[] = {
        VM_Total, VM_Velocity, VM_Acceleration
    };

    iDof->printMultipleOutputAt(stream, atTime, dofchar, EID_MomentumBalance, dofmodes, 3);
}
} // end namespace oofem
