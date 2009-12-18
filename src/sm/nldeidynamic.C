/* $Header: /home/cvs/bp/oofem/sm/src/nldeidynamic.C,v 1.5.4.2 2004/05/14 13:45:45 bp Exp $ */
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
// file NlDEIDynamic.cc
//

#include "nldeidynamic.h"
#include "nlstructuralelement.h"
#include "nummet.h"
#include "ldltfact.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "elementside.h"
#include "dof.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <math.h>
#endif

#include "verbose.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {

#define ZERO_REL_MASS  1.E-6


NlDEIDynamic :: ~NlDEIDynamic()  { }

NumericalMethod *NlDEIDynamic :: giveNumericalMethod(TimeStep *atTime)
// only one has reason for DEIDynamic
//     - SolutionOfLinearEquations

{
    return NULL;  // not necessary here - diagonal matrix is used-simple inversion
}

IRResultType
NlDEIDynamic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralEngngModel :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, dumpingCoef, IFT_NlDEIDynamic_dumpcoef, "dumpcoef"); // C = dumpingCoef * M
    IR_GIVE_FIELD(ir, deltaT, IFT_NlDEIDynamic_deltat, "deltat"); // Macro

    drFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, drFlag, IFT_NlDEIDynamic_drflag, "drflag"); // Macro
    if ( drFlag ) {
        IR_GIVE_FIELD(ir, Tau, IFT_NlDEIDynamic_tau, "tau");
        IR_GIVE_FIELD(ir, pyEstimate, IFT_NlDEIDynamic_py, "py");
    }

    return IRRT_OK;
}





double NlDEIDynamic ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
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

    case VM_Incremental:
        return incrementOfDisplacementVector.at(eq);

    case VM_Velocity:
        return velocityVector.at(eq);

    case VM_Acceleration:
        return accelerationVector.at(eq);

    default:
        _error("giveUnknownComponent: Unknown is of undefined type for this problem");
    }

    return 0.;
}


TimeStep *NlDEIDynamic :: giveNextStep()
{
    int istep = 0;
    double totalTime = 0;
    StateCounterType counter = 1;

    delete previousStep;
    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTime() + deltaT;
        istep     = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, totalTime, deltaT, counter);
    // time and dt variables are set eq to 0 for staics - has no meaning

    return currentStep;
}




void NlDEIDynamic :: solveYourselfAt(TimeStep *tStep) {
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step
    Domain *domain = this->giveDomain(1);

    int nelem = domain->giveNumberOfElements();
    int nman  = domain->giveNumberOfDofManagers();

    IntArray loc;
    Element *element;
    DofManager *node;
    Dof *iDof;
    int nDofs, neq;
    int i, k, n, j, jj;
    double coeff, maxDt, maxOmi, maxOm = 0., maxOmEl, c1, c2;
    FloatArray internalForces;
    FloatArray diagonalStiffMtrx;

    neq = this->giveNumberOfEquations(EID_MomentumBalance);
    if ( initFlag ) {
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling mass matrix\n");
#endif

        //
        // first step  assemble mass Matrix
        //
        FloatMatrix charMtrx, charMtrx2;
	EModelDefaultEquationNumbering en;
        massMatrix.resize(neq);
        massMatrix.zero();
        for ( i = 1; i <= nelem; i++ ) {
            element = domain->giveElement(i);
            element->giveLocationArray(loc, EID_MomentumBalance, en);
            element->giveCharacteristicMatrix(charMtrx, LumpedMassMatrix, tStep);

#ifdef LOCAL_ZERO_MASS_REPLACEMENT
            element->giveCharacteristicMatrix(charMtrx2, StiffnessMatrix, tStep);
#endif
            //
            // assemble it manualy
            //
#ifdef DEBUG
            if ( ( n = loc.giveSize() ) != charMtrx.giveNumberOfRows() ) {
                _error("solveYourselfAt : dimension mismatch");
            }

#endif

            n = loc.giveSize();

#ifdef LOCAL_ZERO_MASS_REPLACEMENT
            maxOmEl = 0.;

            double maxElmass = -1.0;
            for ( j = 1; j <= n; j++ ) {
                maxElmass = max( maxElmass, charMtrx.at(j, j) );
            }

            if ( maxElmass <= 0.0 ) {
                _error("solveYourselfAt: Element with zero (or negative) lumped mass encountered\n");
            }

            for ( j = 1; j <= n; j++ ) {
                if ( charMtrx.at(j, j) > maxElmass * ZERO_REL_MASS ) {
                    maxOmi =  charMtrx2.at(j, j) / charMtrx.at(j, j);
                    maxOmEl = ( maxOmEl > maxOmi ) ? ( maxOmEl ) : ( maxOmi );
                }
            }

            maxOm = ( maxOm > maxOmEl ) ? ( maxOm ) : ( maxOmEl );

            for ( j = 1; j <= n; j++ ) {
                jj = loc.at(j);
                if ( ( jj ) && ( charMtrx.at(j, j) <= maxElmass * ZERO_REL_MASS ) ) {
                    charMtrx.at(j, j) = charMtrx2.at(j, j) / maxOmEl;
                }
            }

#endif

            for ( j = 1; j <= n; j++ ) {
                jj = loc.at(j);
                if ( jj ) {
                    massMatrix.at(jj) += charMtrx.at(j, j);
                }
            }
        }

#ifndef LOCAL_ZERO_MASS_REPLACEMENT
        // if init step - find minimun period of vibration in order to
        // determine maximal admisible time step
        // global variant
        //diagonalStiffMtrx.resize (neq); diagonalStiffMtrx.zero();
        for ( i = 1; i <= nelem; i++ ) {
            element = domain->giveElement(i);
            loc = element->giveLocationArray();
            charMtrx = element->GiveCharacteristicMatrix(StiffnessMatrix, tStep);
            n = loc->giveSize();
            for ( j = 1; j <= n; j++ ) {
                jj = loc->at(j);
                if ( jj ) {
                    diagonalStiffMtrx->at(jj) += charMtrx->at(j, j);
                }
            }

            delete charMtrx;
        }

        // find find minimun period of vibration
        // - global variant
        //

        double maxElmass = -1.0;
        for ( j = 1; j <= n; j++ ) {
            maxElmass = max( maxElmass, charMtrx.at(j, j) );
        }

        if ( maxElmass <= 0.0 ) {
            _error("solveYourselfAt: Element with zero (or negative) lumped mass encountered\n");
        }

        for ( j = 1; j <= neq; j++ ) {
            if ( massMatrix->at(j) > maxElmass * ZERO_REL_MASS ) {
                maxOmi =  diagonalStiffMtrx->at(j) / massMatrix->at(j);
                maxOm = ( maxOm > maxOmi ) ? ( maxOm ) : ( maxOmi );
            }
        }

        // set ZERO MASS members in massMatrix to value which corresponds to
        // maxOm
        // global variant
        //
        for ( i = 1; i <= neq; i++ ) {
            if ( massMatrix->at(i) <= maxElmass * ZERO_REL_MASS ) {
                massMatrix->at(i) = diagonalStiffMtrx->at(i) / maxOm;
            }
        }

        delete diagonalStiffMtrx;
        // end global variant

#endif

        if ( drFlag ) {
            // if Dynamic Relaxation assemble amplitude load vector
            loadRefVector.resize(neq);
            loadRefVector.zero();

            this->assembleVectorFromElements(loadRefVector, tStep, EID_MomentumBalance, ElementForceLoadVector, VM_Total,
					     EModelDefaultEquationNumbering(), domain);
            this->assembleVectorFromDofManagers(loadRefVector, tStep, EID_MomentumBalance, NodalLoadVector, VM_Total,
						EModelDefaultEquationNumbering(), domain);

            // compute the load vector norm pMp
            this->pMp = 0.0;
            for ( i = 1; i <= neq; i++ ) {
                pMp += loadRefVector.at(i) * loadRefVector.at(i) / massMatrix.at(i);
            }

            // solve for rate of loading process (parameter "c") (undamped system assumed)
            if ( dumpingCoef < 1.e-3 ) {
                c = 3.0 * this->pyEstimate / pMp / Tau / Tau;
            } else {
                c = this->pyEstimate * Tau * dumpingCoef * dumpingCoef * dumpingCoef / pMp /
                    ( -3.0 / 2.0 + dumpingCoef * Tau + 2.0 * exp(-dumpingCoef * Tau) - 0.5 * exp(-2.0 * dumpingCoef * Tau) );
            }
        }

        previousIncrementOfDisplacementVector.resize(neq);
        previousIncrementOfDisplacementVector.zero();

        initFlag = 0;
    }


    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        //
        // special init step - compute displacements at tstep 0
        //
        incrementOfDisplacementVector.resize(neq);
        incrementOfDisplacementVector.zero();
        displacementVector.resize(neq);
        displacementVector.zero();
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
                // all dofs are expected to be  DisplacementVector type.
                iDof  =  node->giveDof(k);
                if ( !iDof->isPrimaryDof() ) {
                    continue;
                }

                jj = iDof->__giveEquationNumber();
                if ( jj ) {
                    displacementVector.at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                    velocityVector.at(jj)     = iDof->giveUnknown(EID_MomentumBalance, VM_Velocity, tStep);
                    // accelerationVector = iDof->giveUnknown(AccelerartionVector,tStep) ;
                }
            }
        }

        //
        // set-up numerical model
        //

        // if init - try to determine the best deltaT
        // PI = 3.1415926535897932384626383279; // PI =  3.1415926535897931160E0
        maxDt = 2.0 / sqrt(maxOm);
        if ( deltaT > maxDt ) {
            // print reduced time step increment and minimum period Tmin
            OOFEM_LOG_RELEVANT("deltaT reduced to %e, Tmin is %e\n", maxDt, maxDt * M_PI);
            deltaT = maxDt;
            tStep->setTimeIncrement(deltaT);
        }

        for ( j = 1; j <= neq; j++ ) {
            incrementOfDisplacementVector.at(j) =
                velocityVector.at(j) * ( deltaT ); // becomes previous before used
            displacementVector.at(j) -= incrementOfDisplacementVector.at(j);
        }

        return;
    } // end of init step

    c1 = ( 1. / ( deltaT * deltaT ) );
    c2 = ( 1. / ( 2. * deltaT ) );

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling right hand side\n");
#endif

    for ( i = 1; i <= neq; i++ ) {
        previousIncrementOfDisplacementVector.at(i) = incrementOfDisplacementVector.at(i);
        displacementVector.at(i) += previousIncrementOfDisplacementVector.at(i);
    }

    tStep->incrementStateCounter();            // update solution state counter

    //
    // assembling the element part of load vector
    //
    this->giveInternalForces(internalForces, tStep);
    if ( !drFlag ) {
        loadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
        loadVector.zero();

        this->assembleVectorFromElements(loadVector, tStep, EID_MomentumBalance, ElementForceLoadVector, VM_Total,
					 EModelDefaultEquationNumbering(), domain);
        this->assembleVectorFromDofManagers(loadVector, tStep, EID_MomentumBalance, NodalLoadVector, VM_Total,
					    EModelDefaultEquationNumbering(), domain);
        //
        // assembling additional parts of right hand side
        //
        for ( k = 1; k <= neq; k++ ) {
            loadVector.at(k) -= internalForces.at(k);
        }
    } else {
        // dynamic relaxation
        // compute load factor
        double pt = 0.0;
        for ( k = 1; k <= neq; k++ ) {
            pt += internalForces.at(k) * loadRefVector.at(k) / massMatrix.at(k);
        }

        pt = pt / pMp;
        if ( dumpingCoef < 1.e-3 ) {
            pt += c * ( Tau - tStep->giveTime() ) / Tau;
        } else {
            pt += c * ( 1.0 - exp( dumpingCoef * ( tStep->giveTime() - Tau ) ) ) / dumpingCoef / Tau;
        }

        loadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
        for ( k = 1; k <= neq; k++ ) {
            loadVector.at(k) = pt * loadRefVector.at(k) - internalForces.at(k);
        }


        // compute relative error
        double err = 0.0;
        for ( k = 1; k <= neq; k++ ) {
            err = loadVector.at(k) * loadVector.at(k) / massMatrix.at(k);
        }

        err = err / ( pMp * pt * pt );
        OOFEM_LOG_RELEVANT("Relative error is %e, loadlevel is %e\n", err, pt);
    }

    for ( j = 1; j <= neq; j++ ) {
        coeff =  massMatrix.at(j);
        loadVector.at(j) +=
            coeff * ( c1 - dumpingCoef * c2 ) *
            previousIncrementOfDisplacementVector.at(j);
    }



    //
    // set-up numerical model
    //
    //
    // call numerical model to solve arised problem - done localy here
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTime() );
#endif

    for ( i = 1; i <= neq; i++ ) {
        incrementOfDisplacementVector.at(i) = loadVector.at(i) /
        ( massMatrix.at(i) * ( c1 + dumpingCoef * c2 ) );
        accelerationVector.at(i) = incrementOfDisplacementVector.at(i) -
        previousIncrementOfDisplacementVector.at(i);
        velocityVector.at(i) = incrementOfDisplacementVector.at(i) +
        previousIncrementOfDisplacementVector.at(i);
    }

    accelerationVector.times(c1);
    velocityVector.times(c2);

    // update nodes, elements, etc.
    this->updateYourself( this->giveCurrentStep() );
}


void NlDEIDynamic :: updateYourself(TimeStep *stepN)
{
    // updates internal state to reached one
    // all internal variables are directly updated by
    // numerical method - void function here


    //this->updateInternalState(stepN);
    StructuralEngngModel :: updateYourself(stepN);
}



void
NlDEIDynamic :: giveInternalForces(FloatArray &answer, TimeStep *stepN)
{
    // computes nodal representation of internal forces (real ones)
    // simply assembles contributions from each element in domain

    Element *element;
    IntArray loc;
    FloatArray charVec;
    Domain *domain = this->giveDomain(1);
    int nelems;
    EModelDefaultEquationNumbering en;

    answer.resize( displacementVector.giveSize() );
    answer.zero();

    nelems = domain->giveNumberOfElements();
    for ( int i = 1; i <= nelems; i++ ) {
        element = ( NLStructuralElement * ) domain->giveElement(i);
        element->giveLocationArray(loc, EID_MomentumBalance, en);
        element->giveCharacteristicVector(charVec, NodalInternalForcesVector, VM_Total, stepN);
        if ( charVec.containsOnlyZeroes() ) {
            continue;
        }

        answer.assemble(charVec, loc);
    }

    internalVarUpdateStamp = stepN->giveSolutionStateCounter();
    return;
}


contextIOResultType NlDEIDynamic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = StructuralEngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = incrementOfDisplacementVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = displacementVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->write(& deltaT, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                        // ensure consistent records

    return CIO_OK;
}



contextIOResultType NlDEIDynamic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    int closeFlag = 0;
    int istep, iversion;
    contextIOResultType iores;
    FILE *file;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);
    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    // save element context
    if ( ( iores = StructuralEngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = incrementOfDisplacementVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = displacementVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->read(& deltaT, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                       // ensure consistent records

    return CIO_OK;
}


void
NlDEIDynamic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    static char dofchar[] = "dva";
    static ValueModeType dofmodes[] = {
        VM_Total, VM_Velocity, VM_Acceleration
    };

    iDof->printMultipleOutputAt(stream, atTime, dofchar, EID_MomentumBalance, dofmodes, 3);
}

void
NlDEIDynamic :: terminate(TimeStep *tStep)
{
    StructuralEngngModel :: terminate(tStep);
    this->printReactionForces(tStep, 1);
}

} // end namespace oofem
