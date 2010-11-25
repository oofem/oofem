/* $Header: /home/cvs/bp/oofem/oofemlib/src/stationaryflow.C,v 1.16.4.1 2004/04/05 15:19:44 bp Exp $ */
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
// file stationaryflow.C
//


#include "stationaryflow.h"
#include "nummet.h"
#include "ldltfact.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "dof.h"
#include "datastream.h"
#include "contextioerr.h"
#include "unknownnumberingscheme.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#include "verbose.h"
#include "conTable.h"
#include "skyline.h"

namespace oofem {
NumericalMethod *StationaryFlow :: giveNumericalMethod(TimeStep *tStep)
// only one has reason for StationaryFlow
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



double StationaryFlow ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
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

    if ( chc != EID_ConservationEquation ) {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    switch ( mode ) {
    case VM_Total:
        if ( fluxVector.isNotEmpty() ) {
            return fluxVector.at(eq);
        } else {
            return 0.;
        }

    // return nMethod-> giveUnknownComponent (LinearEquationSolution, eq);

    default:
        _error("giveUnknownComponent: Unknown is of undefined ValueModeType for this problem");
    }

    return 0.;
}


TimeStep *StationaryFlow :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    StateCounterType counter = 1;
    delete previousStep;

    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 0, ( double ) istep, 0., counter);
    // time and dt variables are set eq to 0 for staics - has no meaning
    return currentStep;
}

void StationaryFlow :: solveYourselfAt(TimeStep *tStep) {
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step

    if ( tStep->giveNumber() == 1 ) {
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling conductivity matrix\n");
#endif

        //
        // first step  assemble conductivity Matrix
        //
        /*
         * IntArray* mht = this -> GiveBanWidthVector ();
         * conductivityMatrix = new Skyline ();
         * conductivityMatrix ->  checkSizeTowardsBanWidth (mht) ;
         * delete mht;
         */

        conductivityMatrix = new Skyline();
        conductivityMatrix->buildInternalStructure( this, 1, EID_ConservationEquation, EModelDefaultEquationNumbering() );

        this->assemble( conductivityMatrix, tStep, EID_ConservationEquation, ConductivityMatrix,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        this->assemble( conductivityMatrix, tStep, EID_ConservationEquation, BcLhsDueToConvection,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //
        // alocate space for fluxVector
        //
        //fluxVector = new FloatArray (this->giveNumberOfEquations());
        fluxVector.resize( this->giveNumberOfEquations(EID_ConservationEquation) );
        fluxVector.zero();
    }

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling load\n");
#endif

    //
    // assembling the element part of load vector
    //
    //loadVector = new FloatArray (this->giveNumberOfEquations());
    loadVector.resize( this->giveNumberOfEquations(EID_ConservationEquation) );
    loadVector.zero();
    this->assembleVectorFromElements( loadVector, tStep, EID_ConservationEquation, ElementPPDELoadVector, VM_Total,
                                     EModelDefaultEquationNumbering(), this->giveDomain(1) );

    //
    // assembling the nodal part of load vector
    //
    this->assembleVectorFromDofManagers( loadVector, tStep, EID_ConservationEquation, NodalLoadVector, VM_Total,
                                        EModelDefaultEquationNumbering(), this->giveDomain(1) );

    //
    // set-up numerical model
    //
    /*
     * nMethod -> setSparseMtrxAsComponent ( LinearEquationLhs ,conductivityMatrix) ;
     * nMethod -> setFloatArrayAsComponent ( LinearEquationRhs , &loadVector) ;
     * nMethod -> setFloatArrayAsComponent ( LinearEquationSolution, &fluxVector) ;
     */
    //
    // call numerical model to solve arised problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    //nMethod -> solveYourselfAt(tStep);
    nMethod->solve(conductivityMatrix, & loadVector, & fluxVector);
    //
    // update nodes, elements, etc.
    this->updateYourself( this->giveCurrentStep() );
}



contextIOResultType StationaryFlow :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( ( iores = fluxVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                          // ensure consistent records

    return CIO_OK;
}



contextIOResultType StationaryFlow :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( ( iores = fluxVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                           // ensure consistent records

    return CIO_OK;
}



void
StationaryFlow :: terminate(TimeStep *tStep)
{
    EngngModel :: terminate(tStep);
}

void
StationaryFlow :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    iDof->printSingleOutputAt(stream, atTime, 'f', EID_ConservationEquation, VM_Total);
}
} // end namespace oofem
