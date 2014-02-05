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

// please activate or de-activate next line
//#define LIN_STAB_COMPATIBILITY_MODE

#include "linearstability.h"
#include "timestep.h"
#include "element.h"
#include "contextioerr.h"
#include "floatmatrix.h"
#include "verbose.h"
#include "floatarray.h"
#include "classfactory.h"
#include "datastream.h"
#include "exportmodulemanager.h"
#include "dofmanager.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_EngngModel(LinearStability);

NumericalMethod *LinearStability :: giveNumericalMethod(MetaStep *mStep)
{
    if ( nMethod ) {
        return nMethod;
    }

    nMethod = classFactory.createGeneralizedEigenValueSolver(solverType, this->giveDomain(1), this);
    if ( nMethod == NULL ) {
        _error("giveNumericalMethod:  solver creation failed");
    }

    return nMethod;
}

SparseLinearSystemNM *LinearStability :: giveNumericalMethodForLinStaticProblem(TimeStep *tStep)
{
    if ( nMethodLS ) {
        return nMethodLS;
    }

    nMethodLS = classFactory.createSparseLinSolver(ST_Direct, this->giveDomain(1), this); ///@todo Support other solvers
    if ( nMethodLS == NULL ) {
        _error("giveNumericalMethodForLinStaticProblem:  solver creation failed");
    }

    return nMethodLS;
}

IRResultType
LinearStability :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //StructuralEngngModel::instanciateFrom(ir);
    IR_GIVE_FIELD(ir, numberOfRequiredEigenValues, _IFT_LinearStability_nroot);
    // numberOfSteps set artifficially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    numberOfSteps = numberOfRequiredEigenValues;

    IR_GIVE_FIELD(ir, rtolv, _IFT_LinearStability_rtolv);
    if ( rtolv < 1.e-12 ) {
        rtolv =  1.e-12;
    }

    if ( rtolv > 0.01 ) {
        rtolv =  0.01;
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_LinearStability_stype);
    solverType = ( GenEigvalSolverType ) val;


    nMetaSteps = 0;

    return IRRT_OK;
}


double LinearStability :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, eigen value.
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }
#endif

    int activeVector = ( int ) tStep->giveTargetTime();
    switch ( mode ) {
    case VM_Total: // EigenVector
        if ( activeVector ) {
            return eigVec.at(eq, activeVector);
        }

        return displacementVector.at(eq);

    default:
        _error("giveUnknownComponent: Unknown is of undefined type for this problem");
    }

    return 0.;
}

TimeStep *LinearStability :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( previousStep != NULL ) {
        delete previousStep;
    }

    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, 0., 0., counter);
    // time and dt variables are set eq to 0 for staics - has no meaning

    return currentStep;
}

void LinearStability :: solveYourself()
{
    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    // update state according to new meta step
    this->giveNextStep();
    this->updateAttributes( this->giveCurrentMetaStep() );
    this->solveYourselfAt( this->giveCurrentStep() );
    this->terminate( this->giveCurrentStep() );
}


void LinearStability :: solveYourselfAt(TimeStep *tStep)
{
    //
    // creates system of governing eq's and solves them at given time step
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetatStepumber() ) );
    this->giveNumericalMethodForLinStaticProblem(tStep);

    // first assemble problem at current time step

    if ( tStep->giveNumber() == 1 ) {
        //
        // first step - solve linear static problem
        //
        stiffnessMatrix = classFactory.createSparseMtrx(SMT_Skyline); ///@todo Don't hardcode skyline matrix only
        stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );

        //
        // allocate space for displacement Vector
        //
        displacementVector.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
        //
        // allocate space for load vector
        //
        loadVector.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    }

#ifndef LIN_STAB_COMPATIBILITY_MODE
 #ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling stiffness matrix\n");
 #endif
    stiffnessMatrix->zero();
    this->assemble( stiffnessMatrix, tStep, EID_MomentumBalance, StiffnessMatrix,
                    EModelDefaultEquationNumbering(), this->giveDomain(1) );
#endif


#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling load\n");
#endif

    displacementVector.zero();
    loadVector.zero();

    // Internal forces first, negated;
    this->assembleVector( loadVector, tStep, EID_MomentumBalance, InternalForcesVector, VM_Total,
                          EModelDefaultEquationNumbering(), this->giveDomain(1) );
    loadVector.negated();

    this->assembleVector( loadVector, tStep, EID_MomentumBalance, ExternalForcesVector, VM_Total,
                          EModelDefaultEquationNumbering(), this->giveDomain(1) );
    
#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers(loadVector, EModelDefaultEquationNumbering(), ReactionExchangeTag);
#endif


    //
    // call numerical model to solve problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving linear static problem\n");
#endif

    nMethodLS->solve(stiffnessMatrix, & loadVector, & displacementVector);
    // terminate linear static computation (necessary, in order to compute stresses in elements).
    this->terminateLinStatic( this->giveCurrentStep() );
    /*
     * Normal forces already known, proceed with linear stability
     */

    stiffnessMatrix->zero();
    if ( initialStressMatrix == NULL ) {
        initialStressMatrix = stiffnessMatrix->GiveCopy();
    } else {
        initialStressMatrix->zero();
    }

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling stiffness  matrix\n");
#endif
    this->assemble( stiffnessMatrix, tStep, EID_MomentumBalance, StiffnessMatrix,
                    EModelDefaultEquationNumbering(), this->giveDomain(1) );
#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling  initial stress matrix\n");
#endif
    this->assemble( initialStressMatrix, tStep, EID_MomentumBalance, InitialStressMatrix,
                    EModelDefaultEquationNumbering(), this->giveDomain(1) );
    initialStressMatrix->times(-1.0);

    //  stiffnessMatrix->printYourself();
    //  initialStressMatrix->printYourself();

    //
    // create resulting objects eigVec and eigVal
    //
    eigVec.resize(this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ), numberOfRequiredEigenValues);
    eigVal.resize(numberOfRequiredEigenValues);
    eigVec.zero();
    eigVal.zero();

    //
    // call numerical model to solve arised problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    nMethod->solve(stiffnessMatrix, initialStressMatrix, & eigVal, & eigVec, rtolv, numberOfRequiredEigenValues);
    // compute eigen frequencies
    //for (i = 1; i<= numberOfRequiredEigenValues; i++)
    // eigVal.at(i) = sqrt(eigVal.at(i));
}

void LinearStability :: updateYourself(TimeStep *tStep)
{ }

void
LinearStability :: terminateLinStatic(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    FILE *File;
    int j;
    File = this->giveOutputStream();
    tStep->setTime(0.);

    fprintf( File, "\nOutput for time % .3e \n\n", tStep->giveTargetTime() );
    fprintf(File, "Linear static:\n\n");

    int nnodes = domain->giveNumberOfDofManagers();

    if ( requiresUnknownsDictionaryUpdate() ) {
        for ( j = 1; j <= nnodes; j++ ) {
            this->updateDofUnknownsDictionary(domain->giveDofManager(j), tStep);
        }
    }

    for ( j = 1; j <= nnodes; j++ ) {
        domain->giveDofManager(j)->updateYourself(tStep);
        domain->giveDofManager(j)->printOutputAt(File, tStep);
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Updated nodes & sides ", nnodes)
#  endif


    Element *elem;

    int nelem = domain->giveNumberOfElements();
    for ( j = 1; j <= nelem; j++ ) {
        elem = domain->giveElement(j);
        elem->updateInternalState(tStep);
        elem->updateYourself(tStep);
        elem->printOutputAt(File, tStep);
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Updated Elements ", nelem)
#  endif
    fprintf(File, "\n");
    /*
     * // save context if required
     * // default - save only if ALWAYS is set ( see cltypes.h )
     *
     * if ((domain->giveContextOutputMode() == COM_Always ) ||
     * (domain->giveContextOutputMode() == COM_Required )) {
     * this->saveContext(NULL);
     * }
     * else if (domain->giveContextOutputMode() == COM_UserDefined ) {
     * if (tStep->giveNumber()%domain->giveContextOutputStep() == 0)
     * this->saveContext(NULL);
     * }
     */

    this->printReactionForces(tStep, 1);
}


void LinearStability :: terminate(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    FILE *outputStream = this->giveOutputStream();
    //FloatArray *eigv;
    // Element *elem;
    int i, j;

    // print eigen values on output
    fprintf(outputStream, "\nLinear Stability:");
    fprintf(outputStream, "\nEigen Values are:\n-----------------\n");

    for ( i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        fprintf( outputStream, "%15.8e ", eigVal.at(i) );
        if ( ( i % 5 ) == 0 ) {
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "\n\n");

    int nnodes = domain->giveNumberOfDofManagers();

    for ( i = 1; i <=  numberOfRequiredEigenValues; i++ ) {
        fprintf(outputStream, "\nOutput for eigen value no.  % .3e \n", ( double ) i);
        fprintf( outputStream,
                 "Printing eigen vector no. %d, corresponding eigen value is %15.8e\n\n",
                 i, eigVal.at(i) );
        tStep->setTime( ( double ) i );

        if ( this->requiresUnknownsDictionaryUpdate() ) {
            for ( j = 1; j <= nnodes; j++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(j), tStep);
            }
        }


        for ( j = 1; j <= nnodes; j++ ) {
            domain->giveDofManager(j)->updateYourself(tStep);
            domain->giveDofManager(j)->printOutputAt(outputStream, tStep);
        }

        tStep->setNumber(i);
        exportModuleManager->doOutput(tStep);
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Updated nodes & sides ", nnodes)
#  endif
    fflush( this->giveOutputStream() );
    // save context if required
    this->saveStepContext(tStep);
}


contextIOResultType LinearStability :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    FILE *file = NULL;

    OOFEM_LOG_INFO("Storing context \n");
    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = StructuralEngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = displacementVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = eigVal.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = eigVec.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                        // ensure consistent records

    return CIO_OK;
}



contextIOResultType LinearStability :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    int activeVector, version;
    int istep = 1, iversion = 1;
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file = NULL;

    this->resolveCorrespondingtStepumber(activeVector, version, obj);
    if ( eigVal.isEmpty() ) { // not restored before
        if ( stream == NULL ) {
            if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
                THROW_CIOERR(CIO_IOERR); // override
            }

            stream = new FileDataStream(file);
            closeFlag = 1;
        }

        if ( ( iores = StructuralEngngModel :: restoreContext(stream, mode, ( void * ) & istep) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = displacementVector.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = eigVal.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = eigVec.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( closeFlag ) {
            fclose(file);
            delete stream;
            stream = NULL;
        }                                                    // ensure consistent records
    }

    //  if (istep > numberOfRequiredEigenValues) istep = numberOfRequiredEigenValues ;
    //  printf( "Restoring - corresponding index is %d, EigenValue is %lf\n",
    //     istep,eigVal.at(istep));
    //  setActiveVector (istep);
    if ( activeVector > numberOfRequiredEigenValues ) {
        activeVector = numberOfRequiredEigenValues;
    }

    OOFEM_LOG_INFO( "Restoring - corresponding index is %d, EigenValue is %f\n",
                    activeVector, eigVal.at(activeVector) );
    this->giveCurrentStep()->setTime( ( double ) activeVector );

    return CIO_OK;
}


void
LinearStability :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    iDof->printSingleOutputAt(stream, tStep, 'd', VM_Total);
}
} // end namespace oofem
