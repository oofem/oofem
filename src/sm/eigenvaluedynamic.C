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

#include "eigenvaluedynamic.h"
#include "timestep.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "exportmodulemanager.h"
#include "verbose.h"
#include "usrdefsub.h"
#include "datastream.h"
#include "geneigvalsolvertype.h"
#include "contextioerr.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {

NumericalMethod *EigenValueDynamic :: giveNumericalMethod(MetaStep *mStep)
{
    if ( nMethod ) {
        return nMethod;
    }

    nMethod = CreateUsrDefGeneralizedEigenValueSolver(solverType, 1, this->giveDomain(1), this);
    if ( nMethod == NULL ) {
        _error("giveNumericalMethod:  solver creation failed");
    }

    return nMethod;
}

IRResultType
EigenValueDynamic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    //EngngModel::instanciateFrom (ir);

    IR_GIVE_FIELD(ir, numberOfRequiredEigenValues, IFT_EigenValueDynamic_nroot, "nroot"); // Macro

    // numberOfSteps set artificially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    // numberOfSteps = numberOfRequiredEigenValues;
    numberOfSteps = 1;

    IR_GIVE_FIELD(ir, rtolv, IFT_EigenValueDynamic_rtolv, "rtolv"); // Macro
    if ( rtolv < 1.e-12 ) {
        rtolv =  1.e-12;
    }

    if ( rtolv > 0.01 ) {
        rtolv =  0.01;
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_EigenValueDynamic_stype, "stype"); // Macro
    solverType = ( GenEigvalSolverType ) val;

    val = 0; //Default Skyline
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_EigenValueDynamic_smtype, "smtype");  // Macro
    sparseMtrxType = ( SparseMtrxType ) val;

    return IRRT_OK;
}


double EigenValueDynamic ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
                                                  TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, eigenvalue.
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }

    if ( chc != EID_MomentumBalance ) {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    switch ( mode ) {
    case VM_Total:  // EigenVector
        return eigVec.at( eq, ( int ) tStep->giveTargetTime() );

    default:
        _error("giveUnknownComponent: Unknown is of undefined type for this problem");
    }

    return 0.;
}


double EigenValueDynamic ::  giveUnknownComponent(UnknownType chc, ValueModeType mode,
                                                  TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, eigenvalue.
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }

    if ( chc == EigenValue ) {
        return eigVal.at(eq);
    } else if ( !( chc == EigenVector ) ) {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    switch ( mode ) {
    case VM_Total:  // EigenVector
        return eigVec.at( eq, ( int ) tStep->giveTargetTime() );

    default:
        _error("giveUnknownComponent: Unknown is of undefined type for this problem");
    }

    return 0.;
}


TimeStep *EigenValueDynamic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    StateCounterType counter = 1;

    delete previousStep;
    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, ( double ) istep, 0., counter);
    // time and dt variables are set eq to 0 for statics - has no meaning

    return currentStep;
}


void EigenValueDynamic :: solveYourselfAt(TimeStep *tStep)
{
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling stiffness and mass matrices\n");
#endif

    if ( tStep->giveNumber() == 1 ) {
        //
        // first step  assemble stiffness Matrix
        //

        stiffnessMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
        stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );

        massMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
        massMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );

        this->assemble( stiffnessMatrix, tStep, EID_MomentumBalance, StiffnessMatrix,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        this->assemble( massMatrix, tStep, EID_MomentumBalance, MassMatrix,
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //
        // create resulting objects eigVec and eigVal
        //
        eigVec.resize(this->giveNumberOfEquations(EID_MomentumBalance), numberOfRequiredEigenValues);
        eigVec.zero();
        eigVal.resize(numberOfRequiredEigenValues);
        eigVal.zero();
    }

    //
    // set-up numerical model
    //
    this->giveNumericalMethod(  this->giveMetaStep( tStep->giveMetaStepNumber() ) );

    //
    // call numerical model to solve arised problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    nMethod->solve(stiffnessMatrix, massMatrix, & eigVal, & eigVec, rtolv, numberOfRequiredEigenValues);

    delete stiffnessMatrix;
    delete massMatrix;
    stiffnessMatrix = massMatrix = NULL;
}


void EigenValueDynamic :: updateYourself(TimeStep *stepN)
{ }


void EigenValueDynamic :: terminate(TimeStep *stepN)
{
    Domain *domain = this->giveDomain(1);
    FILE *outputStream = this->giveOutputStream();
    int i, j;

    // print loadcase header
    fprintf(outputStream, "\nOutput for time % .3e \n\n", 1.0);
    // print eigen values on output
    fprintf(outputStream, "\n\nEigen Values (Omega^2) are:\n-----------------\n");

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
        stepN->setTime( ( double ) i ); // we use time as intrinsic eigen value index

        if ( this->requiresUnknownsDictionaryUpdate() ) {
            for ( j = 1; j <= nnodes; j++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(j), stepN);
            }
        }


        for ( j = 1; j <= nnodes; j++ ) {
            domain->giveDofManager(j)->updateYourself(stepN);
            domain->giveDofManager(j)->printOutputAt(outputStream, stepN);
        }
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Updated nodes & sides ", nnodes)
#  endif

    for ( i = 1; i <=  numberOfRequiredEigenValues; i++ ) {
        // export using export manager
        stepN->setTime( ( double ) i ); // we use time as intrinsic eigen value index
        stepN->setNumber(i);
        exportModuleManager->doOutput(stepN);
    }

    this->saveStepContext(stepN);
}


contextIOResultType EigenValueDynamic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = eigVal.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( eigVec.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    } // ensure consistent records

    return CIO_OK;
}


contextIOResultType EigenValueDynamic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    int closeFlag = 0;
    int activeVector = this->resolveCorrespondingEigenStepNumber(obj);
    int istep = 1, iversion = 0;
    contextIOResultType iores;
    FILE *file;

    if ( restoreFlag == 0 ) { // not restored before
        if ( stream == NULL ) {
            if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
                THROW_CIOERR(CIO_IOERR); // override
            }

            stream = new FileDataStream(file);
            closeFlag = 1;
        }

        // save element context

        if ( ( iores = EngngModel :: restoreContext(stream, mode, ( void * ) & istep) ) != CIO_OK ) {
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
        } // ensure consistent records

    }

    if ( activeVector > numberOfRequiredEigenValues ) {
        activeVector = numberOfRequiredEigenValues;
    }

    OOFEM_LOG_INFO( "Restoring - corresponding index is %d, EigenValue is %f\n", activeVector, eigVal.at(activeVector) );
    this->giveCurrentStep()->setTime( ( double ) activeVector );
    this->restoreFlag = 1;

    return CIO_OK;
}


int EigenValueDynamic :: resolveCorrespondingEigenStepNumber(void *obj)
{
    //
    // returns corresponding eigen step number
    //
    if ( obj == NULL ) {
        return 1;
    }

    int *istep = ( int * ) obj;

    if ( * istep > numberOfRequiredEigenValues ) {
        return numberOfRequiredEigenValues;
    }

    if ( * istep <= 0 ) {
        return 1;
    }

    return * istep;
}


void
EigenValueDynamic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    iDof->printSingleOutputAt(stream, atTime, 'd', EID_MomentumBalance, VM_Total);
}

#ifdef __SLEPC_MODULE
void
EigenValueDynamic :: initPetscContexts()
{
    PetscContext *petscContext;

    int i;
    petscContextList->growTo(ndomains);
    for ( i = 0; i < this->ndomains; i++ ) {
        petscContext =  new PetscContext(this, EID_MomentumBalance);
        petscContextList->put(i + 1, petscContext);
    }
}
#endif


} // end namespace oofem
