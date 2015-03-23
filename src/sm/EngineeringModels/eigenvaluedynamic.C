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

#include "../sm/EngineeringModels/eigenvaluedynamic.h"
#include "timestep.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "exportmodulemanager.h"
#include "verbose.h"
#include "classfactory.h"
#include "datastream.h"
#include "geneigvalsolvertype.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dofmanager.h"
#include "dof.h"
#include "domain.h"
#include "unknownnumberingscheme.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_EngngModel(EigenValueDynamic);

NumericalMethod *EigenValueDynamic :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod.reset( classFactory.createGeneralizedEigenValueSolver(solverType, this->giveDomain(1), this) );
        if ( !nMethod ) {
            OOFEM_ERROR("solver creation failed");
        }
    }

    return nMethod.get();
}

IRResultType
EigenValueDynamic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    //EngngModel::instanciateFrom (ir);

    IR_GIVE_FIELD(ir, numberOfRequiredEigenValues, _IFT_EigenValueDynamic_nroot);

    // numberOfSteps set artificially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    // numberOfSteps = numberOfRequiredEigenValues;
    numberOfSteps = 1;

    IR_GIVE_FIELD(ir, rtolv, _IFT_EigenValueDynamic_rtolv);
    if ( rtolv < 1.e-12 ) {
        rtolv =  1.e-12;
    }

    if ( rtolv > 0.01 ) {
        rtolv =  0.01;
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EigenValueDynamic_stype);
    solverType = ( GenEigvalSolverType ) val;

    val = 0; //Default Skyline
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    return IRRT_OK;
}


double EigenValueDynamic :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, eigenvalue.
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }
#endif

    switch ( mode ) {
    case VM_Total:  // EigenVector
    case VM_Incremental:
        return eigVec.at( eq, ( int ) tStep->giveTargetTime() );

    default:
        OOFEM_ERROR("Unknown is of undefined type for this problem");
    }

    return 0.;
}


TimeStep *EigenValueDynamic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(istep, this, 1, ( double ) istep, 0., counter) );

    return currentStep.get();
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

        stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

        massMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        massMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

        this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        this->assemble( *massMatrix, tStep, MassMatrixAssembler(),
                       EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //
        // create resulting objects eigVec and eigVal
        //
        eigVec.resize(this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ), numberOfRequiredEigenValues);
        eigVec.zero();
        eigVal.resize(numberOfRequiredEigenValues);
        eigVal.zero();
    }

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );

    //
    // call numerical model to solve arised problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    nMethod->solve(*stiffnessMatrix, *massMatrix, eigVal, eigVec, rtolv, numberOfRequiredEigenValues);

    stiffnessMatrix.reset(NULL);
    massMatrix.reset(NULL);
}


void EigenValueDynamic :: updateYourself(TimeStep *tStep)
{ }


void EigenValueDynamic :: terminate(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    FILE *outputStream = this->giveOutputStream();

    // print loadcase header
    fprintf(outputStream, "\nOutput for time %.3e \n\n", 1.0);
    // print eigen values on output
    fprintf(outputStream, "\n\nEigen Values (Omega^2) are:\n-----------------\n");

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        fprintf( outputStream, "%15.8e ", eigVal.at(i) );
        if ( ( i % 5 ) == 0 ) {
            fprintf(outputStream, "\n");
        }
    }

    fprintf(outputStream, "\n\n");

    for ( int i = 1; i <=  numberOfRequiredEigenValues; i++ ) {
        fprintf(outputStream, "\nOutput for eigen value no.  %.3e \n", ( double ) i);
        fprintf( outputStream,
                "Printing eigen vector no. %d, corresponding eigen value is %15.8e\n\n",
                i, eigVal.at(i) );
        tStep->setTime( ( double ) i ); // we use time as intrinsic eigen value index

        if ( this->requiresUnknownsDictionaryUpdate() ) {
            for ( auto &dman : domain->giveDofManagers() ) {
                this->updateDofUnknownsDictionary(dman.get(), tStep);
            }
        }


        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself(tStep);
            dman->printOutputAt(outputStream, tStep);
        }
    }

    for ( int i = 1; i <=  numberOfRequiredEigenValues; i++ ) {
        // export using export manager
        tStep->setTime( ( double ) i ); // we use time as intrinsic eigen value index
        tStep->setNumber(i);
        exportModuleManager->doOutput(tStep);
    }
    fflush( this->giveOutputStream() );
    this->saveStepContext(tStep);
}


contextIOResultType EigenValueDynamic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file = NULL;

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

    if ( ( iores = eigVal.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( eigVec.storeYourself(*stream) ) != CIO_OK ) {
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
    FILE *file = NULL;

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

        if ( ( iores = eigVal.restoreYourself(*stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = eigVec.restoreYourself(*stream) ) != CIO_OK ) {
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

} // end namespace oofem
