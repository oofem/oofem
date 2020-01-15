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

#include "sm/EngineeringModels/linearstability.h"
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
#include "dof.h"
#include "unknownnumberingscheme.h"
#include "outputmanager.h"
#include "eigenvectorprimaryfield.h"


#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_EngngModel(LinearStability);

LinearStability :: LinearStability(int i, EngngModel *master) : StructuralEngngModel(i, master),
    numberOfRequiredEigenValues(1),
    rtolv(1e-6),
    solverType(GES_SubspaceIt)
{
    numberOfSteps = 1;
    ndomains = 1;
}


NumericalMethod *LinearStability :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod = classFactory.createGeneralizedEigenValueSolver(solverType, this->giveDomain(1), this);
        if ( !nMethod ) {
            OOFEM_ERROR("solver creation failed");
        }
    }

    return nMethod.get();
}


SparseLinearSystemNM *LinearStability :: giveNumericalMethodForLinStaticProblem(TimeStep *tStep)
{
    if ( !nMethodLS ) {
        nMethodLS = classFactory.createSparseLinSolver(ST_Direct, this->giveDomain(1), this); ///@todo Support other solvers
        if ( !nMethodLS ) {
            OOFEM_ERROR("solver creation failed");
        }
    }

    return nMethodLS.get();
}


void
LinearStability :: initializeFrom(InputRecord &ir)
{
    //StructuralEngngModel::instanciateFrom(ir);
    IR_GIVE_FIELD(ir, numberOfRequiredEigenValues, _IFT_LinearStability_nroot);
    // numberOfSteps set artifficially to numberOfRequiredEigenValues
    // in order to allow
    // use restoreContext function for different eigenValues
    numberOfSteps = numberOfRequiredEigenValues;
    this->field = std::make_unique<EigenVectorPrimaryField>(this, 1, FT_Displacements, numberOfRequiredEigenValues + 1); // +1 for eq. solution

    IR_GIVE_FIELD(ir, rtolv, _IFT_LinearStability_rtolv);
    if ( rtolv < 1.e-12 ) {
        rtolv =  1.e-12;
    } else if ( rtolv > 0.01 ) {
        rtolv =  0.01;
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_LinearStability_stype);
    solverType = ( GenEigvalSolverType ) val;

    nMetaSteps = 0;

    suppressOutput = ir.hasField(_IFT_EngngModel_suppressOutput);

    if (suppressOutput) {
        printf("Suppressing output.\n");
    } else {
        if ( ( outputStream = fopen(this->dataOutputFileName.c_str(), "w") ) == NULL ) {
            OOFEM_ERROR("Can't open output file %s", this->dataOutputFileName.c_str());
        }

        fprintf(outputStream, "%s", PRG_HEADER);
        fprintf(outputStream, "\nStarting analysis on: %s\n", ctime(& this->startTime) );
        fprintf(outputStream, "%s\n", simulationDescription.c_str());
    }
}


int LinearStability :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % (this->numberOfRequiredEigenValues + 1); // +1 for eq. solution 
}


double LinearStability :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return field->giveUnknownValue(dof, mode, tStep);
}


TimeStep *LinearStability :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(istep, this, 1, 0., 0., counter);

    return currentStep.get();
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
    tStep->setNumber(0);
    tStep->setTime(0.0);

    // creates system of governing eq's and solves them at given time step
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );
    this->giveNumericalMethodForLinStaticProblem(tStep);

    // first assemble problem at current time step
    if ( !stiffnessMatrix ) {
        //
        // first step - solve linear static problem
        //
        stiffnessMatrix = classFactory.createSparseMtrx(SMT_Skyline); ///@todo Don't hardcode skyline matrix only
        stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    }

#ifndef LIN_STAB_COMPATIBILITY_MODE
    OOFEM_LOG_INFO("Assembling stiffness matrix\n");
    stiffnessMatrix->zero();
    this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );
#endif

    OOFEM_LOG_INFO("Assembling load\n");
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
    FloatArray displacementVector(neq), loadVector(neq);

    // Internal forces first, negated;
    field->update(VM_Total, tStep, displacementVector, EModelDefaultEquationNumbering());
    this->assembleVector( loadVector, tStep, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    loadVector.negated();

    this->assembleVector( loadVector, tStep, ExternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->updateSharedDofManagers(loadVector, EModelDefaultEquationNumbering(), ReactionExchangeTag);

    OOFEM_LOG_INFO("Solving linear static problem\n");
    nMethodLS->solve(*stiffnessMatrix, loadVector, displacementVector);
    // Initial displacements are stored at position 0; this is a bit of a hack. In the future, a cleaner approach of handling fields could be suitable,
    // but currently, it all converges down to the same giveUnknownComponent, so this is the easisest approach.
    field->update(VM_Total, tStep, displacementVector, EModelDefaultEquationNumbering());
    // terminate linear static computation (necessary, in order to compute stresses in elements).
    // Recompute for updated state:
    this->assembleVector( loadVector, tStep, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->terminateLinStatic( tStep );

    // Normal forces already known, proceed with linear stability
    stiffnessMatrix->zero();
    if ( !initialStressMatrix ) {
        initialStressMatrix = stiffnessMatrix->clone();
    } else {
        initialStressMatrix->zero();
    }

    OOFEM_LOG_INFO("Assembling stiffness  matrix\n");
    this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );

    OOFEM_LOG_INFO("Assembling  initial stress matrix\n");
    this->assemble( *initialStressMatrix, tStep, InitialStressMatrixAssembler(),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );
    initialStressMatrix->times(-1.0);

    //stiffnessMatrix->printYourself();
    //initialStressMatrix->printYourself();

    FloatMatrix eigVec(neq, numberOfRequiredEigenValues);
    eigVal.resize(numberOfRequiredEigenValues);
    eigVal.zero();

    OOFEM_LOG_INFO("Solving ...\n");
    nMethod->solve(*stiffnessMatrix, *initialStressMatrix, eigVal, eigVec, rtolv, numberOfRequiredEigenValues);
    this->field->updateAll(eigVec, EModelDefaultEquationNumbering());
}


void LinearStability :: updateYourself(TimeStep *tStep)
{ }


void
LinearStability :: terminateLinStatic(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);

    for ( auto &dman : domain->giveDofManagers() ) {
        dman->updateYourself(tStep);
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Updated nodes ", domain->giveNumberOfDofManagers())
#  endif

    for ( auto &elem : domain->giveElements() ) {
        elem->updateInternalState(tStep);
        elem->updateYourself(tStep);
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Updated Elements ", domain->giveNumberOfElements())
#  endif

#if 0
    // save context if required
    // default - save only if ALWAYS is set ( see cltypes.h )

    if ((domain->giveContextOutputMode() == COM_Always ) ||
        (domain->giveContextOutputMode() == COM_Required )) {
        this->saveContext(NULL);
    } else if (domain->giveContextOutputMode() == COM_UserDefined ) {
        if (tStep->giveNumber()%domain->giveContextOutputStep() == 0)
            this->saveContext(NULL);
    }
#endif
}


void LinearStability :: doStepOutput(TimeStep *tStep)
{
    if ( !suppressOutput ) {
        this->printOutputAt(this->giveOutputStream(), tStep);
        fflush( this->giveOutputStream() );
    }

    Domain *domain = this->giveDomain(1);
    // i = 0  represents the linear solution, which is followed by the eigen vectors starting at i = 1
    for ( int i = 0; i <= numberOfRequiredEigenValues; i++ ) {
        TimeStep step = *tStep;
        step.setTime( ( double ) i );
        step.setNumber(i);

        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself(&step);
        }

        exportModuleManager.doOutput(&step);
    }
}


void LinearStability :: printOutputAt(FILE *file, TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    if ( !domain->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return;
    }

    fprintf(file, "\nLinear Stability:");
    fprintf(file, "\nEigen Values are:\n-----------------\n");

    for ( int i = 1; i <= numberOfRequiredEigenValues; i++ ) {
        fprintf(file, "%15.8e ", eigVal.at(i) );
        if ( ( i % 5 ) == 0 ) {
            fprintf(file, "\n");
        }
    }

    fprintf(file, "\n\n");

    for ( int i = 0; i <= numberOfRequiredEigenValues; i++ ) {
        TimeStep step = *tStep;
        step.setTime( ( double ) i );
        step.setNumber(i);

        if ( i == 0 ) {
            fprintf(file, "\nLinear solution\n\n");
        } else {
            fprintf(file, "\nEigen vector no. %d, corresponding eigen value is %15.8e\n\n", i, eigVal.at(i));
        }

        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself(&step);
            dman->printOutputAt(file, &step);
        }

        if ( i == 0 ) {
            for ( auto &elem : domain->giveElements() ) {
                elem->printOutputAt(file, &step);
            }
            this->printReactionForces(&step, 1., file);
        }
    }
}


void LinearStability :: setActiveVector(int activeVector)
{
    this->giveCurrentStep()->setTime( ( double ) activeVector );
    this->giveCurrentStep()->setNumber( ( double ) activeVector );
}


void LinearStability :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralEngngModel :: saveContext(stream, mode);

    field->saveContext(stream);

    contextIOResultType iores;
    if ( ( iores = eigVal.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void LinearStability :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralEngngModel :: restoreContext(stream, mode);

    field->restoreContext(stream);

    contextIOResultType iores;
    if ( ( iores = eigVal.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


} // end namespace oofem
