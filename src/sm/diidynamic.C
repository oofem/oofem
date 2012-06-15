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

#include "diidynamic.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "contextioerr.h"
#include "datastream.h"
#include "verbose.h"
#include "structuralelement.h"
#include "structuralelementevaluator.h"
#include "usrdefsub.h"

namespace oofem {
DIIDynamic :: DIIDynamic(int i, EngngModel *_master) : StructuralEngngModel(i, _master),
    loadVector(), previousLoadVector(), rhs(), displacementVector(), velocityVector(), accelerationVector(),
    previousDisplacementVector(), previousVelocityVector(), previousAccelerationVector(), help()
{
    initFlag = true;
    stiffnessMatrix = NULL;
    ndomains = 1;
    nMethod = NULL;

    initialTimeDiscretization = TD_Unspecified;

#ifdef __PARALLEL_MODE
    commMode = ProblemCommMode__NODE_CUT;
    nonlocalExt = 0;
    communicator = nonlocCommunicator = NULL;
    commBuff = NULL;
#endif
}

DIIDynamic :: ~DIIDynamic()
{
    if ( stiffnessMatrix ) {
        delete stiffnessMatrix;
    }

    if ( nMethod ) {
        delete nMethod;
    }
}

NumericalMethod *DIIDynamic :: giveNumericalMethod(MetaStep *mStep)
// Only one has reason for DIIDynamic
// - SolutionOfLinearEquations
{
#ifdef __PARALLEL_MODE

    if ( nMethod ) {
        return nMethod;
    }

    if ( ( solverType == ST_Petsc ) || ( solverType == ST_Feti ) ) {
        if ( nMethod ) {
            return nMethod;
        }

        nMethod = CreateUsrDefSparseLinSolver(solverType, 1, this->giveDomain(1), this);
    }

    if ( nMethod == NULL ) {
        _error("giveNumericalMethod: linear solver creation failed (unknown type or no parallel support)");
    }

    return nMethod;

#endif

    if ( nMethod ) {
        return nMethod;
    }

    nMethod = CreateUsrDefSparseLinSolver(solverType, 1, this->giveDomain(1), this);
    if ( nMethod == NULL ) {
        _error("giveNumericalMethod: linear solver creation failed");
    }

    return nMethod;
}


IRResultType
DIIDynamic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralEngngModel :: initializeFrom(ir);
    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_DIIDynamic_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_DIIDynamic_smtype, "smtype"); // Macro
    sparseMtrxType = ( SparseMtrxType ) val;

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
// Returns unknown quantity displacement, velocity or acceleration of equation eq.
// This function translates this request to numerical method language.
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
    double totalTime = deltaT;
    StateCounterType counter = 1;

    delete previousStep;
    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTargetTime() + deltaT;
        istep     = currentStep->giveNumber() + 1;
        counter   = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep  = new TimeStep(istep, this, 1, totalTime, deltaT, counter);

    return currentStep;
}

void DIIDynamic :: solveYourself()
{
    StructuralEngngModel :: solveYourself();
}

void DIIDynamic :: solveYourselfAt(TimeStep *tStep) {
    // Determine the constants.
    this->determineConstants(tStep);

    Domain *domain = this->giveDomain(1);
    int neq =  this->giveNumberOfEquations(EID_MomentumBalance);

    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        TimeStep *stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0,
                                                 -deltaT, deltaT, 0);
        //
        // Determine initial displacements, velocities and accelerations.
        //
        loadVector.resize(neq);
        loadVector.zero();
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
                //
                // Ask for initial values obtained from boundary conditions and initial conditions.
                //
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

        delete stepWhenIcApply;
    }   // End of initialization.

    if ( initFlag ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling stiffness matrix\n");
#endif
        stiffnessMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
        if ( stiffnessMatrix == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );

        this->assemble(stiffnessMatrix, tStep, EID_MomentumBalance, ModifiedStiffnessMatrix,
                       EModelDefaultEquationNumbering(), domain);

        help.resize(neq);
        help.zero();

        previousDisplacementVector.resize(neq);
        previousVelocityVector.resize(neq);
        previousAccelerationVector.resize(neq);
        previousLoadVector.resize(neq);

        for ( int i = 1; i <= neq; i++ ) {
            previousDisplacementVector.at(i) = displacementVector.at(i);
            previousVelocityVector.at(i)     = velocityVector.at(i);
            previousAccelerationVector.at(i) = accelerationVector.at(i);
            previousLoadVector.at(i)         = loadVector.at(i);
        }

        initFlag = false;
    }

#ifdef VERBOSE
    OOFEM_LOG_DEBUG("Assembling load\n");
#endif

    this->assembleLoadVector(loadVector, domain, VM_Total, tStep);

    //
    // Assemble effective load vector (rhs).
    //

    for ( int i = 1; i <= neq; i++ ) {
        help.at(i) = a1 * previousDisplacementVector.at(i) +
                     a2 * previousVelocityVector.at(i)     +
                     a3 * previousAccelerationVector.at(i);
    }

    this->timesMassMtrx(help, rhs, domain, tStep);
    help.zero();

    for ( int i = 1; i <= neq; i++ ) {
        rhs.at(i) += previousLoadVector.at(i) +
                     Psi * ( loadVector.at(i) - previousLoadVector.at(i) );
    }

    //
    // Call numerical model to solve arised problem.
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    nMethod->solve(stiffnessMatrix, & rhs, & help);

    for ( int i = 1; i <= neq; i++ ) {
        accelerationVector.at(i) = a4 * help.at(i) + 
                                   a5 * previousDisplacementVector.at(i) +
                                   a6 * previousVelocityVector.at(i) +
                                   a7 * previousAccelerationVector.at(i);

        displacementVector.at(i) = previousDisplacementVector.at(i) +
                                   deltaT * previousVelocityVector.at(i) +
                                   a9 * previousAccelerationVector.at(i) +
                                   a10 * accelerationVector.at(i);

        velocityVector.at(i) = previousVelocityVector.at(i) +
                               a8 * ( previousAccelerationVector.at(i) + accelerationVector.at(i) );
    }
}

void
DIIDynamic :: giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                              CharType type, TimeStep *tStep, Domain *domain)
{
    // We don't directly call element ->GiveCharacteristicMatrix() function, because some
    // engngm classes may require special modification of base types supported on
    // element class level.

    if ( type == ModifiedStiffnessMatrix ) {
        Element *element;
        FloatMatrix charMtrx1, charMtrx2;

        element = domain->giveElement(num);
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
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);

    for ( int i = 1; i <= neq; i++ ) {
        previousDisplacementVector.at(i) = displacementVector.at(i);
        previousVelocityVector.at(i)     = velocityVector.at(i);
        previousAccelerationVector.at(i) = accelerationVector.at(i);
        previousLoadVector.at(i)         = loadVector.at(i);
    }

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

void
DIIDynamic :: timesMassMtrx(FloatArray &vec, FloatArray &answer, Domain *domain, TimeStep *tStep)
{
    int nelem = domain->giveNumberOfElements();
    int neq   = this->giveNumberOfEquations(EID_MomentumBalance);
    int i, j, k, jj, kk, n;
    FloatMatrix charMtrx;
    IntArray loc;
    Element *element;
    EModelDefaultEquationNumbering en;

    answer.resize(neq);
    answer.zero();

    for ( i = 1; i <= nelem; i++ ) {
        element = domain->giveElement(i);
#ifdef __PARALLEL_MODE
        // Skip remote elements.
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif

        element->giveLocationArray(loc, EID_MomentumBalance, en);
        element->giveCharacteristicMatrix(charMtrx, MassMatrix, tStep);

#ifdef DEBUG
        if ( ( n = loc.giveSize() ) != charMtrx.giveNumberOfRows() ) {
            _error("solveYourselfAt : dimension mismatch");
        }

#endif

        n = loc.giveSize();
        for ( j = 1; j <= n; j++ ) {
            jj = loc.at(j);
            if ( jj ) {
                for ( k = 1; k <= n; k++ ) {
                    kk = loc.at(k);
                    if ( kk ) {
                        answer.at(jj) += charMtrx.at(j, k) * vec.at(kk);
                    }
                }
            }
        }
    }
}

void
DIIDynamic :: assembleLoadVector(FloatArray &_loadVector, Domain *domain, ValueModeType mode, TimeStep *tStep)
{
    _loadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
    _loadVector.zero();

    this->assembleVector( _loadVector, tStep, EID_MomentumBalance, ExternalForcesVector, mode,
                          EModelDefaultEquationNumbering(), domain);
}

void
DIIDynamic :: determineConstants(TimeStep *tStep)
{
    deltaT = tStep->giveTimeIncrement();
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
    } else { // Wilson Method
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
}

int
DIIDynamic :: checkConsistency()
{
    // Check internal consistency
    int i, nelem;
    Element *ePtr;
    StructuralElement *sePtr;
    StructuralElementEvaluator *see;
    Domain *domain = this->giveDomain(1);

    nelem = domain->giveNumberOfElements();

    // check for proper element type
    for ( i = 1; i <= nelem; i++ ) {
        ePtr = domain->giveElement(i);
        sePtr = dynamic_cast< StructuralElement * >( ePtr );
        see   = dynamic_cast< StructuralElementEvaluator * >( ePtr );

        if ( ( sePtr == NULL ) && ( see == NULL ) ) {
            _warning2("checkConsistency: element %d has no Structural support", i);
            return 0;
        }
    }

    EngngModel :: checkConsistency();

    return 1;
}

contextIOResultType DIIDynamic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR); // Override.
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
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

    if ( ( iores = loadVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    return CIO_OK;
}

contextIOResultType DIIDynamic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    int closeFlag = 0;
    int istep, iversion;
    contextIOResultType iores;
    FILE *file;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);
    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR); // Override.
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
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

    if ( ( iores = loadVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    return CIO_OK;
}
} // end namespace oofem
