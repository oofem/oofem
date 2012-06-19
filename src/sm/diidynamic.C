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
                                                       previousDisplacementVector(), previousVelocityVector(), previousAccelerationVector(), previousIncrementOfDisplacement(), help()
{
    initFlag = true;
    stiffnessMatrix = NULL;
    ndomains = 1;
    nMethod = NULL;

    initialTimeDiscretization = TD_ThreePointBackward;

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

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_DIIDynamic_ddtScheme, "ddtscheme");

    eta = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, eta, IFT_DIIDynamic_eta, "eta");

    delta = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, delta, IFT_DIIDynamic_delta, "delta");

    initialTimeDiscretization = ( TimeDiscretizationType ) val;

    gamma = 0.5; beta = 0.25; // Default Newmark parameters.
    theta = 1.37; // Default Wilson parameter.
    if ( initialTimeDiscretization == TD_Newmark ) {
        OOFEM_LOG_INFO( "Selecting Newmark-beta metod\n" );
        IR_GIVE_OPTIONAL_FIELD(ir, gamma, IFT_DIIDynamic_gamma, "gamma");
        IR_GIVE_OPTIONAL_FIELD(ir, beta, IFT_DIIDynamic_beta, "beta");
    } else if ( initialTimeDiscretization == TD_TwoPointBackward ) {
        OOFEM_LOG_INFO( "Selecting Two-point Backward Euler method\n" );
    } else if ( initialTimeDiscretization == TD_ThreePointBackward ) {
        OOFEM_LOG_INFO( "Selecting Three-point Backward Euler metod\n" );
    } else if ( initialTimeDiscretization == TD_Wilson ) {
        OOFEM_LOG_INFO( "Selecting Wilson-theta metod\n" );
        IR_GIVE_OPTIONAL_FIELD(ir, theta, IFT_DIIDynamic_psi, "theta"); // Macro
        if ( theta < 1.37 ) {
            OOFEM_LOG_WARNING("Found theta < 1.37. Performing correction, theta = 1.37");
            theta = 1.37;
        }
    } else {
        _error("NonLinearDynamic: Time-stepping scheme not found!\n");
    }

    IR_GIVE_FIELD(ir, deltaT, IFT_DIIDynamic_deltat, "deltat"); // Macro

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
    TimeDiscretizationType td = initialTimeDiscretization;

    delete previousStep;
    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTargetTime() + deltaT;
        istep     = currentStep->giveNumber() + 1;
        counter   = currentStep->giveSolutionStateCounter() + 1;
        td        = currentStep->giveTimeDiscretization();
        if ( ( currentStep->giveNumber() == giveNumberOfFirstStep() )
            && ( initialTimeDiscretization == TD_ThreePointBackward ) ) {
            td = TD_ThreePointBackward;
        }
    }

    previousStep = currentStep;

    currentStep  = new TimeStep(istep, this, 1, totalTime, deltaT, counter, td);

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
        previousIncrementOfDisplacement.resize(neq);
        previousIncrementOfDisplacement.zero();

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

        this->assemble(stiffnessMatrix, tStep, EID_MomentumBalance, EffectiveStiffnessMatrix,
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

        initFlag = 0;
    }

    if ( ( previousStep != NULL ) && ( tStep->giveTimeDiscretization() != previousStep->giveTimeDiscretization() ) ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling stiffness matrix\n");
#endif
        stiffnessMatrix->zero();
        this->assemble(stiffnessMatrix, tStep, EID_MomentumBalance, EffectiveStiffnessMatrix,
                       EModelDefaultEquationNumbering(), domain);
    }

#ifdef VERBOSE
    OOFEM_LOG_DEBUG("Assembling load\n");
#endif

    this->assembleLoadVector(loadVector, domain, VM_Total, tStep);

    //
    // Assemble effective load vector (rhs).
    //

    for ( int i = 1; i <= neq; i++ ) {
        help.at(i) = a0 * previousDisplacementVector.at(i)
            + a2 * previousVelocityVector.at(i)
            + a3 * previousAccelerationVector.at(i)
            + eta * ( a1 * previousDisplacementVector.at(i)
                      + a4 * previousVelocityVector.at(i)
                      + a5 * previousAccelerationVector.at(i)
                      + a11 * previousIncrementOfDisplacement.at(i) );
    }

    this->timesMtrx(help, rhs, MassMatrix, domain, tStep);
    help.zero();

    if ( delta != 0 ) {
        for ( int i = 1; i <= neq; i++ ) {
            help.at(i) = delta * ( a1 * previousDisplacementVector.at(i)
                                   + a4 * velocityVector.at(i)
                                   + a5 * accelerationVector.at(i)
                                   + a6 * previousIncrementOfDisplacement.at(i) );
        }
        this->timesMtrx(help, rhs2, StiffnessMatrix, domain, tStep);
        help.zero();
        for ( int i = 1; i <= neq; i++ ) {
            rhs.at(i) += rhs2.at(i);
        }
    }

    if ( tStep->giveTimeDiscretization() == TD_Wilson ) {
        for ( int i = 1; i <= neq; i++ ) {
            rhs.at(i) += previousLoadVector.at(i)
                + theta * ( loadVector.at(i) - previousLoadVector.at(i) );
        }
    } else {
        for ( int i = 1; i <= neq; i++ ) {
            rhs.at(i) += loadVector.at(i);
        }
    }

    //
    // Call numerical model to solve arised problem.
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    nMethod->solve(stiffnessMatrix, & rhs, & displacementVector);

    if ( tStep->giveTimeDiscretization() == TD_Wilson ) {
        OOFEM_LOG_INFO("TD_Wilson: Updating acceleration, velocity and displacement.\n");
        for ( int i = 1; i <= neq; i++ ) {
            accelerationVector.at(i) = a6 * ( displacementVector.at(i)- previousDisplacementVector.at(i) )
                - a7 * previousVelocityVector.at(i)
                + a8 * previousAccelerationVector.at(i);

            velocityVector.at(i) = previousVelocityVector.at(i)
                + a9 * ( previousAccelerationVector.at(i) + accelerationVector.at(i) );

            displacementVector.at(i) = previousDisplacementVector.at(i)
                + deltaT * previousVelocityVector.at(i)
                + a10 * ( accelerationVector.at(i) + 2 * previousAccelerationVector.at(i) );
        }
    } else if ( ( tStep->giveTimeDiscretization() == TD_ThreePointBackward ) 
                || ( tStep->giveTimeDiscretization() == TD_TwoPointBackward ) ) {
        for ( int i = 1; i <= neq; i++ ) {
            accelerationVector.at(i) = a0 * ( displacementVector.at(i) - previousDisplacementVector.at(i) )
                - a2 * previousVelocityVector.at(i);

            velocityVector.at(i) = a1  * ( displacementVector.at(i) - previousDisplacementVector.at(i) )
                - a11 * previousIncrementOfDisplacement.at(i);
        }        
    } else {
        for ( int i = 1; i <= neq; i++ ) {
            accelerationVector.at(i) = a0 * ( displacementVector.at(i) - previousDisplacementVector.at(i) )
                - a2 * previousVelocityVector.at(i)
                - a3 * previousAccelerationVector.at(i);
            
            velocityVector.at(i) = previousVelocityVector.at(i)
                + a6 * previousAccelerationVector.at(i)
                + a7 * accelerationVector.at(i);
        }
    }
}


void
DIIDynamic :: giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                                    CharType type, TimeStep *tStep, Domain *domain)
{
    // We don't directly call element ->GiveCharacteristicMatrix() function, because some
    // engngm classes may require special modification of base types supported on
    // element class level.

    if ( type == EffectiveStiffnessMatrix ) {
        Element *element;
        FloatMatrix charMtrx;

        element = domain->giveElement(num);
        element->giveCharacteristicMatrix(answer, StiffnessMatrix, tStep);
        answer.times(1 + this->delta * a1);

        element->giveCharacteristicMatrix(charMtrx, MassMatrix, tStep);
        charMtrx.times(this->a0 + this->eta * this->a1);

        answer.add(charMtrx);

        return;
    } else {
        StructuralEngngModel :: giveElementCharacteristicMatrix(answer, num, type, tStep, domain);
    }
}


void DIIDynamic :: updateYourself(TimeStep *tStep)
{
    OOFEM_LOG_INFO("updateYourself\n");
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);

    for ( int i = 1; i <= neq; i++ ) {
        previousIncrementOfDisplacement.at(i) = displacementVector.at(i) - previousDisplacementVector.at(i);
        previousDisplacementVector.at(i)      = displacementVector.at(i);
        previousVelocityVector.at(i)          = velocityVector.at(i);
        previousAccelerationVector.at(i)      = accelerationVector.at(i);
        previousLoadVector.at(i)              = loadVector.at(i);
    }

    this->updateInternalState(tStep);
    StructuralEngngModel :: updateYourself(tStep);
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
DIIDynamic :: timesMtrx(FloatArray &vec, FloatArray &answer, CharType type, Domain *domain, TimeStep *tStep)
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
        element->giveCharacteristicMatrix(charMtrx, type, tStep);

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
    if ( ( currentStep->giveNumber() == giveNumberOfFirstStep() )
         && ( initialTimeDiscretization == TD_ThreePointBackward ) ) {
        currentStep->setTimeDiscretization( TD_TwoPointBackward );
    }

    deltaT = tStep->giveTimeIncrement();

    if ( tStep->giveTimeDiscretization() == TD_Newmark ) {
        double dt2 = deltaT * deltaT;
        a0  = 1 / ( beta * dt2 );
        a1  = gamma / ( beta * deltaT );
        a2  = 1 / ( beta * deltaT );
        a3  = 1 / ( 2 *  beta ) - 1;
        a4  = ( gamma / beta ) - 1;
        a5  = deltaT / 2 * ( gamma / beta - 2 );
        a6  = deltaT * ( 1 - gamma );
        a7  = deltaT * gamma;
        a8  = dt2 * ( ( 1 / 2 ) - beta );
        a9  = dt2 * beta;
        a10 = 0;
        a11 = 0;
    } else if ( tStep->giveTimeDiscretization() == TD_TwoPointBackward  ) {
        double dt2 = deltaT * deltaT;
        a0  = 1 / dt2;
        a1  = 1 / deltaT;
        a2  = 1 / deltaT;
        a3  = 0;
        a4  = 0;
        a5  = 0;
        a6  = 0;
        a7  = 0;
        a8  = 0;
        a9  = 0;
        a10 = 0;
        a11 = 0;
    } else if ( tStep->giveTimeDiscretization() == TD_ThreePointBackward ) {
        double dt2 = deltaT * deltaT;
        a0  = 2 / dt2;
        a1  = 3 / ( 2 * deltaT );
        a2  = 2 / deltaT;
        a3  = 0;
        a4  = 0;
        a5  = 0;
        a6  = 0;
        a7  = 0;
        a8  = 0;
        a9  = 0;
        a10 = 0;
        a11 = 1 / ( 2 * deltaT );
    } else if ( tStep->giveTimeDiscretization() ==  TD_Wilson ){
        a0  = 6 / ( ( theta * deltaT ) * ( theta * deltaT ) );
        a1  = 3 / ( theta * deltaT );
        a2  = 2 * a1;
        a3  = 2;
        a4  = 2;
        a5  = theta * deltaT / 2;
        a6  = a0 / theta;
        a7  = a2 / theta;
        a8  = 1 - 3 / theta;
        a9  = deltaT / 2;
        a10 = deltaT * deltaT / 6;
        a11 = 0;
    } else {
        _error("DIIDynamic: Time-stepping scheme not found!\n");
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

    if ( ( iores = previousIncrementOfDisplacement.storeYourself(stream, mode) ) != CIO_OK ) {
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

    if ( ( iores = previousIncrementOfDisplacement.restoreYourself(stream, mode) ) != CIO_OK ) {
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
