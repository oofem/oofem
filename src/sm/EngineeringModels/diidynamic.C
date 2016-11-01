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

#include "../sm/EngineeringModels/diidynamic.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/structuralelementevaluator.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "contextioerr.h"
#include "datastream.h"
#include "verbose.h"
#include "unknownnumberingscheme.h"
#include "classfactory.h"

namespace oofem {
REGISTER_EngngModel(DIIDynamic);

DIIDynamic :: DIIDynamic(int i, EngngModel *_master) : StructuralEngngModel(i, _master),
    loadVector(), previousLoadVector(), rhs(), displacementVector(), velocityVector(), accelerationVector(),
    previousDisplacementVector(), previousVelocityVector(), previousAccelerationVector(), previousIncrementOfDisplacement(), help()
{
    initFlag = true;
    ndomains = 1;

    initialTimeDiscretization = TD_ThreePointBackward;
}

DIIDynamic :: ~DIIDynamic()
{
}

NumericalMethod *DIIDynamic :: giveNumericalMethod(MetaStep *mStep)
// Only one has reason for DIIDynamic
// - SolutionOfLinearEquations
{
    if ( nMethod ) {
        return nMethod.get();
    }

    nMethod.reset( classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this) );
    if ( !nMethod ) {
        OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
    }

    return nMethod.get();
}


IRResultType
DIIDynamic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = StructuralEngngModel :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_DIIDynamic_ddtScheme);

    eta = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, eta, _IFT_DIIDynamic_eta);

    delta = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, delta, _IFT_DIIDynamic_delta);

    initialTimeDiscretization = ( TimeDiscretizationType ) val;

    gamma = 0.5;
    beta = 0.25;              // Default Newmark parameters.
    theta = 1.37; // Default Wilson parameter.
    if ( initialTimeDiscretization == TD_Newmark ) {
        OOFEM_LOG_INFO("Selecting Newmark-beta metod\n");
        IR_GIVE_OPTIONAL_FIELD(ir, gamma, _IFT_DIIDynamic_gamma);
        IR_GIVE_OPTIONAL_FIELD(ir, beta, _IFT_DIIDynamic_beta);
    } else if ( initialTimeDiscretization == TD_TwoPointBackward ) {
        OOFEM_LOG_INFO("Selecting Two-point Backward Euler method\n");
    } else if ( initialTimeDiscretization == TD_ThreePointBackward ) {
        OOFEM_LOG_INFO("Selecting Three-point Backward Euler metod\n");
    } else if ( initialTimeDiscretization == TD_Wilson ) {
        OOFEM_LOG_INFO("Selecting Wilson-theta metod\n");
        IR_GIVE_OPTIONAL_FIELD(ir, theta, _IFT_DIIDynamic_theta);
        if ( theta < 1.37 ) {
            OOFEM_WARNING("Found theta < 1.37. Performing correction, theta = 1.37");
            theta = 1.37;
        }
    } else {
        OOFEM_WARNING("Time-stepping scheme not found!");
        return IRRT_BAD_FORMAT;
    }

    IR_GIVE_FIELD(ir, deltaT, _IFT_DIIDynamic_deltat);

    return IRRT_OK;
}


double DIIDynamic :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// Returns unknown quantity displacement, velocity or acceleration of equation eq.
// This function translates this request to numerical method language.
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }
#endif

    if ( tStep != this->giveCurrentStep() ) {
        OOFEM_ERROR("unknown time step encountered");
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
        OOFEM_ERROR("Unknown is of undefined ValueModeType for this problem");
    }

    return 0.0;
}


TimeStep *DIIDynamic :: giveNextStep()
{
    int istep = giveNumberOfFirstStep();
    double totalTime = deltaT;
    StateCounterType counter = 1;
    TimeDiscretizationType td = initialTimeDiscretization;

    if ( currentStep ) {
        totalTime = currentStep->giveTargetTime() + deltaT;
        istep     = currentStep->giveNumber() + 1;
        counter   = currentStep->giveSolutionStateCounter() + 1;
        td        = currentStep->giveTimeDiscretization();
        if ( currentStep->isTheFirstStep() &&
            ( initialTimeDiscretization == TD_ThreePointBackward ) ) {
            td = TD_ThreePointBackward;
        }
    }

    previousStep = std :: move(currentStep);

    currentStep.reset( new TimeStep(istep, this, 1, totalTime, deltaT, counter, td) );

    return currentStep.get();
}

void DIIDynamic :: initializeYourself(TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
	int neq =  this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());

    if ( tStep->isTheFirstStep() ) {
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

        for ( auto &node : domain->giveDofManagers()) {
            for ( Dof *iDof: *node ) {
                //
                // Ask for initial values obtained from boundary conditions and initial conditions.
                //
                if ( !iDof->isPrimaryDof() ) {
                    continue;
                }

                int jj = iDof->__giveEquationNumber();
                if ( jj ) {
                    displacementVector.at(jj) = iDof->giveUnknown(VM_Total, stepWhenIcApply);
                    velocityVector.at(jj)     = iDof->giveUnknown(VM_Velocity, stepWhenIcApply);
                    accelerationVector.at(jj) = iDof->giveUnknown(VM_Acceleration, stepWhenIcApply);
                }
            }
        }

        delete stepWhenIcApply;
    }   // End of initialization.
}

void DIIDynamic :: solveYourself()
{
    StructuralEngngModel :: solveYourself();
}

void DIIDynamic :: solveYourselfAt(TimeStep *tStep)
{
    // Determine the constants.
    this->determineConstants(tStep);

    Domain *domain = this->giveDomain(1);
    int neq =  this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());


    if ( initFlag ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling stiffness matrix\n");
#endif
        stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        if ( !stiffnessMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

        this->assemble(*stiffnessMatrix, tStep, EffectiveTangentAssembler(false, 1 + this->delta * a1,  this->a0 + this->eta * this->a1),
                       EModelDefaultEquationNumbering(), domain);

        help.resize(neq);
        help.zero();

        previousDisplacementVector = displacementVector;
        previousVelocityVector = velocityVector;
        previousAccelerationVector = accelerationVector;
        previousLoadVector = loadVector;

        initFlag = 0;
    }

    if ( ( tStep->givePreviousStep() != NULL ) && ( tStep->giveTimeDiscretization() != tStep->givePreviousStep()->giveTimeDiscretization() ) ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling stiffness matrix\n");
#endif
        stiffnessMatrix->zero();
        this->assemble(*stiffnessMatrix, tStep, EffectiveTangentAssembler(false, 1 + this->delta * a1,  this->a0 + this->eta * this->a1),
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
        + a2 *previousVelocityVector.at(i)
        + a3 *previousAccelerationVector.at(i)
        + eta * ( a1 * previousDisplacementVector.at(i)
                 + a4 * previousVelocityVector.at(i)
                 + a5 * previousAccelerationVector.at(i)
                 + a11 * previousIncrementOfDisplacement.at(i) );
    }

    this->timesMtrx(help, rhs, MassMatrix, domain, tStep);
    //this->assembleVector(help, tStep, MatrixProductAssembler(MassMatrix(), rhs), VM_Total, 
    //                    EModelDefaultEquationNumbering(), this->giveDomain(1));

    help.zero();

    if ( delta != 0 ) {
        for ( int i = 1; i <= neq; i++ ) {
            help.at(i) = delta * ( a1 * previousDisplacementVector.at(i)
                                  + a4 * previousVelocityVector.at(i)
                                  + a5 * previousAccelerationVector.at(i)
                                  + a6 * previousIncrementOfDisplacement.at(i) );
        }
        this->timesMtrx(help, rhs2, TangentStiffnessMatrix, domain, tStep);
        //this->assembleVector(help, tStep, MatrixProductAssembler(TangentAssembler(), rhs2), VM_Total, 
        //                    EModelDefaultEquationNumbering(), this->giveDomain(1));

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

    nMethod->solve(*stiffnessMatrix, rhs, displacementVector);

    if ( tStep->giveTimeDiscretization() == TD_Wilson ) {
        OOFEM_LOG_INFO("TD_Wilson: Updating acceleration, velocity and displacement.\n");
        for ( int i = 1; i <= neq; i++ ) {
            accelerationVector.at(i) = a6 * ( displacementVector.at(i) - previousDisplacementVector.at(i) )
                                       - a7 *previousVelocityVector.at(i)
            + a8 *previousAccelerationVector.at(i);

            velocityVector.at(i) = previousVelocityVector.at(i)
            + a9 * ( previousAccelerationVector.at(i) + accelerationVector.at(i) );

            displacementVector.at(i) = previousDisplacementVector.at(i)
            + deltaT *previousVelocityVector.at(i)
            + a10 * ( accelerationVector.at(i) + 2 * previousAccelerationVector.at(i) );
        }
    } else if ( ( tStep->giveTimeDiscretization() == TD_ThreePointBackward ) ||
               ( tStep->giveTimeDiscretization() == TD_TwoPointBackward ) ) {
        for ( int i = 1; i <= neq; i++ ) {
            accelerationVector.at(i) = a0 * ( displacementVector.at(i) - previousDisplacementVector.at(i) )
                                       - a2 *previousVelocityVector.at(i);

            velocityVector.at(i) = a1  * ( displacementVector.at(i) - previousDisplacementVector.at(i) )
                                   - a11 *previousIncrementOfDisplacement.at(i);
        }
    } else {
        for ( int i = 1; i <= neq; i++ ) {
            accelerationVector.at(i) = a0 * ( displacementVector.at(i) - previousDisplacementVector.at(i) )
                                       - a2 *previousVelocityVector.at(i)
            - a3 *previousAccelerationVector.at(i);

            velocityVector.at(i) = previousVelocityVector.at(i)
            + a6 *previousAccelerationVector.at(i)
            + a7 *accelerationVector.at(i);
        }
    }
}


void DIIDynamic :: updateYourself(TimeStep *tStep)
{
    previousIncrementOfDisplacement.beDifferenceOf(displacementVector, previousDisplacementVector);
    previousDisplacementVector = displacementVector;
    previousVelocityVector = velocityVector;
    previousAccelerationVector = accelerationVector;
    previousLoadVector = loadVector;

    StructuralEngngModel :: updateYourself(tStep);
}


void
DIIDynamic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    static char dofchar[] = "dva";
    static ValueModeType dofmodes[] = {
        VM_Total, VM_Velocity, VM_Acceleration
    };

    iDof->printMultipleOutputAt(stream, tStep, dofchar, dofmodes, 3);
}


void
DIIDynamic :: timesMtrx(FloatArray &vec, FloatArray &answer, CharType type, Domain *domain, TimeStep *tStep)
{
    int nelem = domain->giveNumberOfElements();
    int neq   = this->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultEquationNumbering() );
    int i, j, k, jj, kk, n;
    FloatMatrix charMtrx, R;
    IntArray loc;
    Element *element;
    EModelDefaultEquationNumbering en;

    answer.resize(neq);
    answer.zero();

    for ( i = 1; i <= nelem; i++ ) {
        element = domain->giveElement(i);
        // Skip remote elements.
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

        element->giveLocationArray(loc, en);
        element->giveCharacteristicMatrix(charMtrx, type, tStep);
        if ( charMtrx.isNotEmpty() ) {
          ///@todo This rotation matrix is not flexible enough.. it can only work with full size matrices and doesn't allow for flexibility in the matrixassembler.
          if ( element->giveRotationMatrix(R) ) {
            charMtrx.rotatedWith(R);
          }
        }

#ifdef DEBUG
        if ( loc.giveSize() != charMtrx.giveNumberOfRows() ) {
            OOFEM_ERROR("dimension mismatch");
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
    _loadVector.resize( this->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultEquationNumbering() ) );
    _loadVector.zero();

    this->assembleVector(_loadVector, tStep, ExternalForceAssembler(), mode,
                         EModelDefaultEquationNumbering(), domain);
    this->updateSharedDofManagers(_loadVector, EModelDefaultEquationNumbering(), LoadExchangeTag);
}

void
DIIDynamic :: determineConstants(TimeStep *tStep)
{
    if ( this->giveCurrentStep()->isTheFirstStep() && initialTimeDiscretization == TD_ThreePointBackward ) {
        this->giveCurrentStep()->setTimeDiscretization( TD_TwoPointBackward );
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
    } else if ( tStep->giveTimeDiscretization() ==  TD_Wilson ) {
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
        OOFEM_ERROR("Time-stepping scheme not found!");
    }
}


contextIOResultType DIIDynamic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file = NULL;

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

    if ( ( iores = displacementVector.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = loadVector.storeYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = previousIncrementOfDisplacement.storeYourself(*stream) ) != CIO_OK ) {
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
    FILE *file = NULL;

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

    if ( ( iores = displacementVector.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = loadVector.restoreYourself(*stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = previousIncrementOfDisplacement.restoreYourself(*stream) ) != CIO_OK ) {
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
