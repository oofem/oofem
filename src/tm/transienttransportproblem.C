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

#include "transienttransportproblem.h"
#include "timestep.h"
#include "dofdistributedprimaryfield.h"
#include "transportelement.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "nrsolver.h"

namespace oofem {
REGISTER_EngngModel(TransientTransportProblem);

TransientTransportProblem :: TransientTransportProblem(int i, EngngModel *_master = NULL) : EngngModel(i, _master)
{
    ndomains = 1;
}

TransientTransportProblem :: ~TransientTransportProblem() {}


NumericalMethod *TransientTransportProblem :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        nMethod.reset( new NRSolver(this->giveDomain(1), this) );
    }
    return nMethod.get();
}


IRResultType
TransientTransportProblem :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int val = SMT_Skyline;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    this->sparseMtrxType = ( SparseMtrxType ) val;

    IR_GIVE_FIELD(ir, this->alpha, _IFT_TransientTransportProblem_alpha);
    
    IR_GIVE_FIELD(ir, this->deltaT, _IFT_TransientTransportProblem_deltaT);
    
    this->keepTangent = ir->hasField(_IFT_TransientTransportProblem_keepTangent);

    this->lumped = ir->hasField(_IFT_TransientTransportProblem_lumped);

    field.reset( new DofDistributedPrimaryField(this, 1, FT_TransportProblemUnknowns, 0) );

    return EngngModel :: initializeFrom(ir);
}


double TransientTransportProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    //return this->field->giveUnknownValue(dof, mode, tStep);
    double val1 = field->giveUnknownValue(dof, VM_Total, tStep);
    double val0 = field->giveUnknownValue(dof, VM_Total, tStep->givePreviousStep());
    if ( mode == VM_Total ) {
        return this->alpha * val1 + (1.-this->alpha) * val0;
    } else if ( mode == VM_Velocity ) {
        return (val1 - val0) / tStep->giveTimeIncrement();
    } else if ( mode == VM_Incremental ) {
        return val1 - val0;
    } else {
        OOFEM_ERROR("Unknown value mode requested");
        return 0;
    }
}


TimeStep *TransientTransportProblem :: giveNextStep()
{
    if ( !currentStep ) {
        // first step -> generate initial step
        currentStep.reset( new TimeStep( *giveSolutionStepWhenIcApply() ) );
    }
    
    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(*previousStep, this->deltaT) );
    return currentStep.get();
}


TimeStep *TransientTransportProblem :: giveSolutionStepWhenIcApply()
{
    if ( !stepWhenIcApply ) {
        ///@todo Should we have this->initT ?
        stepWhenIcApply.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, 0., deltaT, 0) );
        // The initial step goes from [-deltaT, 0], so the intrinsic time is at: -deltaT  + alpha*deltaT
        stepWhenIcApply->setIntrinsicTime(-deltaT + alpha * deltaT);
    }

    return stepWhenIcApply.get();
}


void TransientTransportProblem :: solveYourselfAt(TimeStep *tStep)
{
    Domain *d = this->giveDomain(1);
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    field->advanceSolution(tStep);
    field->applyBoundaryCondition(tStep);
    if ( tStep->isTheFirstStep() ) {
        this->applyIC();
    }
    field->initialize(VM_Total, tStep, solution, EModelDefaultEquationNumbering());


    if ( !effectiveMatrix ) {
        effectiveMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        effectiveMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

        if ( lumped ) {
            capacityDiag.resize(neq);
            this->assembleVector( capacityDiag, tStep, LumpedMassMatrix, VM_Total, EModelDefaultEquationNumbering(), d );
        } else {
            capacityMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
            capacityMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
            this->assemble( *capacityMatrix, tStep, CapacityMatrix, EModelDefaultEquationNumbering(), d );
        }
        
        if ( this->keepTangent ) {
            this->assemble( *effectiveMatrix, tStep, TangentStiffnessMatrix,
                           EModelDefaultEquationNumbering(), d );
            effectiveMatrix->times(alpha);
            if ( lumped ) {
                effectiveMatrix->addDiagonal(1./tStep->giveTimeIncrement(), capacityDiag);
            } else {
                effectiveMatrix->add(1./tStep->giveTimeIncrement(), *capacityMatrix);
            }
        }
    }


    OOFEM_LOG_INFO("Assembling external forces\n");
    FloatArray externalForces(neq);
    externalForces.zero();
    this->assembleVector( externalForces, tStep, ExternalForcesVector, VM_Total, EModelDefaultEquationNumbering(), d );
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);

    // set-up numerical method
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    OOFEM_LOG_INFO("Solving for %d unknowns...\n", neq);

    internalForces.resize(neq);

    FloatArray incrementOfSolution;
    double loadLevel;
    int currentIterations;
    this->nMethod->solve(*this->effectiveMatrix,
                         externalForces,
                         NULL, // ignore
                         this->solution,
                         incrementOfSolution,
                         this->internalForces,
                         this->eNorm,
                         loadLevel, // ignore
                         SparseNonLinearSystemNM :: rlm_total, // ignore
                         currentIterations, // ignore
                         tStep);
}


void
TransientTransportProblem :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    // F(T) + C*dT/dt = Q
    // Linearized:
    // F(T^(k)) + K*a*dT_1 = Q - C * dT/dt^(k) - C/dt * dT_1
    // Rearranged
    // (a*K + C/dt) * dT_1 = Q - (F(T^(k)) + C * dT/dt^(k))
    // K_eff        * dT_1 = Q - F_eff
    // Update:
    // T_1 += dT_1
    
    ///@todo NRSolver should report when the solution changes instead of doing it this way.
    this->field->update(VM_Total, tStep, solution, EModelDefaultEquationNumbering());
    ///@todo Need to reset the boundary conditions properly since some "update" is doing strange 
    /// things such as applying the (wrong) boundary conditions. This call will be removed when that code can be removed.
    this->field->applyBoundaryCondition(tStep);

    if ( cmpn == InternalRhs ) {
        // F_eff = F(T^(k)) + C * dT/dt^(k)
        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, InternalForcesVector, VM_Total,
                             EModelDefaultEquationNumbering(), this->giveDomain(1), & this->eNorm);
        this->updateSharedDofManagers(this->internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);

        if ( lumped ) {
            // Note, inertia contribution cannot be computed on element level when lumped mass matrices are used.
            FloatArray oldSolution, vel;
            this->field->initialize(VM_Total, tStep->givePreviousStep(), oldSolution, EModelDefaultEquationNumbering());
            vel.beDifferenceOf(solution, oldSolution);
            vel.times( 1./tStep->giveTimeIncrement() );
            for ( int i = 0; i < vel.giveSize(); ++i ) {
                this->internalForces[i] += this->capacityDiag[i] * vel[i];
            }
        } else {
            FloatArray tmp;
            this->assembleVector(this->internalForces, tStep, InertiaForcesVector, VM_Total,
                                EModelDefaultEquationNumbering(), this->giveDomain(1), & tmp);
            this->eNorm.add(tmp); ///@todo Fix this, assembleVector shouldn't zero eNorm inside the functions. / Mikael
        }

    } else if ( cmpn == NonLinearLhs ) {
        // K_eff = (a*K + C/dt)
        if ( !this->keepTangent ) {
            this->effectiveMatrix->zero();
            this->assemble( *effectiveMatrix, tStep, TangentStiffnessMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            effectiveMatrix->times(alpha);
            if ( lumped ) {
                effectiveMatrix->addDiagonal(1./tStep->giveTimeIncrement(), capacityDiag);
            } else {
                effectiveMatrix->add(1./tStep->giveTimeIncrement(), *capacityMatrix);
            }
        }
    } else {
        OOFEM_ERROR("Unknown component");
    }
}


void
TransientTransportProblem :: applyIC()
{
    Domain *domain = this->giveDomain(1);
    OOFEM_LOG_INFO("Applying initial conditions\n");

    this->field->applyDefaultInitialCondition();

    ///@todo It's rather strange that the models need the initial values.
    // update element state according to given ic
    TimeStep *s = this->giveSolutionStepWhenIcApply();
    for ( auto &elem : domain->giveElements() ) {
        TransportElement *element = static_cast< TransportElement * >( elem.get() );
        element->updateInternalState(s);
        element->updateYourself(s);
    }
}


int
TransientTransportProblem :: forceEquationNumbering()
{
    this->capacityDiag.clear();
    this->capacityMatrix.reset(NULL);
    this->effectiveMatrix.reset(NULL);
    return EngngModel :: forceEquationNumbering();
}


contextIOResultType
TransientTransportProblem :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    int closeFlag = 0;
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

    if ( ( iores = field->saveContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    // ensure consistent records

    return CIO_OK;
}


contextIOResultType
TransientTransportProblem :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    int closeFlag = 0;
    int istep, iversion;
    FILE *file = NULL;

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

    if ( ( iores = field->restoreContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    // ensure consistent records

    return CIO_OK;
}


int
TransientTransportProblem :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % 2;
}


int
TransientTransportProblem :: requiresUnknownsDictionaryUpdate()
{
    return true;
}

int
TransientTransportProblem :: checkConsistency()
{
    // check for proper element type
    for ( auto &elem : this->giveDomain(1)->giveElements() ) {
        if ( !dynamic_cast< TransportElement * >( elem.get() ) ) {
            OOFEM_WARNING("Element %d has no TransportElement base", elem->giveLabel());
            return 0;
        }
    }

    return EngngModel :: checkConsistency();
}


void
TransientTransportProblem :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}


} // end namespace oofem
