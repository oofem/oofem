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
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "dofdistributedprimaryfield.h"
#include "verbose.h"
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

    field.reset( new DofDistributedPrimaryField(this, 1, FT_TransportProblemUnknowns, 0) );

    return EngngModel :: initializeFrom(ir);
}


double TransientTransportProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
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
    int istep = this->giveNumberOfFirstStep();
    //int mstep = 1;
    StateCounterType counter = 1;
    double tt = 0.;

    delete previousStep;

    if ( currentStep != NULL ) {
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
        tt = currentStep->giveTargetTime();
        previousStep = currentStep;
    } else {
        previousStep = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, 0., this->deltaT, 0);
    }

    currentStep = new TimeStep(istep, this, 1, tt + this->deltaT, this->deltaT, counter);
    // time and dt variables are set eq to 0 for statics - has no meaning
    return currentStep;
}


TimeStep *TransientTransportProblem :: giveSolutionStepWhenIcApply()
{
    if ( stepWhenIcApply == NULL ) {
        ///@todo Should we have this->initT ?
        stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, 0., deltaT, 0);
    }

    return stepWhenIcApply;
}


void TransientTransportProblem :: solveYourselfAt(TimeStep *tStep)
{
    field->advanceSolution(tStep);
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    if ( !effectiveMatrix ) {
        effectiveMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        effectiveMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

        capacityMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        capacityMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    }

    if ( tStep->isTheFirstStep() ) {
        this->field->applyDefaultInitialCondition();
    }

    field->initialize(VM_Total, tStep, solution, EModelDefaultEquationNumbering());

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling external forces\n");
#endif
    FloatArray externalForces(neq);
    externalForces.zero();
    this->assembleVector( externalForces, tStep, ExternalForcesVector, VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);
#if 1
    FloatArray tmp;
    this->capacityMatrix->times(solution, tmp);
    externalForces.add(1.0 / tStep->giveTimeIncrement(), tmp);
#endif

    // set-up numerical method
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    internalForces.resize(neq);

    FloatArray incrementOfSolution;
    double loadLevel;
    int currentIterations;
    this->nMethod->solve(*this->effectiveMatrix,
                         externalForces,
                         NULL,
                         this->solution,
                         incrementOfSolution,
                         this->internalForces,
                         this->eNorm,
                         loadLevel, // Only relevant for incrementalBCLoadVector
                         SparseNonLinearSystemNM :: rlm_total,
                         currentIterations,
                         tStep);
}


void
TransientTransportProblem :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    // F(T) = Q + C*dT/dt
    // Linearized:
    // F(T^(k)) + K*a*dT_1 = Q + C*T_1^(k) + C/dt * dT_1
    // Rearranged
    // (a*K - C/dt) * dT_1 = Q + C*T_1^(k) - F(T^(k))
    // Update:
    // T_1 += dT_1
    
    ///@todo NRSolver should report when the solution changes instead of doing it this way.
    this->field->update(VM_Total, tStep, solution, EModelDefaultEquationNumbering());

    if ( cmpn == InternalRhs ) {
        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, InternalForcesVector, VM_Total,
                             EModelDefaultEquationNumbering(), this->giveDomain(1), & this->eNorm);
        this->updateSharedDofManagers(this->internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
#if 0
        ///@todo We must have this:
        this->assembleVector(this->internalForces, tStep, InertiaForcesVector, VM_Total,
                             EModelDefaultEquationNumbering(), this->giveDomain(1), & this->eNorm);
#endif
        return;
    } else if ( cmpn == NonLinearLhs ) {
        //if ( !this->keepTangent ) {
            this->effectiveMatrix->zero();
            this->assemble( effectiveMatrix.get(), tStep, TangentStiffnessMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            ///@todo Check sign
            effectiveMatrix->times(alpha * tStep->giveTimeIncrement());
            this->assemble( effectiveMatrix.get(), tStep, CapacityMatrix,
                            EModelDefaultEquationNumbering(), this->giveDomain(1) );
            effectiveMatrix->times(1. / tStep->giveTimeIncrement());
        //}
        return;
    } else {
        OOFEM_ERROR("Unknown component");
    }
}


int
TransientTransportProblem :: forceEquationNumbering()
{
    this->capacityMatrix.release();
    this->effectiveMatrix.release();
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
TransientTransportProblem :: checkConsistency()
{
    Domain *domain = this->giveDomain(1);

    // check for proper element type
    for ( auto &elem : domain->giveElements() ) {
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


void
TransientTransportProblem :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    iDof->printSingleOutputAt(stream, tStep, 'd', VM_Total);
}

} // end namespace oofem
