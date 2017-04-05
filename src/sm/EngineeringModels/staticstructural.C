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

#include "../sm/EngineeringModels/staticstructural.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/structuralelementevaluator.h"
#include "timestep.h"
#include "sparsemtrx.h"
#include "nummet.h"
#include "nrsolver.h"
#include "staggeredsolver.h"
#include "dynamicrelaxationsolver.h"
#include "primaryfield.h"
#include "dofdistributedprimaryfield.h"
#include "verbose.h"
#include "error.h"
#include "generalboundarycondition.h"
#include "boundarycondition.h"
#include "activebc.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "communicator.h"
#endif

namespace oofem {
REGISTER_EngngModel(StaticStructural);

StaticStructural :: StaticStructural(int i, EngngModel *_master) : StructuralEngngModel(i, _master),
    internalForces(),
    eNorm(),
    sparseMtrxType(SMT_Skyline)
{
    ndomains = 1;
    solverType = 0;
    mRecomputeStepAfterPropagation = false;
}


StaticStructural :: ~StaticStructural()
{
}


NumericalMethod *StaticStructural :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        if ( solverType == 0 ) {
            nMethod.reset( new NRSolver(this->giveDomain(1), this) );
        } else if ( solverType == 1 ) {
            nMethod.reset( new StaggeredSolver(this->giveDomain(1), this) );
            // Check if sparse matrix is SMT_Skyline
            if ( this->sparseMtrxType != SMT_Skyline ) {
                OOFEM_ERROR("Only Skyline sparse matrix type is currently supported (0) for the staggered solver");
            }
            
        } else if ( solverType == 2 ) {
            nMethod.reset( new DynamicRelaxationSolver(this->giveDomain(1), this) );
        } else {
            OOFEM_ERROR("Unsupported solver (%d). Solvers currently supported are: 0 - NR (default) and 1 - staggered NR, 2 - Dynamic relaxation solver", solverType);
        }
    }
    return nMethod.get();
}

int
StaticStructural :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % 2;
}

IRResultType
StaticStructural :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = StructuralEngngModel :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    int val = SMT_Skyline;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    this->deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, _IFT_StaticStructural_deltat);

    this->solverType = 0; // Default NR
    IR_GIVE_OPTIONAL_FIELD(ir, solverType, _IFT_StaticStructural_solvertype);

    int _val = IG_Tangent;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_EngngModel_initialGuess);
    this->initialGuessType = ( InitialGuess ) _val;
    
    mRecomputeStepAfterPropagation = ir->hasField(_IFT_StaticStructural_recomputeaftercrackpropagation);

#ifdef __PARALLEL_MODE
    ///@todo Where is the best place to create these?
    if ( isParallel() ) {
        delete communicator;
        delete commBuff;
        commBuff = new CommunicatorBuff( this->giveNumberOfProcesses() );
        communicator = new NodeCommunicator(this, commBuff, this->giveRank(),
                                            this->giveNumberOfProcesses());
    }

#endif

    this->field.reset( new DofDistributedPrimaryField(this, 1, FT_Displacements, 0) );

    return IRRT_OK;
}


TimeStep *StaticStructural :: giveNextStep()
{
    if ( !currentStep ) {
        // first step -> generate initial step
        //currentStep.reset( new TimeStep(*giveSolutionStepWhenIcApply()) );
        currentStep.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 1, 0., this->deltaT, 0) );
    }
    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(*previousStep, this->deltaT) );

    return currentStep.get();
}


void StaticStructural :: solveYourself()
{
    ///@todo Generalize this to engngmodel?
#ifdef __PARALLEL_MODE
    if ( this->isParallel() ) {
 #ifdef __VERBOSE_PARALLEL
        // force equation numbering before setting up comm maps
        OOFEM_LOG_INFO( "[process rank %d] neq is %d\n", this->giveRank(), this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()) );
 #endif

        // set up communication patterns
        // needed only for correct shared rection computation
        communicator->setUpCommunicationMaps(this, true);
        if ( nonlocalExt ) {
            nonlocCommunicator->setUpCommunicationMaps(this, true);
        }
    }
#endif

    StructuralEngngModel :: solveYourself();
}


void StaticStructural :: solveYourselfAt(TimeStep *tStep)
{
    int neq;
    int di = 1;

    this->field->advanceSolution(tStep);
    this->field->applyBoundaryCondition(tStep); ///@todo Temporary hack, advanceSolution should apply the boundary conditions directly.

    neq = this->giveNumberOfDomainEquations( di, EModelDefaultEquationNumbering() );
    if (tStep->giveNumber()==1) {
        this->field->initialize(VM_Total, tStep, this->solution, EModelDefaultEquationNumbering() );
    } else {
        this->field->initialize(VM_Total, tStep->givePreviousStep(), this->solution, EModelDefaultEquationNumbering() );
        this->field->update(VM_Total, tStep, this->solution, EModelDefaultEquationNumbering() );
    }
    this->field->applyBoundaryCondition(tStep); ///@todo Temporary hack to override the incorrect values that is set by "update" above. Remove this when that is fixed.

    FloatArray incrementOfSolution(neq), externalForces(neq);

    // Create "stiffness matrix"
    if ( !this->stiffnessMatrix ) {
        this->stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        if ( !this->stiffnessMatrix ) {
            OOFEM_ERROR("Couldn't create requested sparse matrix of type %d", sparseMtrxType);
        }

        this->stiffnessMatrix->buildInternalStructure( this, di, EModelDefaultEquationNumbering() );
    }
    this->internalForces.resize(neq);

    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    this->initMetaStepAttributes( this->giveCurrentMetaStep() );

    if ( this->initialGuessType == IG_Tangent ) {
        OOFEM_LOG_RELEVANT("Computing initial guess\n");
        FloatArray extrapolatedForces(neq);
        this->assembleExtrapolatedForces( extrapolatedForces, tStep, TangentStiffnessMatrix, this->giveDomain(di) );
        extrapolatedForces.negated();
        ///@todo Need to find a general way to support this before enabling it by default.
        //this->assembleVector(extrapolatedForces, tStep, LinearizedDilationForceAssembler(), VM_Incremental, EModelDefaultEquationNumbering(), this->giveDomain(di) );
#if 0
        // Some debug stuff:
        extrapolatedForces.printYourself("extrapolatedForces");
        this->internalForces.zero();
        this->assembleVectorFromElements(this->internalForces, tStep, InternalForceAssembler(), VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(di));
        this->internalForces.printYourself("internal forces");
#endif
        OOFEM_LOG_RELEVANT("Computing old tangent\n");
        this->updateComponent( tStep, NonLinearLhs, this->giveDomain(di) );
        SparseLinearSystemNM *linSolver = nMethod->giveLinearSolver();
        OOFEM_LOG_RELEVANT("Solving for increment\n");
        linSolver->solve(*stiffnessMatrix, extrapolatedForces, incrementOfSolution);
        OOFEM_LOG_RELEVANT("Initial guess found\n");
        this->solution.add(incrementOfSolution);
        
        this->field->update(VM_Total, tStep, this->solution, EModelDefaultEquationNumbering());
        this->field->applyBoundaryCondition(tStep); ///@todo Temporary hack to override the incorrect values that is set by "update" above. Remove this when that is fixed.
    } else if ( this->initialGuessType != IG_None ) {
        OOFEM_ERROR("Initial guess type: %d not supported", initialGuessType);
    } else {
        incrementOfSolution.zero();
    }

    // Build initial/external load
    externalForces.zero();
    this->assembleVector( externalForces, tStep, ExternalForceAssembler(), VM_Total,
                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);

    if ( this->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("\nStaticStructural :: solveYourselfAt - Solving step %d, metastep %d, (neq = %d)\n", tStep->giveNumber(), tStep->giveMetaStepNumber(), neq);
    }

    double loadLevel;
    int currentIterations;
    NM_Status status = this->nMethod->solve(*this->stiffnessMatrix,
                                            externalForces,
                                            NULL,
                                            this->solution,
                                            incrementOfSolution,
                                            this->internalForces,
                                            this->eNorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total,
                                            currentIterations,
                                            tStep);
    if ( !( status & NM_Success ) ) {
        OOFEM_ERROR("No success in solving problem");
    }
}

void StaticStructural :: terminate(TimeStep *tStep)
{
    if ( mRecomputeStepAfterPropagation ) {
        // Propagate cracks and recompute time step
        XfemSolverInterface::propagateXfemInterfaces(tStep, *this, true);
        StructuralEngngModel::terminate(tStep);
    } else {
        // Propagate cracks at the end of the time step
        StructuralEngngModel::terminate(tStep);
        XfemSolverInterface::propagateXfemInterfaces(tStep, *this, false);
    }
    
     for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        this->printReactionForces(tStep, idomain);
    }
}

double StaticStructural :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    double val1 = dof->giveUnknownsDictionaryValue(tStep, VM_Total);
    if ( mode == VM_Total ) {
        return val1;
    } else if ( mode == VM_Incremental ) {
        double val0 = dof->giveUnknownsDictionaryValue(tStep->givePreviousStep(), VM_Total);
        return val1 - val0;
    } else {
        OOFEM_ERROR("Unknown value mode requested");
        return 0;
    }
    return this->field->giveUnknownValue(dof, mode, tStep);
}


void StaticStructural :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    if ( cmpn == InternalRhs ) {
        // Updates the solution in case it has changed 
        ///@todo NRSolver should report when the solution changes instead of doing it this way.
        this->field->update(VM_Total, tStep, this->solution, EModelDefaultEquationNumbering());
        this->field->applyBoundaryCondition(tStep);///@todo Temporary hack to override the incorrect vavues that is set by "update" above. Remove this when that is fixed.

        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, InternalForceAssembler(), VM_Total,
                             EModelDefaultEquationNumbering(), d, & this->eNorm);
        this->updateSharedDofManagers(this->internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);

        internalVarUpdateStamp = tStep->giveSolutionStateCounter(); // Hack for linearstatic
    } else if ( cmpn == NonLinearLhs ) {
        this->stiffnessMatrix->zero();
        this->assemble(*this->stiffnessMatrix, tStep, TangentAssembler(TangentStiffness), EModelDefaultEquationNumbering(), d);
    } else {
        OOFEM_ERROR("Unknown component");
    }
}


contextIOResultType StaticStructural :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( ( iores = StructuralEngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = this->field->saveContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    return CIO_OK;
}


contextIOResultType StaticStructural :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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

    if ( ( iores = StructuralEngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = this->field->restoreContext(*stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }

    return CIO_OK;
}


void
StaticStructural :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}

int
StaticStructural :: forceEquationNumbering()
{
    stiffnessMatrix.reset( NULL );
    return StructuralEngngModel::forceEquationNumbering();
}


void StaticStructural :: setSolution(TimeStep *tStep, const FloatArray &vectorToStore)
{
    this->field->update(VM_Total, tStep, vectorToStore, EModelDefaultEquationNumbering());
}


bool
StaticStructural :: requiresEquationRenumbering(TimeStep *tStep)
{
    if ( tStep->isTheFirstStep() ) {
        return true;
    }
    // Check if Dirichlet b.c.s has changed.
    Domain *d = this->giveDomain(1);
    for ( auto &gbc : d->giveBcs() ) {
        ActiveBoundaryCondition *active_bc = dynamic_cast< ActiveBoundaryCondition * >(gbc.get());
        BoundaryCondition *bc = dynamic_cast< BoundaryCondition * >(gbc.get());
        // We only need to consider Dirichlet b.c.s
        if ( bc || ( active_bc && ( active_bc->requiresActiveDofs() || active_bc->giveNumberOfInternalDofManagers() ) ) ) {
            // Check of the dirichlet b.c. has changed in the last step (if so we need to renumber)
            if ( gbc->isImposed(tStep) != gbc->isImposed(tStep->givePreviousStep()) ) {
                return true;
            }
        }
    }
    return false;
}

int
StaticStructural :: estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType)
{
    int count = 0, pcount = 0;
    Domain *domain = this->giveDomain(1);

    if ( packUnpackType == 0 ) { ///@todo Fix this old ProblemCommMode__NODE_CUT value
        for ( int map: commMap ) {
            DofManager *dman = domain->giveDofManager( map );
            for ( Dof *dof: *dman ) {
                if ( dof->isPrimaryDof() && dof->__giveEquationNumber() > 0 ) {
                    count++;
                } else {
                    pcount++;
                }
            }
        }

        // --------------------------------------------------------------------------------
        // only pcount is relevant here, since only prescribed components are exchanged !!!!
        // --------------------------------------------------------------------------------

        return ( buff.givePackSizeOfDouble(1) * pcount );
    } else if ( packUnpackType == 1 ) {
        for ( int map: commMap ) {
            count += domain->giveElement( map )->estimatePackSize(buff);
        }

        return count;
    }

    return 0;
}

} // end namespace oofem
