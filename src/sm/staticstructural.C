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

#include "staticstructural.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "sparsemtrx.h"
#include "verbose.h"
#include "structuralelement.h"
#include "structuralelementevaluator.h"
#include "primaryfield.h"
#include "nrsolver.h"
#include "classfactory.h"
#include "error.h"
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
    stiffnessMatrix(NULL),
    field(NULL),
    sparseMtrxType(SMT_Skyline),
    nMethod(NULL)
{
    ndomains = 1;

#ifdef __PARALLEL_MODE
    commMode = ProblemCommMode__NODE_CUT;
    nonlocalExt = 0;
    communicator = nonlocCommunicator = NULL;
    commBuff = NULL;
#endif
}


StaticStructural :: ~StaticStructural()
{
    delete field;
    delete stiffnessMatrix;
    delete nMethod;
}


NumericalMethod *StaticStructural :: giveNumericalMethod(MetaStep *mStep)
{
    if ( nMethod ) {
        return nMethod;
    }
    nMethod = new NRSolver(this->giveDomain(1), this);
    return nMethod;
}

IRResultType
StaticStructural :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralEngngModel :: initializeFrom(ir);
    int val = SMT_Skyline;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    this->deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, deltaT, _IFT_StaticStructural_deltat);

#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        delete communicator;
        delete commBuff;
        commBuff = new CommunicatorBuff( this->giveNumberOfProcesses() );
        communicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                               this->giveNumberOfProcesses(),
                                               this->commMode);
    }

#endif

    delete this->field;
    this->field = new PrimaryField(this, 1, FT_Displacements, EID_MomentumBalance, 1);

    return IRRT_OK;
}


TimeStep *StaticStructural :: giveNextStep()
{
    if ( previousStep ) {
        delete previousStep;
    }

    if ( currentStep == NULL ) {
        int istep = this->giveNumberOfFirstStep();
        // first step -> generate initial step
        previousStep = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, -this->deltaT, this->deltaT, 0);
        currentStep = new TimeStep(istep, this, 1, 0.0, this->deltaT, 1);
    } else {
        int istep =  currentStep->giveNumber() + 1;
        StateCounterType counter = currentStep->giveSolutionStateCounter() + 1;
        previousStep = currentStep;
        double dt = currentStep->giveTimeIncrement();
        double totalTime = currentStep->giveTargetTime() + dt;
        currentStep = new TimeStep(istep, this, 1, totalTime, dt, counter);
    }

    return currentStep;
}


void StaticStructural :: solveYourself()
{
    ///@todo Generalize this to engngmodel?
#ifdef __PARALLEL_MODE
    if ( this->isParallel() ) {
 #ifdef __VERBOSE_PARALLEL
        // force equation numbering before setting up comm maps
        OOFEM_LOG_INFO( "[process rank %d] neq is %d\n", this->giveRank(), this->giveNumberOfDomainEquations(EID_MomentumBalance) );
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
    int di = 1;
    int neq = this->giveNumberOfDomainEquations( di, EModelDefaultEquationNumbering() );

    field->advanceSolution(tStep);

    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    this->initMetaStepAttributes( this->giveCurrentMetaStep() );

    // Fetch vector to fill in from primary field.
    this->solution = field->giveSolutionVector(tStep);
    this->solution->resize(neq);
    this->solution->zero();

    // Create "stiffness matrix"
    if ( !this->stiffnessMatrix ) {
        this->stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        if ( !this->stiffnessMatrix ) {
            OOFEM_ERROR("Couldn't create requested sparse matrix of type %d", sparseMtrxType);
        }

        this->stiffnessMatrix->buildInternalStructure( this, di, EID_MomentumBalance, EModelDefaultEquationNumbering() );
    }
    this->internalForces.resize(neq);

    FloatArray incrementOfSolution(neq);


    // Build initial/external load
    FloatArray externalForces(neq);
    externalForces.zero();
    this->assembleVector( externalForces, tStep, EID_MomentumBalance, ExternalForcesVector, VM_Total,
                         EModelDefaultEquationNumbering(), this->giveDomain(1) );
#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);
#endif

    if ( this->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("\nStaticStructural :: solveYourselfAt - Solving step %d, metastep %d, (neq = %d)\n", tStep->giveNumber(), tStep->giveMetaStepNumber(), neq);
    }

    double loadLevel;
    int currentIterations;
    NM_Status status = this->nMethod->solve(this->stiffnessMatrix,
                                            & externalForces,
                                            NULL,
                                            solution,
                                            & incrementOfSolution,
                                            & ( this->internalForces ),
                                            this->eNorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total,
                                            currentIterations,
                                            tStep);

    if ( !( status & NM_Success ) ) {
        OOFEM_ERROR("No success in solving problem");
    }
}


double StaticStructural :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return this->field->giveUnknownValue(dof, mode, tStep);
}


void StaticStructural :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    if ( cmpn == InternalRhs ) {
        this->internalForces.zero();
        this->assembleVector(this->internalForces, tStep, EID_MomentumBalance, InternalForcesVector, VM_Total,
                             EModelDefaultEquationNumbering(), d, & this->eNorm);
#ifdef __PARALLEL_MODE
        this->updateSharedDofManagers(this->internalForces, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
#endif
    } else if ( cmpn == NonLinearLhs ) {
        this->stiffnessMatrix->zero();
        this->assemble(this->stiffnessMatrix, tStep, EID_MomentumBalance, TangentStiffnessMatrix, EModelDefaultEquationNumbering(), d);
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

    if ( ( iores = this->field->saveContext(stream, mode) ) != CIO_OK ) {
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

    if ( ( iores = this->field->restoreContext(stream, mode) ) != CIO_OK ) {
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



void
StaticStructural :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    iDof->printSingleOutputAt(stream, tStep, 'd', VM_Total);
}


#ifdef __PARALLEL_MODE
int
StaticStructural :: estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType)
{
    int count = 0, pcount = 0;
    Domain *domain = this->giveDomain(1);

    if ( packUnpackType == ProblemCommMode__ELEMENT_CUT ) {
        for ( int map: commMap ) {
            count += domain->giveDofManager( map )->giveNumberOfDofs();
        }

        return ( buff.givePackSize(MPI_DOUBLE, 1) * count );
    } else if ( packUnpackType == ProblemCommMode__NODE_CUT ) {
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

        return ( buff.givePackSize(MPI_DOUBLE, 1) * pcount );
    } else if ( packUnpackType == ProblemCommMode__REMOTE_ELEMENT_MODE ) {
        for ( int map: commMap ) {
            count += domain->giveElement( map )->estimatePackSize(buff);
        }

        return count;
    }

    return 0;
}
#endif
} // end namespace oofem
