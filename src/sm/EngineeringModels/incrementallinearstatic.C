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

#include "../sm/EngineeringModels/incrementallinearstatic.h"
#include "timestep.h"
#include "dof.h"
#include "domain.h"
#include "sparsemtrx.h"
#include "dictionary.h"
#include "verbose.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dofmanager.h"
#include "activebc.h"
#include "unknownnumberingscheme.h"

namespace oofem {
REGISTER_EngngModel(IncrementalLinearStatic);

IncrementalLinearStatic :: IncrementalLinearStatic(int i, EngngModel *_master) : StructuralEngngModel(i, _master),
    loadVector(), internalLoadVector(), incrementOfDisplacementVector(), discreteTimes()
{
    ndomains = 1;
    endOfTimeOfInterest = 0;
}


IncrementalLinearStatic :: ~IncrementalLinearStatic()
{}


NumericalMethod *IncrementalLinearStatic :: giveNumericalMethod(MetaStep *mStep)

{
    if ( !nMethod ) {
        nMethod.reset( classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this) );
        if ( !nMethod ) {
            OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
        }
    }

    return nMethod.get();
}


IRResultType IncrementalLinearStatic :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_OPTIONAL_FIELD(ir, discreteTimes, _IFT_IncrementalLinearStatic_prescribedtimes);
    if ( discreteTimes.giveSize() > 0 ) {
        numberOfSteps = discreteTimes.giveSize();
        endOfTimeOfInterest = discreteTimes.at( discreteTimes.giveSize() );
        fixedSteps = false;
    } else {
        deltaT = 1.0;
        IR_GIVE_OPTIONAL_FIELD(ir, deltaT, _IFT_IncrementalLinearStatic_deltat);
        IR_GIVE_FIELD(ir, numberOfSteps, _IFT_EngngModel_nsteps);
        endOfTimeOfInterest = deltaT * numberOfSteps;
        fixedSteps = true;
    }
    IR_GIVE_OPTIONAL_FIELD(ir, endOfTimeOfInterest, _IFT_IncrementalLinearStatic_endoftimeofinterest);

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    //StructuralEngngModel::initializeFrom (ir);
    return IRRT_OK;
}


double IncrementalLinearStatic :: giveDiscreteTime(int iStep)
{
    if ( this->fixedSteps ) {
        return this->deltaT * iStep;
    } else {
        if ( ( iStep > 0 ) && ( iStep <= discreteTimes.giveSize() ) ) {
            return ( discreteTimes.at(iStep) );
        }
    }

    OOFEM_ERROR("invalid iStep");
    return 0.0;
}


TimeStep *IncrementalLinearStatic :: giveNextStep()
{
    if ( !currentStep ) {
        currentStep.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, 0.,  this->giveDiscreteTime(1), 0) );
    }

    previousStep = std :: move(currentStep);
    double dt = this->giveDiscreteTime(previousStep->giveNumber()+1) - previousStep->giveTargetTime();
    currentStep.reset( new TimeStep(*previousStep, dt) );
    return currentStep.get();
}

void IncrementalLinearStatic :: solveYourself()
{
    this->giveDiscreteTime(1); // ensures numberOfSteps defined
    StructuralEngngModel :: solveYourself();
}

void IncrementalLinearStatic :: solveYourselfAt(TimeStep *tStep)
{
    Domain *d = this->giveDomain(1);
    // Creates system of governing eq's and solves them at given time step

    // Initiates the total displacement to zero.
    if ( tStep->isTheFirstStep() ) {
        for ( auto &dofman : d->giveDofManagers() ) {
            for ( Dof *dof: *dofman ) {
                dof->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, 0.);
                dof->updateUnknownsDictionary(tStep, VM_Total, 0.);
            }
        }

        for ( auto &bc : d->giveBcs() ) {
            ActiveBoundaryCondition *abc;

            if ( ( abc = dynamic_cast< ActiveBoundaryCondition * >(bc.get()) ) ) {
                int ndman = abc->giveNumberOfInternalDofManagers();
                for ( int i = 1; i <= ndman; i++ ) {
                    DofManager *dofman = abc->giveInternalDofManager(i);
                    for ( Dof *dof: *dofman ) {
                        dof->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, 0.);
                        dof->updateUnknownsDictionary(tStep, VM_Total, 0.);
                    }
                }
            }
        }
    }

    // Apply dirichlet b.c's on total values
    for ( auto &dofman : d->giveDofManagers() ) {
        for ( Dof *dof: *dofman ) {
            double tot = dof->giveUnknown( VM_Total, tStep->givePreviousStep() );
            if ( dof->hasBc(tStep) ) {
                tot += dof->giveBcValue(VM_Incremental, tStep);
            }

            dof->updateUnknownsDictionary(tStep, VM_Total, tot);
        }
    }


    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

#ifdef VERBOSE
    OOFEM_LOG_RELEVANT("Solving [step number %8d, time %15e, equations %d]\n", tStep->giveNumber(), tStep->giveTargetTime(), neq);
#endif

    if ( neq == 0 ) { // Allows for fully prescribed/empty problems.
        return;
    }

    incrementOfDisplacementVector.resize(neq);
    incrementOfDisplacementVector.zero();

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling load\n");
#endif
    // Assembling the element part of load vector
    internalLoadVector.resize(neq);
    internalLoadVector.zero();
    this->assembleVector( internalLoadVector, tStep, InternalForceAssembler(),
                         VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );

    loadVector.resize(neq);
    loadVector.zero();
    this->assembleVector( loadVector, tStep, ExternalForceAssembler(),
                         VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );

    loadVector.subtract(internalLoadVector);
    this->updateSharedDofManagers(loadVector, EModelDefaultEquationNumbering(), ReactionExchangeTag);


#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling stiffness matrix\n");
#endif
    stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
    if ( !stiffnessMatrix ) {
        OOFEM_ERROR("sparse matrix creation failed");
    }

    stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
    stiffnessMatrix->zero();
    this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );

#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    NM_Status s = nMethod->solve(*stiffnessMatrix, loadVector, incrementOfDisplacementVector);
    if ( !( s & NM_Success ) ) {
        OOFEM_ERROR("No success in solving system.");
    }
}


double IncrementalLinearStatic :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        int hash = this->giveUnknownDictHashIndx(mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            OOFEM_ERROR("Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    } else {
        OOFEM_ERROR("Only the mode requiresUnknownsDictionaryUpdate() is supported");
    }

    return 0.;
}


int IncrementalLinearStatic :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % 2;
}


void IncrementalLinearStatic :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    // update DOF unknowns dictionary, where
    // unknowns are hold instead of keeping them in global unknowns
    // vectors in engng instances
    // this is necessary, because during solution equation numbers for
    // particular DOFs may changed, and it is necessary to keep them
    // in DOF level.

    double val;
    for ( Dof *dof: *inode ) {
        // skip slave DOFs (only master (primary) DOFs have to be updated).
        if ( !dof->isPrimaryDof() ) {
            continue;
        }
        val = dof->giveUnknown(VM_Total, tStep);
        if ( !dof->hasBc(tStep) ) {
            val += this->incrementOfDisplacementVector.at( dof->__giveEquationNumber() );
        }

        dof->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, val);
        dof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}


void IncrementalLinearStatic :: terminate(TimeStep *tStep)
{
    StructuralEngngModel :: terminate(tStep);
    this->printReactionForces(tStep, 1);
    fflush( this->giveOutputStream() );
}


contextIOResultType IncrementalLinearStatic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file = NULL;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = StructuralEngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    } // ensure consistent records

    return CIO_OK;
}


contextIOResultType IncrementalLinearStatic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    int closeFlag = 0, istep, iversion;
    contextIOResultType iores;
    FILE *file = NULL;
    this->resolveCorrespondingStepNumber(istep, iversion, obj);
    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = StructuralEngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    } // ensure consistent records

    return CIO_OK;
}
} // end namespace oofem
