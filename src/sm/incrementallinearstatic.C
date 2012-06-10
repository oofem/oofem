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

#include "incrementallinearstatic.h"
#include "timestep.h"
#include "dof.h"
#include "sparsemtrx.h"
#include "dictionr.h"
#include "verbose.h"
#include "usrdefsub.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
IncrementalLinearStatic :: IncrementalLinearStatic(int i, EngngModel *_master) : StructuralEngngModel(i, _master),
    loadVector(), internalLoadVector(), incrementOfDisplacementVector(), discreteTimes()
{
    stiffnessMatrix = NULL;
    ndomains = 1;
    endOfTimeOfInterest = 0;
    nMethod = NULL;
}


IncrementalLinearStatic :: ~IncrementalLinearStatic()
{
    if ( stiffnessMatrix ) {
        delete stiffnessMatrix;
    }

    if ( nMethod ) {
        delete nMethod;
    }
}


NumericalMethod *IncrementalLinearStatic :: giveNumericalMethod(MetaStep *mStep)

{
    if ( nMethod ) {
        return nMethod;
    }

    nMethod = CreateUsrDefSparseLinSolver(solverType, 1, this->giveDomain(1), this);
    if ( nMethod == NULL ) {
        _error("giveNumericalMethod: linear solver creation failed");
    }

    return nMethod;
}


IRResultType IncrementalLinearStatic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    //StructuralEngngModel::initializeFrom (ir);
    IR_GIVE_OPTIONAL_FIELD(ir, discreteTimes, IFT_IncrementalLinearStatic_prescribedtimes, "prescribedtimes");
    if ( discreteTimes.giveSize() > 0 ) {
        numberOfSteps = discreteTimes.giveSize();
        endOfTimeOfInterest = discreteTimes.at(discreteTimes.giveSize());
        fixedSteps = false;
    } else {
        deltaT = 1.0;
        IR_GIVE_OPTIONAL_FIELD(ir, deltaT, IFT_IncrementalLinearStatic_deltat, "deltat");
        IR_GIVE_FIELD(ir, numberOfSteps, IFT_EngngModel_nsteps, "nsteps");
        endOfTimeOfInterest = deltaT*numberOfSteps;
        fixedSteps = true;
    }
    IR_GIVE_OPTIONAL_FIELD(ir, endOfTimeOfInterest, IFT_IncrementalLinearStatic_endoftimeofinterest, "endoftimeofinterest");

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_IncrementalLinearStatic_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_IncrementalLinearStatic_smtype, "smtype"); // Macro
    sparseMtrxType = ( SparseMtrxType ) val;

    return IRRT_OK;
}


double IncrementalLinearStatic :: giveDiscreteTime(int iStep)
{
    if (this->fixedSteps) {
       return this->deltaT*iStep;
    } else {
        if ( ( iStep > 0 ) && ( iStep <= discreteTimes.giveSize() ) ) {
            return ( discreteTimes.at(iStep) );
        }
    }

    _error("giveDiscreteTime: invalid iStep");
    return 0.0;
}


TimeStep *IncrementalLinearStatic :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    int mstepNum = 1;
    double dt = this->giveDiscreteTime(istep);
    StateCounterType counter = 1;

    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        dt = this->giveDiscreteTime(istep) - this->giveDiscreteTime(istep - 1);
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    if (previousStep != NULL){
        delete previousStep;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, mstepNum, this->giveDiscreteTime(istep), dt, counter);
    return currentStep;
}

void IncrementalLinearStatic :: solveYourself()
{
    this->giveDiscreteTime(1); // ensures numberOfSteps defined
    StructuralEngngModel :: solveYourself();
}

void IncrementalLinearStatic :: solveYourselfAt(TimeStep *tStep)
{
    // Creates system of governing eq's and solves them at given time step

    // Initiates the total displacement to zero.
    if ( tStep->isTheFirstStep() ) {
        Domain *d = this->giveDomain(1);
        for ( int i = 1; i <= d->giveNumberOfDofManagers(); i++ ) {
            DofManager *dofman = d->giveDofManager(i);
            for ( int j = 1; j <= dofman->giveNumberOfDofs(); j++ ) {
                dofman->giveDof(j)->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total_Old, 0);
                dofman->giveDof(j)->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total, 0);
                // This is actually redundant now;
                //dofman->giveDof(j)->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Incremental, 0);
            }
        }
    }

    // Apply dirichlet b.c's on total values
    Domain *d = this->giveDomain(1);
    for ( int i = 1; i <= d->giveNumberOfDofManagers(); i++ ) {
        DofManager *dofman = d->giveDofManager(i);
        for ( int j = 1; j <= dofman->giveNumberOfDofs(); j++ ) {
            Dof *d = dofman->giveDof(j);
            double tot = d->giveUnknown(EID_MomentumBalance, VM_Total_Old, tStep);
            if ( d->hasBc(tStep) ) {
                tot += d->giveBcValue(VM_Incremental, tStep);
            }

            d->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total, tot);
        }
    }


#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    int neq = this->giveNumberOfEquations(EID_MomentumBalance);

    if (neq == 0) { // Allows for fully prescribed/empty problems.
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
    this->assembleVector( internalLoadVector, tStep, EID_MomentumBalance, InternalForcesVector,
                          VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );

    loadVector.resize(neq);
    loadVector.zero();
    this->assembleVector( loadVector, tStep, EID_MomentumBalance, ExternalForcesVector,
                          VM_Total, EModelDefaultEquationNumbering(), this->giveDomain(1) );

    loadVector.subtract(internalLoadVector);

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling stiffness matrix\n");
#endif
    if ( stiffnessMatrix ) {
        delete stiffnessMatrix;
    }

    stiffnessMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
    if ( stiffnessMatrix == NULL ) {
        _error("solveYourselfAt: sparse matrix creation failed");
    }

    stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );
    stiffnessMatrix->zero();
    this->assemble( stiffnessMatrix, tStep, EID_MomentumBalance, StiffnessMatrix,
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );

#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif
    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    NM_Status s = nMethod->solve(stiffnessMatrix, & loadVector, & incrementOfDisplacementVector);
    if ( !(s & NM_Success) ) {
        OOFEM_ERROR("IncrementalLinearStatic :: solverYourselfAt - No success in solving system.");
    }
}


void IncrementalLinearStatic :: updateYourself(TimeStep *stepN)
{
    // updates internal state to reached one
    this->updateInternalState(stepN);
    StructuralEngngModel :: updateYourself(stepN);
}



double IncrementalLinearStatic :: giveUnknownComponent(EquationID type, ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    if ( type != EID_MomentumBalance ) { // heat and mass concetration vector
        OOFEM_ERROR2( "giveUnknownComponent: EquationID %s is undefined for this problem", __EquationIDToString(type) );
        return 0.;
    }

    if ( this->requiresUnknownsDictionaryUpdate() ) {
        int hash = this->giveUnknownDictHashIndx(type, mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    } else {
        OOFEM_ERROR("Only the mode requiresUnknownsDictionaryUpdate() is supported");
    }

    return 0.;
}



void IncrementalLinearStatic :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    // update DOF unknowns dictionary, where
    // unknowns are hold instead of keeping them in global unknowns
    // vectors in engng instances
    // this is necessary, because during solution equation numbers for
    // particular DOFs may changed, and it is necessary to keep them
    // in DOF level.

    int ndofs = inode->giveNumberOfDofs();
    Dof *iDof;
    double val;
    for ( int i = 1; i <= ndofs; i++ ) {
        iDof = inode->giveDof(i);
        val = iDof->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
        if ( !iDof->hasBc(tStep) ) {
            val += this->incrementOfDisplacementVector.at( iDof->__giveEquationNumber() );
        }

        iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total_Old, val);
        iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total, val);
    }
}


void IncrementalLinearStatic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    iDof->printSingleOutputAt(stream, atTime, 'd', EID_MomentumBalance, VM_Total);
}


void IncrementalLinearStatic :: terminate(TimeStep *tStep)
{
    StructuralEngngModel :: terminate(tStep);
    this->printReactionForces(tStep, 1);
}


contextIOResultType IncrementalLinearStatic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file;

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
    FILE *file;
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
