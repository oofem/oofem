/* $Header: /home/cvs/bp/oofem/sm/src/linearstatic.C,v 1.6.4.1 2004/04/05 15:19:47 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//
// file linearstatic.cc
//

#include "linearstatic.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "dofmanager.h"
#include "elementside.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

#include "verbose.h"
#include "conTable.h"
//#include "skyline.h"
#include "structuralelement.h"
#include "usrdefsub.h"
#include "datastream.h"
//#include "compcol.h"
//#include "dyncompcol.h"
#include "contextioerr.h"

#ifdef __PETSC_MODULE
#include "petscsolver.h"
#endif


LinearStatic :: LinearStatic(int i, EngngModel *_master) : StructuralEngngModel(i, _master), loadVector(), displacementVector()
{
    stiffnessMatrix = NULL;
    ndomains = 1;
    nMethod = NULL;
    initFlag = 1;

#ifdef __PARALLEL_MODE
    commMode = ProblemCommMode__NODE_CUT;
    nonlocalExt = 0;
    communicator = nonlocCommunicator = NULL;
    commBuff = NULL;
#endif
}


LinearStatic :: ~LinearStatic()
{
    delete  stiffnessMatrix;
    if ( nMethod ) {
        delete nMethod;
    }
}



NumericalMethod *LinearStatic :: giveNumericalMethod(TimeStep *atTime)
// only one has reason for LinearStatic
//     - SolutionOfLinearEquations

{
    if ( nMethod ) {
        return nMethod;
    }

    nMethod = :: CreateUsrDefSparseLinSolver(solverType, 1, this->giveDomain(1), this);
    if ( nMethod == NULL ) {
        _error("giveNumericalMethod: linear solver creation failed");
    }

    return nMethod;
}

IRResultType
LinearStatic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StructuralEngngModel :: initializeFrom(ir);
    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_LinearStatic_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_LinearStatic_smtype, "smtype"); // Macro
    sparseMtrxType = ( SparseMtrxType ) val;
    /* The following done in updateAttributes
     * if (this->giveNumericalMethod (giveCurrentStep())) nMethod -> instanciateFrom (ir);
     */

#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        commBuff = new CommunicatorBuff( this->giveNumberOfProcesses() );
        communicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                               this->giveNumberOfProcesses(),
                                               this->commMode);
    }

#endif


    return IRRT_OK;
}



double LinearStatic ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
                                             TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displaacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->giveEquationNumber();
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
    case VM_Incremental:
        if ( displacementVector.isNotEmpty() ) {
            return displacementVector.at(eq);
        } else {
            return 0.;
        }

        // return nMethod-> giveUnknownComponent (LinearEquationSolution, eq);

    default:
        _error("giveUnknownComponent: Unknown is of undefined type for this problem");
    }

    return 0.;
}


TimeStep *LinearStatic :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    //int mstep = 1;
    StateCounterType counter = 1;
    delete previousStep;

    if ( currentStep != NULL ) {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, 1, ( double ) istep, 0., counter);
    // time and dt variables are set eq to 0 for staics - has no meaning
    return currentStep;
}


void LinearStatic :: solveYourself()
{
#ifdef __PARALLEL_MODE
    // force equation numbering before setting up comm maps
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);
#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
#endif

    // set up communication patterns
    // needed only for correct shared rection computation
    communicator->setUpCommunicationMaps(this, true);
    if ( nonlocalExt ) {
        nonlocCommunicator->setUpCommunicationMaps(this, true);
    }

#endif

    StructuralEngngModel :: solveYourself();
}



void LinearStatic :: solveYourselfAt(TimeStep *tStep) {
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step

    if ( initFlag ) {
#ifdef VERBOSE
        OOFEM_LOG_INFO("Assembling stiffness matrix\n");
#endif

        //
        // first step  assemble stiffness Matrix
        //
        /*
         * IntArray* mht = this -> GiveBanWidthVector ();
         * stiffnessMatrix = new Skyline ();
         * stiffnessMatrix ->  checkSizeTowardsBanWidth (mht) ;
         * delete mht;
         */
        stiffnessMatrix = :: CreateUsrDefSparseMtrx(sparseMtrxType); // new Skyline ();
        if ( stiffnessMatrix == NULL ) {
            _error("solveYourselfAt: sparse matrix creation failed");
        }

        //stiffnessMatrix = new DynCompCol ();
        //stiffnessMatrix = new CompCol ();

        stiffnessMatrix->buildInternalStructure(this, 1, EID_MomentumBalance);

        this->assemble( stiffnessMatrix, tStep, EID_MomentumBalance, StiffnessMatrix, this->giveDomain(1) );

        //
        // alocate space for displacementVector
        //
        //displacementVector = new FloatArray (this->giveNumberOfEquations());
        displacementVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
        displacementVector.zero();

#ifdef __PETSC_MODULE
        if ( solverType == ST_Petsc ) {
            this->givePetscContext(1, EID_MomentumBalance)->createVecGlobal(& _loadVec);
            this->givePetscContext(1, EID_MomentumBalance)->createVecGlobal(& _dispVec);
        }

#endif
        initFlag = 0;
    }

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling load\n");
#endif

#ifdef __PETSC_MODULE
    // direct interface to PETSC
    if ( solverType == ST_Petsc ) {
      VecSetOption(_loadVec, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
        this->petsc_assembleVectorFromElements( _loadVec, tStep, EID_MomentumBalance, ElementForceLoadVector, VM_Total, this->giveDomain(1) );
        this->petsc_assembleVectorFromElements( _loadVec, tStep, EID_MomentumBalance, ElementNonForceLoadVector, VM_Total, this->giveDomain(1) );
        this->petsc_assembleVectorFromDofManagers( _loadVec, tStep, EID_MomentumBalance, NodalLoadVector, VM_Total, this->giveDomain(1) );
        VecAssemblyBegin(_loadVec);
        VecAssemblyEnd(_loadVec);
        this->giveNumericalMethod(tStep);
#ifdef VERBOSE
        OOFEM_LOG_INFO("Solving ...\n");
#endif

        //nMethod -> solveYourselfAt(tStep);
        PetscSolver *ps = dynamic_cast< PetscSolver * >( nMethod );
        PetscSparseMtrx *psm = dynamic_cast< PetscSparseMtrx * >( stiffnessMatrix );
        ps->petsc_solve(psm, _loadVec, _dispVec);

        this->givePetscContext(1, EID_MomentumBalance)->scatterG2N(_dispVec, & displacementVector, INSERT_VALUES);
    } else
#endif
    {
        //
        // assembling the element part of load vector
        //
        //loadVector = new FloatArray (this->giveNumberOfEquations());
        loadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
        loadVector.zero();

        this->assembleVectorFromElements( loadVector, tStep, EID_MomentumBalance, ElementForceLoadVector, VM_Total, this->giveDomain(1) );
        this->assembleVectorFromElements( loadVector, tStep, EID_MomentumBalance, ElementNonForceLoadVector, VM_Total, this->giveDomain(1) );

        //
        // assembling the nodal part of load vector
        //
        this->assembleVectorFromDofManagers( loadVector, tStep, EID_MomentumBalance, NodalLoadVector, VM_Total, this->giveDomain(1) );

        //
        // set-up numerical model
        //
        this->giveNumericalMethod(tStep);

        /*
         * nMethod -> setSparseMtrxAsComponent ( LinearEquationLhs , stiffnessMatrix) ;
         * nMethod -> setFloatArrayAsComponent ( LinearEquationRhs , &loadVector) ;
         * nMethod -> setFloatArrayAsComponent ( LinearEquationSolution, &displacementVector) ;
         */
        //
        // call numerical model to solve arised problem
        //
#ifdef VERBOSE
        OOFEM_LOG_INFO("Solving ...\n");
#endif

        //nMethod -> solveYourselfAt(tStep);
        nMethod->solve(stiffnessMatrix, & loadVector, & displacementVector);
    }

    tStep->incrementStateCounter();            // update solution state counter
    //
    // update nodes, elements, etc.
    this->updateYourself( this->giveCurrentStep() );
}

void LinearStatic :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    StructuralEngngModel :: updateYourself(stepN);
}


contextIOResultType LinearStatic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    FILE *file;

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

    if ( ( iores = displacementVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                       // ensure consistent records

    return CIO_OK;
}



contextIOResultType LinearStatic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    int istep, iversion;
    FILE *file;

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

    if ( ( iores = displacementVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                        // ensure consistent records

    return CIO_OK;
}





void
LinearStatic :: terminate(TimeStep *tStep)
{
    StructuralEngngModel :: terminate(tStep);
    this->printReactionForces(tStep, 1);
}



int
LinearStatic :: checkConsistency()
{
    // check internal consistency
    // if success returns nonzero
    int i, nelem;
    Element *ePtr;
    StructuralElement *sePtr;
    Domain *domain = this->giveDomain(1);

    nelem = domain->giveNumberOfElements();
    // check for proper element type

    for ( i = 1; i <= nelem; i++ ) {
        ePtr = domain->giveElement(i);
        sePtr = dynamic_cast< StructuralElement * >( ePtr );
        if ( sePtr == NULL ) {
            _warning2("checkConsistency: element %d has no StructuralElement base", i);
            return 0;
        }
    }

    EngngModel :: checkConsistency();

    return 1;
}


void
LinearStatic :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( giveCurrentStep() )->setDomain( this->giveDomain(1) );
}



void
LinearStatic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    iDof->printSingleOutputAt(stream, atTime, 'd', EID_MomentumBalance, VM_Total);
}


#ifdef __PARALLEL_MODE
int
LinearStatic :: estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType)
{
    int mapSize = commMap.giveSize();
    int i, j, ndofs, count = 0, pcount = 0;
    IntArray locationArray;
    Domain *domain = this->giveDomain(1);
    DofManager *dman;
    Dof *jdof;

    if ( packUnpackType == ProblemCommMode__ELEMENT_CUT ) {
        for ( i = 1; i <= mapSize; i++ ) {
            count += domain->giveDofManager( commMap.at(i) )->giveNumberOfDofs();
        }

        return ( buff.givePackSize(MPI_DOUBLE, 1) * count );
    } else if ( packUnpackType == ProblemCommMode__NODE_CUT ) {
        for ( i = 1; i <= mapSize; i++ ) {
            ndofs = ( dman = domain->giveDofManager( commMap.at(i) ) )->giveNumberOfDofs();
            for ( j = 1; j <= ndofs; j++ ) {
                jdof = dman->giveDof(j);
                if ( jdof->isPrimaryDof() && ( jdof->giveEquationNumber() ) ) {
                    count++;
                } else {
                    pcount++;
                }
            }
        }

        // --------------------------------------------------------------------------------
        // only pcount is relevant here, since only prescribed componnts are exchanged !!!!
        // --------------------------------------------------------------------------------

        return ( buff.givePackSize(MPI_DOUBLE, 1) * pcount );
    } else  if ( packUnpackType == ProblemCommMode__REMOTE_ELEMENT_MODE ) {
        for ( i = 1; i <= mapSize; i++ ) {
            count += domain->giveElement( commMap.at(i) )->estimatePackSize(buff);
        }

        return count;
    }

    return 0;
}
#endif
