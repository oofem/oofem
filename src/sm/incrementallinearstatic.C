/* $Header: /home/cvs/bp/oofem/sm/src/incrementallinearstatic.C,v 1.5.4.1 2004/04/05 15:19:47 bp Exp $ */
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
// file incrementallinearstatic.cc
//

#include "incrementallinearstatic.h"
#include "nummet.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "node.h"
#include "dof.h"
#include "boundary.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#include "verbose.h"
#include "conTable.h"
#include "usrdefsub.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
IncrementalLinearStatic :: IncrementalLinearStatic(int i, EngngModel *_master) : StructuralEngngModel(i, _master),
    incrementOfLoadVector(), incrementOfDisplacementVector(), totalDisplacementVector(), discreteTimes()
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

    //if (totalDisplacementVector) delete totalDisplacementVector;
    if ( nMethod ) {
        delete nMethod;
    }
}


NumericalMethod *
IncrementalLinearStatic :: giveNumericalMethod(TimeStep *)
// only one has reason for LinearStatic
//     - SolutionOfLinearEquations

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


/*
 * This function is not necesarry - we use unknowns DOF dictionary instead
 * so particular dofs are not asking emodel, but use their own dictionary
 * - this dictionary is set maintained by emodel->updateDofUnknownsDictionary()
 */
/*
 * double
 * IncrementalLinearStatic ::  giveUnknownComponent (CharType chc, int eq)
 * // returns unknown quantity like displaacement, velocity of equation eq
 * // This function translates this request to numerical method language
 * {
 * switch (chc)
 * {
 * case IncrementOfDisplacementVector:
 * case TotalIncrementOfDisplacementVector:
 *  if (incrementOfDisplacementVector) return incrementOfDisplacementVector -> at(eq);
 *  else return 0.;
 *  break;
 *
 * default:
 *  _error ("giveUnknownComponent: Unknown is of undefined type for this problem");
 * }
 * return 0.;
 * }
 */


IRResultType
IncrementalLinearStatic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //StructuralEngngModel::initializeFrom (ir);
    IR_GIVE_FIELD(ir, endOfTimeOfInterest, IFT_IncrementalLinearStatic_endoftimeofinterest, "endoftimeofinterest"); // Macro
    IR_GIVE_FIELD(ir, discreteTimes, IFT_IncrementalLinearStatic_prescribedtimes, "prescribedtimes"); // Macro
    numberOfSteps = discreteTimes.giveSize();
    //numberOfSteps = readInteger (initString,"nsteps");
    //if (numberOfSteps <= 0) error ("instanciateFrom: nsteps not specified, bad format");

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_IncrementalLinearStatic_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_IncrementalLinearStatic_smtype, "smtype"); // Macro
    sparseMtrxType = ( SparseMtrxType ) val;

    return IRRT_OK;
}

double
IncrementalLinearStatic :: giveDiscreteTime(int iStep)
{
    /*
     * this function returns time valid for iStep timeStep, used in intergation
     * of structure response.
     * This functions, when invoked for the first time, generates table of times.
     * Times in this table are generated according to:
     * .) if there exists some load time function with abrupt change,
     * then time just before and just at abrupt is included in table.
     * .) between these steps under constant loads (or when load changes continuously)
     * we use progressively increasing time step. They are best choosen so that time
     * step be keept constant in the log (t-t') scale
     */
    if ( ( iStep > 0 ) && ( iStep <= discreteTimes.giveSize() ) ) {
        return ( discreteTimes.at(iStep) );
    }

    _error("giveDiscreteTime: invalid iStep");
    /*
     * int i;
     * int nLoadTimeFunct = this->giveDomain()->giveNumberOfLoadTimeFunctions();
     * FloatArray *abruptChanges;
     *
     * dSize = 100;
     * dIndex = 0;
     * discreteLTFTimes = new FloatArray (100);
     *
     * for (i=1; i<= nLoadTimeFunct; i++) {
     * abruptChanges = domain->giveLoadTimeFunction(i)->GiveAbruptChangeTimes();
     * for (j=1; j <= abruptChanges->giveSize(); j++) {
     * if ((abruptChanges->at(j) >0.) && (abruptChanges->at(j) < endOfTimeOfInterest)) {
     *  if (dIndex == dSize-1)
     *   discreteLTFTimes->growTo (dSize+=20);
     *  discreteLTFTimes->at(++dIndex) = abruptChanges->at(j);
     * }
     * }
     * }
     * // sort results
     * discreteLTFTimes -> sortYourself();
     * // generate log time scale
     */
    return 0.0;
}


TimeStep *
IncrementalLinearStatic :: giveNextStep()
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

    delete previousStep;
    previousStep = currentStep;
    currentStep = new TimeStep(istep, this, mstepNum, this->giveDiscreteTime(istep), dt, counter);
    return currentStep;
}

void
IncrementalLinearStatic :: solveYourself()
{
    this->giveNumericalMethod( giveCurrentStep() );
    this->giveDiscreteTime(1); // ensure numberOfSteps defined
    StructuralEngngModel :: solveYourself();
}

/*
 * void
 * IncrementalLinearStatic :: solveYourself ()
 * {
 * int imstep, jstep;
 * int smstep=1, sjstep=1;
 * MetaStep* activeMStep;
 *
 * if (this->currentStep) {
 * smstep = this->currentStep->giveMetaStepNumber();
 * sjstep = this->giveMetaStep(smstep)->giveStepRelativeNumber(this->currentStep->giveNumber()) + 1;
 * }
 *
 * this -> giveNumericalMethod (giveCurrentStep());
 * this -> giveDiscreteTime(1); // ensure numberOfSteps defined
 * for (imstep = smstep; imstep<= nMetaSteps; imstep++) {
 * activeMStep = this->giveMetaStep(imstep);
 * for (jstep = sjstep; jstep <= activeMStep->giveNumberOfSteps(); jstep++)
 * {
 *  this->giveNextStep();
 *  // update state ccording to new meta step
 *  if (jstep == sjstep) this->initMetaStepAttributes (this->giveCurrentStep());
 *  // for current time step force equation renumbering
 *  // to support changes of static system uncomment next line
 *  this->forceEquationNumbering();
 *  this->solveYourselfAt(this->giveCurrentStep());
 *  this->terminate (this->giveCurrentStep());
 * }
 * }
 * }
 */

void
IncrementalLinearStatic :: solveYourselfAt(TimeStep *tStep) {
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step

    // if (tStep->giveNumber() == this->giveNumberOfFirstStep()) {
    //
    // first step  assemble stiffness Matrix
    //
    /*
     * IntArray* mht = this -> GiveBanWidthVector ();
     * stiffnessMatrix = new Skyline ();
     * stiffnessMatrix ->  checkSizeTowardsBanWidth (mht) ;
     * delete mht;
     */

#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTime() );
#endif

    if ( stiffnessMatrix ) {
        delete stiffnessMatrix;
    }

    stiffnessMatrix = CreateUsrDefSparseMtrx(sparseMtrxType); // new Skyline ();
    if ( stiffnessMatrix == NULL ) {
        _error("solveYourselfAt: sparse matrix creation failed");
    }

    stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );

    //if (totalDisplacementVector) delete totalDisplacementVector;
    //totalDisplacementVector = new FloatArray (this->giveNumberOfEquations());


#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling stiffness matrix\n");
#endif
    stiffnessMatrix->zero();     // zero stiffness matrix
    this->assemble( stiffnessMatrix, tStep, EID_MomentumBalance, StiffnessMatrix,
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );
    //
    // alocate space for displacementVector
    //
    incrementOfDisplacementVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
    incrementOfDisplacementVector.zero();

#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling load\n");
#endif

    //
    // assembling the element part of load vector
    //
    //incrementOfLoadVector = new FloatArray (this->giveNumberOfEquations());
    incrementOfLoadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
    incrementOfLoadVector.zero();
    this->assembleVectorFromElements( incrementOfLoadVector, tStep, EID_MomentumBalance, ElementForceLoadVector,
                                     VM_Incremental, EModelDefaultEquationNumbering(), this->giveDomain(1) );
    this->assembleVectorFromElements( incrementOfLoadVector, tStep, EID_MomentumBalance, ElementNonForceLoadVector,
                                     VM_Incremental, EModelDefaultEquationNumbering(), this->giveDomain(1) );

    //
    // assembling the nodal part of load vector
    //
    this->assembleVectorFromDofManagers( incrementOfLoadVector, tStep, EID_MomentumBalance, NodalLoadVector,
                                        VM_Incremental, EModelDefaultEquationNumbering(), this->giveDomain(1) );

    //
    // set-up numerical model
    //
    /*
     * nMethod -> setSparseMtrxAsComponent ( LinearEquationLhs , stiffnessMatrix) ;
     * nMethod -> setFloatArrayAsComponent ( LinearEquationRhs , &incrementOfLoadVector) ;
     * nMethod -> setFloatArrayAsComponent ( LinearEquationSolution, &incrementOfDisplacementVector) ;
     */
    //
    // call numerical model to solve arised problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("Solving ...\n");
#endif

    //nMethod -> solveYourselfAt(tStep);
    nMethod->solve(stiffnessMatrix, & incrementOfLoadVector, & incrementOfDisplacementVector);
    //
    // update nodes, elements, etc.
    this->updateYourself( this->giveCurrentStep() );
}


void
IncrementalLinearStatic :: updateYourself(TimeStep *stepN)
{
    // updates internal state to reached one
    // totalDisplacementVector -> add (incrementOfDisplacementVector);
    this->updateInternalState(stepN);
    StructuralEngngModel :: updateYourself(stepN);
}



contextIOResultType
IncrementalLinearStatic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    int closeFlag = 0;
    contextIOResultType iores;
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

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                        // ensure consistent records

    return CIO_OK;
}



contextIOResultType
IncrementalLinearStatic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    int closeFlag = 0, istep, iversion;
    contextIOResultType iores;
    FILE *file;
    /* if (totalDisplacementVector == NULL) {
     * // allocate some meanigless space - will be overriden by restoreYourself()
     * totalDisplacementVector = new FloatArray (1);
     * }
     */
    this->resolveCorrespondingStepNumber(istep, iversion, obj);
    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR);                                                             // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    // save element context

    if ( ( iores = StructuralEngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                       // ensure consistent records

    return CIO_OK;
}



void
IncrementalLinearStatic :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    // update DOF unknowns dictionary, where
    // unknowns are hold instead of keeping them in global unknowns
    // vectors in engng instancies
    // this is necessary, because during solution equation numbers for
    // particular DOFs may changed, and it is necesary to keep them
    // in DOF level.

    int i, ndofs = inode->giveNumberOfDofs();
    Dof *iDof;
    double val;

    if ( tStep->giveNumber() == this->giveNumberOfFirstStep() ) {
        // special code for first step - increments are also total values
        for ( i = 1; i <= ndofs; i++ ) {
            iDof = inode->giveDof(i);
            if ( iDof->hasBc(tStep) ) { // bound . cond.
                // val = iDof -> giveBcValue() -> give(DisplacementVector,IncrementalMode,tStep) ;
                val = iDof->giveBcValue(VM_Incremental, tStep);
            } else {
                val = this->incrementOfDisplacementVector.at( iDof->__giveEquationNumber() );
            }

            iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Incremental, val);
            iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total, val);
        }
    } else {
        for ( i = 1; i <= ndofs; i++ ) {
            iDof = inode->giveDof(i);
            if ( iDof->hasBc(tStep) ) { // bound . cond.
                //val = iDof -> giveBcValue() -> give(DisplacementVector,IncrementalMode,tStep) ;
                val = iDof->giveBcValue(VM_Incremental, tStep);
            } else {
                val = this->incrementOfDisplacementVector.at( iDof->__giveEquationNumber() );
            }

            iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Incremental, val);
            val += iDof->giveUnknown(EID_MomentumBalance, VM_Total, tStep);

            iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total, val);
        }
    }
}



void
IncrementalLinearStatic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    iDof->printSingleOutputAt(stream, atTime, 'd', EID_MomentumBalance, VM_Total);
}
} // end namespace oofem
