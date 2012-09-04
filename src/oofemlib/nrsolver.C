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

#include "nrsolver.h"
#include "verbose.h"
#include "timestep.h"
#include "mathfem.h"
// includes for ddc - not very clean (NumMethod knows what is "node" and "dof")
#include "node.h"
#include "element.h"
#include "dof.h"
#include "loadtime.h"
#include "linesearch.h"
#include "usrdefsub.h"
#ifdef __PETSC_MODULE
 #include "petscsolver.h"
 #include "petscsparsemtrx.h"
 #include "petscordering.h"
#endif

namespace oofem {
#define nrsolver_ERROR_NORM_SMALL_NUM 1.e-6
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e20
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10


NRSolver :: NRSolver(int i, Domain *d, EngngModel *m, EquationID ut) :
    SparseNonLinearSystemNM(i, d, m, ut), prescribedDofs(), prescribedDofsValues()
{
    //
    // constructor
    //
    nsmax  = 60;     // default maximum number of sweeps allowed
    //Psi    = 0.1;       // displacement control on
    deltaL = 1.0;
    solved = 0;
    NR_Mode = NR_OldMode = nrsolverModifiedNRM;
    NR_ModeTick = -1; // do not switch to calm_NR_OldMode
    MANRMSteps = 0;
    numberOfPrescribedDofs = 0;
    prescribedDofsFlag = false;
    prescribedEqsInitFlag = false;
    prescribedDisplacementLTF = 0;
    linSolver = NULL;
    linesearchSolver = NULL;
    lsFlag = 0; // no line-search
    smConstraintVersion = 0;
#ifdef __PETSC_MODULE
    prescribedEgsIS_defined = false;
#endif
}


NRSolver :: ~NRSolver()
{
    //
    // destructor
    //
    if ( linSolver ) {
        delete linSolver;
    }

    if ( linesearchSolver ) {
        delete linesearchSolver;
    }

#ifdef __PETSC_MODULE
 #ifdef __PARALLEL_MODE
    if ( prescribedEgsIS_defined ) {
        ISDestroy(&prescribedEgsIS);
    }
 #endif
#endif
}


IRResultType
NRSolver :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    nsmax = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nsmax, IFT_NRSolver_maxiter, "maxiter"); // Macro
    if ( nsmax < 30 ) {
        nsmax = 30;
    }

    minIterations = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, minIterations, IFT_NRSolver_miniterations, "miniter"); // Macro



    minStepLength = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, IFT_NRSolver_minsteplength, "minsteplength"); // Macro

    // read if MANRM method is used
    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, IFT_NRSolver_manrmsteps, "manrmsteps"); // Macro
    if ( MANRMSteps > 0 ) {
        NR_Mode = NR_OldMode = nrsolverAccelNRM;
    } else {
        NR_Mode = nrsolverModifiedNRM;
    }

    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NRSolver_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) _val;
    this->giveLinearSolver()->initializeFrom(ir);

    // read relative error tolerances of the solver fo each cc
    IR_GIVE_FIELD(ir, rtol, IFT_NRSolver_rtolv, "rtolv"); // Macro

    prescribedDofs.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDofs, IFT_NRSolver_ddm, "ddm"); // Macro
    prescribedDofsValues.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDofsValues, IFT_NRSolver_ddv, "ddv"); // Macro
    prescribedDisplacementLTF = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDisplacementLTF, IFT_NRSolver_ddltf, "ddltf"); // Macro

    numberOfPrescribedDofs = prescribedDofs.giveSize() / 2;
    if ( numberOfPrescribedDofs != prescribedDofsValues.giveSize() ) {
        _error("instanciateFrom direct displacement mask size mismatch");
    }

    if ( numberOfPrescribedDofs ) {
        prescribedDofsFlag = true;
    } else {
        prescribedDofsFlag = false;
    }

    this->lsFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, lsFlag, IFT_NRSolver_linesearch, "linesearch"); // Macro

    if ( this->lsFlag ) {
        this->giveLineSearchSolver()->initializeFrom(ir);
    }

    /* initialize optional dof groups for convergence criteria evaluation */
    this->nccdg = 0; // default, no dof cc group, all norms evaluated for all dofs
    IR_GIVE_OPTIONAL_FIELD(ir, nccdg, IFT_NRSolver_nccdg, "nccdg"); // Macro
    if ( nccdg >= 1 ) {
        int _i, _j;
        IntArray _val;
        char name [ 12 ];
        // create an empty set
        __DofIDSet _set;
        // resize gof group vector
        this->ccDofGroups.resize(nccdg, _set);
        for ( _i = 0; _i < nccdg; _i++ ) {
            sprintf(name, "ccdg%d", _i + 1);
            // read dof group as int array under ccdg# keyword
            IR_GIVE_FIELD(ir, _val, IFT_NRSolver_ccdg, name); // Macro
            // convert aray into set
            for ( _j = 1; _j <= _val.giveSize(); _j++ ) {
                ccDofGroups.at(_i).insert( (DofIDItem)_val.at(_j) );
            }
        }

        // read relative error tolerances of the solver fo each cc
        // if common rtolv provided, set to this tolerace both rtolf and rtold
        IR_GIVE_OPTIONAL_FIELD(ir, rtolf, IFT_NRSolver_rtolv, "rtolv"); // Macro
        rtold = rtolf;
        // read optional force and displacement tolerances
        IR_GIVE_OPTIONAL_FIELD(ir, rtolf, IFT_NRSolver_rtolf, "rtolf"); // Macro
        IR_GIVE_OPTIONAL_FIELD(ir, rtold, IFT_NRSolver_rtold, "rtold"); // Macro

        if ( ( rtolf.giveSize() != nccdg ) || ( rtold.giveSize() != nccdg ) ) {
            _error2("Incompatible size of rtolf or rtold params, expected size %d (nccdg)", nccdg);
        }
    } else {
        nccdg = 0;
        double _rtol = 1.e-3; // default tolerance
        rtolf.resize(1);
        rtold.resize(1);
        // read relative error tolerances of the solver
        // if common rtolv provided, set to this tolerance both rtolf and rtold
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, IFT_NRSolver_rtolf, "rtolv"); // Macro
        rtolf.at(1) = rtold.at(1) = _rtol;
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, IFT_NRSolver_rtolf, "rtolf"); // Macro
        rtolf.at(1) = _rtol;
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, IFT_NRSolver_rtold, "rtold"); // Macro
        rtold.at(1) = _rtol;
    }

    return IRRT_OK;
}


contextIOResultType
NRSolver :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    return CIO_OK;
}


contextIOResultType
NRSolver :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    return CIO_OK;
}


NM_Status
NRSolver :: solve(SparseMtrx *k, FloatArray *R, FloatArray *R0,
                  FloatArray *X, FloatArray *dX, FloatArray *F,
                  double &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tNow)
//
// this function solve the problem of the unbalanced equilibrium
// using NR scheme
//
//
{
    FloatArray rhs, ddX, RT, XInitial; // residual, iteration increment of solution, total external force, intial value of solution.
    double RRT;
    int neq = X->giveSize();
    int irest = 0;
    NM_Status status;
    bool converged, errorOutOfRangeFlag;
#ifdef __PARALLEL_MODE
 #ifdef __PETSC_MODULE
    PetscContext *parallel_context  = engngModel->givePetscContext(1, ut); ///@todo hard wired domain no 1
 #endif
#endif

    if ( engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("NRSolver:     Iteration       ForceError      DisplError                    \n");
        OOFEM_LOG_INFO("----------------------------------------------------------------------------\n");
    }

    XInitial = * X; // Stored in case of divergence.
    l = 1.0;

    status = NM_None;
    this->giveLinearSolver();

    // compute total load R = R+R0
    RT = * R;
    if ( R0 ) {
        RT.add(*R0);
    }

#ifdef __PARALLEL_MODE
    RRT = parallel_context->norm(RT);
    RRT *= RRT;
#else
    RRT = RT.computeSquaredNorm();
#endif

restart:
    dX->zero();
    ddX.resize(neq);

    // Fetch the matrix before evaluating internal forces.
    // This is intentional, since its a simple way to drastically increase convergence for nonlinear problems.
    // (This old tangent is just used
    engngModel->updateComponent(tNow, NonLinearLhs, domain);
    if ( this->prescribedDofsFlag ) {
        if ( !prescribedEqsInitFlag ) {
            this->initPrescribedEqs();
        }
        applyConstraintsToStiffness(k);
    }

    for (nite = 0; true; ++nite) {
        // Compute the residual
        engngModel->updateComponent(tNow, InternalRhs, domain);
        rhs.beDifferenceOf(RT, *F);

        if ( this->prescribedDofsFlag ) {
            this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tNow);
        }

        // convergence check
        converged = this->checkConvergence(RT, * F, rhs, ddX, * X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag, tNow);

        if (converged && nite >= minIterations ) {
            break;
        } else if ( ( nite >= nsmax ) || errorOutOfRangeFlag ) {
            if ( irest <= NRSOLVER_MAX_RESTARTS ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length

                /*
                 * double time;
                 * time = tNow->giveTime() - tNow->giveTimeIncrement()*(1.-NRSOLVER_RESET_STEP_REDUCE) ;
                 * deltaL =  deltaL * NRSOLVER_RESET_STEP_REDUCE ;
                 * if (deltaL < minStepLength)  deltaL = minStepLength;
                 *
                 * tNow -> setTime(time);
                 * tNow -> setTimeIncrement(tNow->giveTimeIncrement()*NRSOLVER_RESET_STEP_REDUCE);
                 * tNow->incrementStateCounter();              // update solution state counter
                 */

                // restore previous total displacement vector
                *X = XInitial;
                // reset all changes fro previous equilibrium state
                engngModel->initStepIncrements();

                OOFEM_WARNING2("NRSolver:  Iteration Reset: %d\n",irest);

                NR_OldMode  = NR_Mode;
                NR_Mode     = nrsolverFullNRM;
                NR_ModeTick = NRSOLVER_DEFAULT_NRM_TICKS;
                irest++;
                goto restart;
            } else {
                status = NM_NoSuccess;
                OOFEM_WARNING2("NRSolver:  Convergence not reached after %d iterations", nsmax);
                break;
            }
        }

        if ( nite > 0 ) {
            if ( ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
                engngModel->updateComponent(tNow, NonLinearLhs, domain);
                applyConstraintsToStiffness(k);
            }
        }

        if ( ( nite == 0 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equilibrate current state
            rhs.zero();
            R->zero();
            ddX = rhs;
        } else {
            linSolver->solve(k, & rhs, & ddX);
        }

        //
        // update solution
        //
        if ( this->lsFlag && ( nite > 0 ) ) { // Why not nite == 0 ?
            // line search
            LineSearchNM :: LS_status status;
            double eta;
            this->giveLineSearchSolver()->solve(X, & ddX, F, R, R0, prescribedEqs, 1.0, eta, status, tNow);
        }
        X->add(ddX);
        dX->add(ddX);
        tNow->incrementStateCounter(); // update solution state counter
    }
    //
    // end of iteration
    //

    status |= NM_Success;
    solved = 1;

    // Modify Load vector to include "quasi reaction"
    if ( R0 ) {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R->at( prescribedEqs.at(i) ) = F->at( prescribedEqs.at(i) ) - R0->at( prescribedEqs.at(i) ) - R->at( prescribedEqs.at(i) );
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R->at( prescribedEqs.at(i) ) = F->at( prescribedEqs.at(i) ) - R->at( prescribedEqs.at(i) );
        }
    }

    this->lastReactions.resize(numberOfPrescribedDofs);

#ifdef VERBOSE
    if (numberOfPrescribedDofs) {
        // print quasi reactions if direct displacement control used
        OOFEM_LOG_INFO("\n");
        OOFEM_LOG_INFO("NRSolver:     Quasi reaction table                                 \n");
        OOFEM_LOG_INFO("NRSolver:     Node            Dof             Displacement    Force\n");
        double reaction;
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            reaction = R->at( prescribedEqs.at(i) );
            if ( R0 ) {
                reaction += R0->at( prescribedEqs.at(i) );
            }
            lastReactions.at(i) = reaction;
            OOFEM_LOG_INFO("NRSolver:     %-15d %-15d %-+15.5e %-+15.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i),
                           X->at( prescribedEqs.at(i) ), reaction);
        }
        OOFEM_LOG_INFO("\n");
    }
#endif

    return status;
}


SparseLinearSystemNM *
NRSolver :: giveLinearSolver()
{
    if ( linSolver ) {
        if ( linSolver->giveLinSystSolverType() == solverType ) {
            return linSolver;
        } else {
            delete linSolver;
        }
    }

    linSolver = CreateUsrDefSparseLinSolver(solverType, 1, domain, engngModel);
    if ( linSolver == NULL ) {
        _error("giveLinearSolver: linear solver creation failed");
    }

    return linSolver;
}


LineSearchNM *
NRSolver :: giveLineSearchSolver()
{
    if ( linesearchSolver == NULL ) {
        linesearchSolver = new LineSearchNM(1, this->giveDomain(), engngModel);
    }

    return linesearchSolver;
}


void
NRSolver :: initPrescribedEqs()
{
    EModelDefaultEquationNumbering dn;
#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
 #if defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE )
    PetscContext *parallel_context = engngModel->givePetscContext(1, ut); ///@todo hard wired domain no 1
 #endif
    int jglobnum, count = 0, ndofman = domain->giveNumberOfDofManagers();
    int i, j, inode, idof;
    IntArray localPrescribedEqs(numberOfPrescribedDofs);

    for ( j = 1; j <= ndofman; j++ ) {
 #if defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE )
        if ( !parallel_context->isLocal(domain->giveNode(j)) ) {
            continue;
        }
 #endif
        jglobnum = domain->giveNode(j)->giveGlobalNumber();
        for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
            inode = prescribedDofs.at(2 * i - 1);
            idof  = prescribedDofs.at(2 * i);
            if ( inode == jglobnum ) {
                localPrescribedEqs.at(++count) = domain->giveNode(j)->giveDof(idof)->giveEquationNumber(dn);
                continue;
            }
        }
    }

    prescribedEqs.resize(count);
    for ( i = 1; i <= count; i++ ) {
        prescribedEqs.at(i) = localPrescribedEqs.at(i);
    }

    numberOfPrescribedDofs = count;
#else
    int i, inode, idof;
    prescribedEqs.resize(numberOfPrescribedDofs);
    for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
        inode = prescribedDofs.at(2 * i - 1);
        idof  = prescribedDofs.at(2 * i);
        prescribedEqs.at(i) = domain->giveNode(inode)->giveDof(idof)->giveEquationNumber(dn);
    }

#endif
    this->prescribedEqsInitFlag = true;
}


void
NRSolver :: applyConstraintsToStiffness(SparseMtrx *k)
{
    int i;

    if ( this->smConstraintVersion == k->giveVersion() ) {
        return;
    }

#if 0
 #ifdef __PETSC_MODULE
    if ( solverType == ST_Petsc ) {
        PetscScalar diagVal = 1.0;
        if ( k->giveType() != SMT_PetscMtrx ) {
            _error("applyConstraintsToStiffness: PetscSparseMtrx Expected");
        }

        PetscSparseMtrx *lhs = ( PetscSparseMtrx * ) k;

        if ( !prescribedEgsIS_defined ) {
            IntArray eqs;
  #ifdef __PARALLEL_MODE
            PetscNatural2GlobalOrdering *n2lpm = engngModel->givePetscContext(1)->giveN2Gmap();
            int s = prescribedEqs.giveSize();
            eqs.resize(s);
            for ( i = 1; i <= s; i++ ) {
                eqs.at(i) = n2lpm->giveNewEq( prescribedEqs.at(i) );
            }

            ISCreateGeneral(PETSC_COMM_WORLD, s, eqs.givePointer(), & prescribedEgsIS);
            //ISView(prescribedEgsIS,PETSC_VIEWER_STDOUT_WORLD);
  #else
            eqs.resize(numberOfPrescribedDofs);
            for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
                eqs.at(i) = prescribedEqs.at(i) - 1;
            }

            ISCreateGeneral(PETSC_COMM_SELF, numberOfPrescribedDofs, eqs.givePointer(), & prescribedEgsIS);
            //ISView(prescribedEgsIS,PETSC_VIEWER_STDOUT_SELF);
  #endif
            prescribedEgsIS_defined = true;
        }

        //MatView(*(lhs->giveMtrx()),PETSC_VIEWER_STDOUT_WORLD);
        MatZeroRows(* ( lhs->giveMtrx() ), prescribedEgsIS, & diagVal);
        //MatView(*(lhs->giveMtrx()),PETSC_VIEWER_STDOUT_WORLD);
        if ( numberOfPrescribedDofs ) {
            this->smConstraintVersion = k->giveVersion();
        }

        return;
    }

 #endif // __PETSC_MODULE
#else
 #ifdef __PETSC_MODULE
    if ( solverType == ST_Petsc ) {
        if ( k->giveType() != SMT_PetscMtrx ) {
            _error("applyConstraintsToStiffness: PetscSparseMtrx Expected");
        }

        PetscSparseMtrx *lhs = ( PetscSparseMtrx * ) k;

        Vec diag;
        PetscScalar *ptr;
        int eq;

        PetscContext *parallel_context = engngModel->givePetscContext(1, ut);
        parallel_context->createVecGlobal(& diag);
        MatGetDiagonal(* lhs->giveMtrx(), diag);
        VecGetArray(diag, & ptr);
        for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
            eq = prescribedEqs.at(i) - 1;
            MatSetValue(* ( lhs->giveMtrx() ), eq, eq, ptr [ eq ] * 1.e6, INSERT_VALUES);
        }

        MatAssemblyBegin(* lhs->giveMtrx(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(* lhs->giveMtrx(), MAT_FINAL_ASSEMBLY);
        VecRestoreArray(diag, & ptr);
        VecDestroy(&diag);
        if ( numberOfPrescribedDofs ) {
            this->smConstraintVersion = k->giveVersion();
        }

        return;
    }

 #endif // __PETSC_MODULE
#endif
    for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
        k->at( prescribedEqs.at(i), prescribedEqs.at(i) ) *= 1.e6;
    }

    if ( numberOfPrescribedDofs ) {
        this->smConstraintVersion = k->giveVersion();
    }
}


void
NRSolver :: applyConstraintsToLoadIncrement(int nite, const SparseMtrx *k, FloatArray &R,
                                            referenceLoadInputModeType rlm, TimeStep *atTime)
{
    int eq, i;
    double factor = engngModel->giveDomain(1)->giveLoadTimeFunction(prescribedDisplacementLTF)->__at( atTime->giveTargetTime() );
    if ( ( rlm == rlm_total ) && ( !atTime->isTheFirstStep() ) ) {
        //factor -= engngModel->giveDomain(1)->giveLoadTimeFunction(prescribedDisplacementLTF)->
        // at(atTime->givePreviousStep()->giveTime()) ;
        factor -= engngModel->giveDomain(1)->giveLoadTimeFunction(prescribedDisplacementLTF)->
                  __at( atTime->giveTargetTime() - atTime->giveTimeIncrement() );
    }

    if ( nite == 0 ) {
#if 0
 #ifdef __PETSC_MODULE
        if ( solverType == ST_Petsc ) {
  #ifdef __PARALLEL_MODE
            //PetscNatural2LocalOrdering* n2lpm = engngModel->givePetscContext(1)->giveN2Lmap();
            //IntArray* map = n2lpm->giveN2Lmap();
            for ( i = 1; i <= prescribedEqs.giveSize(); i++ ) {
                eq = prescribedEqs.at(i);
                R.at(eq) = prescribedDofsValues.at(i) * factor; // local eq
            }

            return;

  #else // PETSC_SERIAL
            for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
                eq = prescribedEqs.at(i);
                R.at(eq) = prescribedDofsValues.at(i) * factor;
            }

            return;

  #endif
        }

 #endif
#else
 #ifdef __PETSC_MODULE
        if ( solverType == ST_Petsc ) {
            if ( k->giveType() != SMT_PetscMtrx ) {
                _error("applyConstraintsToStiffness: PetscSparseMtrx Expected");
            }

            PetscSparseMtrx *lhs = ( PetscSparseMtrx * ) k;

            Vec diag;
            PetscScalar *ptr;
            int eq;
            engngModel->givePetscContext(1, ut)->createVecGlobal(& diag);
            MatGetDiagonal(* lhs->giveMtrx(), diag);
            VecGetArray(diag, & ptr);

            for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
                eq = prescribedEqs.at(i) - 1;
                R.at(eq + 1) = ptr [ eq ] * prescribedDofsValues.at(i) * factor;
            }

            return;

            VecRestoreArray(diag, & ptr);
            VecDestroy(&diag);
        }

 #endif
#endif
        for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
            eq = prescribedEqs.at(i);
            R.at(eq) = k->at(eq, eq) * prescribedDofsValues.at(i) * factor;
        }
    } else {
        for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = 0.0;
        }
    }
}


void
NRSolver :: printState(FILE *outputStream)
{
#ifdef VERBOSE
    // print quasi reactions if direct displacement control used
    fprintf(outputStream, "\nQuasi reaction table:\n\n");
    fprintf(outputStream, "  node  dof            force\n");
    fprintf(outputStream, "============================\n");
    if ( lastReactions.giveSize() == 0 ) {
        return;
    }

    double reaction;
    for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
        reaction = lastReactions.at(i);
        fprintf(outputStream, "%6d  %3d   %+11.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i), reaction);
    }
    fprintf(outputStream, "============================\n\n");
#endif
}


#if 1
bool
NRSolver :: checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs,  FloatArray &ddX, FloatArray &X,
                             double RRT, double internalForcesEBENorm,
                             int nite, bool &errorOutOfRange, TimeStep *tNow)
{
    int _dg, _idofman, _ielem, _idof, _eq, _ndof, _ng = nccdg;
    int ndofman = domain->giveNumberOfDofManagers();
    int nelem = domain->giveNumberOfElements();
    double forceErr, dispErr, _val;
    DofManager *_idofmanptr;
    Element *_ielemptr;
    Dof *_idofptr;
    FloatArray dg_forceErr(nccdg), dg_dispErr(nccdg), dg_totalLoadLevel(nccdg), dg_totalDisp(nccdg);
    bool answer;
    EModelDefaultEquationNumbering dn;
 #ifdef __PARALLEL_MODE
  #ifdef __PETSC_MODULE
    PetscContext *parallel_context = engngModel->givePetscContext(1, ut); ///@todo Shouldn't be hardwired domain no 1.
    PetscNatural2LocalOrdering *n2l = parallel_context->giveN2Lmap();
  #endif
 #endif

    /*
     * The force errors are (if possible) evaluated as relative errors.
     * If the norm of applied load vector is zero (one may load by temperature, etc)
     * then the norm of reaction forces is used in relative norm evaluation.
     *
     * Note: This is done only when all dofs are included (nccdg = 0). Not implemented if
     * multiple convergence criteria are used.
     *
     */

    answer = true;
    errorOutOfRange = false;

    if ( _ng > 0 ) {
        // zero error norms per group
        dg_forceErr.zero();
        dg_dispErr.zero();
        dg_totalLoadLevel.zero();
        dg_totalDisp.zero();
        // loop over dof managers
        for ( _idofman = 1; _idofman <= ndofman; _idofman++ ) {
            _idofmanptr = domain->giveDofManager(_idofman);
 #if ( defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE ) )
            if ( !parallel_context->isLocal(_idofmanptr) ) {
                continue;
            }

 #endif

            _ndof = _idofmanptr->giveNumberOfDofs();
            // loop over individual dofs
            for ( _idof = 1; _idof <= _ndof; _idof++ ) {
                _idofptr = _idofmanptr->giveDof(_idof);
                // loop over dof groups
                for ( _dg = 1; _dg <= _ng; _dg++ ) {
                    // test if dof ID is in active set
                    if ( ccDofGroups.at(_dg - 1).find( _idofptr->giveDofID() ) != ccDofGroups.at(_dg - 1).end() ) {
                        _eq = _idofptr->giveEquationNumber(dn);

                        if ( _eq ) {
                            _val = rhs.at(_eq);
                            dg_forceErr.at(_dg) += _val * _val;
                            _val = ddX.at(_eq);
                            dg_dispErr.at(_dg)  += _val * _val;
                            // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                            _val = RT.at(_eq);
                            dg_totalLoadLevel.at(_dg) += _val * _val;
                            _val = X.at(_eq);
                            dg_totalDisp.at(_dg) += _val * _val;
                        }
                    }
                } // end loop over dof groups

            } // end loop over DOFs

        } // end loop over dof managers

        // loop over elements and their DOFs
        for ( _ielem = 1; _ielem <= nelem; _ielem++ ) {
            _ielemptr = domain->giveElement(_ielem);
 #ifdef __PARALLEL_MODE
            if ( _ielemptr->giveParallelMode() != Element_local ) {
                continue;
            }

 #endif
            // Aren't internal dofmanagers listed along with all the other dofmanagers, and included in the loop above? / Mikael
            // loop over element internal Dofs
            for (_idofman=1; _idofman<=_ielemptr->giveNumberOfInternalDofManagers(); _idofman++) {
                _ndof = _ielemptr->giveInternalDofManager(_idofman)->giveNumberOfDofs();
                // loop over individual dofs
                for ( _idof = 1; _idof <= _ndof; _idof++ ) {
                    _idofptr = _ielemptr->giveInternalDofManager(_idofman)->giveDof(_idof);
                    // loop over dof groups
                    for ( _dg = 1; _dg <= _ng; _dg++ ) {
                        // test if dof ID is in active set
                        if ( ccDofGroups.at(_dg - 1).find( _idofptr->giveDofID() ) != ccDofGroups.at(_dg - 1).end() ) {
                            _eq = _idofptr->giveEquationNumber(dn);

                            if ( _eq ) {
 #if ( defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE ) )
                                if ( !n2l->giveNewEq(_eq) ) {
                                    continue;
                                }
 #endif

                                _val = rhs.at(_eq);
                                dg_forceErr.at(_dg) += _val * _val;
                                _val = ddX.at(_eq);
                                dg_dispErr.at(_dg)  += _val * _val;
                                // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                                _val = RT.at(_eq);
                                dg_totalLoadLevel.at(_dg) += _val * _val;
                                _val = X.at(_eq);
                                dg_totalDisp.at(_dg) += _val * _val;
                            }
                        }
                    } // end loop over dof groups
                } // end loop over DOFs
            } // end loop over element internal dofmans
        } // end loop over elements

 #ifdef __PARALLEL_MODE
        // exchange individual partition contributions (simultaneously for all groups)
#ifdef __PETSC_MODULE
        FloatArray collectiveErr(_ng);
        parallel_context->accumulate(dg_forceErr,       collectiveErr); dg_forceErr       = collectiveErr;
        parallel_context->accumulate(dg_dispErr,        collectiveErr); dg_dispErr        = collectiveErr;
        parallel_context->accumulate(dg_totalLoadLevel, collectiveErr); dg_totalLoadLevel = collectiveErr;
        parallel_context->accumulate(dg_totalDisp,      collectiveErr); dg_totalDisp      = collectiveErr;
#else
        if (this->engngModel->isParallel()) {
            FloatArray collectiveErr(_ng);
            MPI_Allreduce(dg_forceErr.givePointer(), collectiveErr.givePointer(), _ng, MPI_DOUBLE, MPI_SUM, comm);
            dg_forceErr = collectiveErr;
            MPI_Allreduce(dg_dispErr.givePointer(), collectiveErr.givePointer(), _ng, MPI_DOUBLE, MPI_SUM, comm);
            dg_dispErr = collectiveErr;
            MPI_Allreduce(dg_totalLoadLevel.givePointer(), collectiveErr.givePointer(), _ng, MPI_DOUBLE, MPI_SUM, comm);
            dg_totalLoadLevel = collectiveErr;
            MPI_Allreduce(dg_totalDisp.givePointer(), collectiveErr.givePointer(), _ng, MPI_DOUBLE, MPI_SUM, comm);
            dg_totalDisp = collectiveErr;
            return globalNorm;
        }
#endif
 #endif
        OOFEM_LOG_INFO("NRSolver:     %-10d ", nite);
        // loop over dof groups
        for ( _dg = 1; _dg <= _ng; _dg++ ) {
            //  compute a relative error norm
            if ( ( dg_totalLoadLevel.at(_dg) ) < nrsolver_ERROR_NORM_SMALL_NUM ) {
                dg_forceErr.at(_dg) = sqrt( dg_forceErr.at(_dg) );
            } else {
                dg_forceErr.at(_dg) = sqrt( dg_forceErr.at(_dg) / dg_totalLoadLevel.at(_dg) );
            }

            //
            // compute displacement error
            //
            if ( dg_totalDisp.at(_dg) <  nrsolver_ERROR_NORM_SMALL_NUM ) {
                dg_dispErr.at(_dg) = sqrt( dg_dispErr.at(_dg) );
            } else {
                dg_dispErr.at(_dg) = sqrt( dg_dispErr.at(_dg) / dg_totalDisp.at(_dg) );
            }

            if ( ( fabs( dg_forceErr.at(_dg) ) > rtolf.at(_dg) * NRSOLVER_MAX_REL_ERROR_BOUND ) ||
                ( fabs( dg_dispErr.at(_dg) )  > rtold.at(_dg) * NRSOLVER_MAX_REL_ERROR_BOUND ) ) {
                errorOutOfRange = true;
            }

            if ( ( fabs( dg_forceErr.at(_dg) ) > rtolf.at(_dg) ) || ( fabs( dg_dispErr.at(_dg) ) > rtold.at(_dg) ) ) {
                answer = false;
            }

            OOFEM_LOG_INFO( "%-15e %-15e ", dg_forceErr.at(_dg), dg_dispErr.at(_dg) );
        }
        OOFEM_LOG_INFO("\n");
    } else { // nccdg == 0 -> all dofs included
        double dXX, dXdX;
 #ifdef __PARALLEL_MODE
        forceErr = parallel_context->norm(rhs); forceErr *= forceErr;
        dXX = parallel_context->localNorm(X); dXX *= dXX; // Note: Solutions are always total global values (natural distribution makes little sense for the solution)
        dXdX = parallel_context->norm(ddX); dXdX *= dXdX;
 #else
        forceErr = rhs.computeSquaredNorm();
        dXX = X.computeSquaredNorm();
        dXdX = ddX.computeSquaredNorm();
 #endif
        // we compute a relative error norm
        if ( RRT > nrsolver_ERROR_NORM_SMALL_NUM ) {
            forceErr = sqrt( forceErr / ( RRT ) );
        } else if ( internalForcesEBENorm > nrsolver_ERROR_NORM_SMALL_NUM ) {
            forceErr = sqrt(forceErr / internalForcesEBENorm);
        } else {
            forceErr = sqrt(forceErr); // absolute norm as last resort
        }

 #if 0
        // load vector norm close to zero
        // try to take norm of reactions instead
    } else {
        FloatArray reactions;
        int i, di = 1;   // hard wired domain index =  1
        // ask emodel to evaluate reactions
        ( ( StructuralEngngModel * ) engngModel )->computeReactions(reactions, tNow, di);
        // compute corresponding norm
        double RN;

  #ifdef __PARALLEL_MODE
   #ifdef __PETSC_MODULE
        double myRN = 0.0;
        for ( i = 1; i <= reactions.giveSize(); i++ ) {
            if ( n2l_prescribed->giveNewEq(i) ) {
                myRN += reactions.at(i) * reactions.at(i);
            }
        }

        // account for quasi bc reactions
        for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
            myRN += F.at( prescribedEqs.at(i) ) * F.at( prescribedEqs.at(i) );
        }

        MPI_Allreduce(& myRN, & RN, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif
  #else
        RN = reactions.computeSquaredNorm();
        // account for quasi bc reactions
        for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
            RN += F.at( prescribedEqs.at(i) ) * F.at( prescribedEqs.at(i) );
        }

  #endif // __PARALLEL_MODE
        if ( RN > nrsolver_ERROR_NORM_SMALL_NUM ) {
            forceErr = sqrt( forceErr / ( RN ) );
        } else {
            forceErr = sqrt(forceErr); // absolute norm as last resort
        }
    }

 #endif // if 0

        // compute displacement error
        //
        // err is relative displacement change
        if ( dXX < nrsolver_ERROR_NORM_SMALL_NUM ) {
            dispErr = sqrt(dXdX);
        } else {
            dispErr = sqrt(dXdX / dXX);
        }

        if ( ( fabs(forceErr) > rtolf.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) ||
            ( fabs(dispErr)  > rtold.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) ) {
            errorOutOfRange = true;
        }

        if ( ( fabs(forceErr) > rtolf.at(1) ) || ( fabs(dispErr) > rtold.at(1) ) ) {
            answer = false;
        }
        if ( engngModel->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO("NRSolver:     %-15d %-15e %-15e\n", nite, forceErr, dispErr);
        } else {
            OOFEM_LOG_INFO("  NRSolver:     %-15d %-15e %-15e\n", nite, forceErr, dispErr);
        }
    } // end default case (all dofs conributing)

    return answer;
}

#endif
} // end namespace oofem
