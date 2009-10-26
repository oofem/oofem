/* $Header: /home/cvs/bp/oofem/oofemlib/src/nrsolver.C,v 1.10.4.1 2004/04/05 15:19:43 bp Exp $ */
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
// file nrsolver2.C
//

#include "nrsolver.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <math.h>
#endif

#include "verbose.h"
#include "ldltfact.h"
#include "imlsolver.h"
#include "timestep.h"
#include "flotmtrx.h"
//#include "nlinearstatic.h"
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


#define nrsolver_ERROR_NORM_SMALL_NUM 1.e-6
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e10
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10


NRSolver :: NRSolver(int i, Domain *d, EngngModel *m, EquationID ut) :
    SparseNonLinearSystemNM(i, d, m, ut), prescribedDofs(), prescribedDofsValues()  {
    //
    // constructor
    //
    nsmax  = 60;     // default maximum number of sweeps allowed
    //Psi    = 0.1;       // displacement control on
    deltaL = 1.0;
    solved = 0;
    NR_Mode = NR_OldMode = nrsolverModifiedNRM;
    NR_ModeTick = -1; // do not swith to calm_NR_OldMode
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

NRSolver ::  ~NRSolver() {
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
        ISDestroy(prescribedEgsIS);
    }

 #endif
#endif
}


NM_Status
NRSolver :: solve(SparseMtrx *k, FloatArray *R, FloatArray *R0,
                  FloatArray *Rr, FloatArray *r, FloatArray *DeltaR, FloatArray *F,
                  double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tNow)
//
// this function solve the problem of the unbalanced equilibrium
// using NR scheme
//
//
{
    FloatArray rhs, deltaR, RT;
    FloatArray rInitial;
    //FloatArray F;
    double RRT;
    int neq = r->giveSize();
    int irest = 0;
    NM_Status status;
    bool converged, errorOutOfRangeFlag;
#ifdef __PARALLEL_MODE
 #ifdef __PETSC_MODULE
    // HUHU hard wired domain no 1
    int i;
    PetscNatural2LocalOrdering *n2l = engngModel->givePetscContext(1, ut)->giveN2Lmap();
 #endif
#endif

    OOFEM_LOG_INFO("Time       Iteration       ForceError      DisplError\n__________________________________________________________\n");

    rInitial = * r;
    l = 1.0;

    status = NM_None;
    this->giveLinearSolver();

    if ( this->prescribedDofsFlag ) {
        if ( !prescribedEqsInitFlag ) {
            this->initPrescribedEqs();
        }

        //this->computeBCLoadVector (*Rr, k, tNow);
        applyConstraintsToStiffness(k);
    }

    // compute total load R = R+R0
    RT = * R;
    if ( R0 ) {
        RT.add(R0);
    }

restart:
    DeltaR->zero();

    //linSolver -> setSparseMtrxAsComponent (LinearEquationLhs,k);

    //deltaL = tNow->giveTimeIncrement();

    deltaR.resize(neq);
    // if (tNow ->giveNumber() == 1) {
    rhs =  * R;
    // if (R0) rhs.add(*R0);

    //engngModel->updateComponent (tNow, NonLinearRhs_Total);

#ifdef __PARALLEL_MODE
 #ifdef __PETSC_MODULE
    double myRRT = 0.0;
    for ( i = 1; i <= neq; i++ ) {
        if ( n2l->giveNewEq(i) ) {
            myRRT += RT.at(i) * RT.at(i);
        }
    }

    MPI_Allreduce(& myRRT, & RRT, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 #endif
#else
    RRT = dotProduct(RT.givePointer(), RT.givePointer(), neq);
#endif
    //if (R0) RR0 = dotProduct(R0->givePointer(),R0->givePointer(),neq);
    //else RR0 = 0.0;

    nite = 0;

    do {
        nite++;

        if ( nite > 1 ) {
            if ( ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
                engngModel->updateComponent(tNow, NonLinearLhs, domain);
                //linSolver -> setSparseMtrxAsComponent (LinearEquationLhs,k);
                applyConstraintsToStiffness(k);
            }
        }

        /*
         * linSolver -> setFloatArrayAsComponent (LinearEquationRhs,&rhs);
         * linSolver -> setFloatArrayAsComponent (LinearEquationSolution,&deltaR);
         * linSolver -> solveYourselfAt (tNow);
         * linSolver -> updateYourselfExceptLhs ();
         */
        /*
         * if ((nite == 1) && (numberOfPrescribedDofs)) {
         * // modify t
         * rhs.add (*Rr);
         * }
         */
        if ( this->prescribedDofsFlag ) {
            this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tNow);
        }

        if ( ( nite == 1 ) && ( Rr->giveSize() ) ) {
            rhs.add(* Rr);
        }

        if ( ( nite == 1 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equlibrate current state
            rhs.times(0.0);
            R->times(0.0);
            deltaR = rhs;
        } else {
            linSolver->solve(k, & rhs, & deltaR);
        }

        //
        // update solution
        //
        if ( this->lsFlag && ( nite != 1 ) ) {
            // linesearch
            LineSearchNM :: LS_status status;
            double eta;
            this->giveLineSearchSolver()->solve(r, & deltaR, F, R, R0, prescribedEqs, 1.0, eta, status, tNow);
            DeltaR->add(deltaR);
        } else {
            r->add(deltaR);
            DeltaR->add(deltaR);
            tNow->incrementStateCounter();     // update solution state counter
            //
            // convergency check
            //
            //((NonLinearStatic *)engngModel) -> giveInternalForces(F, *DeltaR, tNow);
            engngModel->updateComponent(tNow, InternalRhs, domain);
            //F->negated();
        }

        //
        // convergency check
        //
        rhs = RT;
        rhs.substract(F);

        //
        // compute forceError
        //
        // err is relative error of unbalanced forces
        // account for quasi BC
        for ( int ii = 1; ii <= numberOfPrescribedDofs; ii++ ) {
            rhs.at( prescribedEqs.at(ii) ) = 0.0;
        }


        converged = this->checkConvergence(RT, rhs, deltaR, * r, RRT, nite, errorOutOfRangeFlag, tNow);


        if ( ( nite >= nsmax ) || errorOutOfRangeFlag ) {
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
                r->times(0.);
                r->add(rInitial);
                // reset all changes fro previous equilibrium state
                engngModel->initStepIncrements();
                DeltaR->zero();
                // restore initial stiffness
                engngModel->updateComponent(tNow, NonLinearLhs, domain);
                // recalculate new Load Vector R
                // engngModel->updateComponent (tNow, NonLinearRhs_Incremental);
                //delete F; F = NULL;
#ifdef VERBOSE
                OOFEM_LOG_INFO("NRSolver iteration Reset ...\n");
#endif
                NR_OldMode  = NR_Mode;
                NR_Mode     = nrsolverFullNRM;
                NR_ModeTick = NRSOLVER_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = NM_NoSuccess;
                _warning2("NRSolver - convergence not reached after %d iterations", nsmax);
                // exit(1);
                break;
            }
        }
    } while ( ( !converged ) || ( nite < minIterations ) );

    //
    // end of iteration
    //
    // ls ->letSolutionBe(deltar);
    // Lambda += DeltaLambda ;      // *
    //
    // update dofs,nodes,Elemms and print result
    //
#ifdef VERBOSE
    // printf ("\nCALM - step iteration finished") ;
#endif

    status |= NM_Success;
    solved = 1;

    // Modify Load vector to include "quasi rection"
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
    // print quasi reactions if direct displacement controll used
    OOFEM_LOG_INFO("\n  Quasi reaction table:\n\n");
    OOFEM_LOG_INFO("  node  dof   displacement         force\n");
    OOFEM_LOG_INFO("========================================\n");
#endif
    double reaction;
    for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
        reaction = R->at( prescribedEqs.at(i) );
        if ( R0 ) {
            reaction += R0->at( prescribedEqs.at(i) );
        }

        lastReactions.at(i) = reaction;
#ifdef VERBOSE
        OOFEM_LOG_INFO("%6d  %3d   %+11.5e  %+11.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i),
                       r->at( prescribedEqs.at(i) ), reaction);
#endif
    }

#ifdef VERBOSE
    OOFEM_LOG_INFO("========================================\n");
#endif

    return status;
}

IRResultType
NRSolver :: initializeFrom(InputRecord *ir)
//
//
//
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

    /** initialize optional dof groups for convergence criteria evaluation */
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
                ccDofGroups.at(_i).insert( _val.at(_j) );
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
        // if common rtolv provided, set to this tolerace both rtolf and rtold
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
NRSolver :: saveContext(DataStream *stream, ContextMode mode, void *obj) {
    return CIO_OK;
}

contextIOResultType
NRSolver :: restoreContext(DataStream *stream, ContextMode mode, void *obj) {
    return CIO_OK;
}


SparseLinearSystemNM *
NRSolver :: giveLinearSolver() {
    if ( linSolver ) {
        if ( linSolver->giveLinSystSolverType() == solverType ) {
            return linSolver;
        } else {
            delete linSolver;
        }
    }

    linSolver = :: CreateUsrDefSparseLinSolver(solverType, 1, domain, engngModel);
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
    PetscNatural2GlobalOrdering *n2lpm = engngModel->givePetscContext(1, ut)->giveN2Gmap();
 #endif
    int jglobnum, count = 0, ndofman = domain->giveNumberOfDofManagers();
    int i, j, inode, idof;
    IntArray localPrescribedEqs(numberOfPrescribedDofs);

    for ( j = 1; j <= ndofman; j++ ) {
        jglobnum = domain->giveNode(j)->giveGlobalNumber();
        for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
            inode = prescribedDofs.at(2 * i - 1);
            idof  = prescribedDofs.at(2 * i);
            if ( inode == jglobnum ) {
 #if defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE )
                // HUHU hard wired domain no 1
                if ( n2lpm->isLocal( domain->giveNode(j) ) ) {
                    localPrescribedEqs.at(++count) = domain->giveNode(j)->giveDof(idof)->giveEquationNumber(dn);
                }

 #else
                localPrescribedEqs.at(++count) = domain->giveNode(j)->giveDof(idof)->giveEquationNumber(dn);
 #endif

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

/*
 * void
 * NRSolver :: computeBCLoadVector (FloatArray& answer, SparseMtrx* k, TimeStep* atTime)
 * {
 * FloatArray rr(k->giveNumberOfRows());
 * int i;
 *
 * double factor = engngModel->giveDomain(1)->giveLoadTimeFunction(prescribedDisplacementLTF)->at(atTime->giveTime());
 *
 * for (i=1; i<=numberOfPrescribedDofs; i++) {
 * rr.at(prescribedEqs.at(i)) = prescribedDofsValues.at(i)*factor;
 * }
 *
 * k->times (rr, answer);
 * answer.times(-1.0);
 * }
 */

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
        int eq, res;

        engngModel->givePetscContext(1, ut)->createVecGlobal(& diag);
        MatGetDiagonal(* lhs->giveMtrx(), diag);
        VecGetArray(diag, & ptr);
        for ( i = 1; i <= numberOfPrescribedDofs; i++ ) {
            eq = prescribedEqs.at(i) - 1;
            res = MatSetValue(* ( lhs->giveMtrx() ), eq, eq, ptr [ eq ] * 1.e6, INSERT_VALUES);
        }

        MatAssemblyBegin(* lhs->giveMtrx(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(* lhs->giveMtrx(), MAT_FINAL_ASSEMBLY);
        VecRestoreArray(diag, & ptr);
        VecDestroy(diag);
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
    double factor = engngModel->giveDomain(1)->giveLoadTimeFunction(prescribedDisplacementLTF)->__at( atTime->giveTime() );
    if ( ( rlm == rlm_total ) && ( !atTime->isTheFirstStep() ) ) {
        //factor -= engngModel->giveDomain(1)->giveLoadTimeFunction(prescribedDisplacementLTF)->
        // at(atTime->givePreviousStep()->giveTime()) ;
        factor -= engngModel->giveDomain(1)->giveLoadTimeFunction(prescribedDisplacementLTF)->
        __at( atTime->giveTime() - atTime->giveTimeIncrement() );
    }

    if ( nite == 1 ) {
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
            VecDestroy(diag);
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
    // print quasi reactions if direct displacement controll used
    fprintf(outputStream, "\nQuasi reaction table:\n\n");
    fprintf(outputStream, "  node  dof            force\n");
    fprintf(outputStream, "============================\n");
#endif
    if ( lastReactions.giveSize() == 0 ) {
        return;
    }

    double reaction;
    for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
        reaction = lastReactions.at(i);
#ifdef VERBOSE
        fprintf(outputStream, "%6d  %3d   %+11.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i), reaction);
#endif
    }

#ifdef VERBOSE
    fprintf(outputStream, "============================\n\n");
#endif
}


#if 1
bool
NRSolver :: checkConvergence(FloatArray &RT, FloatArray &rhs, FloatArray &deltaR, FloatArray &r,
                             double RRT, int nite, bool &errorOutOfRange, TimeStep *tNow)
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
    int i;
    // HUHU hard wired domain no 1
    PetscNatural2LocalOrdering *n2l = engngModel->givePetscContext(1, ut)->giveN2Lmap();
  #endif
 #endif

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
 #ifdef __PARALLEL_MODE
            if ( !_idofmanptr->isLocal() ) {
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
 #if ( defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE ) )
                            if ( !n2l->giveNewEq(_eq) ) {
                                continue;
                            }

 #endif

                            _val = rhs.at(_eq);
                            dg_forceErr.at(_dg) += _val * _val;
                            _val = deltaR.at(_eq);
                            dg_dispErr.at(_dg)  += _val * _val;
                            // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                            _val = RT.at(_eq);
                            dg_totalLoadLevel.at(_dg) += _val * _val;
                            _val = r.at(_eq);
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

            _ndof = _ielemptr->giveNumberOfDofs();
            // loop over individual dofs
            for ( _idof = 1; _idof <= _ndof; _idof++ ) {
                _idofptr = _ielemptr->giveDof(_idof);
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
                            _val = deltaR.at(_eq);
                            dg_dispErr.at(_dg)  += _val * _val;
                            // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                            _val = RT.at(_eq);
                            dg_totalLoadLevel.at(_dg) += _val * _val;
                            _val = r.at(_eq);
                            dg_totalDisp.at(_dg) += _val * _val;
                        }
                    }
                } // end loop over dof groups

            } // end loop over DOFs

        } // end loop over dof managers

 #ifdef __PARALLEL_MODE
        // exchange individual partition contributions (simultaneously for all groups)
        FloatArray collectiveErr(_ng);
        MPI_Allreduce(dg_forceErr.givePointer(), collectiveErr.givePointer(), _ng, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dg_forceErr = collectiveErr;
        MPI_Allreduce(dg_dispErr.givePointer(), collectiveErr.givePointer(), _ng, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dg_dispErr = collectiveErr;
        MPI_Allreduce(dg_totalLoadLevel.givePointer(), collectiveErr.givePointer(), _ng, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dg_totalLoadLevel = collectiveErr;
        MPI_Allreduce(dg_totalDisp.givePointer(), collectiveErr.givePointer(), _ng, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dg_totalDisp = collectiveErr;
 #endif

        OOFEM_LOG_INFO("%-5d %-15e ", ( int ) tNow->giveTime(), nite);
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
        double drr, drdr;
 #ifdef __PARALLEL_MODE
  #ifdef __PETSC_MODULE
        int neq = r.giveSize();
        double myerr [ 3 ] = {
            0., 0., 0.
        }, colerr [ 3 ];
        for ( i = 1; i <= neq; i++ ) {
            if ( n2l->giveNewEq(i) ) {
                myerr [ 0 ] += rhs.at(i) * rhs.at(i);
                myerr [ 1 ] += r.at(i) * r.at(i);
                myerr [ 2 ] += deltaR.at(i) * deltaR.at(i);
            }
        }

        MPI_Allreduce(myerr, colerr, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        forceErr = colerr [ 0 ];
        drr = colerr [ 1 ];
        drdr = colerr [ 2 ];
  #endif
 #else
        forceErr = dotProduct( rhs.givePointer(), rhs.givePointer(), rhs.giveSize() );
        drr = dotProduct( r.givePointer(), r.givePointer(), r.giveSize() );
        drdr = dotProduct( deltaR.givePointer(), deltaR.givePointer(), deltaR.giveSize() );
 #endif

        // we compute a relative error norm
        if ( ( RRT ) > nrsolver_ERROR_NORM_SMALL_NUM ) {
            forceErr = sqrt( forceErr / ( RRT ) );
        } else {
            forceErr = sqrt(forceErr); // absolute norm
        }

        //
        // compute displacement error
        //
        // err is relative displacement change
        if ( drr < nrsolver_ERROR_NORM_SMALL_NUM ) {
            dispErr = drdr;
            dispErr = sqrt(dispErr);
        } else {
            dispErr = drdr / drr;
            dispErr = sqrt(dispErr);
        }

        if ( ( fabs(forceErr) > rtolf.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) ||
            ( fabs(dispErr)  > rtold.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) ) {
            errorOutOfRange = true;
        }

        if ( ( fabs(forceErr) > rtolf.at(1) ) || ( fabs(dispErr) > rtold.at(1) ) ) {
            answer = false;
        }

        OOFEM_LOG_INFO("%-10d %-15d %-15e %-15e\n", ( int ) tNow->giveTime(), nite, forceErr, dispErr);
    } // end default case (all dofs conributing)

    return answer;
}
#endif
