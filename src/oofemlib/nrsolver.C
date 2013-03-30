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
#include "generalbc.h"
#include "dof.h"
#include "loadtime.h"
#include "linesearch.h"
#include "usrdefsub.h"
#ifdef __PETSC_MODULE
 #include "petscsolver.h"
 #include "petscsparsemtrx.h"
 #include "petscordering.h"
#endif

#include <cstdio>

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
    nsmax = 60;     // default maximum number of sweeps allowed
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

    // Choosing a big "enough" number. (Alternative: Force input of maxinter)
    nsmax = 1e8;
    IR_GIVE_OPTIONAL_FIELD(ir, nsmax, _IFT_NRSolver_maxiter);
    if ( nsmax < 0 ) {
        _error("initializeFrom: nsmax < 0");
    }

    minIterations = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, minIterations, _IFT_NRSolver_miniterations);

    minStepLength = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, _IFT_NRSolver_minsteplength);

    // read if MANRM method is used
    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, _IFT_NRSolver_manrmsteps);
    if ( MANRMSteps > 0 ) {
        NR_Mode = NR_OldMode = nrsolverAccelNRM;
    } else {
        NR_Mode = nrsolverModifiedNRM;
    }

    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_NRSolver_lstype);
    solverType = ( LinSystSolverType ) _val;
    this->giveLinearSolver()->initializeFrom(ir);

    // read relative error tolerances of the solver
    // if rtolv provided set to this tolerance both rtolf and rtold
    rtolf.resize(1);
    rtolf.at(1) = 1.e-3; // Default value.
    IR_GIVE_OPTIONAL_FIELD(ir, rtolf.at(1), _IFT_NRSolver_rtolv);
    rtold = rtolf;
    // read optional force and displacement tolerances
    IR_GIVE_OPTIONAL_FIELD(ir, rtolf.at(1), _IFT_NRSolver_rtolf);
    IR_GIVE_OPTIONAL_FIELD(ir, rtold.at(1), _IFT_NRSolver_rtold);

    prescribedDofs.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDofs, _IFT_NRSolver_ddm);
    prescribedDofsValues.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDofsValues, _IFT_NRSolver_ddv);
    prescribedDisplacementLTF = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDisplacementLTF, _IFT_NRSolver_ddltf);

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
    IR_GIVE_OPTIONAL_FIELD(ir, lsFlag, _IFT_NRSolver_linesearch);

    if ( this->lsFlag ) {
        this->giveLineSearchSolver()->initializeFrom(ir);
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
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tNow)
//
// this function solve the problem of the unbalanced equilibrium
// using NR scheme
//
//
{
    // residual, iteration increment of solution, total external force
    FloatArray rhs, ddX, RT;
    double RRT;
    int neq = X->giveSize();
    NM_Status status;
    bool converged, errorOutOfRangeFlag;
#ifdef __PARALLEL_MODE
 #ifdef __PETSC_MODULE
    PetscContext *parallel_context  = engngModel->givePetscContext(this->domain->giveNumber(), ut);
 #endif
#endif

    if ( engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("NRSolver:     Iteration       ForceError      DisplError                    \n");
        OOFEM_LOG_INFO("----------------------------------------------------------------------------\n");
    }

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

    ddX.resize(neq);

    // Fetch the matrix before evaluating internal forces.
    // This is intentional, since its a simple way to drastically increase convergence for nonlinear problems.
    // (This old tangent is just used)

    engngModel->updateComponent(tNow, NonLinearLhs, domain);
    if ( this->prescribedDofsFlag ) {
        if ( !prescribedEqsInitFlag ) {
            this->initPrescribedEqs();
        }
        applyConstraintsToStiffness(k);
    }

    nite = 0;
    do {
        // Compute the residual
        engngModel->updateComponent(tNow, InternalRhs, domain);
        rhs.beDifferenceOf(RT, *F);

        if ( this->prescribedDofsFlag ) {
            this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tNow);
        }

        // convergence check
        converged = this->checkConvergence(RT, * F, rhs, ddX, * X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag, tNow);

        if ( errorOutOfRangeFlag ) {
            status = NM_NoSuccess;
            OOFEM_WARNING2("NRSolver:  Divergence reached after %d iterations", nite);
            break;
        } else if ( converged && ( nite >= minIterations ) ) {
            break;
        } else if ( nite >= nsmax ) {
            OOFEM_LOG_DEBUG("Maximum number of iterations reached\n");
            break;
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

        nite++; // iteration increment
    } while ( true ); // end of iteration

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
 #if defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE )
    PetscContext *parallel_context = engngModel->givePetscContext(this->domain->giveNumber(), ut);
 #endif
    int jglobnum, count = 0, ndofman = domain->giveNumberOfDofManagers();
    int inode, idof;
    IntArray localPrescribedEqs(numberOfPrescribedDofs);

    for ( int j = 1; j <= ndofman; j++ ) {
 #if defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE )
        if ( !parallel_context->isLocal(domain->giveNode(j)) ) {
            continue;
        }
 #endif
        jglobnum = domain->giveNode(j)->giveGlobalNumber();
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            inode = prescribedDofs.at(2 * i - 1);
            idof  = prescribedDofs.at(2 * i);
            if ( inode == jglobnum ) {
                localPrescribedEqs.at(++count) = domain->giveNode(j)->giveDof(idof)->giveEquationNumber(dn);
                continue;
            }
        }
    }

    prescribedEqs.resize(count);
    for ( int i = 1; i <= count; i++ ) {
        prescribedEqs.at(i) = localPrescribedEqs.at(i);
    }

    numberOfPrescribedDofs = count;

    this->prescribedEqsInitFlag = true;
}


void
NRSolver :: applyConstraintsToStiffness(SparseMtrx *k)
{
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
            for ( int i = 1; i <= s; i++ ) {
                eqs.at(i) = n2lpm->giveNewEq( prescribedEqs.at(i) );
            }

            ISCreateGeneral(PETSC_COMM_WORLD, s, eqs.givePointer(), & prescribedEgsIS);
            //ISView(prescribedEgsIS,PETSC_VIEWER_STDOUT_WORLD);
  #else
            eqs.resize(numberOfPrescribedDofs);
            for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
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

        PetscContext *parallel_context = engngModel->givePetscContext(this->domain->giveNumber(), ut);
        parallel_context->createVecGlobal(& diag);
        MatGetDiagonal(* lhs->giveMtrx(), diag);
        VecGetArray(diag, & ptr);
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
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
    for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
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
            engngModel->givePetscContext(this->domain->giveNumber(), ut)->createVecGlobal(& diag);
            MatGetDiagonal(* lhs->giveMtrx(), diag);
            VecGetArray(diag, & ptr);

            for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
                int eq = prescribedEqs.at(i) - 1;
                R.at(eq + 1) = ptr [ eq ] * prescribedDofsValues.at(i) * factor;
            }

            VecRestoreArray(diag, & ptr);
            VecDestroy(&diag);
            return;
        }
 #endif
#endif
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            int eq = prescribedEqs.at(i);
            R.at(eq) = k->at(eq, eq) * prescribedDofsValues.at(i) * factor;
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
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


bool
NRSolver :: checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs,  FloatArray &ddX, FloatArray &X,
                             double RRT, const FloatArray &internalForcesEBENorm,
                             int nite, bool &errorOutOfRange, TimeStep *tNow)
{
    double forceErr, dispErr;
    FloatArray dg_forceErr, dg_dispErr, dg_totalLoadLevel, dg_totalDisp;
    bool answer;
    EModelDefaultEquationNumbering dn;
 #ifdef __PARALLEL_MODE
  #ifdef __PETSC_MODULE
    PetscContext *parallel_context = engngModel->givePetscContext(this->domain->giveNumber(), ut);
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

    if ( internalForcesEBENorm.giveSize() > 1 ) { // Special treatment when just one norm is given; No grouping
        int nccdg = this->domain->giveMaxDofID();
        // Keeps tracks of which dof IDs are actually in use;
        IntArray idsInUse(nccdg);
        idsInUse.zero();
        // zero error norms per group
        dg_forceErr.resize(nccdg); dg_forceErr.zero();
        dg_dispErr.resize(nccdg); dg_dispErr.zero();
        dg_totalLoadLevel.resize(nccdg); dg_totalLoadLevel.zero();
        dg_totalDisp.resize(nccdg); dg_totalDisp.zero();
        // loop over dof managers
        int ndofman = domain->giveNumberOfDofManagers();
        for ( int idofman = 1; idofman <= ndofman; idofman++ ) {
            DofManager *dofman = domain->giveDofManager(idofman);
 #if ( defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE ) )
            if ( !parallel_context->isLocal(dofman) ) {
                continue;
            }

 #endif

            // loop over individual dofs
            int ndof = dofman->giveNumberOfDofs();
            for ( int idof = 1; idof <= ndof; idof++ ) {
                Dof *dof = dofman->giveDof(idof);
                if ( !dof->isPrimaryDof() ) continue;
                int eq = dof->giveEquationNumber(dn);
                int dofid = dof->giveDofID();
                if ( !eq ) continue;
 
                dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                idsInUse.at(dofid) = 1;
            } // end loop over DOFs
        } // end loop over dof managers

        // loop over elements and their DOFs
        int nelem = domain->giveNumberOfElements();
        for ( int ielem = 1; ielem <= nelem; ielem++ ) {
            Element *elem = domain->giveElement(ielem);
 #ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

 #endif
            // loop over element internal Dofs
            for ( int idofman = 1; idofman <= elem->giveNumberOfInternalDofManagers(); idofman++) {
                DofManager *dofman = elem->giveInternalDofManager(idofman);
                int ndof = dofman->giveNumberOfDofs();
                // loop over individual dofs
                for ( int idof = 1; idof <= ndof; idof++ ) {
                    Dof *dof = dofman->giveDof(idof);
                    if ( !dof->isPrimaryDof() ) continue;
                    int eq = dof->giveEquationNumber(dn);
                    int dofid = dof->giveDofID();
                    
                    if ( !eq ) continue;
 #if ( defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE ) )
                    if ( engngModel->isParallel() && !n2l->giveNewEq(eq) ) continue;
 #endif
                    dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                    dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                    dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                    dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                    idsInUse.at(dofid) = 1;
                } // end loop over DOFs
            } // end loop over element internal dofmans
        } // end loop over elements
        
        // loop over boundary conditions and their internal DOFs
        for ( int ibc = 1; ibc <= domain->giveNumberOfBoundaryConditions(); ibc++ ) {
            GeneralBoundaryCondition *bc = domain->giveBc(ibc);

            // loop over element internal Dofs
            for ( int idofman = 1; idofman <= bc->giveNumberOfInternalDofManagers(); idofman++) {
                DofManager *dofman = bc->giveInternalDofManager(idofman);
                int ndof = dofman->giveNumberOfDofs();
                // loop over individual dofs
                for ( int idof = 1; idof <= ndof; idof++ ) {
                    Dof *dof = dofman->giveDof(idof);
                    if ( !dof->isPrimaryDof() ) continue;
                    int eq = dof->giveEquationNumber(dn);
                    int dofid = dof->giveDofID();

                    if ( !eq ) continue;
 #if ( defined ( __PARALLEL_MODE ) && defined ( __PETSC_MODULE ) )
                    if ( engngModel->isParallel() && !n2l->giveNewEq(eq) ) continue;
 #endif
                    dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                    dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                    dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                    dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                    idsInUse.at(dofid) = 1;
                } // end loop over DOFs
            } // end loop over element internal dofmans
        } // end loop over elements

 #ifdef __PARALLEL_MODE
        // exchange individual partition contributions (simultaneously for all groups)
#ifdef __PETSC_MODULE
        FloatArray collectiveErr(nccdg);
        parallel_context->accumulate(dg_forceErr,       collectiveErr); dg_forceErr       = collectiveErr;
        parallel_context->accumulate(dg_dispErr,        collectiveErr); dg_dispErr        = collectiveErr;
        parallel_context->accumulate(dg_totalLoadLevel, collectiveErr); dg_totalLoadLevel = collectiveErr;
        parallel_context->accumulate(dg_totalDisp,      collectiveErr); dg_totalDisp      = collectiveErr;
#else
        if ( this->engngModel->isParallel() ) {
            FloatArray collectiveErr(nccdg);
            MPI_Allreduce(dg_forceErr.givePointer(), collectiveErr.givePointer(), nccdg, MPI_DOUBLE, MPI_SUM, comm);
            dg_forceErr = collectiveErr;
            MPI_Allreduce(dg_dispErr.givePointer(), collectiveErr.givePointer(), nccdg, MPI_DOUBLE, MPI_SUM, comm);
            dg_dispErr = collectiveErr;
            MPI_Allreduce(dg_totalLoadLevel.givePointer(), collectiveErr.givePointer(), nccdg, MPI_DOUBLE, MPI_SUM, comm);
            dg_totalLoadLevel = collectiveErr;
            MPI_Allreduce(dg_totalDisp.givePointer(), collectiveErr.givePointer(), nccdg, MPI_DOUBLE, MPI_SUM, comm);
            dg_totalDisp = collectiveErr;
            return globalNorm;
        }
#endif
 #endif
        OOFEM_LOG_INFO("NRSolver: %-5d", nite);
        // loop over dof groups and check convergence individually
        for ( int dg = 1; dg <= nccdg; dg++ ) {
            // Skips the ones which aren't used in this problem (the residual will be zero for these anyway, but it is annoying to print them all)
            if ( !idsInUse.at(dg) ) {
                continue;
            }
            //  compute a relative error norm
            if ( ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) > nrsolver_ERROR_NORM_SMALL_NUM ) {
                forceErr = sqrt( dg_forceErr.at(dg) / ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) );
            } else {
                 // If both external forces and internal ebe norms are zero, then the residual must be zero.
                forceErr = sqrt( dg_forceErr.at(dg) );
            }

            //
            // compute displacement error
            //
            if ( dg_totalDisp.at(dg) >  nrsolver_ERROR_NORM_SMALL_NUM ) {
                dispErr = sqrt( dg_dispErr.at(dg) / dg_totalDisp.at(dg) );
            } else {
                dispErr = sqrt( dg_dispErr.at(dg) );
            }

            if ( forceErr > rtolf.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ||
                 dispErr  > rtold.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                errorOutOfRange = true;
            }

            if ( forceErr > rtolf.at(1) || dispErr > rtold.at(1) ) {
                answer = false;
            }

            OOFEM_LOG_INFO( "  %s: %.3e %.3e", __DofIDItemToString((DofIDItem)dg), forceErr, dispErr );
        }
        OOFEM_LOG_INFO("\n");
    } else { // No dof grouping
        double dXX, dXdX;
 #ifdef __PARALLEL_MODE
        forceErr = parallel_context->norm(rhs); forceErr *= forceErr;
        dXX = parallel_context->localNorm(X); dXX *= dXX; // Note: Solutions are always total global values (natural distribution makes little sense for the solution)
        dXdX = parallel_context->localNorm(ddX); dXdX *= dXdX;
 #else
        forceErr = rhs.computeSquaredNorm();
        dXX = X.computeSquaredNorm();
        dXdX = ddX.computeSquaredNorm();
 #endif
        // we compute a relative error norm
        if ( ( RRT + internalForcesEBENorm.at(1) ) > nrsolver_ERROR_NORM_SMALL_NUM ) {
            forceErr = sqrt( forceErr / ( RRT + internalForcesEBENorm.at(1) ) );
        } else {
            forceErr = sqrt( forceErr ); // absolute norm as last resort
        }

        // compute displacement error
        // err is relative displacement change
        if ( dXX > nrsolver_ERROR_NORM_SMALL_NUM ) {
            dispErr = sqrt( dXdX / dXX );
        } else {
            dispErr = sqrt( dXdX );
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

} // end namespace oofem
