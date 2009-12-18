/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/petscsolver.C,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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

#include "petscsolver.h"

#ifdef __PETSC_MODULE
#include "petscsparsemtrx.h"
#include "engngm.h"
#include "flotarry.h"
#include "verbose.h"

namespace oofem {

#define TIME_REPORT

#ifdef TIME_REPORT
#ifndef __MAKEDEPEND
#include <time.h>
#endif
#include "clock.h"
#endif

PetscSolver :: PetscSolver(int i, Domain *d, EngngModel *m) : SparseLinearSystemNM(i, d, m)
{
    Lhs = NULL;
    lhsVersion = 0;
    kspInit = false;
}


PetscSolver :: ~PetscSolver() {
    if ( kspInit ) {
        KSPDestroy(ksp);
    }
}

IRResultType
PetscSolver :: initializeFrom(InputRecord *ir)
{
    /*
     * const char *__keyword, *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
     * IRResultType result;                               // Required by IR_GIVE_FIELD macro
     */
    return IRRT_OK;
}


NM_Status
PetscSolver :: solve(SparseMtrx *A, FloatArray *b, FloatArray *x)
{
    int neqs, l_eqs, g_eqs;

    // first check whether Lhs is defined
    if ( !A ) {
        _error("solveYourselfAt: unknown Lhs");
    }

    // and whether Rhs
    if ( !b ) {
        _error("solveYourselfAt: unknown Rhs");
    }

    // and whether previous Solution exist
    if ( !x ) {
        _error("solveYourselfAt: unknown solution array");
    }

    if ( x->giveSize() != ( neqs = b->giveSize() ) ) {
        _error("solveYourselfAt: size mismatch");
    }

    if ( A->giveType() != SMT_PetscMtrx ) {
        _error("solveYourselfAt: PetscSparseMtrx Expected");
    }

    Lhs = ( PetscSparseMtrx * ) A;
    l_eqs = Lhs->giveLeqs();
    g_eqs = Lhs->giveGeqs();

    Vec globRhsVec;
    Vec globSolVec;

#ifdef __PARALLEL_MODE
    /*
     * scatter and gather rhs to global representation
     */
    engngModel->givePetscContext( Lhs->giveDomainIndex(), Lhs->giveEquationID() )->createVecGlobal(& globRhsVec);
    engngModel->givePetscContext( Lhs->giveDomainIndex(), Lhs->giveEquationID() )->scatterL2G(b, globRhsVec, ADD_VALUES);
#else
    VecCreateSeqWithArray(PETSC_COMM_SELF, b->giveSize(), b->givePointer(), & globRhsVec);
#endif

    VecDuplicate(globRhsVec, & globSolVec);
    /* solve */
    //VecView(globRhsVec,PETSC_VIEWER_STDOUT_WORLD);
    this->petsc_solve(Lhs, globRhsVec, globSolVec);

    //VecView(globSolVec,PETSC_VIEWER_STDOUT_WORLD);


#ifdef __PARALLEL_MODE
    engngModel->givePetscContext( Lhs->giveDomainIndex(), Lhs->giveEquationID() )->scatterG2N(globSolVec, x, INSERT_VALUES);
#else
    double *ptr;
    VecGetArray(globSolVec, & ptr);
    x->resize(l_eqs);
    for ( int i = 1; i <= l_eqs; i++ ) {
        x->at(i) = ptr [ i - 1 ];
    }

    VecRestoreArray(globSolVec, & ptr);
#endif

    VecDestroy(globSolVec);
    VecDestroy(globRhsVec);

    return 1;
}

NM_Status
PetscSolver :: petsc_solve(PetscSparseMtrx *Lhs, Vec b, Vec x)
{
    int nite;
    KSPConvergedReason reason;
    if ( Lhs->giveType() != SMT_PetscMtrx ) {
        _error("solveYourselfAt: PetscSparseMtrx Expected");
    }

#ifdef TIME_REPORT
    //clock_t tstart = clock();
    oofem_timeval tstart;
    getUtime(tstart);
#endif

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    *  Create the linear solver and set various options
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    if ( !kspInit ) {
        /*
         * Create linear solver context
         */
#ifdef __PARALLEL_MODE
        KSPCreate(PETSC_COMM_WORLD, & ksp);
#else
        KSPCreate(PETSC_COMM_SELF, & ksp);
#endif
        kspInit = true;
    }

    /*
     * Set operators. Here the matrix that defines the linear system
     * also serves as the preconditioning matrix.
     */
    KSPSetOperators(ksp, * Lhs->giveMtrx(), * Lhs->giveMtrx(), DIFFERENT_NONZERO_PATTERN);

    /*
     * Set linear solver defaults for this problem (optional).
     * - The following two statements are optional; all of these
     * parameters could alternatively be specified at runtime via
     * KSPSetFromOptions().  All of these defaults can be
     * overridden at runtime, as indicated below.
     */

    //KSPSetTolerances(ksp,1.e-2/((m+1)*(n+1)),1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    /*
     * Set runtime options, e.g.,
     * -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
     * These options will override those specified above as long as
     * KSPSetFromOptions() is called _after_ any other customization
     * routines.
     */
    KSPSetFromOptions(ksp);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    *  Solve the linear system
    *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    //MatView(*Lhs->giveMtrx(),PETSC_VIEWER_STDOUT_SELF);

    //KSPSetRhs(ksp,b);
    //KSPSetSolution(ksp,x);
    KSPSolve(ksp, b, x);
    KSPGetConvergedReason(ksp, & reason);
    KSPGetIterationNumber(ksp, & nite);
    OOFEM_LOG_INFO("PetscSolver::petsc_solve KSPConvergedReason: %d, number of iterations: %d\n", reason, nite);

#ifdef TIME_REPORT
    oofem_timeval ut;
    getRelativeUtime(ut, tstart);
    OOFEM_LOG_INFO( "PetscSolver info: user time consumed by solution: %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif

    return 1;
}

void
PetscSolver :: reinitialize()
{
    if ( kspInit ) {
        KSPDestroy(ksp);
    }

    kspInit = false;
}

} // end namespace oofem
#endif //ifdef __PETSC_MODULE

#ifndef __PETSC_MODULE
namespace oofem {

PetscSolver :: PetscSolver(int i, Domain *d, EngngModel *m) : SparseLinearSystemNM(i, d, m)
{
    _error("PetscSolver: can't create, PETSc support not compiled");
}

PetscSolver :: ~PetscSolver() { }

IRResultType
PetscSolver :: initializeFrom(InputRecord *ir) { return IRRT_OK; }

NM_Status
PetscSolver :: solve(SparseMtrx *A, FloatArray *b, FloatArray *x) { return NM_NoSuccess; }

void
PetscSolver :: reinitialize()
{ }

} // end namespace oofem
#endif
