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

#include "slepcsolver.h"

#define TIME_REPORT
#include "petscsparsemtrx.h"
#include "engngm.h"
#include "floatarray.h"
#include "verbose.h"

#ifdef TIME_REPORT
 #include "timer.h"
#endif

namespace oofem {
SLEPcSolver :: SLEPcSolver(Domain *d, EngngModel *m) : SparseGeneralEigenValueSystemNM(d, m)
{
    A = B = NULL;
    epsInit = false;
}


SLEPcSolver :: ~SLEPcSolver()
{
    if ( epsInit ) {
        EPSDestroy(& eps);
    }
}

NM_Status
SLEPcSolver :: solve(SparseMtrx &a, SparseMtrx &b, FloatArray &_eigv, FloatMatrix &_r, double rtol, int nroot)
{
    FILE *outStream;
    PetscErrorCode ierr;
    int size;
    ST st;

    outStream = domain->giveEngngModel()->giveOutputStream();

    // first check whether Lhs is defined

    if ( a->giveNumberOfRows() != a->giveNumberOfColumns() ||
        b->giveNumberOfRows() != b->giveNumberOfRows() ||
        a->giveNumberOfColumns() != b->giveNumberOfColumns() ) {
        OOFEM_ERROR("matrices size mismatch");
    }

    A = dynamic_cast< PetscSparseMtrx * >(&a);
    B = dynamic_cast< PetscSparseMtrx * >(&b);

    if ( !A || !B ) {
        OOFEM_ERROR("PetscSparseMtrx Expected");
    }

    size = engngModel->giveParallelContext( A->giveDomainIndex() )->giveNumberOfNaturalEqs(); // A->giveLeqs();

    _r.resize(size, nroot);
    _eigv.resize(nroot);


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *             Create the eigensolver and set various options
     *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    int nconv, nite;
    EPSConvergedReason reason;

#ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
#endif

    if ( !epsInit ) {
        /*
         * Create eigensolver context
         */
#ifdef __PARALLEL_MODE
        MPI_Comm comm = engngModel->giveParallelComm();
#else
        MPI_Comm comm = PETSC_COMM_SELF;
#endif
        ierr = EPSCreate(comm, & eps);
        CHKERRQ(ierr);
        epsInit = true;
    }

    /*
     * Set operators. In this case, it is a generalized eigenvalue problem
     */

    ierr = EPSSetOperators( eps, * A->giveMtrx(), * B->giveMtrx() );
    CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_GHEP);
    CHKERRQ(ierr);
    ierr = EPSGetST(eps, & st);
    CHKERRQ(ierr);
    ierr = STSetType(st, STSINVERT);
    CHKERRQ(ierr);
    ierr = STSetMatStructure(st, SAME_NONZERO_PATTERN);
    CHKERRQ(ierr);
    ierr = EPSSetTolerances(eps, ( PetscReal ) rtol, PETSC_DECIDE);
    CHKERRQ(ierr);
    ierr = EPSSetDimensions(eps, ( PetscInt ) nroot, PETSC_DECIDE, PETSC_DECIDE);
    CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps, EPS_SMALLEST_MAGNITUDE);
    CHKERRQ(ierr);

    /*
     * Set solver parameters at runtime
     */

    ierr = EPSSetFromOptions(eps);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *                   Solve the eigensystem
     *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = EPSSolve(eps);
    CHKERRQ(ierr);

    ierr = EPSGetConvergedReason(eps, & reason);
    CHKERRQ(ierr);
    ierr = EPSGetIterationNumber(eps, & nite);
    CHKERRQ(ierr);
    OOFEM_LOG_INFO("SLEPcSolver::solve EPSConvergedReason: %d, number of iterations: %d\n", reason, nite);

    ierr = EPSGetConverged(eps, & nconv);
    CHKERRQ(ierr);

    if ( nconv > 0 ) {
        fprintf(outStream, "SLEPcSolver :: solveYourselfAt: Convergence reached for RTOL=%20.15f", rtol);
        PetscScalar kr;
        Vec Vr;

        ierr = MatGetVecs(* B->giveMtrx(), PETSC_NULL, & Vr);
        CHKERRQ(ierr);

        FloatArray Vr_loc;

        for ( int i = 0; i < nconv && i < nroot; i++ ) {
            ierr = EPSGetEigenpair(eps, nconv - i - 1, & kr, PETSC_NULL, Vr, PETSC_NULL);
            CHKERRQ(ierr);

            //Store the eigenvalue
            _eigv->at(i + 1) = kr;

            //Store the eigenvector
            A->scatterG2L(Vr, Vr_loc);
            for ( int j = 0; j < size; j++ ) {
                _r->at(j + 1, i + 1) = Vr_loc.at(j + 1);
            }
        }

        ierr = VecDestroy(& Vr);
        CHKERRQ(ierr);
    } else {
        OOFEM_ERROR("No converged eigenpairs");
    }

#ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_INFO( "SLEPcSolver info: user time consumed by solution: %.2fs\n", timer.getUtime() );
#endif

    return NM_Success;
}
} // end namespace oofem
