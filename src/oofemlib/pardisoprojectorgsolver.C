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

#include "pardisoprojectorgsolver.h"

#include "compcol.h"
#include "symcompcol.h"
#include "engngm.h"
#include "floatarray.h"
#include "verbose.h"
#include "timer.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
REGISTER_SparseLinSolver(PardisoProjectOrgSolver, ST_PardisoProjectOrg);

PardisoProjectOrgSolver :: PardisoProjectOrgSolver(Domain *d, EngngModel *m) : SparseLinearSystemNM(d, m) { }

PardisoProjectOrgSolver :: ~PardisoProjectOrgSolver() { }

NM_Status PardisoProjectOrgSolver :: solve(SparseMtrx &A, FloatArray &b, FloatArray &x)
{
    int neqs = b.giveSize();
    x.resize(neqs);

    int mtype = -2;        // Real symmetric positive definite matrix
    CompCol *mat = dynamic_cast< SymCompCol * >(& A);
    if ( !mat ) {
        mtype = 11;        // Real unsymmetric matrix
        mat = dynamic_cast< CompCol * >(& A);
        if ( !mat ) {
            OOFEM_ERROR("CompCol matrix needed for Pardiso solver");
        }
    }

    // Pardiso's CGS-implementation can't handle b = 0.
    if ( b.computeSquaredNorm() == 0 ) {
        return NM_Success;
    }

    int *ia = mat->giveColPtr().givePointer();
    int *ja = mat->giveRowIndex().givePointer();
    double *a = mat->giveValues().givePointer();


    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
    /*     notation.                                                        */
    /* -------------------------------------------------------------------- */
    int n =  mat->giveNumberOfRows();
    for ( int i = 0; i < n + 1; i++ ) {
        ia [ i ] += 1;
    }
    int nnz = ia [ n ];
    for ( int i = 0; i < nnz; i++ ) {
        ja [ i ] += 1;
    }




    Timer timer;
    timer.startTimer();

    // RHS and solution vectors.
    int nrhs = 1;          // Number of right hand sides.

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void *pt [ 64 ];

    // Pardiso control parameters.
    IntArray iparm(64); ///@todo pardisoinit seems to write outside this array
    FloatArray dparm(64);
    int maxfct, mnum, phase, error, msglvl;

    double ddum = 0.;           // Double dummy
    int idum = 0;              // Integer dummy.

    // Setup Pardiso control parameters
    /* -------------------------------------------------------------------- */
    error = 0;
    int solver = 0; // sparse direct
    pardisoinit(pt, & mtype, & solver, iparm.givePointer(), dparm.givePointer(), & error);  // INITIALIZATION!
    if ( error != 0 ) {
        OOFEM_WARNING("Error during pardiso init, error code: %d", error);
        return NM_NoSuccess;
    }
    // Settings are here:
    // https://software.intel.com/en-us/articles/pardiso-parameter-table#table2

    iparm [ 0 ] = 1;
    ///@todo I might be misunderstanding something, but this iterative solver still does a full factorization. No options for incomplete factorizations.
    //iparm[4-1] = 32; // 10*L + K. K = 1 implies CGS (instead of LU), K = 2 implies CG. L specifies exponent tolerance.
    iparm [ 8 - 1 ] = 2; /* Max numbers of iterative refinement steps. */  ///@todo I have no idea if this is suitable value. Examples use 2. / Mikael
    iparm [ 12 - 1 ] = 2; // Transpose (we have a CSC matrix representation here instead of the expected CSR)
    iparm [ 35 - 1 ] = 1; // 1 implies 0-indexing
    //iparm[27-1] = 1; // Checks the matrix (only in MKL)
    ///@todo This is not included in the table of options for some reason!

    maxfct = 1;         // Maximum number of numerical factorizations.
    mnum   = 1;         // Which factorization to use.
    msglvl = 0;         // Print statistical information
    error  = 0;         // Initialize error flag

    /* -------------------------------------------------------------------- */
    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    /*     all memory that is necessary for the factorization.              */
    /* -------------------------------------------------------------------- */
    phase = 11;

    pardiso( pt, & maxfct, & mnum, & mtype, & phase, & neqs,
             a, ( int * ) ia, ( int * ) ja,
             & idum, & nrhs, iparm.givePointer(), & msglvl, & ddum, & ddum, & error, dparm.givePointer() ); // FACTORIZATION!

    if ( error != 0 ) {
        OOFEM_WARNING("Error during symbolic factorization: %d", error);
        return NM_NoSuccess;
    }
    OOFEM_LOG_DEBUG("Reordering completed: %d nonzero factors, %d factorization MFLOPS\n", iparm [ 17 - 1 ], iparm [ 18 - 1 ]);

    /* -------------------------------------------------------------------- */
    /* ..  Numerical factorization.                                         */
    /* -------------------------------------------------------------------- */
    phase = 22;

    pardiso( pt, & maxfct, & mnum, & mtype, & phase, & neqs,
             a, ( int * ) ia, ( int * ) ja,
             & idum, & nrhs, iparm.givePointer(), & msglvl, & ddum, & ddum, & error, dparm.givePointer() );

    if ( error != 0 ) {
        OOFEM_WARNING("ERROR during numerical factorization: %d", error);
        return NM_NoSuccess;
    }
    OOFEM_LOG_DEBUG("Factorization completed ...\n");

    /* -------------------------------------------------------------------- */
    /* ..  Back substitution and iterative refinement.                      */
    /* -------------------------------------------------------------------- */
    phase = 33;

    pardiso( pt, & maxfct, & mnum, & mtype, & phase, & neqs,
             a, ( int * ) ia, ( int * ) ja,
             & idum, & nrhs, iparm.givePointer(), & msglvl, b.givePointer(), x.givePointer(), & error, dparm.givePointer() );

    printf("iparm(20) = %d\n", iparm [ 20 ]);
    if ( error != 0 ) {
        OOFEM_WARNING("ERROR during solution: %d, iparm(20) = %d", error, iparm [ 20 - 1 ]);
        return NM_NoSuccess;
    }

    OOFEM_LOG_DEBUG("Solve completed ... \n");


    /* -------------------------------------------------------------------- */
    /* ..  Convert matrix back to 0-based C-notation.                       */
    /* -------------------------------------------------------------------- */
    for ( int i = 0; i < n + 1; i++ ) {
        ia [ i ] -= 1;
    }
    for ( int i = 0; i < nnz; i++ ) {
        ja [ i ] -= 1;
    }


    /* -------------------------------------------------------------------- */
    /* ..  Termination and release of memory.                               */
    /* -------------------------------------------------------------------- */


    phase = -1;                 /* Release internal memory. */

    pardiso( pt, & maxfct, & mnum, & mtype, & phase,
             & neqs, & ddum, & idum, & idum, & idum, & nrhs,
             iparm.givePointer(), & msglvl, & ddum, & ddum, & error, dparm.givePointer() );

    timer.stopTimer();
    OOFEM_LOG_INFO( "MKLPardisoSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );

    NM_Status s = NM_Success;
    return s;
}

#if 0
///@todo Parallel mode of this.
NM_Status PardisoProjectOrgSolver :: solve(SparseMtrx &A, FloatMatrix &B, FloatMatrix &X)
{
    ///@todo Write subfunction for this as to not repeat everything. Should be easy to add. It is very important to add this support(!!!!!)
}
#endif
} // end namespace oofem
