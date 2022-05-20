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

#include "sparsemtrx.h"
#include "floatarray.h"
#include "compcol.h"
#include "linsystsolvertype.h"
#include "classfactory.h"

#include "SUPERLU_MT/include/slu_mt_ddefs.h"
#include "superlusolver.h"
#include <stdlib.h>
#include <math.h>
#include "verbose.h"
//#include "globals.h"

#ifdef TIME_REPORT
 #include "timer.h"
#endif

namespace oofem {
REGISTER_SparseLinSolver(SuperLUSolver, ST_SuperLU_MT)


SuperLUSolver :: SuperLUSolver(Domain *d, EngngModel *m) : SparseLinearSystemNM(d, m)
{}


SuperLUSolver :: ~SuperLUSolver() {
  if (this->AAllocated) {
    Destroy_SuperMatrix_Store(& this->A);
    Destroy_SuperNode_SCP(& this->L);
    Destroy_CompCol_NCP(& this->U);
  }
  if (this->permAllocated) {
    SUPERLU_FREE(this->perm_r);
    SUPERLU_FREE(this->perm_c);
  }
}


void
SuperLUSolver :: initializeFrom(InputRecord &ir)
{
  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = 0: natural ordering
   *   permc_spec = 1: minimum degree ordering on structure of A'*A
   *   permc_spec = 2: minimum degree ordering on structure of A'+A
   *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
   */

  this->permc_spec = 0; // default
  IR_GIVE_OPTIONAL_FIELD(ir, this->permc_spec, _IFT_SuperLUSolver_Permcspec);
}

NM_Status
SuperLUSolver :: solve(SparseMtrx &Lhs, FloatArray &b, FloatArray &x)
{
    //1. Step: Transform SparseMtrx *A to SuperMatrix
    //2. Step: Transfrom FloatArray *b to SuperVector
    //3. Step: Transfrom FLoatArray *x to SuperVector
    //3. Step: call superLU Solve with tranformed input parameters
    //4. Step: Transfrom result back to FloatArray *x
    //5. Step: return SUCCESS!
#ifdef TIME_REPORT
  Timer timer;
  timer.startTimer();
#endif
	  

    CompCol *CC = dynamic_cast< CompCol * >( & Lhs );

    if ( CC ) {
        SuperMatrix B, X;
        int_t nprocs;
        fact_t fact;
        trans_t trans;
        yes_no_t refact, usepr;
        equed_t equed;
        double *a;
        int_t *asub, *xa;
        void *work;
        superlumt_options_t superlumt_options;
        int_t info, lwork, nrhs, /*ldx,*/ panel_size, relax;
        int_t m, n, nnz;
        double *rhsb, *rhsx /*, *xact*/;
        double *R, *C;
        double *ferr, *berr;
        double /*u,*/ drop_tol, rpg, rcond;
        superlu_memusage_t superlu_memusage;
        void parse_command_line();


        /* Default parameters to control factorization. */
        nprocs = omp_get_max_threads(); //omp_get_num_threads() does not work as we are not in parallel region!!;
        printf("Use number of LU threads: %u\n", nprocs);
        fact  = DOFACT; //EQUILIBRATE;
        trans = NOTRANS;
        equed = NOEQUIL;
        refact = NO;
        panel_size = sp_ienv(1);
        relax = sp_ienv(2);
        //u     = 1.0;
        usepr = NO;
        drop_tol = 0.0;
        lwork = 0;
        nrhs  = 1;

        m =  CC->giveNumberOfRows();
        n = CC->giveNumberOfColumns();
        nnz = CC->giveNumberOfNonzeros();
        a =  CC->giveValues().givePointer();
        asub =  CC->giveRowIndex().givePointer();
        xa = CC->giveColPtr().givePointer();

	if (0) {
	  dCreate_CompCol_Matrix(& this->A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
	  if ( !( this->perm_r = intMalloc(m) ) ) {
	    SUPERLU_ABORT("Malloc fails for perm_r[].");
	  }
	  if ( !( this->perm_c = intMalloc(n) ) ) {
	    SUPERLU_ABORT("Malloc fails for perm_c[].");
	  }
	  dCreate_Dense_Matrix(& B, m, nrhs, b.givePointer(), m, SLU_DN, SLU_D, SLU_GE);
	  
	  get_perm_c(this->permc_spec, &this->A, this->perm_c);
	  //dPrint_Dense_Matrix(&B);
	  pdgssv(nprocs, &this->A, this->perm_c, this->perm_r, &this->L, &this->U, &B, &info);
	  //dPrint_Dense_Matrix(&B);
	  x.resize( m );
	  this->convertRhs(& B, x);
	  Destroy_SuperMatrix_Store(& this->A);
	  Destroy_SuperMatrix_Store(& B);
	  Destroy_SuperNode_SCP(& this->L);
	  Destroy_CompCol_NCP(& this->U);
	  SUPERLU_FREE(this->perm_r);
	  SUPERLU_FREE(this->perm_c);
	} else { 

	  /* Command line options to modify default behavior. */
	  //parse_command_line(argc, argv, &nprocs, &lwork, &panel_size, &relax,
	  //	       &u, &fact, &trans, &refact, &equed);
	  
	  if ( lwork > 0 ) {
            work = SUPERLU_MALLOC(lwork);
            OOFEM_LOG_DEBUG("Use work space of size LWORK = " IFMT " bytes\n", lwork);
            if ( !work ) {
	      SUPERLU_ABORT("cannot allocate work[]");
            }
	  }
	  
	  //printf("Use work space of size LWORK = " IFMT " bytes\n", lwork);
	  
#if ( PRNTlevel == 1 )
	  cpp_defs();
	  printf( "int_t %d bytes\n", sizeof( int_t ) );
#endif
	  
	  if (CC->giveVersion() == this->lhsVersion) {
	    // reuse existing factorization
	    fact = FACTORED;
	    OOFEM_LOG_DEBUG ("SuperLU_MT:LHS already factored\n");
	  } else {
	    fact  = DOFACT; //EQUILIBRATE;
	    // solve for new (updated) lhs
	    if (this->AAllocated) Destroy_SuperMatrix_Store(& this->A);
	    dCreate_CompCol_Matrix(& this->A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
	    this->AAllocated = true;
	    this->lhsVersion = CC->giveVersion();
	    OOFEM_LOG_DEBUG ("SuperLU_MT:Factoring.....\n");
	    
	    if (!this->permAllocated) {
	      if ( !( this->perm_r = intMalloc(m) ) ) {
		SUPERLU_ABORT("Malloc fails for perm_r[].");
	      }
	      if ( !( this->perm_c = intMalloc(n) ) ) {
		SUPERLU_ABORT("Malloc fails for perm_c[].");
	      }
	      this->permAllocated = true;
	    }
	  }
	  
	  if ( !( rhsb = doubleMalloc(m * nrhs) ) ) {
            SUPERLU_ABORT("Malloc fails for rhsb[].");
	  }
	  if ( !( rhsx = doubleMalloc(m * nrhs) ) ) {
            SUPERLU_ABORT("Malloc fails for rhsx[].");
	  }
	  dCreate_Dense_Matrix(& B, m, nrhs, b.givePointer(), m, SLU_DN, SLU_D, SLU_GE);
	  dCreate_Dense_Matrix(& X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
	  //dPrint_Dense_Matrix(&B);
	  //dPrint_Dense_Matrix(&X);
	  
	  if ( !( R = ( double * ) SUPERLU_MALLOC( this->A.nrow * sizeof( double ) ) ) ) {
            SUPERLU_ABORT("SUPERLU_MALLOC fails for R[].");
	  }
	  if ( !( C = ( double * ) SUPERLU_MALLOC( this->A.ncol * sizeof( double ) ) ) ) {
            SUPERLU_ABORT("SUPERLU_MALLOC fails for C[].");
	  }
	  if ( !( ferr = ( double * ) SUPERLU_MALLOC( nrhs * sizeof( double ) ) ) ) {
            SUPERLU_ABORT("SUPERLU_MALLOC fails for ferr[].");
	  }
	  if ( !( berr = ( double * ) SUPERLU_MALLOC( nrhs * sizeof( double ) ) ) ) {
            SUPERLU_ABORT("SUPERLU_MALLOC fails for berr[].");
	  }
	  
	  get_perm_c(this->permc_spec, & this->A, this->perm_c);
	  
	  superlumt_options.SymmetricMode = YES;
	  superlumt_options.diag_pivot_thresh = 0.0;
	  
	  superlumt_options.nprocs = nprocs;
	  superlumt_options.fact = fact;
	  superlumt_options.trans = trans;
	  superlumt_options.refact = refact;
	  superlumt_options.panel_size = panel_size;
	  superlumt_options.relax = relax;
	  superlumt_options.usepr = usepr;
	  superlumt_options.drop_tol = drop_tol;
	  superlumt_options.PrintStat = NO;
	  superlumt_options.perm_c = perm_c;
	  superlumt_options.perm_r = perm_r;
	  //superlumt_options.work = work;
	  superlumt_options.lwork = lwork;
	  if ( !( superlumt_options.etree = intMalloc(n) ) ) {
	    SUPERLU_ABORT("Malloc fails for etree[].");
	  }
	  if ( !( superlumt_options.colcnt_h = intMalloc(n) ) ) {
	    SUPERLU_ABORT("Malloc fails for colcnt_h[].");
	  }
	  if ( !( superlumt_options.part_super_h = intMalloc(n) ) ) {
	    SUPERLU_ABORT("Malloc fails for colcnt_h[].");
	  }
	  
	  OOFEM_LOG_DEBUG ("sym_mode %d\tdiag_pivot_thresh %.4e\n",
		 superlumt_options.SymmetricMode,
		 superlumt_options.diag_pivot_thresh);
	  
	  /*
	   * Solve the system and compute the condition number
	   * and error bounds using pdgssvx.
	   */
	  pdgssvx(nprocs, & superlumt_options, & this->A, this->perm_c, this->perm_r, & equed, R, C, & this->L, & this->U, & B, & X, & rpg, & rcond, ferr, berr, & superlu_memusage, & info);
	  x.resize( b.giveSize() );
	  this->convertRhs(& X, x);
	  
#if 0
	  printf("psgssvx(): info " IFMT "\n", info);
	  
	  if ( info == 0 || info == n + 1 ) {
            SCPformat *Lstore;
            NCPformat *Ustore;
	    
            printf("Recip. pivot growth = %e\n", rpg);
            printf("Recip. condition number = %e\n", rcond);
            printf("%8s%16s%16s\n", "rhs", "FERR", "BERR");
            for ( int i = 0; i < nrhs; ++i ) {
	      printf(IFMT "%16e%16e\n", i + 1, ferr [ i ], berr [ i ]);
            }
	    
            Lstore = ( SCPformat * ) L.Store;
            Ustore = ( NCPformat * ) U.Store;
            printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
            printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
            printf("No of nonzeros in L+U = " IFMT "\n", Lstore->nnz + Ustore->nnz - n);
            printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions " IFMT "\n",
                   superlu_memusage.for_lu / 1e6, superlu_memusage.total_needed / 1e6,
                   superlu_memusage.expansions);
	    
            fflush(stdout);
	  } else if ( info > 0 && lwork == -1 ) {
            printf("** Estimated memory: " IFMT " bytes\n", info - n);
	  }
#else
	  if ( info > 0 && lwork == -1 ) {
            printf("** Estimated memory: " IFMT " bytes\n", info - n);
	  }
#endif
	  
	  SUPERLU_FREE(rhsb);
	  SUPERLU_FREE(rhsx);
	  //SUPERLU_FREE (xact);
	  SUPERLU_FREE(R);
	  SUPERLU_FREE(C);
	  SUPERLU_FREE(ferr);
	  SUPERLU_FREE(berr);
	  //Destroy_SuperMatrix_Store(& this->A);
	  Destroy_SuperMatrix_Store(& B);
	  Destroy_SuperMatrix_Store(& X);
	  //    SUPERLU_FREE (superlumt_options.etree);
	  //    SUPERLU_FREE (superlumt_options.colcnt_h);
	  ///   SUPERLU_FREE (superlumt_options.part_super_h);
	  if ( lwork == 0 ) {
	    //Destroy_SuperNode_SCP(& this->L);
	    //Destroy_CompCol_NCP(& this->U);
	  } else if ( lwork > 0 ) {
            SUPERLU_FREE(work);
	  }
	}
#ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_INFO( "SuperLU_MT info: user time consumed by solution: %.2fs\n", timer.getUtime() );
    OOFEM_LOG_INFO( "SuperLU_MT info: real time consumed by solution: %.2fs\n", timer.getWtime() );
#endif
    
    } else {
      OOFEM_ERROR("Incompatible sparse storage encountered");
    }
    
    //dPrint_Dense_Matrix(&B);
    /*
     * Lstore = (SCPformat *) L.Store;
     * Ustore = (NCPformat *) U.Store;
     * printf("No of nonzeros in factor L = " IFMT "\n", Lstore->nnz);
     * printf("No of nonzeros in factor U = " IFMT "\n", Ustore->nnz);
     *
     * fflush(stdout);
     * dPrint_CompCol_Matrix(&L);
     * dPrint_CompCol_Matrix(&U);
     */
    //solved = 1;
    return NM_Success;
}

void
SuperLUSolver :: convertRhs(SuperMatrix *X, FloatArray &x)
{
    x.zero();

    if ( ( X->Stype == SLU_DN ) && ( X->Dtype == SLU_D ) ) {
        // process only first column                                                    //int r;

        DNformat *data = ( ( DNformat * ) ( X->Store ) );

#pragma omp parallel
        {
#pragma omp parallel for
            for ( int r = 0; r <= data->lda - 1; r++ ) {
                //r = data->nzval[i];
                x.at(r + 1) = ( ( double * ) data->nzval ) [ r ];
            }
        }
    }   else {
        OOFEM_ERROR("convertRhs: unsupported matrix storage type or data type");
    }
}

//dPrint - following code use to print SuperLU matrix
int_t SuperLUSolver :: dPrint_CompCol_Matrix(SuperMatrix *A)
{
    NCformat *Astore;
    double *dp;

    printf("\nCompCol matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype, A->Dtype, A->Mtype);
    Astore = ( NCformat * ) A->Store;
    dp = ( double * ) Astore->nzval;
    printf("nrow " IFMT ", ncol " IFMT ", nnz " IFMT "\n", A->nrow, A->ncol, Astore->nnz);
    printf("\nnzval: ");
    for ( int_t i = 0; i < Astore->nnz; ++i ) {
        printf("%f  ", dp [ i ]);
    }
    printf("\nrowind: ");
    for ( int_t i = 0; i < Astore->nnz; ++i ) {
        printf(IFMT, Astore->rowind [ i ]);
    }
    printf("\ncolptr: ");
    for ( int_t i = 0; i <= A->ncol; ++i ) {
        printf(IFMT, Astore->colptr [ i ]);
    }
    printf("\nend CompCol matrix.\n");

    return 0;
}

int_t SuperLUSolver :: dPrint_Dense_Matrix(SuperMatrix *A)
{
    DNformat *Astore;
    double *dp;

    printf("\nDense matrix: ");
    printf("Stype %d , Dtype %d , Mtype %d\n", A->Stype, A->Dtype, A->Mtype);
    Astore = ( DNformat * ) A->Store;
    dp = ( double * ) Astore->nzval;
    printf("nrow " IFMT ", ncol " IFMT ", lda " IFMT "\n",
           A->nrow, A->ncol, Astore->lda);
    printf("\nnzval: ");
    for ( int_t i = 0; i < A->nrow; ++i ) {
        printf("%f  ", dp [ i ]);
    }
    printf("\nend Dense matrix.\n");

    return 0;
}

int_t
dCheckZeroDiagonal(int_t n, int_t *rowind, int_t *colbeg,
                   int_t *colend, int_t *perm)
{
    int_t nd = 0;

    for ( int_t j = 0; j < n; ++j ) {
        int_t nzd = 0;
        for ( int_t i = colbeg [ j ]; i < colend [ j ]; ++i ) {
            if ( perm [ rowind [ i ] ] == j ) {
                nzd = 1;
                ++nd;
                break;
            }
        }
        if ( nzd == 0 ) {
            printf("Zero diagonal at column " IFMT "\n", j);
        }
    }

    printf(".. dCheckZeroDiagonal() -- # diagonals " IFMT "\n", nd);

    return 0;
}
} // end namespace oofem
