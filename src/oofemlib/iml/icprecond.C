/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***               v. 1.5c        */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// adapted & optimized by Borek Patzak

#include "icprecond.h"
#include "symcompcol.h"
#include "mathfem.h"

namespace oofem {
CompCol_ICPreconditioner :: CompCol_ICPreconditioner(const SparseMtrx &A, InputRecord &attributes) :
    Preconditioner(A, attributes)
{ }

IRResultType
CompCol_ICPreconditioner :: initializeFrom(InputRecord *ir)
{
    return Preconditioner :: initializeFrom(ir);
}

void
CompCol_ICPreconditioner :: init(const SparseMtrx &A)
{
    if ( dynamic_cast< const SymCompCol * >(&A) ) {
        this->initialize( * ( ( SymCompCol * ) & A ) );
    } else if ( dynamic_cast< const CompCol * >(&A) ) {
        this->initialize( * ( ( CompCol * ) & A ) );
    } else {
        OOFEM_ERROR("unsupported sparse matrix type");
    }
}




void
CompCol_ICPreconditioner :: initialize(const CompCol &A)
{
    dim_ [ 0 ] = A.dim(0);
    dim_ [ 1 ] = A.dim(1);

    pntr_.resize(A.dim(1) + 1);

    nz_ = 0;
    for ( int k = 0; k < dim_ [ 1 ]; k++ ) {
        for ( int j = A.col_ptr(k); j < A.col_ptr(k + 1); j++ ) {
            if ( A.row_ind(j) >= k ) {
                nz_++;
            }
        }
    }

    val_.resize(nz_);
    indx_.resize(nz_);

    // Copy just triangular part
    pntr_[0] = 0;
    for ( int k = 0; k < dim_ [ 1 ]; k++ ) {
        pntr_[k + 1] = pntr_[k];
        for ( int j = A.col_ptr(k); j < A.col_ptr(k + 1); j++ ) {
            if ( A.row_ind(j) >= k ) {
                int i = pntr_[k + 1]++;
                val_[i] = A.val(j);
                indx_[i] = A.row_ind(j);
            }
        }
    }

    for ( int i = 0; i < dim_ [ 1 ]; i++ ) {
        if ( indx_[ pntr_[i] ] != i ) {
            OOFEM_ERROR("diagonal not found!");
        }
    }

    this->ICFactor();
}


void
CompCol_ICPreconditioner :: solve(const FloatArray &x, FloatArray &y) const
{
    y = x;
    this->ICSolve(y);
}


void
CompCol_ICPreconditioner :: trans_solve(const FloatArray &x, FloatArray &y) const
{
    y = x;
    this->ICSolve(y);
}


void
CompCol_ICPreconditioner :: ICFactor()
{
    int n = pntr_.giveSize() - 1;
 
    for ( int k = 0; k < n - 1; k++ ) {
        int d = pntr_[k];
        double z = val_[d] = sqrt( val_[d] );

        for ( int i = d + 1; i < pntr_[k + 1]; i++ ) {
            val_[i] /= z;
        }

        for ( int i = d + 1; i < pntr_[k + 1]; i++ ) {
            double z = val_[i];
            int h = indx_[i];
            int g = i;

            for ( int j = pntr_[h]; j < pntr_[h + 1]; j++ ) {
                for ( ; g < pntr_[k + 1] && indx_[g + 1] <= indx_[j]; g++ ) {
                    if ( indx_[g] == indx_[j] ) {
                        val_[j] -= z * val_[g];
                    }
                }
            }
        }
    }

    int d = pntr_[n - 1];
    val_[d] = sqrt( val_[d] );
}


void
CompCol_ICPreconditioner :: ICSolve(FloatArray &dest) const
{
    int M = dest.giveSize();
    FloatArray work(M);
    // lower diag

    // solve Uw=x
    for ( int i = 0; i < M; i++ ) {
        double val;
        work[i] = val = ( dest[i] + work[i] ) / val_( pntr_[i] );
        for ( int t = pntr_[i] + 1; t < pntr_[i + 1]; t++ ) {
            work( indx_[t] ) -= val_[t] * val;
        }
    }

    dest.zero();
    // solve U^Ty=w
    for ( int i = M - 1; i >= 0; i-- ) {
        for ( int t = pntr_[i] + 1; t < pntr_[i + 1]; t++ ) {
            dest[i] -= val_[t] * dest( indx_[t] );
        }

        dest[i] = ( work[i] + dest[i] ) / val_( pntr_[i] );
    }
}


void
CompCol_ICPreconditioner :: qsortRow(IntArray &src, FloatArray &val, int l, int r)
{
    if ( r <= l ) {
        return;
    }

    int i = qsortRowPartition(src, val, l, r);
    qsortRow(src, val, l, i - 1);
    qsortRow(src, val, i + 1, r);
}


int
CompCol_ICPreconditioner :: qsortRowPartition(IntArray &src, FloatArray &val, int l, int r)
{
    int i = l - 1, j = r;
    int v = src[r];
    int swap;
    double dswap;

    for ( ; ; ) {
        while ( ( src[++i] <  v ) ) {
            ;
        }

        while ( ( v < src[--j] ) ) {
            if ( j == l ) {
                break;
            }
        }

        if ( i >= j ) {
            break;
        }

        swap = src[i];
        src[i] = src[j];
        src[j] = swap;
        dswap = val[i];
        val[i] = val[j];
        val[j] = dswap;
    }

    swap = src[i];
    src[i] = src[r];
    src[r] = swap;
    dswap = val[i];
    val[i] = val[r];
    val[r] = dswap;

    return i;
}
} // end namespace oofem
