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
    if ( A.giveType() == SMT_SymCompCol ) {
        this->initialize( * ( ( SymCompCol * ) & A ) );
    } else if ( A.giveType() == SMT_CompCol ) {
        this->initialize( * ( ( CompCol * ) & A ) );
    } else {
        OOFEM_ERROR("CompCol_ICPreconditioner::init : unsupported sparse matrix type");
    }
}




void
CompCol_ICPreconditioner :: initialize(const CompCol &A)
{
    dim_ [ 0 ] = A.dim(0);
    dim_ [ 1 ] = A.dim(1);

    pntr_.resize(A.dim(1) + 1);


    int i, j, k;

    nz_ = 0;
    for ( k = 0; k < dim_ [ 1 ]; k++ ) {
        for ( j = A.col_ptr(k); j < A.col_ptr(k + 1); j++ ) {
            if ( A.row_ind(j) >= k ) {
                nz_++;
            }
        }
    }

    val_.resize(nz_);
    indx_.resize(nz_);

    // Copy just triangular part
    pntr_(0) = 0;
    for ( k = 0; k < dim_ [ 1 ]; k++ ) {
        pntr_(k + 1) = pntr_(k);
        for ( j = A.col_ptr(k); j < A.col_ptr(k + 1); j++ ) {
            if ( A.row_ind(j) >= k ) {
                i = pntr_(k + 1)++;
                val_(i) = A.val(j);
                indx_(i) = A.row_ind(j);
            }
        }
    }

    //  for (i = 0; i < dim_[1]; i++)
    //    qsortRow(indx_, val_, pntr_(i), pntr_(i+1) -1);

    for ( i = 0; i < dim_ [ 1 ]; i++ ) {
        if ( indx_( pntr_(i) ) != i ) {
            OOFEM_ERROR("CompCol_ICPreconditioner::initialize -  diagonal not found!");
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
    int d, g, h, i, j, k, n = pntr_.giveSize() - 1;
    double z;

    for ( k = 0; k < n - 1; k++ ) {
        d = pntr_(k);
        z = val_(d) = sqrt( val_(d) );

        for ( i = d + 1; i < pntr_(k + 1); i++ ) {
            val_(i) /= z;
        }

        for ( i = d + 1; i < pntr_(k + 1); i++ ) {
            z = val_(i);
            h = indx_(i);
            g = i;

            for ( j = pntr_(h); j < pntr_(h + 1); j++ ) {
                for ( ; g < pntr_(k + 1) && indx_(g + 1) <= indx_(j); g++ ) {
                    if ( indx_(g) == indx_(j) ) {
                        val_(j) -= z * val_(g);
                    }
                }
            }
        }
    }

    d = pntr_(n - 1);
    val_(d) = sqrt( val_(d) );
}


void
CompCol_ICPreconditioner :: ICSolve(FloatArray &dest) const
{
    int i, t, M = dest.giveSize();
    FloatArray work(M);
    double val;
    // lower diag

    work.zero();
    // solve Uw=x
    for ( i = 0; i < M; i++ ) {
        work(i) = val = ( dest(i) + work(i) ) / val_( pntr_(i) );
        for ( t = pntr_(i) + 1; t < pntr_(i + 1); t++ ) {
            work( indx_(t) ) -= val_(t) * val;
        }
    }

    dest.zero();
    // solve U^Ty=w
    for ( i = M - 1; i >= 0; i-- ) {
        for ( t = pntr_(i) + 1; t < pntr_(i + 1); t++ ) {
            dest(i) -= val_(t) * dest( indx_(t) );
        }

        dest(i) = ( work(i) + dest(i) ) / val_( pntr_(i) );
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
    int v = src(r);
    int swap;
    double dswap;

    for ( ; ; ) {
        while ( ( src(++i) <  v ) ) {
            ;
        }

        while ( ( v < src(--j) ) ) {
            if ( j == l ) {
                break;
            }
        }

        if ( i >= j ) {
            break;
        }

        swap = src(i);
        src(i) = src(j);
        src(j) = swap;
        dswap = val(i);
        val(i) = val(j);
        val(j) = dswap;
    }

    swap = src(i);
    src(i) = src(r);
    src(r) = swap;
    dswap = val(i);
    val(i) = val(r);
    val(r) = dswap;

    return i;
}
} // end namespace oofem
