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


void
CompCol_ICPreconditioner :: initializeFrom(InputRecord &ir)
{
    Preconditioner :: initializeFrom(ir);
}


void
CompCol_ICPreconditioner :: init(const SparseMtrx &A)
{
    if ( dynamic_cast< const SymCompCol * >(&A) ) {
        this->initialize( static_cast< const SymCompCol & >(A) );
    } else if ( dynamic_cast< const CompCol * >(&A) ) {
        this->initialize( static_cast< const CompCol & >(A) );
    } else {
        OOFEM_ERROR("unsupported sparse matrix type");
    }
}


void
CompCol_ICPreconditioner :: initialize(const CompCol &A)
{
    dim [ 0 ] = A.giveNumberOfRows();
    dim [ 1 ] = A.giveNumberOfColumns();

    pntr.resize(dim [ 1 ] + 1);

    nz = 0;
    for ( int k = 0; k < dim [ 1 ]; k++ ) {
        for ( int j = A.col_ptr(k); j < A.col_ptr(k + 1); j++ ) {
            if ( A.row_ind(j) >= k ) {
                nz++;
            }
        }
    }

    val.resize(nz);
    indx.resize(nz);

    // Copy just triangular part
    pntr[0] = 0;
    for ( int k = 0; k < dim [ 1 ]; k++ ) {
        pntr[k + 1] = pntr[k];
        for ( int j = A.col_ptr(k); j < A.col_ptr(k + 1); j++ ) {
            if ( A.row_ind(j) >= k ) {
                int i = pntr[k + 1]++;
                val[i] = A.values(j);
                indx[i] = A.row_ind(j);
            }
        }
    }

    for ( int i = 0; i < dim [ 1 ]; i++ ) {
        if ( indx[ pntr[i] ] != i ) {
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
    int n = pntr.giveSize() - 1;
 
    for ( int k = 0; k < n - 1; k++ ) {
        int d = pntr[k];
        double z = val[d] = sqrt( val[d] );

        for ( int i = d + 1; i < pntr[k + 1]; i++ ) {
            val[i] /= z;
        }

        for ( int i = d + 1; i < pntr[k + 1]; i++ ) {
            double z = val[i];
            int h = indx[i];
            int g = i;

            for ( int j = pntr[h]; j < pntr[h + 1]; j++ ) {
                for ( ; g < pntr[k + 1] && indx[g + 1] <= indx[j]; g++ ) {
                    if ( indx[g] == indx[j] ) {
                        val[j] -= z * val[g];
                    }
                }
            }
        }
    }

    int d = pntr[n - 1];
    val[d] = sqrt( val[d] );
}


void
CompCol_ICPreconditioner :: ICSolve(FloatArray &dest) const
{
    int M = dest.giveSize();
    FloatArray work(M);
    // lower diag

    // solve Uw=x
    for ( int i = 0; i < M; i++ ) {
        double temp;
        work[i] = temp = ( dest[i] + work[i] ) / val[ pntr[i] ];
        for ( int t = pntr[i] + 1; t < pntr[i + 1]; t++ ) {
            work[ indx[t] ] -= val[t] * temp;
        }
    }

    dest.zero();
    // solve U^Ty=w
    for ( int i = M - 1; i >= 0; i-- ) {
        for ( int t = pntr[i] + 1; t < pntr[i + 1]; t++ ) {
            dest[i] -= val[t] * dest( indx[t] );
        }

        dest[i] = ( work[i] + dest[i] ) / val[ pntr[i] ];
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

        std::swap(src[i], src[j]);
        std::swap(val[i], val[j]);
    }

    std::swap(src[i], src[r]);
    std::swap(val[i], val[r]);

    return i;
}
} // end namespace oofem
