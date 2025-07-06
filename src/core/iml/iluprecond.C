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

#include "dyncompcol.h"
#include "iluprecond.h"
#include "verbose.h"

#ifdef TIME_REPORT
 #include "timer.h"
#endif

namespace oofem {
CompCol_ILUPreconditioner ::
CompCol_ILUPreconditioner(const SparseMtrx &A, InputRecord &attributes) : Preconditioner(A, attributes)
{ }

void
CompCol_ILUPreconditioner :: initializeFrom(InputRecord &ir)
{
    Preconditioner :: initializeFrom(ir);
}


void
CompCol_ILUPreconditioner :: init(const SparseMtrx &A)
{
#ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
#endif

    if ( dynamic_cast< const CompCol * >(&A) ) {
        this->initialize( static_cast< const CompCol & >(A) );
    } else if ( dynamic_cast< const DynCompCol * >(&A) ) {
        this->initialize( static_cast< const DynCompCol & >(A) );
    } else {
        OOFEM_ERROR("unsupported sparse matrix type");
    }

#ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_INFO( "ILUP: user time consumed by factorization: %.2fs\n", timer.getUtime() );
#endif
}


void
CompCol_ILUPreconditioner :: initialize(const CompCol &A)
{
    int pn, qn, rn;

    // Copy
    dim [ 0 ] = A.giveNumberOfRows();
    dim [ 1 ] = A.giveNumberOfColumns();

    l_colptr.resize(dim [ 1 ] + 1);
    u_colptr.resize(dim [ 1 ] + 1);

    l_nz = 0;
    u_nz = 0;

    // Get size of l and u
    for ( int i = 0; i < dim [ 1 ]; i++ ) {
        for ( int j = A.col_ptr(i); j < A.col_ptr(i + 1); j++ ) {
            if ( A.row_ind(j) > i ) {
                l_nz++;
            } else {
                u_nz++;
            }
        }
    }

    l_val.resize(l_nz);
    u_val.resize(u_nz);
    l_rowind.resize(l_nz);
    u_rowind.resize(u_nz);

    l_colptr[0] = u_colptr[0] = 0;

    // Split up A into l and u
    for ( int i = 0; i < dim [ 1 ]; i++ ) {
        l_colptr[i + 1] = l_colptr[i];
        u_colptr[i + 1] = u_colptr[i];

        for ( int j = A.col_ptr(i); j < A.col_ptr(i + 1); j++ ) {
            if ( A.row_ind(j) > i ) {
                int k = l_colptr[i + 1]++;
                l_val[k] = A.values(j);
                l_rowind[k] = A.row_ind(j);
            } else if ( A.row_ind(j) <= i ) {
                int k = u_colptr[i + 1]++;
                u_val[k] = A.values(j);
                u_rowind[k] = A.row_ind(j);
            }
        }
    }

    /* sort entries to assure entries stored with increasing row index */
    /*
     * for (i = 0; i < dim[1]; i++) {
     *  QSort(l_rowind_, l_val_, l_colptr[i], l_colptr[i+1] - l_colptr[i]);
     *  QSort(u_rowind_, u_val_, u_colptr[i], u_colptr[i+1] - u_colptr[i]);
     * }
     */
    // Factor matrix
    for ( int i = 0; i < dim [ 0 ] - 1; i++ ) {
        double multiplier = u_val[u_colptr[i + 1] - 1];

        for ( int j = l_colptr[i]; j < l_colptr[i + 1]; j++ ) {
            l_val[j] /= multiplier;
        }

        for ( int j = u_colptr[i + 1]; j < u_colptr[i + 2] - 1; j++ ) {
            multiplier = u_val[j];
            qn = j + 1;
            rn = l_colptr[i + 1];
            for ( pn = l_colptr[ u_rowind[j] ]; pn < l_colptr[u_rowind[j] + 1] && l_rowind[pn] <= i + 1; pn++ ) {
                while ( qn < u_colptr[i + 2] && u_rowind[qn] < l_rowind[pn] ) {
                    qn++;
                }

                if ( qn < u_colptr[i + 2] && l_rowind[pn] == u_rowind[qn] ) {
                    u_val[qn] -= multiplier * l_val[pn];
                }
            }

            for ( ; pn < l_colptr[u_rowind[j] + 1]; pn++ ) {
                while ( rn < l_colptr[i + 2] && l_rowind[rn] < l_rowind[pn] ) {
                    rn++;
                }

                if ( rn < l_colptr[i + 2] && l_rowind[pn] == l_rowind[rn] ) {
                    l_val[rn] -= multiplier * l_val[pn];
                }
            }
        }
    }
}


void
CompCol_ILUPreconditioner :: initialize(const DynCompCol &A)
{
    int pn, qn, rn;

    // Copy
    dim [ 0 ] = A.giveNumberOfRows();
    dim [ 1 ] = A.giveNumberOfColumns();

    l_colptr.resize(A.giveNumberOfColumns() + 1);
    u_colptr.resize(A.giveNumberOfColumns() + 1);

    l_nz = 0;
    u_nz = 0;

    // Get size of l and u
    for ( int i = 0; i < dim [ 1 ]; i++ ) {
        for ( int j = 0; j < A.row_ind(i).giveSize(); j++ ) {
            if ( A.row_ind(i).at(j + 1) > i ) {
                l_nz++;
            } else {
                u_nz++;
            }
        }
    }

    l_val.resize(l_nz);
    u_val.resize(u_nz);
    l_rowind.resize(l_nz);
    u_rowind.resize(u_nz);

    l_colptr[0] = u_colptr[0] = 0;

    // Split up A into l and u
    for ( int i = 0; i < dim [ 1 ]; i++ ) {
        l_colptr[i + 1] = l_colptr[i];
        u_colptr[i + 1] = u_colptr[i];

        for ( int j = 0; j < A.row_ind(i).giveSize(); j++ ) {
            if ( A.row_ind(i).at(j + 1) > i ) {
                int k = l_colptr[i + 1]++;
                l_val[k] = A.column(i).at(j + 1);
                l_rowind[k] = A.row_ind(i).at(j + 1);
            } else if ( A.row_ind(i).at(j + 1) <= i ) {
                int k = u_colptr[i + 1]++;
                u_val[k] = A.column(i).at(j + 1);
                u_rowind[k] = A.row_ind(i).at(j + 1);
            }
        }
    }

    /* sort entries to assure entries stored with increasing row index */

    for ( int i = 0; i < dim [ 1 ]; i++ ) {
        qsortRow(l_rowind, l_val, l_colptr[i], l_colptr[i + 1] - 1);
        qsortRow(u_rowind, u_val, u_colptr[i], u_colptr[i + 1] - 1);
    }

    // Factor matrix
    for ( int i = 0; i < dim [ 0 ] - 1; i++ ) {
        double multiplier = u_val[u_colptr[i + 1] - 1];

        for ( int j = l_colptr[i]; j < l_colptr[i + 1]; j++ ) {
            l_val[j] /= multiplier;
        }

        for ( int j = u_colptr[i + 1]; j < u_colptr[i + 2] - 1; j++ ) {
            multiplier = u_val[j];
            qn = j + 1;
            rn = l_colptr[i + 1];
            for ( pn = l_colptr[ u_rowind[j] ]; pn < l_colptr[u_rowind[j] + 1] && l_rowind[pn] <= i + 1; pn++ ) {
                while ( qn < u_colptr[i + 2] && u_rowind[qn] < l_rowind[pn] ) {
                    qn++;
                }

                if ( qn < u_colptr[i + 2] && l_rowind[pn] == u_rowind[qn] ) {
                    u_val[qn] -= multiplier * l_val[pn];
                }
            }

            for ( ; pn < l_colptr[u_rowind[j] + 1]; pn++ ) {
                while ( rn < l_colptr[i + 2] && l_rowind[rn] < l_rowind[pn] ) {
                    rn++;
                }

                if ( rn < l_colptr[i + 2] && l_rowind[pn] == l_rowind[rn] ) {
                    l_val[rn] -= multiplier * l_val[pn];
                }
            }
        }
    }
}


void
CompCol_ILUPreconditioner :: solve(const FloatArray &x, FloatArray &y) const
{
    int M = x.giveSize();
    FloatArray work(M);

    y.resize(M);

    // solve Lw=x
    for ( int i = 0; i < M; i++ ) {
        work[i] += x[i];
        for ( int t = l_colptr[i]; t < l_colptr[i + 1]; t++ ) {
            work[ l_rowind[t] ] -= l_val[t] * work[i];
        }
    }

    y.zero();
    // solve Uy=w
    for ( int i = M - 1; i >= 0; i-- ) {
        y[i] = ( work[i] ) / u_val[u_colptr[i + 1] - 1];
        for ( int t = u_colptr[i]; t < u_colptr[i + 1] - 1; t++ ) {
            work[ u_rowind[t] ] -= u_val[t] * y[i];
        }
    }
}


void
CompCol_ILUPreconditioner :: trans_solve(const FloatArray &x, FloatArray &y) const
{
    int M = x.giveSize();
    FloatArray work(M);

    y.resize(M);

    // solve for U^Tw = x
    for ( int i = 0; i < M; i++ ) {
        double val = 0.0;
        for ( int t = u_colptr[i]; t < u_colptr[i + 1] - 1; t++ ) {
            val += u_val[t] * x[ u_rowind[t] ];
        }

        work[i] = ( x[i] - val ) / u_val[u_colptr[i + 1] - 1];
    }


    // solve for L^T y = w
    for ( int i = M - 1; i >= 0; i-- ) {
        double val = 0.0;
        for ( int t = l_colptr[i]; t < l_colptr[i + 1]; t++ ) {
            val += l_val[t] * work( l_rowind[t] );
        }

        y[i] = work[i] - val;
    }
}


void
CompCol_ILUPreconditioner :: qsortRow(IntArray &src, FloatArray &val, int l, int r)
{
    if ( r <= l ) {
        return;
    }

    int i = qsortRowPartition(src, val, l, r);
    qsortRow(src, val, l, i - 1);
    qsortRow(src, val, i + 1, r);
}


int
CompCol_ILUPreconditioner :: qsortRowPartition(IntArray &src, FloatArray &val, int l, int r)
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
