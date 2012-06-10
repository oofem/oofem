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
 #include "clock.h"
#endif

#ifdef DynCompCol_USE_STL_SETS
 #ifndef __MAKEDEPEND
  #include <map>
 #endif
#endif

namespace oofem {
CompCol_ILUPreconditioner ::
CompCol_ILUPreconditioner(const SparseMtrx &A, InputRecord &attributes) : Preconditioner(A, attributes)
{ }

IRResultType
CompCol_ILUPreconditioner :: initializeFrom(InputRecord *ir)
{
    Preconditioner :: initializeFrom(ir);
    return IRRT_OK;
}


void
CompCol_ILUPreconditioner :: init(const SparseMtrx &A)
{
#ifdef TIME_REPORT
    //clock_t tstart = clock();
    oofem_timeval tstart;
    getUtime(tstart);
#endif

    if ( A.giveType() == SMT_CompCol ) {
        this->initialize( * ( ( CompCol * ) & A ) );
    } else if ( A.giveType() == SMT_DynCompCol ) {
        this->initialize( * ( ( DynCompCol * ) & A ) );
    } else {
        OOFEM_ERROR("CompCol_ILUPreconditioner::init : unsupported sparse matrix type");
    }

#ifdef TIME_REPORT
    oofem_timeval ut;
    getRelativeUtime(ut, tstart);
    OOFEM_LOG_INFO( "ILUP: user time consumed by factorization: %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif
}


void
CompCol_ILUPreconditioner :: initialize(const CompCol &A)
{
    int i, j, k, pn, qn, rn;
    double multiplier;

    // Copy
    dim_ [ 0 ] = A.giveNumberOfRows();
    dim_ [ 1 ] = A.giveNumberOfColumns();

    l_colptr_.resize(A.dim(1) + 1);
    u_colptr_.resize(A.dim(1) + 1);

    l_nz_ = 0;
    u_nz_ = 0;

    // Get size of l and u
    for ( i = 0; i < dim_ [ 1 ]; i++ ) {
        for ( j = A.col_ptr(i); j < A.col_ptr(i + 1); j++ ) {
            if ( A.row_ind(j) > i ) {
                l_nz_++;
            } else {
                u_nz_++;
            }
        }
    }

    l_val_.resize(l_nz_);
    u_val_.resize(u_nz_);
    l_rowind_.resize(l_nz_);
    u_rowind_.resize(u_nz_);

    l_colptr_(0) = u_colptr_(0) = 0;

    // Split up A into l and u
    for ( i = 0; i < dim_ [ 1 ]; i++ ) {
        l_colptr_(i + 1) = l_colptr_(i);
        u_colptr_(i + 1) = u_colptr_(i);

        for ( j = A.col_ptr(i); j < A.col_ptr(i + 1); j++ ) {
            if ( A.row_ind(j) > i ) {
                k = l_colptr_(i + 1)++;
                l_val_(k) = A.val(j);
                l_rowind_(k) = A.row_ind(j);
            } else if ( A.row_ind(j) <= i ) {
                k = u_colptr_(i + 1)++;
                u_val_(k) = A.val(j);
                u_rowind_(k) = A.row_ind(j);
            }
        }
    }

    /* sort entries to assure entries stored with increasing row index */
    /*
     * for (i = 0; i < dim_[1]; i++) {
     *  QSort(l_rowind_, l_val_, l_colptr_[i], l_colptr_[i+1] - l_colptr_[i]);
     *  QSort(u_rowind_, u_val_, u_colptr_[i], u_colptr_[i+1] - u_colptr_[i]);
     * }
     */
    // Factor matrix
    for ( i = 0; i < dim_ [ 0 ] - 1; i++ ) {
        multiplier = u_val_(u_colptr_(i + 1) - 1);

        for ( j = l_colptr_(i); j < l_colptr_(i + 1); j++ ) {
            l_val_(j) /= multiplier;
        }

        for ( j = u_colptr_(i + 1); j < u_colptr_(i + 2) - 1; j++ ) {
            multiplier = u_val_(j);
            qn = j + 1;
            rn = l_colptr_(i + 1);
            for ( pn = l_colptr_( u_rowind_(j) );
                  pn < l_colptr_(u_rowind_(j) + 1) && l_rowind_(pn) <= i + 1; pn++ ) {
                while ( qn < u_colptr_(i + 2) && u_rowind_(qn) < l_rowind_(pn) ) {
                    qn++;
                }

                if ( qn < u_colptr_(i + 2) && l_rowind_(pn) == u_rowind_(qn) ) {
                    u_val_(qn) -= multiplier * l_val_(pn);
                }
            }

            for ( ; pn < l_colptr_(u_rowind_(j) + 1); pn++ ) {
                while ( rn < l_colptr_(i + 2) && l_rowind_(rn) < l_rowind_(pn) ) {
                    rn++;
                }

                if ( rn < l_colptr_(i + 2) && l_rowind_(pn) == l_rowind_(rn) ) {
                    l_val_(rn) -= multiplier * l_val_(pn);
                }
            }
        }
    }
}


void
CompCol_ILUPreconditioner :: initialize(const DynCompCol &A)
{
    int i, j, k, pn, qn, rn;
    double multiplier;

    // Copy
    dim_ [ 0 ] = A.giveNumberOfRows();
    dim_ [ 1 ] = A.giveNumberOfColumns();

    l_colptr_.resize(A.giveNumberOfColumns() + 1);
    u_colptr_.resize(A.giveNumberOfColumns() + 1);

    l_nz_ = 0;
    u_nz_ = 0;

#ifndef DynCompCol_USE_STL_SETS
    // Get size of l and u
    for ( i = 0; i < dim_ [ 1 ]; i++ ) {
        for ( j = 0; j < A.row_ind(i)->giveSize(); j++ ) {
            if ( A.row_ind(i)->at(j + 1) > i ) {
                l_nz_++;
            } else {
                u_nz_++;
            }
        }
    }

    l_val_.resize(l_nz_);
    u_val_.resize(u_nz_);
    l_rowind_.resize(l_nz_);
    u_rowind_.resize(u_nz_);

    l_colptr_(0) = u_colptr_(0) = 0;

    // Split up A into l and u
    for ( i = 0; i < dim_ [ 1 ]; i++ ) {
        l_colptr_(i + 1) = l_colptr_(i);
        u_colptr_(i + 1) = u_colptr_(i);

        for ( j = 0; j < A.row_ind(i)->giveSize(); j++ ) {
            if ( A.row_ind(i)->at(j + 1) > i ) {
                k = l_colptr_(i + 1)++;
                l_val_(k) = A.column(i)->at(j + 1);
                l_rowind_(k) = A.row_ind(i)->at(j + 1);
            } else if ( A.row_ind(i)->at(j + 1) <= i ) {
                k = u_colptr_(i + 1)++;
                u_val_(k) = A.column(i)->at(j + 1);
                u_rowind_(k) = A.row_ind(i)->at(j + 1);
            }
        }
    }

#else
    std :: map< int, double > :: iterator pos;
    // Get size of l and u
    for ( i = 0; i < dim_ [ 1 ]; i++ ) {
        for ( pos = A.column(i)->begin(); pos != A.column(i)->end(); ++pos ) {
            if ( pos->first > i ) {
                l_nz_++;
            } else {
                u_nz_++;
            }
        }
    }

    l_val_.resize(l_nz_);
    u_val_.resize(u_nz_);
    l_rowind_.resize(l_nz_);
    u_rowind_.resize(u_nz_);

    l_colptr_(0) = u_colptr_(0) = 0;

    // Split up A into l and u
    for ( i = 0; i < dim_ [ 1 ]; i++ ) {
        l_colptr_(i + 1) = l_colptr_(i);
        u_colptr_(i + 1) = u_colptr_(i);

        for ( pos = A.column(i)->begin(); pos != A.column(i)->end(); ++pos ) {
            if ( pos->first > i ) {
                k = l_colptr_(i + 1)++;
                l_val_(k) = pos->second;
                l_rowind_(k) = pos->first;
            } else if ( pos->first <= i ) {
                k = u_colptr_(i + 1)++;
                u_val_(k) = pos->second;
                u_rowind_(k) = pos->first;
            }
        }
    }

#endif
    /* sort entries to assure entries stored with increasing row index */

    for ( i = 0; i < dim_ [ 1 ]; i++ ) {
        qsortRow(l_rowind_, l_val_, l_colptr_(i), l_colptr_(i + 1) - 1);
        qsortRow(u_rowind_, u_val_, u_colptr_(i), u_colptr_(i + 1) - 1);
    }

    // Factor matrix
    for ( i = 0; i < dim_ [ 0 ] - 1; i++ ) {
        multiplier = u_val_(u_colptr_(i + 1) - 1);

        for ( j = l_colptr_(i); j < l_colptr_(i + 1); j++ ) {
            l_val_(j) /= multiplier;
        }

        for ( j = u_colptr_(i + 1); j < u_colptr_(i + 2) - 1; j++ ) {
            multiplier = u_val_(j);
            qn = j + 1;
            rn = l_colptr_(i + 1);
            for ( pn = l_colptr_( u_rowind_(j) );
                  pn < l_colptr_(u_rowind_(j) + 1) && l_rowind_(pn) <= i + 1; pn++ ) {
                while ( qn < u_colptr_(i + 2) && u_rowind_(qn) < l_rowind_(pn) ) {
                    qn++;
                }

                if ( qn < u_colptr_(i + 2) && l_rowind_(pn) == u_rowind_(qn) ) {
                    u_val_(qn) -= multiplier * l_val_(pn);
                }
            }

            for ( ; pn < l_colptr_(u_rowind_(j) + 1); pn++ ) {
                while ( rn < l_colptr_(i + 2) && l_rowind_(rn) < l_rowind_(pn) ) {
                    rn++;
                }

                if ( rn < l_colptr_(i + 2) && l_rowind_(pn) == l_rowind_(rn) ) {
                    l_val_(rn) -= multiplier * l_val_(pn);
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

    int i, t;

    y.resize(M);
    work.zero();
    // solve Lw=x
    for ( i = 0; i < M; i++ ) {
        work(i) = x(i) + work(i);
        for ( t = l_colptr_(i); t < l_colptr_(i + 1); t++ ) {
            work( l_rowind_(t) ) -= l_val_(t) * work(i);
        }
    }

    y.zero();
    // solve Uy=w
    for ( i = M - 1; i >= 0; i-- ) {
        y(i) = ( work(i) ) / u_val_(u_colptr_(i + 1) - 1);
        for ( t = u_colptr_(i); t < u_colptr_(i + 1) - 1; t++ ) {
            work( u_rowind_(t) ) -= u_val_(t) * y(i);
        }
    }

    //return y;
}


void
CompCol_ILUPreconditioner :: trans_solve(const FloatArray &x, FloatArray &y) const
{
    int M = x.giveSize();
    FloatArray work(M);

    int i, t;
    double val;

    y.resize(M);
    //work.zero();
    // solve for U^Tw = x
    for ( i = 0; i < M; i++ ) {
        val = 0.0;
        for ( t = u_colptr_(i); t < u_colptr_(i + 1) - 1; t++ ) {
            val += u_val_(t) * x( u_rowind_(t) );
        }

        work(i) = ( x(i) - val ) / u_val_(u_colptr_(i + 1) - 1);
    }


    // solve for L^T y = w
    for ( i = M - 1; i >= 0; i-- ) {
        val = 0.0;
        for ( t = l_colptr_(i); t < l_colptr_(i + 1); t++ ) {
            val += l_val_(t) * work( l_rowind_(t) );
        }

        y.at(i) = work(i) - val;
    }

    //return y;
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
