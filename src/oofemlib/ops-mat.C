#include"floatmatrix.h"
#include"floatarray.h"
#include"intarray.h"
#ifndef __MATH_INTERNAL
    #define __MATH_INTERNAL
#endif
#include"ops-mat.h"

#include <numeric>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>

namespace oofem::mat {


// Some forward declarations for LAPACK. Remember to append the underscore to the function name.
#ifdef __LAPACK_MODULE
extern "C" {
/// Computes the reciprocal condition number for a LU decomposed function.
extern void dgecon_(const char *norm, const int *n, const double *a, const int *lda,
                    const double *anorm, double *rcond, double *work, int *iwork, int *info, int norm_len);
/// Replaces a with the LU-decomposition.
extern int dgetrf_(const int *m, const int *n, double *a, const int *lda, int *lpiv, int *info);
/// Replaces a with its inverse.
extern int dgetri_(const int *n, double *a, const int *lda, int *ipiv, double *work, const int *lwork, int *info);
/// Solves a system of equations.
extern int dgesv_(const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, const double *b, const int *ldb, int *info);
/// Computes the norm.
extern double dlange_(const char *norm, const int *m, const int *n, const double *a, const int *lda, double *work, int norm_len);
/// Computes eigenvalues and vectors.
extern int dsyevx_(const char *jobz,  const char *range, const char *uplo, const int *n, double *a, const int *lda,
                   const double *vl, const double *vu, const int *il, const int *iu,
                   const double *abstol, int *m, double *w, double *z, const int *ldz,
                   double *work, int *lwork, int *iwork, int *ifail, int *info,
                   int jobz_len, int range_len, int uplo_len);
/// Solves system which has been LU-factorized.
extern void dgetrs_(const char *trans, const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, const double *b, const int *ldb, int *info);
/// General matrix multiplication
extern void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha,
                   const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc,
                   int a_columns, int b_columns, int c_columns);
/// General dyad product of vectors
extern void dger_(const int *m, const int *n, const double *alpha, const double *x, const int *incx,
                  const double *y, const int *incy, double *a, const int *lda,
                  int x_len, int y_len, int a_columns);
/// Symmetric dyad product of vector
extern void dsyr_(const char *uplo, const int *n, const double *alpha, const double *x, const int *incx,
                  double *a, const int *lda, int x_len, int a_columns);
/// Y = Y + alpha * X
extern void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy, int xsize, int ysize);
/// X = alpha * X
extern void dscal_(const int *n, const double *alpha, const double *x, const int *incx, int size);
}
#endif


void assemble(FloatMatrix& dst, const FloatMatrix &src, const IntArray &loc)
{
    std::size_t ii, jj, size = src.giveRowSize();

#ifndef NDEBUG
    if ( size != loc.size() ) {
        OOFEM_ERROR("dimensions of 'src' and 'loc' mismatch");
    }

    if ( !src.isSquare() ) {
        OOFEM_ERROR("'src' is not sqaure matrix");
    }
#endif

    for ( std::size_t i = 1; i <= size; i++ ) {
        if ( ( ii = loc.at(i) ) ) {
            for ( std::size_t j = 1; j <= size; j++ ) {
                if ( ( jj = loc.at(j) ) ) {
                    dst.at(ii, jj) += src.at(i, j);
                }
            }
        }
    }
};

void assemble(FloatMatrix& dst, const FloatMatrix &src, const IntArray &rowind, const IntArray &colind)
{
    int ii, jj;
    int nr = src.giveNumberOfRows();
    int nc = src.giveNumberOfColumns();

#ifndef NDEBUG
    if ( nr != rowind.giveSize() ) {
        OOFEM_ERROR("row dimensions of 'src' and 'rowind' mismatch");
    }

    if ( nc != colind.giveSize() ) {
        OOFEM_ERROR("column dimensions of 'src' and 'colind' mismatch");
    }
#endif

    for ( int i = 1; i <= nr; i++ ) {
        if ( ( ii = rowind.at(i) ) ) {
            for ( int j = 1; j <= nc; j++ ) {
                if ( ( jj = colind.at(j) ) ) {
                    dst.at(ii, jj) += src.at(i, j);
                }
            }
        }
    }
}


void assembleT(FloatMatrix& dst, const FloatMatrix &src, const IntArray &rowind, const IntArray &colind)
{
    int ii, jj;
    int nr = src.giveNumberOfRows();
    int nc = src.giveNumberOfColumns();

#ifndef NDEBUG
    if ( nr != rowind.giveSize() ) {
        OOFEM_ERROR("row dimensions of 'src' and 'rowind' mismatch");
    }

    if ( nc != colind.giveSize() ) {
        OOFEM_ERROR("column dimensions of 'src' and 'colind' mismatch");
    }
#endif

    for ( int i = 1; i <= nr; i++ ) {
        if ( ( ii = rowind.at(i) ) ) {
            for ( int j = 1; j <= nc; j++ ) {
                if ( ( jj = colind.at(j) ) ) {
                    dst.at(jj, ii) += src.at(i, j);
                }
            }
        }
    }
}


void assemble(FloatMatrix& dst, const FloatMatrix &src, const int *rowind, const int *colind)
{
    int ii, jj;
    int nr = src.giveNumberOfRows();
    int nc = src.giveNumberOfColumns();

    for ( int i = 1; i <= nr; i++ ) {
        if ( ( ii = rowind [ i - 1 ] ) ) {
            for ( int j = 1; j <= nc; j++ ) {
                if ( ( jj = colind [ j - 1 ] ) ) {
                    dst.at(ii, jj) += src.at(i, j);
                }
            }
        }
    }
}



double computeNorm(const FloatMatrix& m, char p)
{
#  ifdef __LAPACK_MODULE
    FloatArray work( m.giveNumberOfRows() );
    int lda = max(m.nRows, 1);
    double norm = dlange_(& p, & m.nRows, & m.nColumns, m.givePointer(), & lda, work.givePointer(), 1);
    return norm;

#  else
    if ( p == '1' ) { // Maximum absolute column sum.
        double col_sum, max_col = 0.0;
        for (std::size_t j = 1; j <= (size_t) m.giveNumberOfColumns(); j++ ) {
            col_sum  = 0.0;
            for (std::size_t i = 1; i <= (size_t) m.giveNumberOfRows(); i++ ) {
                col_sum += std::fabs( m.at(i, j) );
            }
            if ( col_sum > max_col ) {
                max_col = col_sum;
            }
        }
        return max_col;
    }
    ///@todo Use m when obtaining eigen values is implemented.
    /*else if (p == '2') {
     *  double lambda_max;
     *  FloatMatrix AtA;
     *  FloatArray eigs;
     *  AtA.beTProductOf(*m,m);
     *  Ata.eigenValues(eigs, 1);
     *  return sqrt(eigs(0));
     * } */else {
        OOFEM_ERROR("p == %d not implemented.\n", p);
        // return 0.0;
    }
#  endif
}



double computeReciprocalCondition(const FloatMatrix& m, char p)
{
#  ifndef NDEBUG
    if ( !m.isSquare() ) {
        OOFEM_ERROR("receiver must be square (is %d by %d)", m.giveNRows(), m.giveNColumns());
    }
#  endif
    double anorm = ::oofem::mat::computeNorm(m,p);

#  ifdef __LAPACK_MODULE
    int n = m.giveNumberOfRows();
    FloatArray work(4 *n);
    IntArray iwork(n);
    int info;
    double rcond;

    if ( n > 3 ) { // Only use m routine for larger matrices. (not sure if it's suitable for n = 3) // Mikael
        FloatMatrix a_cpy = m;
        dgetrf_(& n, & n, a_cpy.givePointer(), & n, iwork.givePointer(), & info);
        if ( info < 0 ) {
            OOFEM_ERROR("dgetfr error %d\n", info);
        }
        dgecon_(& p, & ( m.nRows ), a_cpy.givePointer(), & m.nRows, & anorm, & rcond, work.givePointer(), iwork.givePointer(), & info, 1);
        if ( info < 0 ) {
            OOFEM_ERROR("dgecon error %d\n", info);
        }
        return rcond;
    }
#  endif
    if ( m.giveDeterminant() <= 1e-6 * anorm ) {
        return 0.0;
    }
    FloatMatrix inv;
    inv=beInverseOf(m);
    return 1.0 / ( ::oofem::mat::computeNorm(inv, p) * anorm );
}




FloatMatrix beLocalCoordSys(const FloatArray &normal)
{ //normal should be at the first position, easier for interface material models
    FloatMatrix ret;
    if ( normal.giveSize() == 1 ) {
        ret.resize(1, 1);
        ret.at(1, 1) = normal(0);
    } else if ( normal.giveSize() == 2 ) {
        ret.resize(2, 2);
        ret.at(1, 1) = normal(0);
        ret.at(1, 2) = normal(1);
        
        ret.at(2, 1) = normal(1);
        ret.at(2, 2) = -normal(0);
        
    } else if ( normal.giveSize() == 3 ) {
        // Create a permutated vector of n, *always* length 1 and significantly different from n.
        FloatArray b, t = Vec3(
            normal(1), -normal(2), normal(0)
        );                                                    // binormal and tangent

        // Construct orthogonal vector
        double npn = t.dotProduct(normal);
        t.add(-npn, normal);
        t.normalize();
        b.beVectorProductOf(t, normal);

        ret.resize(3, 3);
        ret.at(1, 1) = normal.at(1);
        ret.at(1, 2) = normal.at(2);
        ret.at(1, 3) = normal.at(3);
        
        ret.at(2, 1) = b.at(1);
        ret.at(2, 2) = b.at(2);
        ret.at(2, 3) = b.at(3);
        
        ret.at(3, 1) = t.at(1);
        ret.at(3, 2) = t.at(2);
        ret.at(3, 3) = t.at(3);
    } else {
        OOFEM_ERROR("Normal needs 1 to 3 components.");
    }
    return ret;
}



FloatMatrix beInverseOf(const FloatMatrix& src){ FloatMatrix inv; beInverseOf(inv,src); return inv; }


bool beInverseOf(FloatMatrix& inv, const FloatMatrix &src)
// Receiver becomes inverse of given parameter src. If necessary, size is adjusted.
{
    double det;

#  ifndef NDEBUG
    if ( !src.isSquare() ) {
        OOFEM_ERROR("cannot inverse a %d by %d matrix", src.giveNRows(), src.giveNColumns());
    }
#  endif

    inv.resize_funny(src.giveNumberOfRows(), src.giveNumberOfColumns());
    size_t nRows = (size_t) src.giveNumberOfRows();

    if ( nRows == 1 ) {
        if (fabs(src.at(1,1)) > 1.e-30) {
            inv.at(1, 1) = 1. / src.at(1, 1);
            return true;
        } else {
            return false;
        }
    } else if ( nRows == 2 ) {
        det = src.at(1, 1) * src.at(2, 2) - src.at(1, 2) * src.at(2, 1);
        if (fabs(det)>1.e-30) {
            inv.at(1, 1) =  src.at(2, 2) / det;
            inv.at(2, 1) = -src.at(2, 1) / det;
            inv.at(1, 2) = -src.at(1, 2) / det;
            inv.at(2, 2) =  src.at(1, 1) / det;
            return true;
        } else {
            return false;
        }
    } else if ( nRows == 3 ) {
        det = src.at(1, 1) * src.at(2, 2) * src.at(3, 3) + src.at(1, 2) * src.at(2, 3) * src.at(3, 1) +
              src.at(1, 3) * src.at(2, 1) * src.at(3, 2) - src.at(1, 3) * src.at(2, 2) * src.at(3, 1) -
              src.at(2, 3) * src.at(3, 2) * src.at(1, 1) - src.at(3, 3) * src.at(1, 2) * src.at(2, 1);
        if (fabs(det)>1.e-30) {
            inv.at(1, 1) = ( src.at(2, 2) * src.at(3, 3) - src.at(2, 3) * src.at(3, 2) ) / det;
            inv.at(2, 1) = ( src.at(2, 3) * src.at(3, 1) - src.at(2, 1) * src.at(3, 3) ) / det;
            inv.at(3, 1) = ( src.at(2, 1) * src.at(3, 2) - src.at(2, 2) * src.at(3, 1) ) / det;
            inv.at(1, 2) = ( src.at(1, 3) * src.at(3, 2) - src.at(1, 2) * src.at(3, 3) ) / det;
            inv.at(2, 2) = ( src.at(1, 1) * src.at(3, 3) - src.at(1, 3) * src.at(3, 1) ) / det;
            inv.at(3, 2) = ( src.at(1, 2) * src.at(3, 1) - src.at(1, 1) * src.at(3, 2) ) / det;
            inv.at(1, 3) = ( src.at(1, 2) * src.at(2, 3) - src.at(1, 3) * src.at(2, 2) ) / det;
            inv.at(2, 3) = ( src.at(1, 3) * src.at(2, 1) - src.at(1, 1) * src.at(2, 3) ) / det;
            inv.at(3, 3) = ( src.at(1, 1) * src.at(2, 2) - src.at(1, 2) * src.at(2, 1) ) / det;
            return true;
        } else {
            return false;
        }
    } else {
#ifdef __LAPACK_MODULE
        int n = inv.nRows;
        IntArray ipiv(n);
        int lwork, info;
        * inv = src;

        // LU-factorization
        dgetrf_(& n, & n, inv.givePointer(), & n, ipiv.givePointer(), & info);
        if ( info != 0 ) {
            OOFEM_WARNING("dgetrf error %d", info);
            return false;
        }

        // Inverse
        lwork = n * n;
        FloatArray work(lwork);
        dgetri_(& inv.nRows, inv.givePointer(), & inv.nRows, ipiv.givePointer(), work.givePointer(), & lwork, & info);
        if ( info > 0 ) {
            OOFEM_WARNING("Singular at %d", info);
            return false;
        } else if ( info < 0 ) {
            OOFEM_ERROR("Error on input %d", info);
        }
#else
        // size >3 ... gaussian elimination - slow but safe
        //
        double piv, linkomb;
        FloatMatrix tmp = src;
        // initialize answer to be unity matrix;
        inv.zero();
        for (std::size_t i = 1; i <= nRows; i++ ) {
            inv.at(i, i) = 1.0;
        }

        // lower triangle elimination by columns
        for (std::size_t i = 1; i < nRows; i++ ) {
            piv = tmp.at(i, i);
            if ( fabs(piv) < 1.e-30 ) {
                OOFEM_WARNING("pivot (%d,%d) to close to small (< 1.e-20)", i, i);
                return false;
            }

            for (std::size_t j = i + 1; j <= nRows; j++ ) {
                linkomb = tmp.at(j, i) / tmp.at(i, i);
                for (std::size_t k = i; k <= nRows; k++ ) {
                    tmp.at(j, k) -= tmp.at(i, k) * linkomb;
                }

                for (std::size_t k = 1; k <= nRows; k++ ) {
                    inv.at(j, k) -= inv.at(i, k) * linkomb;
                }
            }
        }

        // upper triangle elimination by columns
        for ( std::size_t i = nRows; i > 1; i-- ) {
            piv = tmp.at(i, i);
            for ( std::size_t j = i - 1; j > 0; j-- ) {
                linkomb = tmp.at(j, i) / piv;
                for ( std::size_t k = i; k > 0; k-- ) {
                    tmp.at(j, k) -= tmp.at(i, k) * linkomb;
                }

                for ( std::size_t k = nRows; k > 0; k-- ) {
                    // tmp -> at(j,k)-= tmp  ->at(i,k)*linkomb;
                    inv.at(j, k) -= inv.at(i, k) * linkomb;
                }
            }
        }

        // diagonal scaling
        for ( std::size_t i = 1; i <= nRows; i++ ) {
            for ( std::size_t j = 1; j <= nRows; j++ ) {
                inv.at(i, j) /= tmp.at(i, i);
            }
        }
        return true;
#endif
    }
}



FloatMatrix beMatrixFormOfStress(const FloatArray &aArray)
{
    FloatMatrix ret;
    // Revrites the  matrix on vector form (symmetrized matrix used), order: 11, 22, 33, 23, 13, 12
#  ifndef NDEBUG
    if ( aArray.giveSize() != 6 && aArray.giveSize() != 9 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }
#  endif
    ret.resize(3, 3);
    if ( aArray.giveSize() == 9 ) {
        ret.at(1, 1) = aArray.at(1);
        ret.at(2, 2) = aArray.at(2);
        ret.at(3, 3) = aArray.at(3);
        ret.at(2, 3) = aArray.at(4);
        ret.at(1, 3) = aArray.at(5);
        ret.at(1, 2) = aArray.at(6);
        ret.at(3, 2) = aArray.at(7);
        ret.at(3, 1) = aArray.at(8);
        ret.at(2, 1) = aArray.at(9);
    } else if ( aArray.giveSize() == 6 ) {
        ret.at(1, 1) = aArray.at(1);
        ret.at(2, 2) = aArray.at(2);
        ret.at(3, 3) = aArray.at(3);
        ret.at(2, 3) = aArray.at(4);
        ret.at(1, 3) = aArray.at(5);
        ret.at(1, 2) = aArray.at(6);
        ret.at(3, 2) = aArray.at(4);
        ret.at(3, 1) = aArray.at(5);
        ret.at(2, 1) = aArray.at(6);
    }
    return ret;
}

bool solveForRhs(FloatMatrix& m, const FloatArray &b, FloatArray &answer, bool transpose)
// solves equation b = this * x
{
#  ifndef NDEBUG
    if ( !m.isSquare() ) {
        OOFEM_ERROR("cannot solve a %d by %d matrix", m.giveNRows(), m.giveNColumns());
    }

    if ( m.giveNRows() != b.size() ) {
        OOFEM_ERROR("dimension mismatch");
    }
#  endif

#ifdef __LAPACK_MODULE
    int info, nrhs = 1;
    IntArray ipiv(m.nRows);
    answer = b;
    //dgesv_( &m.nRows, &nrhs, m.givePointer(), &m.nRows, ipiv.givePointer(), answer.givePointer(), &m.nRows, &info );
    dgetrf_(& m.nRows, & m.nRows, m.givePointer(), & m.nRows, ipiv.givePointer(), & info);
    if ( info == 0 ) {
        dgetrs_(transpose ? "Transpose" : "No transpose", & m.nRows, & nrhs, m.givePointer(), & m.nRows, ipiv.givePointer(), answer.givePointer(), & m.nRows, & info);
    }
    if ( info != 0 ) {
        return false;
    }
#else
    std::size_t pivRow;
    double piv, linkomb, help;
    FloatMatrix *mtrx, trans;
    if ( transpose ) {
        trans.beTranspositionOf(m);
        mtrx = & trans;
    } else {
        mtrx = &m;
    }

    size_t nRows = mtrx->giveNumberOfRows();

    answer = b;

    // initialize answer to be unity matrix;
    // lower triangle elimination by columns
    for ( std::size_t i = 1; i < nRows; i++ ) {
        // find the suitable row and pivot
        piv = fabs( mtrx->at(i, i) );
        pivRow = i;
        for (std::size_t j = i + 1; j <= nRows; j++ ) {
            if ( fabs( mtrx->at(j, i) ) > piv ) {
                pivRow = j;
                piv = fabs( mtrx->at(j, i) );
            }
        }

        if ( piv < 1.e-20 ) {
            return false;
        }

        // exchange rows
        if ( pivRow != i ) {
            for (std::size_t j = i; j <= nRows; j++ ) {
                help = mtrx->at(i, j);
                mtrx->at(i, j) = mtrx->at(pivRow, j);
                mtrx->at(pivRow, j) = help;
            }
            help = answer.at(i);
            answer.at(i) = answer.at(pivRow);
            answer.at(pivRow) = help;
        }

        for (std::size_t j = i + 1; j <= nRows; j++ ) {
            linkomb = mtrx->at(j, i) / mtrx->at(i, i);
            for (std::size_t k = i; k <= nRows; k++ ) {
                mtrx->at(j, k) -= mtrx->at(i, k) * linkomb;
            }

            answer.at(j) -= answer.at(i) * linkomb;
        }
    }

    // back substitution
    for ( std::size_t i = nRows; i >= 1; i-- ) {
        help = 0.;
        for ( std::size_t j = i + 1; j <= nRows; j++ ) {
            help += mtrx->at(i, j) * answer.at(j);
        }

        answer.at(i) = ( answer.at(i) - help ) / mtrx->at(i, i);
    }
#endif
    return true;
}

bool solveForRhs(FloatMatrix& m, const FloatMatrix &b, FloatMatrix &answer, bool transpose)
// solves equation b = m * x
// returns x. m and b are kept untouched
//
// gaussian elimination - slow but safe
//
{
#  ifndef NDEBUG
    if ( !m.isSquare() ) {
        OOFEM_ERROR("cannot solve a %d by %d matrix", m.giveNRows(), m.giveNColumns());
    }

    if ( m.giveNRows() != b.giveRowSize() ) {
        OOFEM_ERROR("dimension mismatch");
    }
#  endif

#ifdef __LAPACK_MODULE
    int info;
    IntArray ipiv(m.nRows);
    answer = b;
    dgetrf_(& m.nRows, & m.nRows, m.givePointer(), & m.nRows, ipiv.givePointer(), & info);
    if ( info == 0 ) {
        dgetrs_(transpose ? "t" : "n", & m.nRows, & answer.nColumns, m.givePointer(), & m.nRows, ipiv.givePointer(), answer.givePointer(), & m.nRows, & info);
    }
    if ( info != 0 ) {
        return false; //OOFEM_ERROR("error %d", info);
    }
#else
    std::size_t pivRow, nPs;
    double piv, linkomb, help;
    FloatMatrix *mtrx, trans;
    if ( transpose ) {
        trans.beTranspositionOf(m);
        mtrx = & trans;
    } else {
        mtrx = & m;
    }

    size_t nRows = mtrx->giveNumberOfRows();

    nPs = b.giveNumberOfColumns();
    answer = b;
    // initialize answer to be unity matrix;
    // lower triangle elimination by columns
    for (std::size_t i = 1; i < nRows; i++ ) {
        // find the suitable row and pivot
        piv = fabs( mtrx->at(i, i) );
        pivRow = i;
        for (std::size_t j = i + 1; j <= nRows; j++ ) {
            if ( fabs( mtrx->at(j, i) ) > piv ) {
                pivRow = j;
                piv = fabs( mtrx->at(j, i) );
            }
        }

        if ( fabs(piv) < 1.e-20 ) {
            return false; //OOFEM_ERROR("pivot too small, cannot solve %d by %d matrix", nRows, nColumns);
        }

        // exchange rows
        if ( pivRow != i ) {
            for (std::size_t j = i; j <= nRows; j++ ) {
                help = mtrx->at(i, j);
                mtrx->at(i, j) = mtrx->at(pivRow, j);
                mtrx->at(pivRow, j) = help;
            }

            for (std::size_t j = 1; j <= nPs; j++ ) {
                help = answer.at(i, j);
                answer.at(i, j) = answer.at(pivRow, j);
                answer.at(pivRow, j) = help;
            }
        }

        for (std::size_t j = i + 1; j <= nRows; j++ ) {
            linkomb = mtrx->at(j, i) / mtrx->at(i, i);
            for (std::size_t k = i; k <= nRows; k++ ) {
                mtrx->at(j, k) -= mtrx->at(i, k) * linkomb;
            }

            for (std::size_t k = 1; k <= nPs; k++ ) {
                answer.at(j, k) -= answer.at(i, k) * linkomb;
            }
        }
    }

    // back substitution
    for ( std::size_t i = nRows; i >= 1; i-- ) {
        for (std::size_t k = 1; k <= nPs; k++ ) {
            help = 0.;
            for ( std::size_t j = i + 1; j <= nRows; j++ ) {
                help += mtrx->at(i, j) * answer.at(j, k);
            }

            answer.at(i, k) = ( answer.at(i, k) - help ) / mtrx->at(i, i);
        }
    }
    return true;
#endif
}


void plusProductSymmUpper(FloatMatrix& recv, const FloatMatrix &a, const FloatMatrix &b, double dV)
// Adds to the receiver the product  a(transposed).b dV .
// The receiver size is adjusted, if necessary.
// This method assumes that both the receiver and the product above are
// symmetric matrices, and therefore computes only the upper half of the
// receiver ; the lower half is not modified. Other advantage : it does
// not compute the transposition of matrix a.
{
    if ( !recv.isNotEmpty() ) {
        recv.resize(a.giveNColumns(),b.giveNColumns());
    }
    size_t recv_nRows=recv.giveNRows();
    size_t recv_nColumns=recv.giveNColumns();

#ifdef __LAPACK_MODULE
    double beta = 1.;
    ///@todo We should determine which is the best choice overall. For large systems more block matrix operations is necessary.
    /// For smaller systems the overhead from function calls might be larger, but the overhead might be tiny, or using symmetry at all might be undesireable.
    if ( recv_nRows < 20 ) {
        // Split the matrix into 2 columns, s1 + s2 = n ( = nRows = nColumns ).
        int s1 = recv_nRows / 2;
        int s2 = recv_nRows - s1;
        // First column block, we only take the first s rows by only taking the first s columns in the matrix a.
        dgemm_("t", "n", & s1, & s1, & a.nRows, & dV, a.givePointer(), & a.nRows, b.givePointer(), & b.nRows, & beta, recv.givePointer(), & recv_nRows, a.nColumns, b.nColumns, recv_nColumns);
        // Second column block starting a memory position c * nRows
        dgemm_("t", "n", & recv_nRows, & s2, & a.nRows, & dV, a.givePointer(), & a.nRows, & b.givePointer() [ s1 * b.nRows ], & b.nRows, & beta, & recv.givePointer() [ s1 * recv_nRows ], & recv_nRows, a.nColumns, b.nColumns, recv_nColumns);
    } else {
        // Get suitable blocksize. Around 10 rows should be suitable (slightly adjusted to minimize number of blocks):
        // Smaller blocks than ~10 didn't show any performance gains in my benchmarks. / Mikael
        int block = ( recv_nRows - 1 ) / ( recv_nRows / 10 ) + 1;
        int start = 0;
        int end = block;
        while ( start < recv_nRows ) {
            int s = end - start;
            dgemm_("t", "n", & end, & s, & a.nRows, & dV, a.givePointer(), & a.nRows, & b.givePointer() [ start * b.nRows ], & b.nRows, & beta, & recv.givePointer() [ start * recv_nRows ],
                   & recv_nRows, a.nColumns, b.nColumns, recv_nColumns);
            start = end;
            end += block;
            if ( end > recv_nRows ) {
                end = recv_nRows;
            }
        }
    }
#else
    for (std::size_t i = 1; i <= recv_nRows; i++ ) {
        for (std::size_t j = i; j <= recv_nColumns; j++ ) {
            double summ = 0.;
            for (std::size_t k = 1; k <= a.giveNRows(); k++ ) {
                summ += a.at(k, i) * b.at(k, j);
            }

            recv.at(i, j) += summ * dV;
        }
    }
#endif
}


void plusDyadSymmUpper(FloatMatrix& recv, const FloatArray &a, double dV)
{
    if ( !recv.isNotEmpty() ) {
        recv.resize(a.giveSize(),a.giveSize());
    }
    size_t recv_nRows=recv.giveNumberOfRows();
    size_t recv_nColumns=recv.giveNumberOfColumns();
#ifdef __LAPACK_MODULE
    int inc = 1;
    int sizeA = a.giveSize();
    dsyr_("u", & sizeA, & dV, a.givePointer(), & inc,
          recv.givePointer(), & recv_nRows,
          sizeA, recv_nColumns);
#else
    for (std::size_t i = 1; i <= recv_nRows; i++ ) {
        for (std::size_t j = i; j <= recv_nColumns; j++ ) {
            recv.at(i, j) += a.at(i) * a.at(j) * dV;
        }
    }
#endif
}



void plusProductUnsym(FloatMatrix& recv, const FloatMatrix &a, const FloatMatrix &b, double dV)
// Adds to the receiver the product  a(transposed).b dV .
// If the receiver has a null size, it is expanded.
// Advantage : does not compute the transposition of matrix a.
{
    if ( !recv.isNotEmpty() ) {
        recv.resize(a.giveNumberOfColumns(),b.giveNumberOfColumns());
    }
    size_t recv_nRows=recv.giveNRows();
    size_t recv_nColumns=recv.giveNColumns();
#ifdef __LAPACK_MODULE
    double beta = 1.;
    dgemm_("t", "n", & recv_nRows, & recv_nColumns, & a.nRows,
           & dV, a.givePointer(), & a.nRows, b.givePointer(), & b.nRows,
           & beta, recv.givePointer(), & recv_nRows,
           a.nColumns, b.nColumns, recv_nColumns);
#else
    for (std::size_t i = 1; i <= recv_nRows; i++ ) {
        for (std::size_t j = 1; j <= recv_nColumns; j++ ) {
            double summ = 0.;
            for (std::size_t k = 1; k <= a.giveNRows(); k++ ) {
                summ += a.at(k, i) * b.at(k, j);
            }

            recv.at(i, j) += summ * dV;
        }
    }
#endif
}


void plusDyadUnsym(FloatMatrix& recv, const FloatArray &a, const FloatArray &b, double dV)
{
    if ( !recv.isNotEmpty() ) {
        recv.resize(a.giveSize(),b.giveSize());
    }
    size_t recv_nRows=recv.giveNRows();
    size_t recv_nColumns=recv.giveNColumns();
#ifdef __LAPACK_MODULE
    int inc = 1;
    int sizeA = a.giveSize();
    int sizeB = b.giveSize();
    dger_(& sizeA, & sizeB, & dV, a.givePointer(), & inc,
          b.givePointer(), & inc, recv.givePointer(), & recv_nRows,
          sizeA, sizeB, recv_nColumns);
#else
    for (std::size_t i = 1; i <= recv_nRows; i++ ) {
        for (std::size_t j = 1; j <= recv_nColumns; j++ ) {
            recv.at(i, j) += a.at(i) * b.at(j) * dV;
        }
    }
#endif
}


FloatMatrix bePinvID()
// this matrix is the product of the 6x6 deviatoric projection matrix ID
// and the inverse scaling matrix Pinv
{
    FloatMatrix ret;
    ret.resize(6, 6);
    ret(0,0)=ret(1,1)=ret(2,2)=2./3.;
    ret(1,0)=ret(2,0)=ret(0,1)=ret(0,2)=ret(1,2)=ret(2,1)=-1./3.;
    ret(3,3)=ret(4,4)=ret(5,5)=1/2.;
    #if 0
        values [ 0 ] = values [ 7 ] = values [ 14 ] = 2. / 3.;
        values [ 1 ] = values [ 2 ] = values [ 6 ] = values [ 8 ] = values [ 12 ] = values [ 13 ] = -1. / 3.;
        values [ 21 ] = values [ 28 ] = values [ 35 ] = 0.5;
    #endif
    return ret;
}



FloatMatrix beMatrixForm(const FloatArray &aArray)
{
    // Revrites the vector on matrix form (symmetrized matrix used if size is 6),
    // order: 11, 22, 33, 23, 13, 12
    // order: 11, 22, 33, 23, 13, 12, 32, 31, 21
    FloatMatrix ret;
#  ifndef NDEBUG
    if ( aArray.giveSize() != 6 && aArray.giveSize() != 9 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }
#  endif
    ret.resize(3, 3);
    if ( aArray.giveSize() == 9 ) {
        ret.at(1, 1) = aArray.at(1);
        ret.at(2, 2) = aArray.at(2);
        ret.at(3, 3) = aArray.at(3);
        ret.at(2, 3) = aArray.at(4);
        ret.at(1, 3) = aArray.at(5);
        ret.at(1, 2) = aArray.at(6);
        ret.at(3, 2) = aArray.at(7);
        ret.at(3, 1) = aArray.at(8);
        ret.at(2, 1) = aArray.at(9);
    } else if ( aArray.giveSize() == 6 ) {
        ret.at(1, 1) = aArray.at(1);
        ret.at(2, 2) = aArray.at(2);
        ret.at(3, 3) = aArray.at(3);
        ret.at(2, 3) = aArray.at(4);
        ret.at(1, 3) = aArray.at(5);
        ret.at(1, 2) = aArray.at(6);
        ret.at(3, 2) = aArray.at(4);
        ret.at(3, 1) = aArray.at(5);
        ret.at(2, 1) = aArray.at(6);
    }
    return ret;
}

void changeComponentOrder(FloatMatrix& m)
{
    // Changes index order between abaqus <-> OOFEM
//#  ifndef NDEBUG
//    if ( nRows != 6 || nColumns != 6 ) {
//        OOFEM_ERROR("matrix dimension is not 6x6");
//    }
//#  endif
   int nRows = m.giveNumberOfRows();
   int nColumns = m.giveNumberOfColumns();

    if ( nRows == 6 && nColumns == 6 ) {
        // This could probably be done more beautifully + efficiently.

        std :: swap( m.at(4, 1), m.at(6, 1) );

        std :: swap( m.at(4, 2), m.at(6, 2) );

        std :: swap( m.at(4, 3), m.at(6, 3) );

        std :: swap( m.at(1, 4), m.at(1, 6) );
        std :: swap( m.at(2, 4), m.at(2, 6) );
        std :: swap( m.at(3, 4), m.at(3, 6) );
        std :: swap( m.at(4, 4), m.at(6, 6) );
        std :: swap( m.at(5, 4), m.at(5, 6) );
        std :: swap( m.at(6, 4), m.at(4, 6) );

        std :: swap( m.at(4, 5), m.at(6, 5) );
    } else if ( nRows == 9 && nColumns == 9 ) {
        // OOFEM:           11, 22, 33, 23, 13, 12, 32, 31, 21
        // UMAT:            11, 22, 33, 12, 13, 23, 32, 21, 31
        const int abq2oo [ 9 ] = {
            1,  2,  3,  6,  5,  4,  7,  9,  8
        };

        FloatMatrix tmp(9, 9);
        for ( int i = 1; i <= 9; i++ ) {
            for ( int j = 1; j <= 9; j++ ) {
                tmp.at(i, j) = m.at(abq2oo [ i - 1 ], abq2oo [ j - 1 ]);
            }
        }

        m = tmp;
    }
}


bool jaco_(FloatMatrix& M, FloatArray &eval, FloatMatrix &v, int nf)
{
    /*
     * Solves the eigenvalues and eigenvectors of real
     * symmetric matrix by jacobi method.
     *  Written by bp. Inspired by ED WILSON jaco_ procedure.
     *
     * Parameters (input):
     * nf - number of significant figures
     *
     * Output params:
     * eval - eigen values (not sorted)
     * v    - eigenvectors (stored columvise)
     */


    /* Local variables */
    int neq = M.giveNumberOfRows();

    double c_b2 = .10;
    //double c_b27 = .01;

    /* Function Body */
#ifndef NDEBUG
    if ( !M.isSquare() ) {
        OOFEM_ERROR("Not square matrix");
    }
    // check for symmetry
    for ( int i = 1; i <= neq; i++ ) {
        for ( int j = i + 1; j <= neq; j++ ) {
            //if ( M.at(i, j) != M.at(j, i) ) {
            if ( fabs( M.at(i, j) - M.at(j, i) ) > 1.0e-6 ) {
                OOFEM_ERROR("Not Symmetric matrix");
            }
        }
    }

#endif

    eval.resize(neq);
    v.resize(neq, neq);

    for ( int i = 1; i <= neq; i++ ) {
        eval.at(i) = M.at(i, i);
    }

    double tol = pow(c_b2, nf);
    double sum = 0.0;
    for ( int i = 1; i <= neq; ++i ) {
        for ( int j = 1; j <= neq; ++j ) {
            sum += fabs( M.at(i, j) );
            v.at(i, j) = 0.0;
        }

        v.at(i, i) = 1.0;
    }

    if ( sum <= 0.0 ) {
        return 0;
    }


    /* ---- REDUCE MATRIX TO DIAGONAL ---------------- */
    int ite = 0;
    double ssum;
    do {
        ssum = 0.0;
        for ( int j = 2; j <= neq; ++j ) {
            int ih = j - 1;
            for ( int i = 1; i <= ih; ++i ) {
                if ( ( fabs( M.at(i, j) ) / sum ) > tol ) {
                    ssum += fabs( M.at(i, j) );
                    /* ---- CALCULATE ROTATION ANGLE ----------------- */
                    double aa = atan2( M.at(i, j) * 2.0, eval.at(i) - eval.at(j) ) /  2.0;
                    double si = sin(aa);
                    double co = cos(aa);
                    /*
                     *   // ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
                     *   for (k = 1; k <= neq; ++k) {
                     *    tt = M.at(k, i);
                     *    M.at(k, i) = co * tt + si * M.at(k, j);
                     *    M.at(k, j) = -si * tt + co * M.at(k, j);
                     *    tt = v.at(k, i);
                     *    v.at(k, i) = co * tt + si * v.at(k, j);
                     *    // L500:
                     *    v.at(k, j) = -si * tt + co * v.at(k, j);
                     *   }
                     *   // ---- MODIFY DIAGONAL TERMS --------------------
                     *   M.at(i, i) = co * M.at(i, i) + si * M.at(j, i);
                     *   M.at(j, j) = -si * M.at(i, j) + co * M.at(j, j);
                     *   M.at(i, j) = 0.0;
                     *   // ---- MAKE "A" MATRIX SYMMETRICAL --------------
                     *   for (k = 1; k <= neq; ++k) {
                     *    M.at(i, k) = M.at(k, i);
                     *    M.at(j, k) = M.at(k, j);
                     *    // L600:
                     *   }
                     */
                    // ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
                    for ( int k = 1; k < i; ++k ) {
                        double tt = M.at(k, i);
                        M.at(k, i) = co * tt + si * M.at(k, j);
                        M.at(k, j) = -si * tt + co * M.at(k, j);
                        tt = v.at(k, i);
                        v.at(k, i) = co * tt + si *v.at(k, j);
                        v.at(k, j) = -si * tt + co *v.at(k, j);
                    }

                    // diagonal term (i,i)
                    double tt = eval.at(i);
                    eval.at(i) = co * tt + si * M.at(i, j);
                    double aij = -si * tt + co *M.at(i, j);
                    tt = v.at(i, i);
                    v.at(i, i) = co * tt + si *v.at(i, j);
                    v.at(i, j) = -si * tt + co *v.at(i, j);

                    for ( int k = i + 1; k < j; ++k ) {
                        double tt = M.at(i, k);
                        M.at(i, k) = co * tt + si * M.at(k, j);
                        M.at(k, j) = -si * tt + co * M.at(k, j);
                        tt = v.at(k, i);
                        v.at(k, i) = co * tt + si *v.at(k, j);
                        v.at(k, j) = -si * tt + co *v.at(k, j);
                    }

                    // diagonal term (j,j)
                    tt = M.at(i, j);
                    double aji = co * tt + si *eval.at(j);
                    eval.at(j) = -si * tt + co *eval.at(j);

                    tt = v.at(j, i);
                    v.at(j, i) = co * tt + si *v.at(j, j);
                    v.at(j, j) = -si * tt + co *v.at(j, j);
                    //
                    for ( int k = j + 1; k <= neq; ++k ) {
                        double tt = M.at(i, k);
                        M.at(i, k) = co * tt + si *M.at(j, k);
                        M.at(j, k) = -si * tt + co *M.at(j, k);
                        tt = v.at(k, i);
                        v.at(k, i) = co * tt + si *v.at(k, j);
                        v.at(k, j) = -si * tt + co *v.at(k, j);
                    }

                    // ---- MODIFY DIAGONAL TERMS --------------------
                    eval.at(i) = co * eval.at(i) + si * aji;
                    eval.at(j) = -si * aij + co *eval.at(j);
                    M.at(i, j) = 0.0;
                } else {
                    /* ---- A(I,J) MADE ZERO BY ROTATION ------------- */
                    ;
                }
            }
        }

        /* ---- CHECK FOR CONVERGENCE -------------------- */
        if ( ++ite > 50 ) {
            OOFEM_ERROR("too many iterations");
        }
    } while ( fabs(ssum) / sum > tol );

    // restore original matrix
    for ( int i = 1; i <= neq; i++ ) {
        for ( int j = i; j <= neq; j++ ) {
            M.at(i, j) = M.at(j, i);
        }
    }

    return 0;
} /* jaco_ */



void rotatedWith(FloatMatrix& a, const FloatMatrix &r, char mode)
// Returns the receiver 'a' rotated according the change-of-base matrix r.
// The method performs the operation  a = r^T . a . r . or the inverse
{
    FloatMatrix rta;

    if ( mode == 'n' ) {
        rta.beTProductOf(r, a);     //  r^T . a
        a.beProductOf(rta, r);      //  r^T . a . r
    } else if ( mode == 't' ) {
        rta.beProductOf(r, a);      //  r . a
        a.beProductTOf(rta, r);     //  r . a . r^T
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}




void printYourself(const FloatMatrix& m)
// Prints the receiver on screen.
{
    size_t nRows=m.giveNRows(), nColumns=m.giveNColumns();
    printf("FloatMatrix with dimensions : %zu %zu\n",
           nRows, nColumns);
    if ( nRows <= 250 && nColumns <= 250 ) {
        for (std::size_t i = 1; i <= nRows; ++i ) {
            for (std::size_t j = 1; j <= nColumns && j <= 100; ++j ) {
                printf( "%10.3e  ", m.at(i, j) );
            }

            printf("\n");
        }
    } else {
        printf("   large matrix : coefficients not printed \n");
    }
}


void printYourselfToFile(const FloatMatrix& m, const std::string filename, const bool showDimensions)
// Prints the receiver to file.
{
    size_t nRows=m.giveNRows(), nColumns=m.giveNColumns();
    std :: ofstream matrixfile (filename);
    if (matrixfile.is_open()) {
        if (showDimensions)
            matrixfile << "FloatMatrix with dimensions : " << nRows << ", " << nColumns << "\n";
        matrixfile << std::scientific << std::right << std::setprecision(3);
        for (std::size_t i = 1; i <= nRows; ++i ) {
            for (std::size_t j = 1; j <= nColumns; ++j ) {
                matrixfile << std::setw(10) << m.at(i, j) << "\t";
            }

            matrixfile << "\n";
        }
        matrixfile.close();
    } else {
        OOFEM_ERROR("Failed to write to file");
    }
}


void printYourself(const FloatMatrix& m, const std::string &name)
// Prints the receiver on screen.
{
    size_t nRows=m.giveNRows(), nColumns=m.giveNColumns();
    printf("%s (%zu x %zu): \n", name.c_str(), nRows, nColumns);
    if ( nRows <= 250 && nColumns <= 250 ) {
        for (std::size_t i = 1; i <= nRows; ++i ) {
            for (std::size_t j = 1; j <= nColumns && j <= 100; ++j ) {
                printf( "%10.3e  ", m.at(i, j) );
            }

            printf("\n");
        }
    } else {
        for (std::size_t i = 1; i <= nRows && i <= 20; ++i ) {
            for (std::size_t j = 1; j <= nColumns && j <= 10; ++j ) {
                printf( "%10.3e  ", m.at(i, j) );
            }
            if ( nColumns > 10 ) printf(" ...");
            printf("\n");
        }
        if ( nRows > 20 )  printf(" ...\n");
    }
}


void pY(const FloatMatrix& m)
// Prints the receiver on screen with higher accuracy than printYourself.
{
    size_t nRows=m.giveNRows(), nColumns=m.giveNColumns();
    printf("[");
    for (std::size_t i = 1; i <= nRows; ++i ) {
        for (std::size_t j = 1; j <= nColumns; ++j ) {
            printf( "%20.15e", m.at(i, j) );
            if ( j < nColumns ) {
                printf(",");
            } else {
                printf(";");
            }
        }
    }

    printf("];\n");
}


void writeCSV(const FloatMatrix& m, const std :: string &name)
{
    size_t nRows=m.giveNRows(), nColumns=m.giveNColumns();
    FILE *file = fopen(name.c_str(), "w");
    for (std::size_t i = 1; i <= nRows; ++i ) {
        for (std::size_t j = 1; j <= nColumns; ++j ) {
            fprintf(file, "%10.3e, ", m.at(i, j) );
        }

        fprintf(file, "\n");
    }
    fclose(file);
}

} /* namespace oofem::mat */

