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
/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "error.h"
#include "datastream.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ostream>
#include <numeric>
#define RESIZE(nr, nc) \
    { \
        this->nRows = nr; this->nColumns = nc; \
        int nsize = this->nRows * this->nColumns; \
        if ( nsize < ( int ) this->values.size() ) { \
            this->values.resize(nsize); \
        } else if ( nsize > ( int ) this->values.size() ) { \
            this->values.assign(nsize, 0.); \
        } \
    }

#ifdef BOOST_PYTHON
 #include <boost/python.hpp>
 #include <boost/python/extract.hpp>
#endif

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


namespace oofem {
FloatMatrix :: FloatMatrix(const FloatArray &vector, bool transpose)
//
// constructor : creates (vector->giveSize(),1) FloatMatrix
// if transpose = 1 creates (1,vector->giveSize()) FloatMatrix
//
{
    if ( transpose ) {
        nRows = 1; // column vector
        nColumns = vector.giveSize();
    } else {
        nRows = vector.giveSize(); // row vector- default
        nColumns = 1;
    }

    values = vector.values;
}


FloatMatrix :: FloatMatrix(std :: initializer_list< std :: initializer_list< double > >mat)
{
    RESIZE( mat.size(), mat.begin()->size() )
    auto p = this->values.begin();
    for ( auto col : mat ) {
#if DEBUG
        if ( this->nRows != ( int ) col.size() ) {
            OOFEM_ERROR("Initializer list has inconsistent column sizes.");
        }
#endif
        for ( auto x : col ) {
            * p = x;
            p++;
        }
    }
}


FloatMatrix &FloatMatrix :: operator = ( std :: initializer_list< std :: initializer_list< double > >mat )
{
    RESIZE( ( int ) mat.begin()->size(), ( int ) mat.size() );
    auto p = this->values.begin();
    for ( auto col : mat ) {
#if DEBUG
        if ( this->nRows != ( int ) col.size() ) {
                OOFEM_ERROR("Initializer list has inconsistent column sizes.");
        }
#endif
        for ( auto x : col ) {
            * p = x;
            p++;
        }
    }

    return * this;
}


FloatMatrix &FloatMatrix :: operator = ( std :: initializer_list< FloatArray >mat )
{
    RESIZE( mat.begin()->giveSize(), ( int ) mat.size() );
    auto p = this->values.begin();
    for ( auto col : mat ) {
#if DEBUG
        if ( this->nRows != col.giveSize() ) {
                OOFEM_ERROR("Initializer list has inconsistent column sizes.");
        }
#endif
        for ( auto x : col ) {
            * p = x;
            p++;
        }
    }

    return * this;
}


void FloatMatrix :: checkBounds(int i, int j) const
// Checks that the receiver includes a position (i,j).
{
    if ( i <= 0 ) {
        OOFEM_ERROR("matrix error on rows : %d < 0", i);
    }
    if ( j <= 0 ) {
        OOFEM_ERROR("matrix error on columns : %d < 0", j);
    }

    if ( i > nRows ) {
        OOFEM_ERROR("matrix error on rows : %d > %d", i, nRows);
    }

    if ( j > nColumns ) {
        OOFEM_ERROR("matrix error on columns : %d > %d", j, nColumns);
    }
}

bool FloatMatrix :: isFinite() const
{
    for (double val : values) {
        if ( !std::isfinite(val) ) {
            return false;
        }
    }

    return true;
}

#ifdef DEBUG
double &FloatMatrix :: at(int i, int j)
// Returns the coefficient (i,j) of the receiver. Safer but slower than
// the inline version of method 'at'.
{
    this->checkBounds(i, j);
    return values [ ( j - 1 ) * nRows + i - 1 ];
}

double FloatMatrix :: at(int i, int j) const
// Returns the coefficient (i,j) of the receiver. Safer but slower than
// the inline version of method 'at'.
{
    this->checkBounds(i, j);
    return values [ ( j - 1 ) * nRows + i - 1 ];
}

double &FloatMatrix :: operator() (int i, int j)
{
    this->checkBounds(i + 1, j + 1);
    return values [ j * nRows + i ];
}

double FloatMatrix :: operator() (int i, int j) const
{
    this->checkBounds(i + 1, j + 1);
    return values [ j * nRows + i ];
}
#endif


void FloatMatrix :: assemble(const FloatMatrix &src, const IntArray &loc)
{
    int ii, jj, size = src.giveNumberOfRows();

#ifdef DEBUG
    if ( size != loc.giveSize() ) {
        OOFEM_ERROR("dimensions of 'src' and 'loc' mismatch");
    }

    if ( !src.isSquare() ) {
        OOFEM_ERROR("'src' is not sqaure matrix");
    }
#endif

    for ( int i = 1; i <= size; i++ ) {
        if ( ( ii = loc.at(i) ) ) {
            for ( int j = 1; j <= size; j++ ) {
                if ( ( jj = loc.at(j) ) ) {
                    this->at(ii, jj) += src.at(i, j);
                }
            }
        }
    }
}


void FloatMatrix :: assemble(const FloatMatrix &src, const IntArray &rowind, const IntArray &colind)
{
    int ii, jj;
    int nr = src.giveNumberOfRows();
    int nc = src.giveNumberOfColumns();

#ifdef DEBUG
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
                    this->at(ii, jj) += src.at(i, j);
                }
            }
        }
    }
}


void FloatMatrix :: assemble(const FloatMatrix &src, const int *rowind, const int *colind)
{
    int ii, jj;
    int nr = src.giveNumberOfRows();
    int nc = src.giveNumberOfColumns();

    for ( int i = 1; i <= nr; i++ ) {
        if ( ( ii = rowind [ i - 1 ] ) ) {
            for ( int j = 1; j <= nc; j++ ) {
                if ( ( jj = colind [ j - 1 ] ) ) {
                    this->at(ii, jj) += src.at(i, j);
                }
            }
        }
    }
}


void FloatMatrix :: beTranspositionOf(const FloatMatrix &src)
{
    // receiver becomes a transposition of src
    int nrows = src.giveNumberOfColumns(), ncols = src.giveNumberOfRows();
    RESIZE(nrows, ncols);

    for ( int i = 1; i <= nrows; i++ ) {
        for ( int j = 1; j <= ncols; j++ ) {
            this->at(i, j) = src.at(j, i);
        }
    }
}


void FloatMatrix :: beProductOf(const FloatMatrix &aMatrix, const FloatMatrix &bMatrix)
// Receiver = aMatrix * bMatrix
{
#  ifdef DEBUG
    if ( aMatrix.nColumns != bMatrix.nRows ) {
        OOFEM_ERROR("error in product A*B : dimensions do not match");
    }
#  endif
    RESIZE(aMatrix.nRows, bMatrix.nColumns);
#  ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    dgemm_("n", "n", & this->nRows, & this->nColumns, & aMatrix.nColumns,
           & alpha, aMatrix.givePointer(), & aMatrix.nRows, bMatrix.givePointer(), & bMatrix.nRows,
           & beta, this->givePointer(), & this->nRows,
           aMatrix.nColumns, bMatrix.nColumns, this->nColumns);
#  else
    for ( int i = 1; i <= aMatrix.nRows; i++ ) {
        for ( int j = 1; j <= bMatrix.nColumns; j++ ) {
            double coeff = 0.;
            for ( int k = 1; k <= aMatrix.nColumns; k++ ) {
                coeff += aMatrix.at(i, k) * bMatrix.at(k, j);
            }

            this->at(i, j) = coeff;
        }
    }
#  endif
}


void FloatMatrix :: beTProductOf(const FloatMatrix &aMatrix, const FloatMatrix &bMatrix)
// Receiver = aMatrix^T * bMatrix
{
#  ifdef DEBUG
    if ( aMatrix.nRows != bMatrix.nRows ) {
        OOFEM_ERROR("error in product A*B : dimensions do not match");
    }
#  endif
    RESIZE(aMatrix.nColumns, bMatrix.nColumns);
#  ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    dgemm_("t", "n", & this->nRows, & this->nColumns, & aMatrix.nRows,
           & alpha, aMatrix.givePointer(), & aMatrix.nRows, bMatrix.givePointer(), & bMatrix.nRows,
           & beta, this->givePointer(), & this->nRows,
           aMatrix.nColumns, bMatrix.nColumns, this->nColumns);
#  else
    for ( int i = 1; i <= aMatrix.nColumns; i++ ) {
        for ( int j = 1; j <= bMatrix.nColumns; j++ ) {
            double coeff = 0.;
            for ( int k = 1; k <= aMatrix.nRows; k++ ) {
                coeff += aMatrix.at(k, i) * bMatrix.at(k, j);
            }

            this->at(i, j) = coeff;
        }
    }
#endif
}


void FloatMatrix :: beProductTOf(const FloatMatrix &aMatrix, const FloatMatrix &bMatrix)
// Receiver = aMatrix * bMatrix^T
{
#  ifdef DEBUG
    if ( aMatrix.nColumns != bMatrix.nColumns ) {
        OOFEM_ERROR("error in product A*B : dimensions do not match");
    }
#  endif
    RESIZE(aMatrix.nRows, bMatrix.nRows);
#  ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    dgemm_("n", "t", & this->nRows, & this->nColumns, & aMatrix.nColumns,
           & alpha, aMatrix.givePointer(), & aMatrix.nRows, bMatrix.givePointer(), & bMatrix.nRows,
           & beta, this->givePointer(), & this->nRows,
           aMatrix.nColumns, bMatrix.nColumns, this->nColumns);
#  else
    for ( int i = 1; i <= aMatrix.nRows; i++ ) {
        for ( int j = 1; j <= bMatrix.nRows; j++ ) {
            double coeff = 0.;
            for ( int k = 1; k <= aMatrix.nColumns; k++ ) {
                coeff += aMatrix.at(i, k) * bMatrix.at(j, k);
            }

            this->at(i, j) = coeff;
        }
    }
#  endif
}


void FloatMatrix :: addProductOf(const FloatMatrix &aMatrix, const FloatMatrix &bMatrix)
// Receiver = aMatrix * bMatrix
{
#  ifdef DEBUG
    if ( aMatrix.nColumns != bMatrix.nRows ) {
        OOFEM_ERROR("error in product A*B : dimensions do not match");
    }
    if ( aMatrix.nRows != this->nRows || bMatrix.nColumns != this->nColumns ) {
        OOFEM_ERROR("error in product receiver : dimensions do not match");
    }
#  endif

#  ifdef __LAPACK_MODULE
    double alpha = 1., beta = 1.;
    dgemm_("n", "n", & this->nRows, & this->nColumns, & aMatrix.nColumns,
           & alpha, aMatrix.givePointer(), & aMatrix.nRows, bMatrix.givePointer(), & bMatrix.nRows,
           & beta, this->givePointer(), & this->nRows,
           aMatrix.nColumns, bMatrix.nColumns, this->nColumns);
#  else
    for ( int i = 1; i <= aMatrix.nRows; i++ ) {
        for ( int j = 1; j <= bMatrix.nColumns; j++ ) {
            double coeff = 0.;
            for ( int k = 1; k <= aMatrix.nColumns; k++ ) {
                coeff += aMatrix.at(i, k) * bMatrix.at(k, j);
            }

            this->at(i, j) += coeff;
        }
    }
#  endif
}

void FloatMatrix :: addTProductOf(const FloatMatrix &aMatrix, const FloatMatrix &bMatrix)
// Receiver += aMatrix^T * bMatrix
{
#  ifdef DEBUG
    if ( aMatrix.nRows != bMatrix.nRows ) {
        OOFEM_ERROR("error in product A*B : dimensions do not match");
    }
    if ( aMatrix.nColumns != this->nColumns || bMatrix.nColumns != this->nRows ) {
        OOFEM_ERROR("error in product receiver : dimensions do not match");
    }
#  endif

#  ifdef __LAPACK_MODULE
    double alpha = 1., beta = 1.;
    dgemm_("t", "n",  & this->nRows, & this->nColumns, & aMatrix.nRows,
           & alpha, aMatrix.givePointer(), & aMatrix.nRows, bMatrix.givePointer(), & bMatrix.nRows,
           & beta, this->givePointer(), & this->nRows,
           aMatrix.nColumns, bMatrix.nColumns, this->nColumns);
#  else
    for ( int i = 1; i <= aMatrix.nColumns; i++ ) {
        for ( int j = 1; j <= bMatrix.nColumns; j++ ) {
            double coeff = 0.;
            for ( int k = 1; k <= aMatrix.nRows; k++ ) {
                coeff += aMatrix.at(k, i) * bMatrix.at(k, j);
            }

            this->at(i, j) += coeff;
        }
    }
#  endif
}


void FloatMatrix :: beDyadicProductOf(const FloatArray &vec1, const FloatArray &vec2)
// Receiver = vec1 * vec2^T
{
    int n1 = vec1.giveSize();
    int n2 = vec2.giveSize();
    RESIZE(n1, n2);
    for ( int j = 1; j <= n2; j++ ) {
        for ( int i = 1; i <= n1; i++ ) {
            this->at(i, j) = vec1.at(i) * vec2.at(j);
        }
    }
}

void FloatMatrix :: beNMatrixOf(const FloatArray &n, int nsd)
{
    this->resize(nsd, n.giveSize() * nsd);
    for ( int i = 0; i < n.giveSize(); ++i ) {
        for ( int j = 0; j < nsd; ++j ) {
            ( * this )( j, i * nsd + j ) = n(i);
        }
    }
}

void FloatMatrix :: beLocalCoordSys(const FloatArray &normal)
{
    if ( normal.giveSize() == 1 ) {
        this->resize(1, 1);
        this->at(1, 1) = normal(0);
    } else if ( normal.giveSize() == 2 ) {
        this->resize(2, 2);
        this->at(1, 1) = normal(1);
        this->at(1, 2) = -normal(0);

        this->at(2, 1) = normal(0);
        this->at(2, 2) = normal(1);
    } else if ( normal.giveSize() == 3 ) {
        // Create a permutated vector of n, *always* length 1 and significantly different from n.
        FloatArray b, t = {
            normal(1), -normal(2), normal(0)
        };                                                    // binormal and tangent

        // Construct orthogonal vector
        double npn = t.dotProduct(normal);
        t.add(-npn, normal);
        t.normalize();
        b.beVectorProductOf(t, normal);

        this->resize(3, 3);

        this->at(1, 1) = t(0);
        this->at(1, 2) = t(1);
        this->at(1, 3) = t(2);

        this->at(2, 1) = b(0);
        this->at(2, 2) = b(1);
        this->at(2, 3) = b(2);

        this->at(3, 1) = normal(0);
        this->at(3, 2) = normal(1);
        this->at(3, 3) = normal(2);
    } else {
        OOFEM_ERROR("Normal needs 1 to 3 components.");
    }
}

void FloatMatrix :: setSubMatrix(const FloatMatrix &src, int sr, int sc)
{
    sr--;
    sc--;

    int srcRows = src.giveNumberOfRows(), srcCols = src.giveNumberOfColumns();
#ifdef DEBUG
    int nr = sr + srcRows;
    int nc = sc + srcCols;

    if ( ( this->giveNumberOfRows() < nr ) || ( this->giveNumberOfColumns() < nc ) ) {
        OOFEM_ERROR("Sub matrix doesn't fit inside allocated space.");
    }
#endif

    // add sub-matrix
    for ( int j = 0; j < srcCols; j++ ) {
        for ( int i = 0; i < srcRows; i++ ) {
            ( * this )( sr + i, sc + j ) = src(i, j);
        }
    }
    //memcpy( &(*this)(sr + 0, sc + j), &src(0,j), srcRows * sizeof(int));
}


void FloatMatrix :: setTSubMatrix(const FloatMatrix &src, int sr, int sc)
{
    sr--;
    sc--;

    int srcRows = src.giveNumberOfRows(), srcCols = src.giveNumberOfColumns();
#ifdef DEBUG
    int nr = sr + srcCols;
    int nc = sc + srcRows;

    if ( ( this->giveNumberOfRows() < nr ) || ( this->giveNumberOfColumns() < nc ) ) {
        OOFEM_ERROR("Sub matrix doesn't fit inside allocated space");
    }
#endif

    // add sub-matrix
    for ( int i = 0; i < srcCols; i++ ) {
        for ( int j = 0; j < srcRows; j++ ) {
            ( * this )( sr + i, sc + j ) = src(j, i);
        }
    }
}



void FloatMatrix :: addSubVectorRow(const FloatArray &src, int sr, int sc)
{
    sc--;

    int srcCols = src.giveSize();

    int nr = sr;
    int nc = sc + srcCols;

    if ( ( this->giveNumberOfRows() < nr ) || ( this->giveNumberOfColumns() < nc ) ) {
        this->resizeWithData( max(this->giveNumberOfRows(), nr), max(this->giveNumberOfColumns(), nc) );
    }

    // add sub-matrix
    for ( int j = 1; j <= srcCols; j++ ) {
        this->at(sr, sc + j) += src.at(j);
    }
}


void FloatMatrix :: addSubVectorCol(const FloatArray &src, int sr, int sc)
{
    sr--;

    int srcRows = src.giveSize();

    int nr = sr + srcRows;
    int nc = sc;

    if ( ( this->giveNumberOfRows() < nr ) || ( this->giveNumberOfColumns() < nc ) ) {
        this->resizeWithData( max(this->giveNumberOfRows(), nr), max(this->giveNumberOfColumns(), nc) );
    }

    // add sub-matrix
    for ( int j = 1; j <= srcRows; j++ ) {
        this->at(sr + j, sc) += src.at(j);
    }
}



void FloatMatrix :: setColumn(const FloatArray &src, int c)
{
    int nr = src.giveSize();
#ifdef DEBUG
    if ( this->giveNumberOfRows() != nr || c < 1 || c > this->giveNumberOfColumns() ) {
        OOFEM_ERROR("Size mismatch");
    }
#endif

    auto P = this->values.begin() + ( c - 1 ) * nr;
    std :: copy(src.begin(), src.end(), P);
}


void FloatMatrix :: copyColumn(FloatArray &dest, int c) const
{
    int nr = this->giveNumberOfRows();
#ifdef DEBUG
    if ( c < 1 || c > this->giveNumberOfColumns() ) {
        OOFEM_ERROR("Column outside range (%d)", c);
    }
#endif

    dest.resize(nr);
    auto P = this->values.begin() + ( c - 1 ) * nr;
    std :: copy( P, P + nr, dest.begin() );
}


void FloatMatrix :: copySubVectorRow(const FloatArray &src, int sr, int sc)
{
    sc--;

    int srcCols = src.giveSize();

    int nr = sr;
    int nc = sc + srcCols;

    if ( ( this->giveNumberOfRows() < nr ) || ( this->giveNumberOfColumns() < nc ) ) {
        this->resizeWithData( max(this->giveNumberOfRows(), nr), max(this->giveNumberOfColumns(), nc) );
    }

    // add sub-matrix
    for ( int j = 1; j <= srcCols; j++ ) {
        this->at(sr, sc + j) = src.at(j);
    }
}



void FloatMatrix :: plusProductSymmUpper(const FloatMatrix &a, const FloatMatrix &b, double dV)
// Adds to the receiver the product  a(transposed).b dV .
// The receiver size is adjusted, if necessary.
// This method assumes that both the receiver and the product above are
// symmetric matrices, and therefore computes only the upper half of the
// receiver ; the lower half is not modified. Other advantage : it does
// not compute the transposition of matrix a.
{
    if ( !this->isNotEmpty() ) {
        this->nRows = a.nColumns;
        this->nColumns = b.nColumns;
        this->values.assign(a.nColumns * b.nColumns, 0.);
    }

#ifdef __LAPACK_MODULE
    double beta = 1.;
    ///@todo We should determine which is the best choice overall. For large systems more block matrix operations is necessary.
    /// For smaller systems the overhead from function calls might be larger, but the overhead might be tiny, or using symmetry at all might be undesireable.
    if ( this->nRows < 20 ) {
        // Split the matrix into 2 columns, s1 + s2 = n ( = nRows = nColumns ).
        int s1 = this->nRows / 2;
        int s2 = this->nRows - s1;
        // First column block, we only take the first s rows by only taking the first s columns in the matrix a.
        dgemm_("t", "n", & s1, & s1, & a.nRows, & dV, a.givePointer(), & a.nRows, b.givePointer(), & b.nRows, & beta, this->givePointer(), & this->nRows, a.nColumns, b.nColumns, this->nColumns);
        // Second column block starting a memory position c * nRows
        dgemm_("t", "n", & this->nRows, & s2, & a.nRows, & dV, a.givePointer(), & a.nRows, & b.givePointer() [ s1 * b.nRows ], & b.nRows, & beta, & this->givePointer() [ s1 * this->nRows ], & this->nRows, a.nColumns, b.nColumns, this->nColumns);
    } else {
        // Get suitable blocksize. Around 10 rows should be suitable (slightly adjusted to minimize number of blocks):
        // Smaller blocks than ~10 didn't show any performance gains in my benchmarks. / Mikael
        int block = ( this->nRows - 1 ) / ( this->nRows / 10 ) + 1;
        int start = 0;
        int end = block;
        while ( start < this->nRows ) {
            int s = end - start;
            dgemm_("t", "n", & end, & s, & a.nRows, & dV, a.givePointer(), & a.nRows, & b.givePointer() [ start * b.nRows ], & b.nRows, & beta, & this->givePointer() [ start * this->nRows ],
                   & this->nRows, a.nColumns, b.nColumns, this->nColumns);
            start = end;
            end += block;
            if ( end > this->nRows ) {
                end = this->nRows;
            }
        }
    }
#else
    for ( int i = 1; i <= nRows; i++ ) {
        for ( int j = i; j <= nColumns; j++ ) {
            double summ = 0.;
            for ( int k = 1; k <= a.nRows; k++ ) {
                summ += a.at(k, i) * b.at(k, j);
            }

            this->at(i, j) += summ * dV;
        }
    }
#endif
}


void FloatMatrix :: plusDyadSymmUpper(const FloatArray &a, double dV)
{
    if ( !this->isNotEmpty() ) {
        this->nRows = a.giveSize();
        this->nColumns = a.giveSize();
        this->values.assign(this->nRows * this->nColumns, 0.);
    }
#ifdef __LAPACK_MODULE
    int inc = 1;
    int sizeA = a.giveSize();
    dsyr_("u", & sizeA, & dV, a.givePointer(), & inc,
          this->givePointer(), & this->nRows,
          sizeA, this->nColumns);
#else
    for ( int i = 1; i <= nRows; i++ ) {
        for ( int j = i; j <= nColumns; j++ ) {
            this->at(i, j) += a.at(i) * a.at(j) * dV;
        }
    }
#endif
}



void FloatMatrix :: plusProductUnsym(const FloatMatrix &a, const FloatMatrix &b, double dV)
// Adds to the receiver the product  a(transposed).b dV .
// If the receiver has a null size, it is expanded.
// Advantage : does not compute the transposition of matrix a.
{
    if ( !this->isNotEmpty() ) {
        this->nRows = a.nColumns;
        this->nColumns = b.nColumns;
        this->values.assign(this->nRows * this->nColumns, 0.);
    }
#ifdef __LAPACK_MODULE
    double beta = 1.;
    dgemm_("t", "n", & this->nRows, & this->nColumns, & a.nRows,
           & dV, a.givePointer(), & a.nRows, b.givePointer(), & b.nRows,
           & beta, this->givePointer(), & this->nRows,
           a.nColumns, b.nColumns, this->nColumns);
#else
    for ( int i = 1; i <= nRows; i++ ) {
        for ( int j = 1; j <= nColumns; j++ ) {
            double summ = 0.;
            for ( int k = 1; k <= a.nRows; k++ ) {
                summ += a.at(k, i) * b.at(k, j);
            }

            this->at(i, j) += summ * dV;
        }
    }
#endif
}


void FloatMatrix :: plusDyadUnsym(const FloatArray &a, const FloatArray &b, double dV)
{
    if ( !this->isNotEmpty() ) {
        this->nRows = a.giveSize();
        this->nColumns = b.giveSize();
        this->values.assign(this->nRows * this->nColumns, 0.);
    }
#ifdef __LAPACK_MODULE
    int inc = 1;
    int sizeA = a.giveSize();
    int sizeB = b.giveSize();
    dger_(& sizeA, & sizeB, & dV, a.givePointer(), & inc,
          b.givePointer(), & inc, this->givePointer(), & this->nRows,
          sizeA, sizeB, this->nColumns);
#else
    for ( int i = 1; i <= nRows; i++ ) {
        for ( int j = 1; j <= nColumns; j++ ) {
            this->at(i, j) += a.at(i) * b.at(j) * dV;
        }
    }
#endif
}


void FloatMatrix :: beInverseOf(const FloatMatrix &src)
// Receiver becomes inverse of given parameter src. If necessary, size is adjusted.
{
    double det;

#  ifdef DEBUG
    if ( !src.isSquare() ) {
        OOFEM_ERROR("cannot inverse a %d by %d matrix", src.nRows, src.nColumns);
    }
#  endif

    RESIZE(src.nRows, src.nColumns);

    if ( nRows == 1 ) {
        this->at(1, 1) = 1. / src.at(1, 1);
        return;
    } else if ( nRows == 2 ) {
        det = src.at(1, 1) * src.at(2, 2) - src.at(1, 2) * src.at(2, 1);
        this->at(1, 1) =  src.at(2, 2) / det;
        this->at(2, 1) = -src.at(2, 1) / det;
        this->at(1, 2) = -src.at(1, 2) / det;
        this->at(2, 2) =  src.at(1, 1) / det;
        return;
    } else if ( nRows == 3 ) {
        det = src.at(1, 1) * src.at(2, 2) * src.at(3, 3) + src.at(1, 2) * src.at(2, 3) * src.at(3, 1) +
              src.at(1, 3) * src.at(2, 1) * src.at(3, 2) - src.at(1, 3) * src.at(2, 2) * src.at(3, 1) -
              src.at(2, 3) * src.at(3, 2) * src.at(1, 1) - src.at(3, 3) * src.at(1, 2) * src.at(2, 1);

        this->at(1, 1) = ( src.at(2, 2) * src.at(3, 3) - src.at(2, 3) * src.at(3, 2) ) / det;
        this->at(2, 1) = ( src.at(2, 3) * src.at(3, 1) - src.at(2, 1) * src.at(3, 3) ) / det;
        this->at(3, 1) = ( src.at(2, 1) * src.at(3, 2) - src.at(2, 2) * src.at(3, 1) ) / det;
        this->at(1, 2) = ( src.at(1, 3) * src.at(3, 2) - src.at(1, 2) * src.at(3, 3) ) / det;
        this->at(2, 2) = ( src.at(1, 1) * src.at(3, 3) - src.at(1, 3) * src.at(3, 1) ) / det;
        this->at(3, 2) = ( src.at(1, 2) * src.at(3, 1) - src.at(1, 1) * src.at(3, 2) ) / det;
        this->at(1, 3) = ( src.at(1, 2) * src.at(2, 3) - src.at(1, 3) * src.at(2, 2) ) / det;
        this->at(2, 3) = ( src.at(1, 3) * src.at(2, 1) - src.at(1, 1) * src.at(2, 3) ) / det;
        this->at(3, 3) = ( src.at(1, 1) * src.at(2, 2) - src.at(1, 2) * src.at(2, 1) ) / det;

        //p[0]= (values[4]*values[8]-values[7]*values[5])/det ;
        //p[1]= (values[7]*values[2]-values[1]*values[8])/det ;
        //p[2]= (values[1]*values[5]-values[4]*values[2])/det ;
        //p[3]= (values[6]*values[5]-values[3]*values[8])/det ;
        //p[4]= (values[0]*values[8]-values[6]*values[2])/det ;
        //p[5]= (values[3]*values[2]-values[0]*values[5])/det ;
        //p[6]= (values[3]*values[7]-values[6]*values[4])/det ;
        //p[7]= (values[6]*values[1]-values[0]*values[7])/det ;
        //p[8]= (values[0]*values[4]-values[3]*values[1])/det ;
        return;
    } else {
#ifdef __LAPACK_MODULE
        int n = this->nRows;
        IntArray ipiv(n);
        int lwork, info;
        * this = src;

        // LU-factorization
        dgetrf_(& n, & n, this->givePointer(), & n, ipiv.givePointer(), & info);
        if ( info != 0 ) {
            OOFEM_ERROR("dgetrf error %d", info);
        }

        // Inverse
        lwork = n * n;
        FloatArray work(lwork);
        dgetri_(& this->nRows, this->givePointer(), & this->nRows, ipiv.givePointer(), work.givePointer(), & lwork, & info);
        if ( info > 0 ) {
            OOFEM_ERROR("Singular at %d", info);
        } else if ( info < 0 ) {
            OOFEM_ERROR("Error on input %d", info);
        }
#else
        // size >3 ... gaussian elimination - slow but safe
        //
        double piv, linkomb;
        FloatMatrix tmp = src;
        // initialize answer to be unity matrix;
        this->zero();
        for ( int i = 1; i <= nRows; i++ ) {
            this->at(i, i) = 1.0;
        }

        // lower triangle elimination by columns
        for ( int i = 1; i < nRows; i++ ) {
            piv = tmp.at(i, i);
            if ( fabs(piv) < 1.e-20 ) {
                OOFEM_ERROR("cannot inverse a %d by %d matrix", nRows, nColumns);
            }

            for ( int j = i + 1; j <= nRows; j++ ) {
                linkomb = tmp.at(j, i) / tmp.at(i, i);
                for ( int k = i; k <= nRows; k++ ) {
                    tmp.at(j, k) -= tmp.at(i, k) * linkomb;
                }

                for ( int k = 1; k <= nRows; k++ ) {
                    this->at(j, k) -= this->at(i, k) * linkomb;
                }
            }
        }

        // upper triangle elimination by columns
        for ( int i = nRows; i > 1; i-- ) {
            piv = tmp.at(i, i);
            for ( int j = i - 1; j > 0; j-- ) {
                linkomb = tmp.at(j, i) / piv;
                for ( int k = i; k > 0; k-- ) {
                    tmp.at(j, k) -= tmp.at(i, k) * linkomb;
                }

                for ( int k = nRows; k > 0; k-- ) {
                    // tmp -> at(j,k)-= tmp  ->at(i,k)*linkomb;
                    this->at(j, k) -= this->at(i, k) * linkomb;
                }
            }
        }

        // diagonal scaling
        for ( int i = 1; i <= nRows; i++ ) {
            for ( int j = 1; j <= nRows; j++ ) {
                this->at(i, j) /= tmp.at(i, i);
            }
        }
#endif
    }
}


void FloatMatrix :: beSubMatrixOf(const FloatMatrix &src,
                                  int topRow, int bottomRow, int topCol, int bottomCol)
/*
 * modifies receiver to be  submatrix of the src matrix
 * size of receiver  submatrix is determined from
 * input parameters
 */
{
#ifdef DEBUG
    if ( ( topRow < 1 ) || ( bottomRow < 1 ) || ( topCol < 1 ) || ( bottomCol < 1 ) ) {
        OOFEM_ERROR("subindexes size mismatch");
    }

    if ( ( src.nRows < bottomRow ) || ( src.nColumns < bottomCol ) || ( ( bottomRow - topRow ) > src.nRows ) ||
         ( ( bottomCol - topCol ) > src.nColumns ) ) {
        OOFEM_ERROR("subindexes size mismatch");
    }
#endif


    int topRm1, topCm1;
    topRm1 = topRow - 1;
    topCm1 = topCol - 1;

    // allocate return value
    this->resize(bottomRow - topRm1, bottomCol - topCm1);
    for ( int i = topRow; i <= bottomRow; i++ ) {
        for ( int j = topCol; j <= bottomCol; j++ ) {
            this->at(i - topRm1, j - topCm1) = src.at(i, j);
        }
    }
}



void
FloatMatrix :: beSubMatrixOf(const FloatMatrix &src, const IntArray &indxRow, const IntArray &indxCol)
/*
 * Modifies receiver to be a sub-matrix of the src matrix.
 * sub-matrix has size(indxRow) x size(indxCol) with values given as
 * this(i,j) = src( indxRow(i), indxCol(j) )
 */
{
#  ifdef DEBUG
    if ( indxRow.maximum() > src.giveNumberOfRows()  ||  indxCol.maximum() > src.giveNumberOfColumns()  ||
         indxRow.minimum() < 1  ||  indxCol.minimum() < 1 ) {
        OOFEM_ERROR("index exceeds source dimensions");
    }
# endif

    int szRow = indxRow.giveSize();
    int szCol = indxCol.giveSize();
    this->resize(szRow, szCol);

    for ( int i = 1; i <= szRow; i++ ) {
        for ( int j = 1; j <= szCol; j++ ) {
            this->at(i, j) = src.at( indxRow.at(i), indxCol.at(j) );
        }
    }
}

void FloatMatrix :: add(const FloatMatrix &aMatrix)
// Adds aMatrix to the receiver. If the receiver has a null size,
// adjusts its size to that of aMatrix. Returns the modified receiver.
{
    if ( aMatrix.nRows == 0 || aMatrix.nColumns == 0 ) {
        return;
    }

    if ( !this->isNotEmpty() ) {
        this->operator = ( aMatrix );
        return;
    }
#     ifdef DEBUG
    if ( ( aMatrix.nRows != nRows || aMatrix.nColumns != nColumns ) && aMatrix.isNotEmpty() ) {
        OOFEM_ERROR("dimensions mismatch : (r1,c1)+(r2,c2) : (%d,%d)+(%d,%d)", nRows, nColumns, aMatrix.nRows, aMatrix.nColumns);
    }
#     endif

#ifdef __LAPACK_MODULE
    int aSize = aMatrix.nRows * aMatrix.nColumns;
    int inc = 1;
    double s = 1.;
    daxpy_(& aSize, & s, aMatrix.givePointer(), & inc, this->givePointer(), & inc, aSize, aSize);
#else
    for ( size_t i = 0; i < this->values.size(); i++ ) {
        this->values [ i ] += aMatrix.values [ i ];
    }
#endif
}


void FloatMatrix :: add(double s, const FloatMatrix &aMatrix)
// Adds aMatrix to the receiver. If the receiver has a null size,
// adjusts its size to that of aMatrix. Returns the modified receiver.
{
    if ( aMatrix.nRows == 0 || aMatrix.nColumns == 0 ) {
        return;
    }

    if ( !this->isNotEmpty() ) {
        this->operator = ( aMatrix );
        this->times(s);
        return;
    }
#     ifdef DEBUG
    if ( ( aMatrix.nRows != nRows || aMatrix.nColumns != nColumns ) && aMatrix.isNotEmpty() ) {
        OOFEM_ERROR("dimensions mismatch : (r1,c1)+(r2,c2) : (%d,%d)+(%d,%d)", nRows, nColumns, aMatrix.nRows, aMatrix.nColumns);
    }
#     endif

#ifdef __LAPACK_MODULE
    int aSize = aMatrix.nRows * aMatrix.nColumns;
    int inc = 1;
    daxpy_(& aSize, & s, aMatrix.givePointer(), & inc, this->givePointer(), & inc, aSize, aSize);
#else
    for ( size_t i = 0; i < this->values.size(); i++ ) {
        this->values [ i ] += s * aMatrix.values [ i ];
    }
#endif
}

void FloatMatrix :: subtract(const FloatMatrix &aMatrix)
// Adds aMatrix to the receiver. If the receiver has a null size,
// adjusts its size to that of aMatrix. Returns the modified receiver.
{
    if ( !this->isNotEmpty() ) {
        this->operator = ( aMatrix );
        this->negated();
        return;
    }
#     ifdef DEBUG
    if ( ( aMatrix.nRows != nRows || aMatrix.nColumns != nColumns ) && aMatrix.isNotEmpty() ) {
        OOFEM_ERROR("dimensions mismatch : (r1,c1)-(r2,c2) : (%d,%d)-(%d,%d)", nRows, nColumns, aMatrix.nRows, aMatrix.nColumns);
    }
#     endif

#ifdef __LAPACK_MODULE
    int aSize = aMatrix.nRows * aMatrix.nColumns;
    int inc = 1;
    double s = -1.;
    daxpy_(& aSize, & s, aMatrix.givePointer(), & inc, this->givePointer(), & inc, aSize, aSize);
#else
    for ( size_t i = 0; i < this->values.size(); i++ ) {
        this->values [ i ] -= aMatrix.values [ i ];
    }
#endif
}


bool FloatMatrix :: solveForRhs(const FloatArray &b, FloatArray &answer, bool transpose)
// solves equation b = this * x
{
#  ifdef DEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR("cannot solve a %d by %d matrix", nRows, nColumns);
    }

    if ( nRows != b.giveSize() ) {
        OOFEM_ERROR("dimension mismatch");
    }
#  endif

#ifdef __LAPACK_MODULE
    int info, nrhs = 1;
    IntArray ipiv(this->nRows);
    answer = b;
    //dgesv_( &this->nRows, &nrhs, this->givePointer(), &this->nRows, ipiv.givePointer(), answer.givePointer(), &this->nRows, &info );
    dgetrf_(& this->nRows, & this->nRows, this->givePointer(), & this->nRows, ipiv.givePointer(), & info);
    if ( info == 0 ) {
        dgetrs_(transpose ? "Transpose" : "No transpose", & this->nRows, & nrhs, this->givePointer(), & this->nRows, ipiv.givePointer(), answer.givePointer(), & this->nRows, & info);
    }
    if ( info != 0 ) {
        return false;
    }
#else
    int pivRow;
    double piv, linkomb, help;
    FloatMatrix *mtrx, trans;
    if ( transpose ) {
        trans.beTranspositionOf(* this);
        mtrx = & trans;
    } else {
        mtrx = this;
    }

    answer = b;

    // initialize answer to be unity matrix;
    // lower triangle elimination by columns
    for ( int i = 1; i < nRows; i++ ) {
        // find the suitable row and pivot
        piv = fabs( mtrx->at(i, i) );
        pivRow = i;
        for ( int j = i + 1; j <= nRows; j++ ) {
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
            for ( int j = i; j <= nRows; j++ ) {
                help = mtrx->at(i, j);
                mtrx->at(i, j) = mtrx->at(pivRow, j);
                mtrx->at(pivRow, j) = help;
            }
            help = answer.at(i);
            answer.at(i) = answer.at(pivRow);
            answer.at(pivRow) = help;
        }

        for ( int j = i + 1; j <= nRows; j++ ) {
            linkomb = mtrx->at(j, i) / mtrx->at(i, i);
            for ( int k = i; k <= nRows; k++ ) {
                mtrx->at(j, k) -= mtrx->at(i, k) * linkomb;
            }

            answer.at(j) -= answer.at(i) * linkomb;
        }
    }

    // back substitution
    for ( int i = nRows; i >= 1; i-- ) {
        help = 0.;
        for ( int j = i + 1; j <= nRows; j++ ) {
            help += mtrx->at(i, j) * answer.at(j);
        }

        answer.at(i) = ( answer.at(i) - help ) / mtrx->at(i, i);
    }
#endif

    return true;
}


void FloatMatrix :: solveForRhs(const FloatMatrix &b, FloatMatrix &answer, bool transpose)
// solves equation b = this * x
// returns x. this and b are kept untouched
//
// gaussian elimination - slow but safe
//
{
#  ifdef DEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR("cannot solve a %d by %d matrix", nRows, nColumns);
    }

    if ( nRows != b.giveNumberOfRows() ) {
        OOFEM_ERROR("dimension mismatch");
    }
#  endif

#ifdef __LAPACK_MODULE
    int info;
    IntArray ipiv(this->nRows);
    answer = b;
    dgetrf_(& this->nRows, & this->nRows, this->givePointer(), & this->nRows, ipiv.givePointer(), & info);
    if ( info == 0 ) {
        dgetrs_(transpose ? "t" : "n", & this->nRows, & answer.nColumns, this->givePointer(), & this->nRows, ipiv.givePointer(), answer.givePointer(), & this->nRows, & info);
    }
    if ( info != 0 ) {
        OOFEM_ERROR("error %d", info);
    }
#else
    int pivRow, nPs;
    double piv, linkomb, help;
    FloatMatrix *mtrx, trans;
    if ( transpose ) {
        trans.beTranspositionOf(* this);
        mtrx = & trans;
    } else {
        mtrx = this;
    }

    nPs = b.giveNumberOfColumns();
    answer = b;
    // initialize answer to be unity matrix;
    // lower triangle elimination by columns
    for ( int i = 1; i < nRows; i++ ) {
        // find the suitable row and pivot
        piv = fabs( mtrx->at(i, i) );
        pivRow = i;
        for ( int j = i + 1; j <= nRows; j++ ) {
            if ( fabs( mtrx->at(j, i) ) > piv ) {
                pivRow = j;
                piv = fabs( mtrx->at(j, i) );
            }
        }

        if ( fabs(piv) < 1.e-20 ) {
            OOFEM_ERROR("pivot too small, cannot solve %d by %d matrix", nRows, nColumns);
        }

        // exchange rows
        if ( pivRow != i ) {
            for ( int j = i; j <= nRows; j++ ) {
                help = mtrx->at(i, j);
                mtrx->at(i, j) = mtrx->at(pivRow, j);
                mtrx->at(pivRow, j) = help;
            }

            for ( int j = 1; j <= nPs; j++ ) {
                help = answer.at(i, j);
                answer.at(i, j) = answer.at(pivRow, j);
                answer.at(pivRow, j) = help;
            }
        }

        if ( fabs(piv) < 1.e-20 ) {
            OOFEM_ERROR("cannot solve, zero pivot encountered");
        }

        for ( int j = i + 1; j <= nRows; j++ ) {
            linkomb = mtrx->at(j, i) / mtrx->at(i, i);
            for ( int k = i; k <= nRows; k++ ) {
                mtrx->at(j, k) -= mtrx->at(i, k) * linkomb;
            }

            for ( int k = 1; k <= nPs; k++ ) {
                answer.at(j, k) -= answer.at(i, k) * linkomb;
            }
        }
    }

    // back substitution
    for ( int i = nRows; i >= 1; i-- ) {
        for ( int k = 1; k <= nPs; k++ ) {
            help = 0.;
            for ( int j = i + 1; j <= nRows; j++ ) {
                help += mtrx->at(i, j) * answer.at(j, k);
            }

            answer.at(i, k) = ( answer.at(i, k) - help ) / mtrx->at(i, i);
        }
    }
#endif
}


void FloatMatrix :: initFromVector(const FloatArray &vector, bool transposed)
//
// constructor : creates (vector->giveSize(),1) FloatMatrix
// if transpose = 1 creates (1,vector->giveSize()) FloatMatrix
//
{
    if ( transposed ) {
        this->nRows = 1;
        this->nColumns = vector.giveSize();
    } else {
        this->nRows = vector.giveSize();
        this->nColumns = 1;
    }

    this->values = vector.values;
}


void FloatMatrix :: zero()
{
    std :: fill(this->values.begin(), this->values.end(), 0.);
}


void FloatMatrix :: beUnitMatrix()
{
#ifdef DEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR("cannot make unit matrix of %d by %d matrix", nRows, nColumns);
    }
#endif

    this->zero();
    for ( int i = 1; i <= nRows; i++ ) {
        this->at(i, i) = 1.0;
    }
}


void FloatMatrix :: bePinvID()
// this matrix is the product of the 6x6 deviatoric projection matrix ID
// and the inverse scaling matrix Pinv
{
    this->resize(6, 6);
    values [ 0 ] = values [ 7 ] = values [ 14 ] = 2. / 3.;
    values [ 1 ] = values [ 2 ] = values [ 6 ] = values [ 8 ] = values [ 12 ] = values [ 13 ] = -1. / 3.;
    values [ 21 ] = values [ 28 ] = values [ 35 ] = 0.5;
}


void FloatMatrix :: resize(int rows, int columns)
//
// resizes receiver, all data will be lost
//
{
    this->nRows = rows;
    this->nColumns = columns;
    this->values.assign(rows * columns, 0.);
}


void FloatMatrix :: resizeWithData(int rows, int columns)
//
// resizes receiver, all data kept
//
{
    // Check of resize if necessary at all.
    if ( rows == this->nRows && columns == this->nColumns ) {
        return;
    }

    FloatMatrix old( std :: move(* this) );

    this->nRows = rows;
    this->nColumns = columns;
    this->values.resize(rows * columns);

    int ii = min( rows, old.giveNumberOfRows() );
    int jj = min( columns, old.giveNumberOfColumns() );
    // copy old values if possible
    for ( int i = 1; i <= ii; i++ ) {
        for ( int j = 1; j <= jj; j++ ) {
            this->at(i, j) = old.at(i, j);
        }
    }
}


void FloatMatrix :: hardResize(int rows, int columns)
//
// resizes receiver, all data will be lost
//
{
    this->nRows = rows;
    this->nColumns = columns;
    values.assign(rows * columns, 0.);
    this->values.shrink_to_fit();
}


double FloatMatrix :: giveDeterminant() const
// Returns the determinant of the receiver.
{
#  ifdef DEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR("cannot compute the determinant of a non-square %d by %d matrix", nRows, nColumns);
    }
#  endif

    if ( nRows == 1 ) {
        return values [ 0 ];
    } else if ( nRows == 2 ) {
        return ( values [ 0 ] * values [ 3 ] - values [ 1 ] * values [ 2 ] );
    } else if ( nRows == 3 ) {
        return ( values [ 0 ] * values [ 4 ] * values [ 8 ] + values [ 3 ] * values [ 7 ] * values [ 2 ] +
                 values [ 6 ] * values [ 1 ] * values [ 5 ] - values [ 6 ] * values [ 4 ] * values [ 2 ] -
                 values [ 7 ] * values [ 5 ] * values [ 0 ] - values [ 8 ] * values [ 3 ] * values [ 1 ] );
    } else {
        OOFEM_ERROR("sorry, cannot compute the determinant of a matrix larger than 3x3");
    }

    return 0.;
}


void FloatMatrix :: beDiagonal(const FloatArray &diag)
{
    int n = diag.giveSize();
    this->resize(n, n);
    for ( int i = 0; i < n; ++i ) {
        (*this)(i, i) = diag[i];
    }
}


double FloatMatrix :: giveTrace() const
{
#  ifdef DEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR("cannot compute the trace of a non-square %d by %d matrix", nRows, nColumns);
    }
#  endif
    double answer = 0.;
    for ( int k = 0; k < nRows; k++ ) {
        answer += values [ k * ( nRows + 1 ) ];
    }
    return answer;
}


void FloatMatrix :: printYourself() const
// Prints the receiver on screen.
{
    printf("FloatMatrix with dimensions : %d %d\n",
           nRows, nColumns);
    if ( nRows <= 250 && nColumns <= 250 ) {
        for ( int i = 1; i <= nRows; ++i ) {
            for ( int j = 1; j <= nColumns && j <= 100; ++j ) {
                printf( "%10.3e  ", this->at(i, j) );
            }

            printf("\n");
        }
    } else {
        printf("   large matrix : coefficients not printed \n");
    }
}


void FloatMatrix :: printYourself(const std::string &name) const
// Prints the receiver on screen.
{
    printf("%s (%d x %d): \n", name.c_str(), nRows, nColumns);
    if ( nRows <= 250 && nColumns <= 250 ) {
        for ( int i = 1; i <= nRows; ++i ) {
            for ( int j = 1; j <= nColumns && j <= 100; ++j ) {
                printf( "%10.3e  ", this->at(i, j) );
            }

            printf("\n");
        }
    } else {
        printf("   large matrix : coefficients not printed \n");
    }
}


void FloatMatrix :: pY() const
// Prints the receiver on screen with higher accuracy than printYourself.
{
    printf("[");
    for ( int i = 1; i <= nRows; ++i ) {
        for ( int j = 1; j <= nColumns; ++j ) {
            printf( "%20.15e", this->at(i, j) );
            if ( j < nColumns ) {
                printf(",");
            } else {
                printf(";");
            }
        }
    }

    printf("];\n");
}


void FloatMatrix :: writeCSV(const std :: string &name) const
{
    FILE *file = fopen(name.c_str(), "w");
    for ( int i = 1; i <= nRows; ++i ) {
        for ( int j = 1; j <= nColumns; ++j ) {
            fprintf(file, "%10.3e, ", this->at(i, j) );
        }

        fprintf(file, "\n");
    }
    fclose(file);
}


void FloatMatrix :: rotatedWith(const FloatMatrix &r, char mode)
// Returns the receiver 'a' rotated according the change-of-base matrix r.
// The method performs the operation  a = r^T . a . r . or the inverse
{
    FloatMatrix rta;

    if ( mode == 'n' ) {
        rta.beTProductOf(r, * this);     //  r^T . a
        this->beProductOf(rta, r);       //  r^T . a . r
    } else if ( mode == 't' ) {
        rta.beProductOf(r, * this);      //  r . a
        this->beProductTOf(rta, r);      //  r . a . r^T
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}



void FloatMatrix :: symmetrized()
// Initializes the lower half of the receiver to the upper half.
{
#  ifdef DEBUG
    if ( nRows != nColumns ) {
        OOFEM_ERROR("cannot symmetrize a non-square matrix");
    }

#   endif

    for ( int i = 2; i <= nRows; i++ ) {
        for ( int j = 1; j < i; j++ ) {
            this->at(i, j) = this->at(j, i);
        }
    }
}


void FloatMatrix :: times(double factor)
// Multiplies every coefficient of the receiver by factor. Answers the
// modified receiver.
{
    for ( double &x : this->values ) {
        x *= factor;
    }
    // dscal_ seemed to be slower for typical usage of this function.
}


void FloatMatrix :: negated()
{
    for ( double &x : this->values ) {
        x = -x;
    }
}


double FloatMatrix :: computeFrobeniusNorm() const
{
    return sqrt( std :: inner_product(this->values.begin(), this->values.end(), this->values.begin(), 0.) );
}

double FloatMatrix :: computeNorm(char p) const
{
#  ifdef __LAPACK_MODULE
    FloatArray work( this->giveNumberOfRows() );
    int lda = max(this->nRows, 1);
    double norm = dlange_(& p, & this->nRows, & this->nColumns, this->givePointer(), & lda, work.givePointer(), 1);
    return norm;

#  else
    if ( p == '1' ) { // Maximum absolute column sum.
        double col_sum, max_col = 0.0;
        for ( int j = 1; j <= this->nColumns; j++ ) {
            col_sum  = 0.0;
            for ( int i = 1; i <= this->nRows; i++ ) {
                col_sum += fabs( this->at(i, j) );
            }
            if ( col_sum > max_col ) {
                max_col = col_sum;
            }
        }
        return max_col;
    }
    ///@todo Use this when obtaining eigen values is implemented.
    /*else if (p == '2') {
     *  double lambda_max;
     *  FloatMatrix AtA;
     *  FloatArray eigs;
     *  AtA.beTProductOf(*this,this);
     *  Ata.eigenValues(eigs, 1);
     *  return sqrt(eigs(0));
     * } */else {
        OOFEM_ERROR("p == %d not implemented.\n", p);
        return 0.0;
    }
#  endif
}



void FloatMatrix :: beMatrixForm(const FloatArray &aArray)
{
    // Revrites the vector on matrix form (symmetrized matrix used if size is 6),
    // order: 11, 22, 33, 23, 13, 12
    // order: 11, 22, 33, 23, 13, 12, 32, 31, 21
#  ifdef DEBUG
    if ( aArray.giveSize() != 6 && aArray.giveSize() != 9 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }
#  endif
    this->resize(3, 3);
    if ( aArray.giveSize() == 9 ) {
        this->at(1, 1) = aArray.at(1);
        this->at(2, 2) = aArray.at(2);
        this->at(3, 3) = aArray.at(3);
        this->at(2, 3) = aArray.at(4);
        this->at(1, 3) = aArray.at(5);
        this->at(1, 2) = aArray.at(6);
        this->at(3, 2) = aArray.at(7);
        this->at(3, 1) = aArray.at(8);
        this->at(2, 1) = aArray.at(9);
    } else if ( aArray.giveSize() == 6 ) {
        this->at(1, 1) = aArray.at(1);
        this->at(2, 2) = aArray.at(2);
        this->at(3, 3) = aArray.at(3);
        this->at(2, 3) = aArray.at(4);
        this->at(1, 3) = aArray.at(5);
        this->at(1, 2) = aArray.at(6);
        this->at(3, 2) = aArray.at(4);
        this->at(3, 1) = aArray.at(5);
        this->at(2, 1) = aArray.at(6);
    }
}

void FloatMatrix :: changeComponentOrder()
{
    // Changes index order between abaqus <-> OOFEM
//#  ifdef DEBUG
//    if ( nRows != 6 || nColumns != 6 ) {
//        OOFEM_ERROR("matrix dimension is not 6x6");
//    }
//#  endif

    if ( nRows == 6 && nColumns == 6 ) {
        // This could probably be done more beautifully + efficiently.

        std :: swap( this->at(4, 1), this->at(6, 1) );

        std :: swap( this->at(4, 2), this->at(6, 2) );

        std :: swap( this->at(4, 3), this->at(6, 3) );

        std :: swap( this->at(1, 4), this->at(1, 6) );
        std :: swap( this->at(2, 4), this->at(2, 6) );
        std :: swap( this->at(3, 4), this->at(3, 6) );
        std :: swap( this->at(4, 4), this->at(6, 6) );
        std :: swap( this->at(5, 4), this->at(5, 6) );
        std :: swap( this->at(6, 4), this->at(4, 6) );

        std :: swap( this->at(4, 5), this->at(6, 5) );
    } else if ( nRows == 9 && nColumns == 9 ) {
        // OOFEM:           11, 22, 33, 23, 13, 12, 32, 31, 21
        // UMAT:            11, 22, 33, 12, 13, 23, 32, 21, 31
        const int abq2oo [ 9 ] = {
            1,  2,  3,  6,  5,  4,  7,  9,  8
        };

        FloatMatrix tmp(9, 9);
        for ( int i = 1; i <= 9; i++ ) {
            for ( int j = 1; j <= 9; j++ ) {
                tmp.at(i, j) = this->at(abq2oo [ i - 1 ], abq2oo [ j - 1 ]);
            }
        }

        * this = tmp;
    }
}



double FloatMatrix :: computeReciprocalCondition(char p) const
{
#  ifdef DEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR("receiver must be square (is %d by %d)", this->nRows, this->nColumns);
    }
#  endif
    double anorm = this->computeNorm(p);

#  ifdef __LAPACK_MODULE
    int n = this->nRows;
    FloatArray work(4 *n);
    IntArray iwork(n);
    int info;
    double rcond;

    if ( n > 3 ) { // Only use this routine for larger matrices. (not sure if it's suitable for n = 3) // Mikael
        FloatMatrix a_cpy = * this;
        dgetrf_(& n, & n, a_cpy.givePointer(), & n, iwork.givePointer(), & info);
        if ( info < 0 ) {
            OOFEM_ERROR("dgetfr error %d\n", info);
        }
        dgecon_(& p, & ( this->nRows ), a_cpy.givePointer(), & this->nRows, & anorm, & rcond, work.givePointer(), iwork.givePointer(), & info, 1);
        if ( info < 0 ) {
            OOFEM_ERROR("dgecon error %d\n", info);
        }
        return rcond;
    }
#  endif
    if ( this->giveDeterminant() <= 1e-6 * anorm ) {
        return 0.0;
    }
    FloatMatrix inv;
    inv.beInverseOf(* this);
    return 1.0 / ( inv.computeNorm(p) * anorm );
}

void FloatMatrix :: beMatrixFormOfStress(const FloatArray &aArray)
{
    // Revrites the  matrix on vector form (symmetrized matrix used), order: 11, 22, 33, 23, 13, 12
#  ifdef DEBUG
    if ( aArray.giveSize() != 6 && aArray.giveSize() != 9 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }
#  endif
    this->resize(3, 3);
    if ( aArray.giveSize() == 9 ) {
        this->at(1, 1) = aArray.at(1);
        this->at(2, 2) = aArray.at(2);
        this->at(3, 3) = aArray.at(3);
        this->at(2, 3) = aArray.at(4);
        this->at(1, 3) = aArray.at(5);
        this->at(1, 2) = aArray.at(6);
        this->at(3, 2) = aArray.at(7);
        this->at(3, 1) = aArray.at(8);
        this->at(2, 1) = aArray.at(9);
    } else if ( aArray.giveSize() == 6 ) {
        this->at(1, 1) = aArray.at(1);
        this->at(2, 2) = aArray.at(2);
        this->at(3, 3) = aArray.at(3);
        this->at(2, 3) = aArray.at(4);
        this->at(1, 3) = aArray.at(5);
        this->at(1, 2) = aArray.at(6);
        this->at(3, 2) = aArray.at(4);
        this->at(3, 1) = aArray.at(5);
        this->at(2, 1) = aArray.at(6);
    }
}

#if 0
bool FloatMatrix :: computeEigenValuesSymmetric(FloatArray &lambda, FloatMatrix &v, int neigs) const
{
 #ifdef __LAPACK_MODULE
    double abstol = 1.0; ///@todo Find suitable tolerance.
    int lda, n, ldz, info, found, lwork;
    n = this->nRows;
    lda = n;
    ldz = n;
    FloatMatrix a;
    a = * this;
    if ( neigs == 0 ) {
        neigs = n;
    }

    lambda.resize(neigs);
    v.resize(n, neigs);

    IntArray ifail(n), iwork(5 * n);
    FloatArray work( ( n + 3 ) *n ); // maximum block size should be less than n (?)
    lwork = ( n + 3 ) * n;
    if ( neigs > 0 ) {
        int one = 1;
        dsyevx_("N", "I", "U",
                & n, a.givePointer(), & n,
                NULL, NULL, & one, & neigs,
                & abstol, & found, lambda.givePointer(), v.givePointer(), & ldz,
                work.givePointer(), & lwork, iwork.givePointer(), ifail.givePointer(), & info, 1, 1, 1);
    } else {
        dsyevx_("N", "A", "U",
                & n, a.givePointer(), & n,
                NULL, NULL, NULL, NULL,
                & abstol, & found, lambda.givePointer(), v.givePointer(), & ldz,
                work.givePointer(), & lwork, iwork.givePointer(), ifail.givePointer(), & info, 1, 1, 1);
    }

    return info == 0;

 #else
    OOFEM_ERROR("Requires LAPACK");
    return false;

 #endif
}
#endif

contextIOResultType FloatMatrix :: storeYourself(DataStream &stream) const
// writes receiver's binary image into stream
// use id to distinguish some instances
// return value >0 success
//              =0 file i/o error
{
    // write size
    if ( !stream.write(nRows) ) {
        return ( CIO_IOERR );
    }

    if ( !stream.write(nColumns) ) {
        return ( CIO_IOERR );
    }

    // write raw data
    if ( !stream.write(this->givePointer(), nRows * nColumns) ) {
        return ( CIO_IOERR );
    }

    // return result back
    return CIO_OK;
}


contextIOResultType FloatMatrix :: restoreYourself(DataStream &stream)
// reads receiver from stream
// warning - overwrites existing data!
// returns 0 if file i/o error
//        -1 if id of class id is not correct
{
    // read size
    if ( !stream.read(nRows) ) {
        return ( CIO_IOERR );
    }

    if ( !stream.read(nColumns) ) {
        return ( CIO_IOERR );
    }

    this->values.resize(nRows * nColumns);

    // read raw data
    if ( !stream.read(this->givePointer(), nRows * nColumns) ) {
        return ( CIO_IOERR );
    }

    // return result back
    return CIO_OK;
}


int
FloatMatrix :: givePackSize(DataStream &buff) const
{
    return buff.givePackSizeOfInt(1) + buff.givePackSizeOfInt(1) +
           buff.givePackSizeOfDouble(nRows * nColumns);
}


bool FloatMatrix :: jaco_(FloatArray &eval, FloatMatrix &v, int nf)
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
    double ssum, aa, co, si, tt, tol, sum, aij, aji;
    int ite, i, j, k, ih;
    int neq = this->giveNumberOfRows();

    double c_b2 = .10;
    //double c_b27 = .01;

    /* Function Body */
#ifdef DEBUG
    if ( !isSquare() ) {
        OOFEM_ERROR("Not square matrix");
    }
    // check for symmetry
    for ( i = 1; i <= neq; i++ ) {
        for ( j = i + 1; j <= neq; j++ ) {
            //if ( this->at(i, j) != this->at(j, i) ) {
            if ( fabs( this->at(i, j) - this->at(j, i) ) > 1.0e-6 ) {
                OOFEM_ERROR("Not Symmetric matrix");
            }
        }
    }

#endif

    eval.resize(neq);
    v.resize(neq, neq);

    for ( i = 1; i <= neq; i++ ) {
        eval.at(i) = this->at(i, i);
    }

    tol = pow(c_b2, nf);
    sum = 0.0;
    for ( i = 1; i <= neq; ++i ) {
        for ( j = 1; j <= neq; ++j ) {
            sum += fabs( this->at(i, j) );
            v.at(i, j) = 0.0;
        }

        v.at(i, i) = 1.0;
    }

    if ( sum <= 0.0 ) {
        return 0;
    }


    /* ---- REDUCE MATRIX TO DIAGONAL ---------------- */
    ite = 0;
    do {
        ssum = 0.0;
        for ( j = 2; j <= neq; ++j ) {
            ih = j - 1;
            for ( i = 1; i <= ih; ++i ) {
                if ( ( fabs( this->at(i, j) ) / sum ) > tol ) {
                    ssum += fabs( this->at(i, j) );
                    /* ---- CALCULATE ROTATION ANGLE ----------------- */
                    aa = atan2( this->at(i, j) * 2.0, eval.at(i) - eval.at(j) ) /  2.0;
                    si = sin(aa);
                    co = cos(aa);
                    /*
                     *   // ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
                     *   for (k = 1; k <= neq; ++k) {
                     *    tt = this->at(k, i);
                     *    this->at(k, i) = co * tt + si * this->at(k, j);
                     *    this->at(k, j) = -si * tt + co * this->at(k, j);
                     *    tt = v.at(k, i);
                     *    v.at(k, i) = co * tt + si * v.at(k, j);
                     *    // L500:
                     *    v.at(k, j) = -si * tt + co * v.at(k, j);
                     *   }
                     *   // ---- MODIFY DIAGONAL TERMS --------------------
                     *   this->at(i, i) = co * this->at(i, i) + si * this->at(j, i);
                     *   this->at(j, j) = -si * this->at(i, j) + co * this->at(j, j);
                     *   this->at(i, j) = 0.0;
                     *   // ---- MAKE "A" MATRIX SYMMETRICAL --------------
                     *   for (k = 1; k <= neq; ++k) {
                     *    this->at(i, k) = this->at(k, i);
                     *    this->at(j, k) = this->at(k, j);
                     *    // L600:
                     *   }
                     */
                    // ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
                    for ( k = 1; k < i; ++k ) {
                        tt = this->at(k, i);
                        this->at(k, i) = co * tt + si *this->at(k, j);
                        this->at(k, j) = -si * tt + co *this->at(k, j);
                        tt = v.at(k, i);
                        v.at(k, i) = co * tt + si *v.at(k, j);
                        v.at(k, j) = -si * tt + co *v.at(k, j);
                    }

                    // diagonal term (i,i)
                    tt = eval.at(i);
                    eval.at(i) = co * tt + si *this->at(i, j);
                    aij = -si * tt + co *this->at(i, j);
                    tt = v.at(i, i);
                    v.at(i, i) = co * tt + si *v.at(i, j);
                    v.at(i, j) = -si * tt + co *v.at(i, j);

                    for ( k = i + 1; k < j; ++k ) {
                        tt = this->at(i, k);
                        this->at(i, k) = co * tt + si *this->at(k, j);
                        this->at(k, j) = -si * tt + co *this->at(k, j);
                        tt = v.at(k, i);
                        v.at(k, i) = co * tt + si *v.at(k, j);
                        v.at(k, j) = -si * tt + co *v.at(k, j);
                    }

                    // diagonal term (j,j)
                    tt = this->at(i, j);
                    aji = co * tt + si *eval.at(j);
                    eval.at(j) = -si * tt + co *eval.at(j);

                    tt = v.at(j, i);
                    v.at(j, i) = co * tt + si *v.at(j, j);
                    v.at(j, j) = -si * tt + co *v.at(j, j);
                    //
                    for ( k = j + 1; k <= neq; ++k ) {
                        tt = this->at(i, k);
                        this->at(i, k) = co * tt + si *this->at(j, k);
                        this->at(j, k) = -si * tt + co *this->at(j, k);
                        tt = v.at(k, i);
                        v.at(k, i) = co * tt + si *v.at(k, j);
                        v.at(k, j) = -si * tt + co *v.at(k, j);
                    }

                    // ---- MODIFY DIAGONAL TERMS --------------------
                    eval.at(i) = co * eval.at(i) + si * aji;
                    eval.at(j) = -si * aij + co *eval.at(j);
                    this->at(i, j) = 0.0;
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
    for ( i = 1; i <= neq; i++ ) {
        for ( j = i; j <= neq; j++ ) {
            this->at(i, j) = this->at(j, i);
        }
    }

    return 0;
} /* jaco_ */


#ifdef BOOST_PYTHON
void
FloatMatrix :: __setitem__(boost :: python :: api :: object t, double val)
{
    this->at(boost :: python :: extract< int >(t [ 0 ]) + 1, boost :: python :: extract< int >(t [ 1 ]) + 1) = val;
}

double
FloatMatrix :: __getitem__(boost :: python :: api :: object t)
{
    return this->at(boost :: python :: extract< int >(t [ 0 ]) + 1, boost :: python :: extract< int >(t [ 1 ]) + 1);
}
#endif

std :: ostream &operator << ( std :: ostream & out, const FloatMatrix & x )
{
    out << x.nRows << " " << x.nColumns << " {";
    for ( int i = 0; i < x.nRows; ++i ) {
        for ( int j = 0; j < x.nColumns; ++j ) {
            out << " " << x(i, j);
        }
        out << ";";
    }
    out << "}";
    return out;
}
} // end namespace oofem
