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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

#define __MATH_INTERNAL
#include "ops-mat.h"


namespace oofem{
    void FloatMatrix::resize_funny(int nr, int nc){
        this->nRows = nr; this->nColumns = nc; \
        std::size_t nsize = this->nRows * this->nColumns; \
        if ( nsize < this->values.size() ) { \
            this->values.resize(nsize); \
        } else if ( nsize > this->values.size() ) { \
            this->values.assign(nsize, 0.); \
        } \
    }

}

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
    resize_funny( mat.size(), mat.begin()->size() );
    auto p = this->values.begin();
    for ( auto col : mat ) {
#if DEBUG
        if ( this->nRows != col.size() ) {
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
    resize_funny( ( int ) mat.begin()->size(), ( int ) mat.size() );
    auto p = this->values.begin();
    for ( auto col : mat ) {
#if DEBUG
        if ( this->nRows != col.size() ) {
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
    resize_funny( mat.begin()->giveSize(), ( int ) mat.size() );
    auto p = this->values.begin();
    for ( auto col : mat ) {
#if DEBUG
        if ( this->nRows != col.size() ) {
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


void FloatMatrix :: checkBounds(std::size_t i, std::size_t j) const
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
		printf("APA \n");
        OOFEM_ERROR("matrix error on columns : %d > %d", j, nColumns);
    }
}

bool FloatMatrix :: isAllFinite() const
{
    for (double val : values) {
        if ( !std::isfinite(val) ) {
            return false;
        }
    }

    return true;
}

void FloatMatrix :: assemble(const FloatMatrix &src, const IntArray &loc)
{
    mat::assemble(*this,src,loc);
}


void FloatMatrix :: assemble(const FloatMatrix &src, const IntArray &rowind, const IntArray &colind)
{
    mat::assemble(*this,src,rowind,colind);
}

void FloatMatrix :: assembleT(const FloatMatrix &src, const IntArray &rowind, const IntArray &colind)
{
    mat::assembleT(*this,src,rowind,colind);
}


void FloatMatrix :: assemble(const FloatMatrix &src, const int *rowind, const int *colind)
{
    mat::assemble(*this,src,rowind,colind);
}


void FloatMatrix :: beTranspositionOf(const FloatMatrix &src)
{
    // receiver becomes a transposition of src
    int nrows = src.giveNumberOfColumns(), ncols = src.giveNumberOfRows();
    resize_funny(nrows, ncols);

    for ( int i = 1; i <= nrows; i++ ) {
        for ( int j = 1; j <= ncols; j++ ) {
            this->at(i, j) = src.at(j, i);
        }
    }
}


void FloatMatrix :: beProductOf(const FloatMatrix &aMatrix, const FloatMatrix &bMatrix)
// Receiver = aMatrix * bMatrix
{
#  ifndef NDEBUG
    if ( aMatrix.nColumns != bMatrix.nRows ) {
        OOFEM_ERROR("error in product A*B : dimensions do not match, A(*,%d), B(%d,*)", aMatrix.nColumns, bMatrix.nRows);
    }
#  endif
    resize_funny(aMatrix.nRows, bMatrix.nColumns);
#  ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    dgemm_("n", "n", & this->nRows, & this->nColumns, & aMatrix.nColumns,
           & alpha, aMatrix.givePointer(), & aMatrix.nRows, bMatrix.givePointer(), & bMatrix.nRows,
           & beta, this->givePointer(), & this->nRows,
           aMatrix.nColumns, bMatrix.nColumns, this->nColumns);
#  else
    for ( std::size_t i = 1; i <= aMatrix.nRows; i++ ) {
        for ( std::size_t j = 1; j <= bMatrix.nColumns; j++ ) {
            double coeff = 0.;
            for ( std::size_t k = 1; k <= aMatrix.nColumns; k++ ) {
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
#  ifndef NDEBUG
    if ( aMatrix.nRows != bMatrix.nRows ) {
        OOFEM_ERROR("error in product A*B : dimensions do not match");
    }
#  endif
    resize_funny(aMatrix.nColumns, bMatrix.nColumns);
#  ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    dgemm_("t", "n", & this->nRows, & this->nColumns, & aMatrix.nRows,
           & alpha, aMatrix.givePointer(), & aMatrix.nRows, bMatrix.givePointer(), & bMatrix.nRows,
           & beta, this->givePointer(), & this->nRows,
           aMatrix.nColumns, bMatrix.nColumns, this->nColumns);
#  else
    for ( std::size_t i = 1; i <= aMatrix.nColumns; i++ ) {
        for (std::size_t j = 1; j <= bMatrix.nColumns; j++ ) {
            double coeff = 0.;
            for (std::size_t k = 1; k <= aMatrix.nRows; k++ ) {
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
#  ifndef NDEBUG
    if ( aMatrix.nColumns != bMatrix.nColumns ) {
        OOFEM_ERROR("error in product A*B : dimensions do not match");
    }
#  endif
    resize_funny(aMatrix.nRows, bMatrix.nRows);
#  ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    dgemm_("n", "t", & this->nRows, & this->nColumns, & aMatrix.nColumns,
           & alpha, aMatrix.givePointer(), & aMatrix.nRows, bMatrix.givePointer(), & bMatrix.nRows,
           & beta, this->givePointer(), & this->nRows,
           aMatrix.nColumns, bMatrix.nColumns, this->nColumns);
#  else
    for (std::size_t i = 1; i <= aMatrix.nRows; i++ ) {
        for (std::size_t j = 1; j <= bMatrix.nRows; j++ ) {
            double coeff = 0.;
            for (std::size_t k = 1; k <= aMatrix.nColumns; k++ ) {
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
#  ifndef NDEBUG
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
    for (std::size_t i = 1; i <= aMatrix.nRows; i++ ) {
        for (std::size_t j = 1; j <= bMatrix.nColumns; j++ ) {
            double coeff = 0.;
            for (std::size_t k = 1; k <= aMatrix.nColumns; k++ ) {
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
#  ifndef NDEBUG
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
    for (std::size_t i = 1; i <= aMatrix.nColumns; i++ ) {
        for (std::size_t j = 1; j <= bMatrix.nColumns; j++ ) {
            double coeff = 0.;
            for (std::size_t k = 1; k <= aMatrix.nRows; k++ ) {
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
    resize_funny(n1, n2);
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
    *this=mat::beLocalCoordSys(normal);
}

void FloatMatrix :: setSubMatrix(const FloatMatrix &src, int sr, int sc)
{
    sr--;
    sc--;

    int srcRows = src.giveNumberOfRows(), srcCols = src.giveNumberOfColumns();
#ifndef NDEBUG
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
#ifndef NDEBUG
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
#ifndef NDEBUG
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
#ifndef NDEBUG
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
{
    mat::plusProductSymmUpper(*this,a,b,dV);
}


void FloatMatrix :: plusDyadSymmUpper(const FloatArray &a, double dV)
{
    mat::plusDyadSymmUpper(*this,a,dV);
}



void FloatMatrix :: plusProductUnsym(const FloatMatrix &a, const FloatMatrix &b, double dV)
{
    mat::plusProductUnsym(*this,a,b,dV);
}


void FloatMatrix :: plusDyadUnsym(const FloatArray &a, const FloatArray &b, double dV)
{
    mat::plusDyadUnsym(*this,a,b,dV);
}


bool FloatMatrix :: beInverseOf(const FloatMatrix &src)
{
    return mat::beInverseOf(*this,src);
}


void FloatMatrix :: beSubMatrixOf(const FloatMatrix &src,
                                  std::size_t topRow, std::size_t bottomRow, std::size_t topCol, std::size_t bottomCol)
/*
 * modifies receiver to be  submatrix of the src matrix
 * size of receiver  submatrix is determined from
 * input parameters
 */
{
#ifndef NDEBUG
    if ( ( topRow < 1 ) || ( bottomRow < 1 ) || ( topCol < 1 ) || ( bottomCol < 1 ) ) {
        OOFEM_ERROR("subindexes size mismatch");
    }

    if ( ( src.nRows < bottomRow ) || ( src.nColumns < bottomCol ) || ( ( bottomRow - topRow ) > src.nRows ) ||
         ( ( bottomCol - topCol ) > src.nColumns ) ) {
        OOFEM_ERROR("subindexes size mismatch");
    }
#endif


    std::size_t topRm1, topCm1;
    topRm1 = topRow - 1;
    topCm1 = topCol - 1;

    // allocate return value
    this->resize(bottomRow - topRm1, bottomCol - topCm1);
    for (std::size_t i = topRow; i <= bottomRow; i++ ) {
        for (std::size_t j = topCol; j <= bottomCol; j++ ) {
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
#  ifndef NDEBUG
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
#     ifndef NDEBUG
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
#     ifndef NDEBUG
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
#     ifndef NDEBUG
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
{
    return mat::solveForRhs(*this,b,answer,transpose);
}


bool FloatMatrix :: solveForRhs(const FloatMatrix &b, FloatMatrix &answer, bool transpose)
{
    return mat::solveForRhs(*this,b,answer,transpose);
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
#ifndef NDEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR("cannot make unit matrix of %d by %d matrix", nRows, nColumns);
    }
#endif

    this->zero();
    for (std::size_t i = 1; i <= nRows; i++ ) {
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


void FloatMatrix :: resize(std::size_t rows, std::size_t columns)
//
// resizes receiver, all data will be lost
//
{
    this->nRows = rows;
    this->nColumns = columns;
    this->values.assign(rows * columns, 0.);
}


void FloatMatrix :: resizeWithData(std::size_t rows, std::size_t columns)
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

    std::size_t ii = min( rows, old.giveRowSize() );
    std::size_t jj = min( columns, old.giveColSize() );
    // copy old values if possible
    for (std::size_t i = 1; i <= ii; i++ ) {
        for (std::size_t j = 1; j <= jj; j++ ) {
            this->at(i, j) = old.at(i, j);
        }
    }
}



double FloatMatrix :: giveDeterminant() const
// Returns the determinant of the receiver.
{
#  ifndef NDEBUG
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

    // return 0.;
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
#  ifndef NDEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR("cannot compute the trace of a non-square %d by %d matrix", nRows, nColumns);
    }
#  endif
    double answer = 0.;
    for (std::size_t k = 0; k < nRows; k++ ) {
        answer += values [ k * ( nRows + 1 ) ];
    }
    return answer;
}


void FloatMatrix :: printYourself() const
{
    mat::printYourself(*this);
}


void FloatMatrix :: printYourselfToFile(const std::string filename, const bool showDimensions) const
{
    mat::printYourselfToFile(*this,filename,showDimensions);
}


void FloatMatrix :: printYourself(const std::string &name) const
{
    mat::printYourself(*this,name);
}


void FloatMatrix :: pY() const
{
    mat::pY(*this);
}


void FloatMatrix :: writeCSV(const std :: string &name) const
{
    mat::writeCSV(*this,name);
}


void FloatMatrix :: rotatedWith(const FloatMatrix &r, char mode)
{
    mat::rotatedWith(*this,r,mode);
}



void FloatMatrix :: symmetrized()
// Initializes the lower half of the receiver to the upper half.
{
#  ifndef NDEBUG
    if ( nRows != nColumns ) {
        OOFEM_ERROR("cannot symmetrize a non-square matrix");
    }

#   endif

    for (std::size_t i = 2; i <= nRows; i++ ) {
        for (std::size_t j = 1; j < i; j++ ) {
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
    return mat::computeNorm(*this,p);
}



void FloatMatrix :: beMatrixForm(const FloatArray &aArray)
{
    *this=mat::beMatrixForm(aArray);
}

void FloatMatrix :: changeComponentOrder()
{
    mat::changeComponentOrder(*this);
}



double FloatMatrix :: computeReciprocalCondition(char p) const
{
    return mat::computeReciprocalCondition(*this,p);
}

void FloatMatrix :: beMatrixFormOfStress(const FloatArray &aArray)
{
    // Revrites the  matrix on vector form (symmetrized matrix used), order: 11, 22, 33, 23, 13, 12
#  ifndef NDEBUG
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
    return buff.givePackSizeOfSizet(1) + buff.givePackSizeOfSizet(1) +
           buff.givePackSizeOfDouble(nRows * nColumns);
}


bool FloatMatrix :: jaco_(FloatArray &eval, FloatMatrix &v, int nf)
{
    return mat::jaco_(*this,eval,v,nf);
}



std :: ostream &operator << ( std :: ostream & out, const FloatMatrix & x )
{
    out << x.nRows << " " << x.nColumns << " {";
    for (std::size_t i = 0; i < x.nRows; ++i ) {
        for (std::size_t j = 0; j < x.nColumns; ++j ) {
            out << " " << x(i, j);
        }
        out << ";";
    }
    out << "}";
    return out;
}

FloatMatrix &operator *= ( FloatMatrix & x, const double & a ) {x.times(a); return x;}
FloatMatrix operator *( const FloatMatrix & a, const FloatMatrix & b ) {FloatMatrix ans; ans.beProductOf (a,b); return ans;}
FloatArray operator *( const FloatMatrix & a, const FloatArray & b ) {FloatArray ans; ans.beProductOf (a,b); return ans;}
FloatMatrix operator +( const FloatMatrix & a, const FloatMatrix & b ) {FloatMatrix ans(a); ans.add(b); return ans;}
FloatMatrix operator -( const FloatMatrix & a, const FloatMatrix & b ) {FloatMatrix ans(a); ans.subtract(b); return ans;}
FloatMatrix &operator += ( FloatMatrix & a, const FloatMatrix & b ) {a.add(b); return a;}
FloatMatrix &operator -= ( FloatMatrix & a, const FloatMatrix & b ) {a.subtract(b); return a;}






} // end namespace oofem
