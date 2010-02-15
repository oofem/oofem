/* $Header: /home/cvs/bp/oofem/oofemlib/src/flotmtrx.C,v 1.22.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */


//   file FLOTMTRX.CC

#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"
#include "error.h"
#include "datastream.h"
#include "classtype.h"

#ifndef __MAKEDEPEND
#include <math.h>
#endif

#ifdef __PARALLEL_MODE
#include "combuff.h"
#endif

namespace oofem {

FloatMatrix :: FloatMatrix(FloatArray *vector, int transpose)
//
// constructor : creates (vector->giveSize(),1) FloatMatrix
// if transpose = 1 creates (1,vector->giveSize()) FloatMatrix
//
{
    if ( transpose ) {
        nRows = 1; // column vector
        nColumns = vector->giveSize();
    } else {
        nRows = vector->giveSize(); // row vector- default
        nColumns = 1;
    }

    allocatedSize = nRows * nColumns;
    values = allocDouble(allocatedSize);

    if ( transpose ) {
        for ( int i = 1; i <= nColumns; i++ ) {
            this->at(1, i) = vector->at(i);
        }
    } else {
        for ( int i = 1; i <= nRows; i++ ) {
            this->at(i, 1) = vector->at(i);
        }
    }
}


FloatMatrix ::
FloatMatrix(const FloatMatrix &src) : Matrix(src.nRows, src.nColumns)
{
    // copy constructor
    double *P1, *P2;
    int i;

    //this->nRows = src.nRows;
    //this->nColumns = src.nColumns;

    allocatedSize = nRows * nColumns;
    values = allocDouble(allocatedSize);

    P1 = values;
    P2 = src.values;
    for ( i = 0; i < nRows * nColumns; i++ ) {
        P1 [ i ] = P2 [ i ];
    }
}

FloatMatrix &
FloatMatrix :: operator=(const FloatMatrix &src)
{
    // assignment: cleanup and copy
    double *P1, *P2;
    int i;
    this->resize(src.nRows, src.nColumns);

    P1 = values;
    P2 = src.values;
    for ( i = 0; i < nRows * nColumns; i++ ) {
        P1 [ i ] = P2 [ i ];
    }

    return * this;
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

#endif


double &FloatMatrix :: operator()(int i, int j)
{
#ifdef DEBUG
    assert(0 <= i && i < nRows);
    assert(0 <= j && j < nColumns);
#endif
    return values [ ( j ) * nRows + i ];
}

double FloatMatrix :: operator()(int i, int j) const
{
#ifdef DEBUG
    assert(0 <= i && i < nRows);
    assert(0 <= j && j < nColumns);
#endif
    return values [ ( j ) * nRows + i ];
}




void
FloatMatrix :: assemble(const FloatMatrix &src, const IntArray &loc)
{
    int i, j, ii, jj, size;

    if ( ( size = src.giveNumberOfRows() ) != loc.giveSize() ) {
        OOFEM_ERROR("FloatMatrix :: assemble : dimensions of 'src' and 'loc' mismatch");
    }

    if ( !src.isSquare() ) {
        OOFEM_ERROR("FloatMatrix :: assemble : 'src' is not sqaure matrix");
    }

    for ( i = 1; i <= size; i++ ) {
        if ( ( ii = loc.at(i) ) ) {
            for ( j = 1; j <= size; j++ ) {
                if ( ( jj = loc.at(j) ) ) {
                    this->at(ii, jj) += src.at(i, j);
                }
            }
        }
    }
}

void
FloatMatrix :: assemble(const FloatMatrix &src, const IntArray &rowind, const IntArray &colind)
{
    int i, j, ii, jj;
    int nr = src.giveNumberOfRows();
    int nc = src.giveNumberOfColumns();

    if ( nr != rowind.giveSize() ) {
        OOFEM_ERROR("FloatMatrix :: assemble : row dimensions of 'src' and 'rowind' mismatch");
    }

    if ( nc != colind.giveSize() ) {
        OOFEM_ERROR("FloatMatrix :: assemble : column dimensions of 'src' and 'colind' mismatch");
    }

    for ( i = 1; i <= nr; i++ ) {
        if ( ( ii = rowind.at(i) ) ) {
            for ( j = 1; j <= nc; j++ ) {
                if ( ( jj = colind.at(j) ) ) {
                    this->at(ii, jj) += src.at(i, j);
                }
            }
        }
    }
}

void
FloatMatrix :: assemble(const FloatMatrix &src, const int *rowind, const int *colind)
{
    int i, j, ii, jj;
    int nr = src.giveNumberOfRows();
    int nc = src.giveNumberOfColumns();

    for ( i = 1; i <= nr; i++ ) {
        if ( ( ii = rowind [ i - 1 ] ) ) {
            for ( j = 1; j <= nc; j++ ) {
                if ( ( jj = colind [ j - 1 ] ) ) {
                    this->at(ii, jj) += src.at(i, j);
                }
            }
        }
    }
}


void
FloatMatrix :: beTranspositionOf(const FloatMatrix &src)
{
    // receiver becomes a transposition of src
    int i, j, nrows = src.giveNumberOfColumns(), ncols = src.giveNumberOfRows();
    this->resize(nrows, ncols);

    for ( i = 1; i <= nrows; i++ ) {
        for ( j = 1; j <= ncols; j++ ) {
            this->at(i, j) = src.at(j, i);
        }
    }
}


void
FloatMatrix :: beProductOf(const FloatMatrix &aMatrix, const FloatMatrix &bMatrix)
// Receiver = aMatrix * bMatrix
{
    int i, j, k, p;
    double coeff;
    // FloatMatrix* answer ;

#  ifdef DEBUG
    if ( aMatrix.nColumns != bMatrix.nRows ) {
        OOFEM_ERROR("FloatMatrix::beProductOf : error in product A*B : dimensions do not match");
    }

#  endif

    p      = bMatrix.nColumns;
    this->resize(aMatrix.nRows, p);
    // answer = new FloatMatrix(nRows,p) ;
    for ( i = 1; i <= aMatrix.nRows; i++ ) {
        for ( j = 1; j <= p; j++ ) {
            coeff = 0.;
            for ( k = 1; k <= aMatrix.nColumns; k++ ) {
                coeff += aMatrix.at(i, k) * bMatrix.at(k, j);
            }

            this->at(i, j) = coeff;
        }
    }

    return;
}

void
FloatMatrix :: beDyadicProductOf(const FloatArray &vec1, const FloatArray &vec2)
// Receiver = vec1 * vec2^T
{
  int i,j;
  int n1 = vec1.giveSize();
  int n2 = vec2.giveSize();
  this->resize(n1,n2);
  for ( i = 1; i <= n1; i++ )
    for ( j = 1; j <= n2; j++ )
      this->at(i,j) = vec1.at(i) * vec2.at(j);
  return;
}


void
FloatMatrix :: beTProductOf(const FloatMatrix &aMatrix, const FloatMatrix &bMatrix)
// Receiver = aMatrix^T * bMatrix
{
    int i, j, k, p;
    double coeff;
    // FloatMatrix* answer ;

#  ifdef DEBUG
    if ( aMatrix.nRows != bMatrix.nRows ) {
        OOFEM_ERROR("FloatMatrix::beTProductOf : error in product A*B : dimensions do not match");
    }

#  endif

    p      = bMatrix.nColumns;
    this->resize(aMatrix.nColumns, p);
    // answer = new FloatMatrix(nColumns,p) ;
    for ( i = 1; i <= aMatrix.nColumns; i++ ) {
        for ( j = 1; j <= p; j++ ) {
            coeff = 0.;
            for ( k = 1; k <= aMatrix.nRows; k++ ) {
                coeff += aMatrix.at(k, i) * bMatrix.at(k, j);
            }

            this->at(i, j) = coeff;
        }
    }

    return;
}



void
FloatMatrix :: beProductTOf(const FloatMatrix &aMatrix, const FloatMatrix &bMatrix)
// Receiver = aMatrix * bMatrix^T
{
    int i, j, k, p;
    double coeff;
    // FloatMatrix* answer ;

#  ifdef DEBUG
    if ( aMatrix.nColumns != bMatrix.nColumns ) {
        OOFEM_ERROR("FloatMatrix::beProductTOf : error in product A*B : dimensions do not match");
    }

#  endif

    p      = bMatrix.nRows;
    this->resize(aMatrix.nRows, p);
    // answer = new FloatMatrix(nRows,p) ;
    for ( i = 1; i <= aMatrix.nRows; i++ ) {
        for ( j = 1; j <= p; j++ ) {
            coeff = 0.;
            for ( k = 1; k <= aMatrix.nColumns; k++ ) {
                coeff += aMatrix.at(i, k) * bMatrix.at(j, k);
            }

            this->at(i, j) = coeff;
        }
    }

    return;
}



void
FloatMatrix :: addSubMatrix(const FloatMatrix &src, int sr, int sc)
{
    sr--;
    sc--;

    int srcRows = src.giveNumberOfRows(), srcCols = src.giveNumberOfColumns();

    int nr = sr + srcRows;
    int nc = sc + srcCols;

    if ( ( this->giveNumberOfRows() < nr ) || ( this->giveNumberOfColumns() < nc ) ) {
        this->resizeWithData( max(this->giveNumberOfRows(), nr), max(this->giveNumberOfColumns(), nc) );
    }

    // add sub-matrix
    int i, j;
    for ( i = 1; i <= srcRows; i++ ) {
        for ( j = 1; j <= srcCols; j++ ) {
            this->at(sr + i, sc + j) += src.at(i, j);
        }
    }
}



void
FloatMatrix :: beSubMatrixOf(const FloatMatrix &src,
                             int topRow, int bottomRow, int topCol, int bottomCol)
/*
 * modifies receiver to be  submatrix of the src matrix
 * size of receiver  submatrix is determined from
 * input parametrs
 */
{
#ifdef DEBUG
    // check input params
    if ( ( topRow < 1 ) || ( bottomRow < 1 ) || ( topCol < 1 ) || ( bottomCol < 1 ) ) {
        OOFEM_ERROR("FloatMatrix::beSubMatrixOf : subindexes size mismatch");
    }

    if ( ( src.nRows < bottomRow ) || ( src.nColumns < bottomCol ) || ( ( bottomRow - topRow ) > src.nRows ) ||
        ( ( bottomCol - topCol ) > src.nColumns ) ) {
        OOFEM_ERROR("FloatMatrix::beSubMatrixOf : subindexes size mismatch");
    }

#endif


    int i, j, topRm1, topCm1;
    topRm1 = topRow - 1;
    topCm1 = topCol - 1;

    // allocate return value
    this->resize(bottomRow - topRm1, bottomCol - topCm1);
    // FloatMatrix * subM = new FloatMatrix (bottomRow-topRm1,bottomCol-topCm1);
    for ( i = topRow; i <= bottomRow; i++ ) {
        for ( j = topCol; j <= bottomCol; j++ ) {
            this->at(i - topRm1, j - topCm1) = src.at(i, j);
        }
    }

    return;
}

void
FloatMatrix :: addSubVectorRow(const FloatArray &src, int sr, int sc)
{
    sc--;

    int srcCols = src.giveSize();

    int nr = sr;
    int nc = sc + srcCols;

    if ( ( this->giveNumberOfRows() < nr ) || ( this->giveNumberOfColumns() < nc ) ) {
        this->resizeWithData( max(this->giveNumberOfRows(), nr), max(this->giveNumberOfColumns(), nc) );
    }

    // add sub-matrix
    int j;
    for ( j = 1; j <= srcCols; j++ ) {
        this->at(sr, sc + j) += src.at(j);
    }
}

void
FloatMatrix :: copySubVectorRow(const FloatArray &src, int sr, int sc)
{
    sc--;

    int srcCols = src.giveSize();

    int nr = sr;
    int nc = sc + srcCols;

    if ( ( this->giveNumberOfRows() < nr ) || ( this->giveNumberOfColumns() < nc ) ) {
        this->resizeWithData( max(this->giveNumberOfRows(), nr), max(this->giveNumberOfColumns(), nc) );
    }

    // add sub-matrix
    int j;
    for ( j = 1; j <= srcCols; j++ ) {
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
    int i, j, k;
    double summ;

    if ( !this->isNotEmpty() ) {
        resize(a.nColumns, b.nColumns);
        zero();
    }

    for ( i = 1; i <= nRows; i++ ) {
        for ( j = i; j <= nColumns; j++ ) {
            summ = 0.;
            for ( k = 1; k <= a.nRows; k++ ) {
                summ += a.at(k, i) * b.at(k, j);
            }

            this->at(i, j) += summ * dV;
        }
    }
}



void FloatMatrix :: plusProductUnsym(const FloatMatrix &a, const FloatMatrix &b, double dV)
// Adds to the receiver the product  a(transposed).b dV .
// If the receiver has a null size, it is expanded.
// Advantage : does not compute the transposition of matrix a.
{
    int i, j, k;
    double summ;

    if ( !this->isNotEmpty() ) {
        resize(a.nColumns, b.nColumns);
        zero();
    }

    for ( i = 1; i <= nRows; i++ ) {
        for ( j = 1; j <= nColumns; j++ ) {
            summ = 0.;
            for ( k = 1; k <= a.nRows; k++ ) {
                summ += a.at(k, i) * b.at(k, j);
            }

            this->at(i, j) += summ * dV;
        }
    }
}


void
FloatMatrix :: beInverseOf(const FloatMatrix &src)
// Receiver becomes inverse of given parameter src. If necessary, size is adjusted.

{
    //FloatMatrix* answer ;
    double det;

#  ifdef DEBUG
    if ( !src.isSquare() ) {
        OOFEM_ERROR3("FloatMatrix::beInverseOf : cannot inverse a %d by %d matrix", src.nRows, src.nColumns);
    }

#  endif

    this->resize(src.nRows, src.nColumns);

    if ( nRows == 1 ) {
        this->at(1, 1) = 1. / src.at(1, 1);
        return;
    } else if ( nRows == 2 ) {
        det  = src.at(1, 1) * src.at(2, 2) - src.at(1, 2) * src.at(2, 1);
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
    } else    {
        // size >3 ... gaussian elimination - slow but safe
        //
        int i, j, k;
        double piv, linkomb;
        FloatMatrix tmp = src;
        this->zero();
        // initialize answer to be unity matrix;
        for ( i = 1; i <= nRows; i++ ) {
            this->at(i, i) = 1.0;
        }

        // lower triangle elimination by columns
        for ( i = 1; i < nRows; i++ ) {
            piv = tmp.at(i, i);
            if ( fabs(piv) < 1.e-20 ) {
                OOFEM_ERROR3("FloatMatrix::beInverseOf : cannot inverse a %d by %d matrix", nRows, nColumns);
            }

            for ( j = i + 1; j <= nRows; j++ ) {
                linkomb = tmp.at(j, i) / tmp.at(i, i);
                for ( k = i; k <= nRows; k++ ) {
                    tmp.at(j, k) -= tmp.at(i, k) * linkomb;
                }

                for ( k = 1; k <= nRows; k++ ) {
                    this->at(j, k) -= this->at(i, k) * linkomb;
                }
            }
        }

        // upper triangle elimination by columns
        for ( i = nRows; i > 1; i-- ) {
            piv = tmp.at(i, i);
            for ( j = i - 1; j > 0; j-- ) {
                linkomb = tmp.at(j, i) / tmp.at(i, i);
                for ( k = i; k > 0; k-- ) {
                    tmp.at(j, k) -= tmp.at(i, k) * linkomb;
                }

                for ( k = nRows; k > 0; k-- ) {
                    // tmp -> at(j,k)-= tmp  ->at(i,k)*linkomb;
                    this->at(j, k) -= this->at(i, k) * linkomb;
                }
            }
        }

        // diagonal scaling
        for ( i = 1; i <= nRows; i++ ) {
            for ( j = 1; j <= nRows; j++ ) {
                this->at(i, j) /= tmp.at(i, i);
            }
        }

        //delete tmp;
        return;
    }
}


void
FloatMatrix :: beSubMatrixOf(const FloatMatrix &src, const IntArray &indx)
/*
 * modifies receiver to be  (sub)matrix of the src.
 * (sub)matrix has size of max value in indx
 * and on its position (indx->at(i),indx->at(j)) are values from src at (i, j).
 * if indx->at(i) or indx->at(j) are <= 0 then at position i,j is zero.
 *
 * Warning:
 * This method should produce also bigger matrix than src
 * Works only for square matrices.
 */
{
    int size, n, i, j, ii, jj;
    // FloatMatrix *answer;

    if ( ( n = indx.giveSize() ) == 0 ) {
        this->resize(0, 0);
        return;
    }

#  ifdef DEBUG
    if ( !src.isSquare() ) {
        OOFEM_ERROR("FloatMatrix::beSubMatrixOf : cannot construct required submatrix");
    }

# endif

    if ( n != src.nRows ) {
        OOFEM_ERROR("FloatMatrix::beSubMatrixOf : giveSubMatrix size mismatch");
    }

    for ( size = 0, i = 1; i <= n; i++ ) {
        if ( indx.at(i) > size ) {
            size = indx.at(i);
        }
    }

    //answer = new FloatMatrix (size,size);
    this->resize(size, size);

    for ( i = 1; i <= n; i++ ) {
        for ( j = 1; j <= n; j++ ) {
            if ( ( ( ii = indx.at(i) ) != 0 ) && ( ( jj = indx.at(j) ) != 0 ) ) {
                this->at(ii, jj) = src.at(i, j);
            }
        }
    }

    return;
}



void
FloatMatrix :: beSubMatrixOfSizeOf(const FloatMatrix &src, const IntArray &indx, int size)
/*
 * modifies receiver to be  (sub)matrix of the src matrix.
 * (sub)matrix has size size
 * and on its position (indx->at(i),indx->at(j)) are values from src at (i, j).
 * if indx->at(i) or indx->at(j) are <= 0 then at position i,j is zero.
 *
 * Warning:
 * This method should produce also bigger matrix than src
 * Works only for square matrices.
 */
{
    int tsize, n, i, j, ii, jj;
    //FloatMatrix *answer;

    if ( ( n = indx.giveSize() ) == 0 ) {
        this->resize(0, 0);
        return;
    }

#  ifdef DEBUG
    if ( !src.isSquare() ) {
        OOFEM_ERROR("FloatMatrix::beSubMatrixOfSizeOf : cannot construct submatrix");
    }

# endif

    if ( n != src.nRows ) {
        OOFEM_ERROR("FloatMatrix::beSubMatrixOfSizeOf : giveSubMatrix size mismatch");
    }

    for ( tsize = 0, i = 1; i <= n; i++ ) {
        if ( indx.at(i) > tsize ) {
            tsize = indx.at(i);
        }
    }

    if ( tsize > size ) {
        OOFEM_ERROR("FloatMatrix::beSubMatrixOfSizeOf : index in mask exceed size");
    }

    //answer = new FloatMatrix (size,size);
    this->resize(size, size);

    for ( i = 1; i <= n; i++ ) {
        for ( j = 1; j <= n; j++ ) {
            if ( ( ( ii = indx.at(i) ) != 0 ) && ( ( jj = indx.at(j) ) != 0 ) ) {
                this->at(ii, jj) = src.at(i, j);
            }
        }
    }

    return;
}


void FloatMatrix :: plus(const FloatMatrix &aMatrix)
// Adds aMatrix to the receiver. If the receiver has a null size,
// adjusts its size to that of aMatrix. Returns the modified receiver.
{
    register int i;
    int n, m;
    double *P1, *P2;

    n = aMatrix.nRows;
    m = aMatrix.nColumns;
    if ( nRows * nColumns == 0 ) {
        this->operator=(aMatrix);
    } else {
#     ifdef DEBUG
        if ( n - nRows || m - nColumns ) {
            fprintf(stderr, "FloatMatrix::plus : dimensions mismatch : (r1,c1)+(r2,c2) : (%d,%d)+(%d,%d)", nRows, nColumns, n, m);
            exit(1);
        }

#     endif

        P1 = values;
        P2 = aMatrix.values;
        i  = n * m;
        while ( i-- ) {
            * P1++ += * P2++;
        }
    }

    return;
}

/*
 * void
 * FloatMatrix :: solveForRhs (const FloatArray& b, FloatArray& answer) const
 * // solves equation b = this * x
 * // returns x. this and b are kept untouched
 * //
 * // gaussian elimination - slow but safe
 * //
 * {
 * int i,j,k;
 * double piv,linkomb, help;
 *
 #  ifdef DEBUG
 * if (! this->isSquare()) {
 * printf ("error : cannot solve a %d by %d matrix ! \nfile %s, line %d\n",
 *    nRows,nColumns,__FILE__,__LINE__);
 * exit(0) ;}
 * if (nRows != b.giveSize()) {
 * printf ("error : dimension mismatch \nfile %s, line %d\n",__FILE__,__LINE__);
 * exit(0) ;}
 *
 #  endif
 *
 * FloatMatrix tmp = *this;
 * //FloatArray *answer = b -> GiveCopy();
 * answer = b ;
 * // initialize answer to be unity matrix;
 * // lower triangle elimination by columns
 * for (i=1; i< nRows; i++) {
 * piv = tmp.at(i,i);
 * if (fabs(piv) < 1.e-20) {
 * printf ("error : cannot solve  a %d by %d matrix ! \nfile %s, line %d\n",
 *     nRows,nColumns,__FILE__,__LINE__);
 * exit(0) ;
 * }
 * for (j=i+1; j<=nRows; j++) {
 * linkomb = tmp.at(j,i) / tmp.at(i,i);
 * for (k=i; k<= nRows; k++)
 *  tmp.at(j,k)-= tmp.at(i,k)*linkomb;
 * answer.at(j)-= answer.at(i)*linkomb;
 * }
 * }
 * // back substitution
 * for (i=nRows; i>=1; i--) {
 * help = 0.;
 * for (j=i+1; j<=nRows; j++)
 * help += tmp.at(i,j) * answer.at(j);
 * answer.at(i) = (answer.at(i) - help) / tmp.at(i,i);
 * }
 * // delete tmp;
 * return ;
 * }
 *
 * void
 * FloatMatrix :: solveForRhs (const FloatMatrix& b, FloatMatrix& answer) const
 * // solves equation b = this * x
 * // returns x. this and b are kept untouched
 * //
 * // gaussian elimination - slow but safe
 * //
 * {
 * int i,j,k,nPs;
 * double piv,linkomb, help;
 *
 #  ifdef DEBUG
 * if (! this->isSquare()) {
 * printf ("error : cannot solve a %d by %d matrix ! \nfile %s, line %d\n",
 *    nRows,nColumns,__FILE__,__LINE__);
 * exit(0) ;}
 * if (nRows != b.giveNumberOfRows()) {
 * printf ("error : dimension mismatch \nfile %s, line %d\n",__FILE__,__LINE__);
 * exit(0) ;}
 *
 #  endif
 *
 * nPs = b.giveNumberOfColumns();
 * FloatMatrix tmp = *this;
 * //FloatArray *answer = b -> GiveCopy();
 * answer = b ;
 * // initialize answer to be unity matrix;
 * // lower triangle elimination by columns
 * for (i=1; i< nRows; i++) {
 * piv = tmp.at(i,i);
 * if (fabs(piv) < 1.e-20) {
 * printf ("error : cannot solve  a %d by %d matrix ! \nfile %s, line %d\n",
 *     nRows,nColumns,__FILE__,__LINE__);
 * exit(0) ;
 * }
 * for (j=i+1; j<=nRows; j++) {
 * linkomb = tmp.at(j,i) / tmp.at(i,i);
 * for (k=i; k<= nRows; k++)
 *  tmp.at(j,k)-= tmp.at(i,k)*linkomb;
 * for (k=1; k<= nPs; k++)
 *  answer.at(j,k)-= answer.at(i,k)*linkomb;
 * }
 * }
 * // back substitution
 * for (i=nRows; i>=1; i--) {
 * for (k=1; k<= nPs; k++) {
 * help = 0.;
 * for (j=i+1; j<=nRows; j++)
 *  help += tmp.at(i,j) * answer.at(j,k);
 * answer.at(i,k) = (answer.at(i,k) - help) / tmp.at(i,i);
 * }
 * }
 * // delete tmp;
 * return ;
 * }
 */

void
FloatMatrix :: solveForRhs(const FloatArray &b, FloatArray &answer)
// solves equation b = this * x
// returns x. this and b are kept untouched
//
// gaussian elimination - slow but safe
// with row pivoting
//
{
    int i, j, k, pivRow;
    double piv, linkomb, help;

#  ifdef DEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR3("FloatMatrix::solveForRhs : cannot solve a %d by %d matrix", nRows, nColumns);
    }

    if ( nRows != b.giveSize() ) {
        OOFEM_ERROR("FloatMatrix::solveForRhs : dimension mismatch");
    }

#  endif

    //FloatArray *answer = b -> GiveCopy();
    answer = b;

    // initialize answer to be unity matrix;
    // lower triangle elimination by columns
    for ( i = 1; i < nRows; i++ ) {
        // find the suitable row and pivot
        piv = fabs( this->at(i, i) );
        pivRow = i;
        for ( j = i + 1; j <= nRows; j++ ) {
            if ( fabs( this->at(j, i) ) > piv ) {
                pivRow = j;
                piv = fabs( this->at(j, i) );
            }
        }

        if ( fabs(piv) < 1.e-20 ) {
            OOFEM_ERROR3("FloatMatrix::solveForRhs : cannot solve  a %d by %d matrix", nRows, nColumns);
        }

        // exchange rows
        if ( pivRow != i ) {
            for ( j = i; j <= nRows; j++ ) {
                help = this->at(i, j);
                this->at(i, j) = this->at(pivRow, j);
                this->at(pivRow, j) = help;
            }

            help = answer.at(i);
            answer.at(i) = answer.at(pivRow);
            answer.at(pivRow) = help;
        }

        for ( j = i + 1; j <= nRows; j++ ) {
            linkomb = this->at(j, i) / this->at(i, i);
            for ( k = i; k <= nRows; k++ ) {
                this->at(j, k) -= this->at(i, k) * linkomb;
            }

            answer.at(j) -= answer.at(i) * linkomb;
        }
    }

    // back substitution
    for ( i = nRows; i >= 1; i-- ) {
        help = 0.;
        for ( j = i + 1; j <= nRows; j++ ) {
            help += this->at(i, j) * answer.at(j);
        }

        answer.at(i) = ( answer.at(i) - help ) / this->at(i, i);
    }

    // delete tmp;
    return;
}

void
FloatMatrix :: solveForRhs(const FloatMatrix &b, FloatMatrix &answer)
// solves equation b = this * x
// returns x. this and b are kept untouched
//
// gaussian elimination - slow but safe
//
{
    int i, j, k, pivRow, nPs;
    double piv, linkomb, help;

#  ifdef DEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR3("FloatMatrix::solveForRhs : cannot solve a %d by %d matrix", nRows, nColumns);
    }

    if ( nRows != b.giveNumberOfRows() ) {
        OOFEM_ERROR("FloatMatrix::solveForRhs : dimension mismatch");
    }

#  endif

    nPs = b.giveNumberOfColumns();
    //FloatArray *answer = b -> GiveCopy();
    answer = b;
    // initialize answer to be unity matrix;
    // lower triangle elimination by columns
    for ( i = 1; i < nRows; i++ ) {
        // find the suitable row and pivot
        piv = fabs( this->at(i, i) );
        pivRow = i;
        for ( j = i + 1; j <= nRows; j++ ) {
            if ( fabs( this->at(j, i) ) > piv ) {
                pivRow = j;
                piv = fabs( this->at(j, i) );
            }
        }

        if ( fabs(piv) < 1.e-20 ) {
            OOFEM_ERROR3("FloatMatrix::solveForRhs : cannot solve  a %d by %d matrix", nRows, nColumns);
        }

        // exchange rows
        if ( pivRow != i ) {
            for ( j = i; j <= nRows; j++ ) {
                help = this->at(i, j);
                this->at(i, j) = this->at(pivRow, j);
                this->at(pivRow, j) = help;
            }

            for ( j = 1; j <= nPs; j++ ) {
                help = answer.at(i, j);
                answer.at(i, j) = answer.at(pivRow, j);
                answer.at(pivRow, j) = help;
            }
        }

        if ( fabs(piv) < 1.e-20 ) {
            OOFEM_ERROR("FloatMatrix::solveForRhs : cannot solve, zero pivot encountered");
        }

        for ( j = i + 1; j <= nRows; j++ ) {
            linkomb = this->at(j, i) / this->at(i, i);
            for ( k = i; k <= nRows; k++ ) {
                this->at(j, k) -= this->at(i, k) * linkomb;
            }

            for ( k = 1; k <= nPs; k++ ) {
                answer.at(j, k) -= answer.at(i, k) * linkomb;
            }
        }
    }

    // back substitution
    for ( i = nRows; i >= 1; i-- ) {
        for ( k = 1; k <= nPs; k++ ) {
            help = 0.;
            for ( j = i + 1; j <= nRows; j++ ) {
                help += this->at(i, j) * answer.at(j, k);
            }

            answer.at(i, k) = ( answer.at(i, k) - help ) / this->at(i, i);
        }
    }

    // delete tmp;
    return;
}



void
FloatMatrix :: initFromVector(const FloatArray &vector, int transposed)
//
// constructor : creates (vector->giveSize(),1) FloatMatrix
// if transpose = 1 creates (1,vector->giveSize()) FloatMatrix
//
{
    if ( transposed ) {
        resize( 1, vector.giveSize() );
    } else {
        resize(vector.giveSize(), 1); // row vector- default
    }

    if ( transposed ) {
        for ( int i = 1; i <= nColumns; i++ ) {
            this->at(1, i) = vector.at(i);
        }
    } else {
        for ( int i = 1; i <= nRows; i++ ) {
            this->at(i, 1) = vector.at(i);
        }
    }
}


//#########################################################################################//


FloatMatrix *FloatMatrix :: GiveCopy()
// Creates and returns a copy of the receiver.
{
    FloatMatrix *answer;
    double *P1, *P2;
    int i;

    answer = new FloatMatrix(nRows, nColumns);
    P1 = answer->values;
    P2 = values;
    i  = nRows * nColumns;
    while ( i-- ) {
        * P1++ = * P2++;
    }

    return answer;
}

void
FloatMatrix :: zero() const
{
    // zeroing the receiver - fast implementation
    double *P1;
    int i;

    P1 = this->values;
    i  = nRows * nColumns;
    while ( i-- ) {
        * P1++ = 0.;
    }
}

void FloatMatrix :: beUnitMatrix()
{
    int i;

    if ( !this->isSquare() ) {
        OOFEM_ERROR3("FloatMatrix::beUnitMatrix : cannot make unit matrix of %d by %d matrix", nRows, nColumns);
    }

    this->zero();
    for ( i = 1; i <= nRows; i++ ) {
        this->at(i, i) = 1.0;
    }
}

// this matrix is the product of the 6x6 deviatoric projection matrix ID
// and the inverse scaling matrix Pinv
void FloatMatrix :: bePinvID()
{
    this->resize(6,6);
    this->zero();
    values[0] = values[7] = values[14] = 2./3.;
    values[1] = values[2] = values[6] = values[8] = values[12] = values[13] = -1./3.;
    values[21] = values[28] = values[35] = 0.5;
}

void
FloatMatrix :: resize(int rows, int columns, int allocChunk)
//
// resizes receiver, all data will be lost
//
{
    if ( rows * columns > allocatedSize ) {
        // memory realocation necessary
        if ( values ) {
            freeDouble(values);
        }

        if ( allocChunk < 0 ) {
            allocChunk = 0;
        }

        allocatedSize = rows * columns + allocChunk; // REMEMBER NEW ALLOCATED SIZE
        values = allocDouble(allocatedSize);
    } else {
        // reuse previously allocated space
    }

    this->nRows = rows;
    this->nColumns = columns;
}

void
FloatMatrix :: resizeWithData(int rows, int columns)
//
// resizes receiver, all data kept
//
{
    FloatMatrix old(* this);

    if ( rows * columns > allocatedSize ) {
        // memory realocation necessary
        if ( values ) {
            freeDouble(values);
        }

        allocatedSize = rows * columns; // REMEMBER NEW ALLOCATED SIZE
        values = allocDouble(allocatedSize);
    } else {
        // reuse previously allocated space
    }

    this->nRows = rows;
    this->nColumns = columns;

    int ii, jj, i, j;

    ii = min( rows, old.giveNumberOfRows() );
    jj = min( columns, old.giveNumberOfColumns() );
    // copy old values if possible
    for ( i = 1; i <= ii; i++ ) {
        for ( j = 1; j <= jj; j++ ) {
            this->at(i, j) = old.at(i, j);
        }
    }
}


FloatMatrix *
FloatMatrix :: hardResize(int rows, int columns)
//
// resizes receiver, all data will be lost
//
{
    // memory realocation necessary
    if ( values ) {
        freeDouble(values);
    }

    allocatedSize = rows * columns; // REMEMBER NEW ALLOCATED SIZE
    values = allocDouble(allocatedSize);

    this->nRows = rows;
    this->nColumns = columns;

    return this;
}


double FloatMatrix :: giveDeterminant()
// Returns the determinant of the receiver.
{
#  ifdef DEBUG
    if ( !this->isSquare() ) {
        OOFEM_ERROR3("FloatMatrix::giveDeterminant : cannot compute determinant of a %d by %d matrix", nRows, nColumns);
    }

#  endif

    if ( nRows == 1 ) {
        return values [ 0 ];
    } else if ( nRows == 2 ) {
        return  ( values [ 0 ] * values [ 3 ] - values [ 1 ] * values [ 2 ] );
    } else if ( nRows == 3 ) {
        return ( values [ 0 ] * values [ 4 ] * values [ 8 ] + values [ 3 ] * values [ 7 ] * values [ 2 ] +
                 values [ 6 ] * values [ 1 ] * values [ 5 ] - values [ 6 ] * values [ 4 ] * values [ 2 ] -
                 values [ 7 ] * values [ 5 ] * values [ 0 ] - values [ 8 ] * values [ 3 ] * values [ 1 ] );
    } else {
        OOFEM_ERROR3("FloatMatrix::giveDeterminant : sorry, cannot inverse %d by %d matrices", nRows, nColumns);
    }

    return 0.;
}

/*
 * FloatArray*  FloatMatrix :: SolveForRhs (FloatArray *b)
 * // solves equation b = this * x
 * // returns x. this and b are kept untouched
 * //
 * // gaussian elimination - slow but safe
 * //
 * {
 * int i,j,k;
 * double piv,linkomb, help;
 *
 #  ifdef DEBUG
 * if (! this->isSquare()) {
 * printf ("error : cannot solve a %d by %d matrix ! \nfile %s, line %d\n",
 *    nRows,nColumns,__FILE__,__LINE__);
 * exit(0) ;}
 * if (nRows != b->giveSize()) {
 * printf ("error : dimension mismatch \nfile %s, line %d\n",__FILE__,__LINE__);
 * exit(0) ;}
 *
 #  endif
 *
 * FloatMatrix *tmp = this->GiveCopy();
 * FloatArray *answer = b -> GiveCopy();
 * // initialize answer to be unity matrix;
 * // lower triangle elimination by columns
 * for (i=1; i< nRows; i++) {
 * piv = tmp->at(i,i);
 * if (fabs(piv) < 1.e-20) {
 * printf ("error : cannot solve  a %d by %d matrix ! \nfile %s, line %d\n",
 *     nRows,nColumns,__FILE__,__LINE__);
 * exit(0) ;
 * }
 * for (j=i+1; j<=nRows; j++) {
 * linkomb = tmp->at(j,i) / tmp->at(i,i);
 * for (k=i; k<= nRows; k++)
 *  tmp -> at(j,k)-= tmp  ->at(i,k)*linkomb;
 * answer->at(j)-= answer->at(i)*linkomb;
 * }
 * }
 * // back substitution
 * for (i=nRows; i>=1; i--) {
 * help = 0.;
 * for (j=i+1; j<=nRows; j++)
 * help += tmp -> at(i,j) * answer->at(j);
 * answer->at(i) = (answer->at(i) - help) / tmp->at(i,i);
 * }
 * delete tmp;
 * return answer;
 * }
 */
/*
 * FloatMatrix*  FloatMatrix :: GiveInverse ()
 * // Returns a new matrix, the inverse of the receiver. (implemented only
 * // for 1x1 and 2x2 matrices)
 * {
 * FloatMatrix* answer ;
 * double       det ;
 * double*      p ;
 *
 #  ifdef DEBUG
 *    if (! this->isSquare()) {
 * printf ("error : cannot inverse a %d by %d matrix ! \nfile %s, line %d\n",
 *  nRows,nColumns,__FILE__,__LINE__);
 * exit(0) ;}
 #  endif
 *
 * answer = new FloatMatrix(nRows,nRows) ;
 * p      = answer->values ;
 *
 * if (nRows == 1) {
 *    p[0] = 1. / values[0] ;
 * return answer;
 * }
 * else if (nRows == 2) {
 *    det  = values[0]*values[3] - values[1]*values[2] ;
 *    p[0] =  values[3] / det ;
 *    p[1] = -values[1] / det ;
 *    p[2] = -values[2] / det ;
 *    p[3] =  values[0] / det ;
 * return answer;
 * }
 * else if (nRows == 3) {
 *    det = values[0]*values[4]*values[8]+values[3]*values[7]*values[2] +
 *   values[6]*values[1]*values[5]-values[6]*values[4]*values[2] -
 *   values[7]*values[5]*values[0]-values[8]*values[3]*values[1] ;
 *    p[0]= (values[4]*values[8]-values[7]*values[5])/det ;
 *    p[1]= (values[7]*values[2]-values[1]*values[8])/det ;
 *    p[2]= (values[1]*values[5]-values[4]*values[2])/det ;
 *    p[3]= (values[6]*values[5]-values[3]*values[8])/det ;
 *    p[4]= (values[0]*values[8]-values[6]*values[2])/det ;
 *    p[5]= (values[3]*values[2]-values[0]*values[5])/det ;
 *    p[6]= (values[3]*values[7]-values[6]*values[4])/det ;
 *    p[7]= (values[6]*values[1]-values[0]*values[7])/det ;
 *    p[8]= (values[0]*values[4]-values[3]*values[1])/det ;
 * return answer;
 *
 *  }
 * else {
 * // size >3 ... gaussian elimination - slow but safe
 * //
 * int i,j,k;
 * double piv,linkomb;
 * FloatMatrix *tmp = this->GiveCopy();
 * // initialize answer to be unity matrix;
 * for (i=1; i<= nRows; i++) answer->at(i,i) = 1.0;
 * // lower triangle elimination by columns
 * for (i=1; i< nRows; i++) {
 * piv = tmp->at(i,i);
 * if (fabs(piv) < 1.e-20) {
 *  printf ("error : cannot inverse a %d by %d matrix ! \nfile %s, line %d\n",
 *      nRows,nColumns,__FILE__,__LINE__);
 *  exit(0) ;
 * }
 * for (j=i+1; j<=nRows; j++) {
 *  linkomb = tmp->at(j,i) / tmp->at(i,i);
 *  for (k=i; k<= nRows; k++)
 *   tmp -> at(j,k)-= tmp  ->at(i,k)*linkomb;
 *  for (k=1; k<= nRows; k++)
 *   answer->at(j,k)-= answer->at(i,k)*linkomb;
 *
 * }
 * }
 * // upper triangle elimination by columns
 * for (i=nRows; i>1 ; i--) {
 * piv = this->at(i,i);
 * for (j=i-1; j > 0; j--) {
 *  linkomb = tmp->at(j,i) / tmp->at(i,i);
 *  for (k=nRows ; k > 0 ; k--) {
 *   // tmp -> at(j,k)-= tmp  ->at(i,k)*linkomb;
 *   answer->at(j,k)-= answer->at(i,k)*linkomb;
 *  }
 * }
 * }
 * // diagonal scaling
 * for (i=1; i<=nRows; i++)
 * for (j=1; j<=nRows; j++)
 *  answer->at(i,j)/= tmp->at(i,i);
 * delete tmp;
 * return answer ;
 * }
 * }
 */

FloatMatrix *FloatMatrix :: GiveTransposition()
// Returns a new matrix, the transposition of the receiver.
{
    int i, j;
    FloatMatrix *answer;

    answer = new FloatMatrix(nColumns, nRows);
    for ( i = 1; i <= nRows; i++ ) {
        for ( j = 1; j <= nColumns; j++ ) {
            answer->at(j, i) = this->at(i, j);
        }
    }

    return answer;
}


/*
 * FloatMatrix*  FloatMatrix :: Lumped ()
 * // Returns a new diagonal matrix, which is the lumped receiver : all
 * // coefficients on a column are concentrated on the diagonal.
 * // for creating of lumped mass matrix - good, but only for linear elements.
 * {
 * DiagonalMatrix *answer ;
 * double         sum,*p ;
 * int            j ;
 * register int   i ;
 *
 * answer = new DiagonalMatrix(nRows) ;
 * p      = values ;
 * for (j=1 ; j<=nColumns ; j++) {
 *    sum = 0. ;
 *    i   = nRows ;
 *    while (i--)
 * sum += *p++ ;
 *    answer -> at(j,j) = sum ;}
 *
 * return answer ;
 * }
 */
/*
 * void
 * FloatMatrix :: beLumpedOf (const FloatMatrix& src)
 * // Modifies a receiver to become  a new  matrix, which is the lumped src : all
 * // coefficients of src  on a column are concentrated on the diagonal.
 * // for creating of lumped mass matrix - good, but only for linear elements.
 *
 * // REMOVED DUE To UNREALISTIC IMPLEMENTATION
 * {
 * double         sum,*p ;
 * int            j ;
 * register int   i ;
 *
 * this->resize (src.nRows, src.nColumns) ;
 * p      = src.values ;
 * for (j=1 ; j<=nColumns ; j++) {
 *    sum = 0. ;
 *    i   = nRows ;
 *    while (i--)
 * sum += *p++ ;
 *    this -> at(j,j) = sum ;}
 *
 * return  ;
 * }
 */
/*
 * FloatMatrix*  FloatMatrix :: plus (FloatMatrix* aMatrix)
 * // Adds aMatrix to the receiver. If the receiver has a null size,
 * // adjusts its size to that of aMatrix. Returns the modified receiver.
 * {
 * register int i ;
 * int      n,m ;
 * double   *P1,*P2 ;
 *
 * if (aMatrix -> isDiagonal())
 *    return  this->plusDiagonalMatrix((DiagonalMatrix*)aMatrix) ;
 *
 * n = aMatrix -> nRows ;
 * m = aMatrix -> nColumns ;
 * if (nRows*nColumns == 0) {
 *    if (values)
 * freeDouble (values) ;
 *    nRows    = n ;
 *    nColumns = m ;
 *    i        = n * m ;
 *    values   = allocDouble(i) ;
 *    P1       = values ;
 *    P2       = aMatrix->values ;
 *    while (i--)
 *P1++ = *P2++ ;}
 * else {
 #     ifdef DEBUG
 * if (n-nRows || m-nColumns) {
 *    printf ("dimensions mismatch : r1,c1,r2,c2 : %d %d %d %d\nfile %s, line %d\n",
 *      nRows,n,nColumns,m,__FILE__,__LINE__) ;
 *    exit(0) ;}
 #     endif
 *
 *    P1 = values ;
 *    P2 = aMatrix->values ;
 *    i  = n * m ;
 *    while (i--)
 *P1++ += *P2++ ;}
 *
 * return this ;
 * }
 */
/*
 * FloatMatrix*  FloatMatrix :: plusDiagonalMatrix (DiagonalMatrix* aMatrix)
 * {
 * register int i ;
 * int      n ;
 *
 * n = aMatrix -> giveNumberOfRows() ;
 * if (nRows*nColumns == 0) {
 *    if (values)
 * freeDouble (values) ;
 *    nRows    = n ;
 *    nColumns = n ;
 *    values   = allocDouble(n*n) ;}
 *
 #  ifdef DEBUG
 *    if (n-nRows) {
 * printf ("dimensions mismatch in FloatMatrix+DiagonalMatrix\nfile %s,line %d\n",
 *    __FILE__,__LINE__) ;
 * exit(0) ;}
 #  endif
 *
 * for (i=1 ; i<=nRows ; i++)
 *    this->at(i,i) += aMatrix->at(i,i) ;
 *
 * return this ;
 * }
 */
/*
 * void  FloatMatrix :: plusProductSymmUpper (FloatMatrix* a, FloatMatrix* b, double dV)
 * // Adds to the receiver the product  a(transposed).b dV .
 * // If the receiver has a null size, it is expanded.
 * // This method assumes that both the receiver and the product above are
 * // symmetric matrices, and therefore computes only the upper half of the
 * // receiver ; the lower half is not modified. Other advantage : it does
 * // not compute the transposition of matrix a.
 * {
 * this->plusProductSymmUpper (*a, *b, dV);
 * }
 */

void FloatMatrix :: printYourself() const
// Prints the receiver on screen.
{
    int i, j;

    printf("FloatMatrix with dimensions : %d %d\n",
           nRows, nColumns);
    if ( nRows <= 50 && nColumns <= 50 ) {
        for ( i = 1; i <= nRows; ++i ) {
            for ( j = 1; j <= nColumns && j <= 50; ++j ) {
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
    int i, j;

    printf("[");
        for ( i = 1; i <= nRows; ++i ) {
            for ( j = 1; j <= nColumns; ++j ) {
                printf( "%20.15e", this->at(i, j) );
		if (j<nColumns)
		  printf(",");
		else
		  printf(";");
	    }
        }
	printf("];\n");
}

void
FloatMatrix :: rotatedWith(const FloatMatrix &r)
// Returns the receiver 'a' rotated according the change-of-base matrix r.
// The method performs the operation  a = r(transp) * a * r .
// Warning : this works only for square matrices (see copying to receiver)
// when the receiver and r has the same size!
{
    FloatMatrix rt, rta;

    rt.beTranspositionOf(r);         //  r(transp)
    rta.beProductOf(rt, * this);     //  r(transp) . a
    this->beProductOf(rta, r);       //  r(transp) . a . r

    return;
}


/*
 * FloatMatrix*  FloatMatrix :: RotatedWith (FloatMatrix* r)
 * // Returns the receiver 'a' rotated according the change-of-base matrix r.
 * // The method performs the operation  a = r(transp) * a * r .
 * // This method works for arbitrary matrices
 * {
 * // double        *p1,*p2 ;
 * FloatMatrix   *rt,*rta,*rtar ;
 *
 * rt   = r -> GiveTransposition() ;          //  r(transp)
 * rta  = rt -> Times(this) ;                 //  r(transp) . a
 * rtar = rta -> Times(r) ;                   //  r(transp) . a . r
 *
 * delete rt ;
 * delete rta ;
 * return rtar ;
 * }
 */



void FloatMatrix :: symmetrized()
// Initializes the lower half of the receiver to the upper half.
{
    int i, j;

#  ifdef DEBUG
    if ( nRows != nColumns ) {
        OOFEM_ERROR("FloatMatrix::symmetrized : cannot symmetrize a non-square matrix");
    }

#   endif

    for ( i = 2; i <= nRows; i++ ) {
        for ( j = 1; j < i; j++ ) {
            this->at(i, j) = this->at(j, i);
        }
    }

    return;
}


FloatMatrix *FloatMatrix :: times(double factor)
// Multiplies every coefficient of the receiver by factor. Answers the
// modified receiver.
{
    register int i;
    double *p;

    p = values;
    i = nRows * nColumns;
    while ( i-- ) {
        * p++ *= factor;
    }

    return this;
}


FloatArray *FloatMatrix :: Times(FloatArray *anArray)
// Returns the product of the receiver and anArray.
{
    int i, j;
    double sum;
    FloatArray *answer;

#  ifdef DEBUG
    if ( nColumns - anArray->giveSize() ) {
        OOFEM_ERROR4( "FloatMatrix::Times : in product A(%d,%d)*x(%d) : dimension mismatch", nRows, nColumns, anArray->giveSize() );
    }

#  endif

    answer = new FloatArray(nRows);
    for ( i = 1; i <= nRows; i++ ) {
        sum = 0.;
        for ( j = 1; j <= nColumns; j++ ) {
            sum += this->at(i, j) * anArray->at(j);
        }

        answer->at(i) = sum;
    }

    return answer;
}


FloatMatrix *FloatMatrix :: Times(FloatMatrix *aMatrix)
// Returns the product of the receiver and aMatrix. Easier to use than
// operator * .
{
    int i, j, k, p;
    double coeff;
    FloatMatrix *answer;

#  ifdef DEBUG
    if ( nColumns != aMatrix->nRows ) {
        fprintf(stderr, "FloatMatrix::Times : in product A(%d,%d)*B(%d,%d) : dimensions do not match",
                nRows, nColumns, aMatrix->nRows, aMatrix->nColumns);
        exit(0);
    }

#  endif

    p      = aMatrix->nColumns;
    answer = new FloatMatrix(nRows, p);
    for ( i = 1; i <= nRows; i++ ) {
        for ( j = 1; j <= p; j++ ) {
            coeff = 0.;
            for ( k = 1; k <= nColumns; k++ ) {
                coeff += this->at(i, k) * aMatrix->at(k, j);
            }

            answer->at(i, j) = coeff;
        }
    }

    return answer;
}

double
FloatMatrix :: computeFrobeniusNorm()
{
    int i, j;
    double answer = 0.0;
    for ( i = 1; i <= nRows; i++ ) {
        for ( j = 1; j <= nColumns; j++ ) {
            answer += this->at(i, j) * this->at(i, j);
        }
    }

    return sqrt(answer);
}



/*
 * FloatArray*  FloatMatrix :: timesTo (FloatArray* anArray, FloatArray* toArray)
 * // Stores the product of the receiver and anArray to toArray and returns
 * // pointer to toArray.
 * // toArray must exists and have proper dimensions
 * {
 * int         i,j ;
 * double      sum ;
 * // FloatArray* answer ;
 *
 * if (toArray==NULL) {
 * printf("TimesTo : toArray is NULL pointer\nfile %s,line %d\n",
 *    __FILE__, __LINE__);
 * exit(0);
 * }
 #  ifdef DEBUG
 *    if ((nColumns != anArray->giveSize()) || (nRows != toArray->giveSize())) {
 * printf ("error in product A*x : dimension mismatch \nfile %s, line %d\n",
 *    __FILE__,__LINE__) ;
 * exit(0) ;}
 #  endif
 *
 * for (i=1 ; i<=nRows ; i++) {
 *    sum = 0. ;
 *    for (j=1 ; j<=nColumns ; j++)
 * sum += this->at(i,j) * anArray->at(j) ;
 *    toArray->at(i) = sum ;}
 * return toArray ;
 * }
 */

/*
 * FloatMatrix*  FloatMatrix :: timesTo (FloatMatrix* aMatrix, FloatMatrix* toMatrix)
 * // stores the product of the receiver and aMatrix to toMatrix. Easier to use than
 * // operator * .
 * // toMatrix must have proper dimensions
 * {
 * int          i,j,k,p ;
 * double       coeff ;
 * // FloatMatrix* answer ;
 *
 #  ifdef DEBUG
 *    if (nColumns != aMatrix->nRows) {
 *  printf ("error in product A*B : dimensions do not match \nfile %s, line %d\n",
 *     __FILE__,__LINE__) ;
 *  exit(0) ;}
 * if ((toMatrix->nRows != this->nRows)||(toMatrix->nColumns!=aMatrix->nColumns)){
 *  printf ("error in product A*B : dimensions do not match \nfile %s, line %d\n",
 *     __FILE__,__LINE__) ;
 * exit(0) ;}
 #  endif
 *
 * p      = aMatrix -> nColumns ;
 * for (i=1 ; i<=nRows ; i++)
 *    for (j=1 ; j<=p ; j++) {
 *       coeff = 0. ;
 * for (k=1 ; k<=nColumns ; k++)
 *   coeff += this->at(i,k) * aMatrix->at(k,j) ;
 * toMatrix->at(i,j) = coeff ;}
 * return toMatrix ;
 * }
 */


contextIOResultType FloatMatrix :: storeYourself(DataStream *stream, ContextMode mode)
// writes receiver's binary image into stream
// use id to distinguish some instances
// return value >0 succes
//              =0 file i/o error
{
    int type_id = FloatMatrixClass;
    int size = nRows * nColumns;
    // write class header
    if ( !stream->write(& type_id, 1) ) {
        return ( CIO_IOERR );
    }

    // write size
    if ( !stream->write(& nRows, 1) ) {
        return ( CIO_IOERR );
    }

    if ( !stream->write(& nColumns, 1) ) {
        return ( CIO_IOERR );
    }

    // write raw data
    if ( !stream->write(values, size) ) {
        return ( CIO_IOERR );
    }

    // return result back
    return CIO_OK;
}

contextIOResultType FloatMatrix :: restoreYourself(DataStream *stream, ContextMode mode)
// reads receiver from stream
// warning - overwrites existing data!
// returns 0 if file i/o error
//        -1 if id od class id is not correct
{
    int class_id;
    // read class header
    if ( !stream->read(& class_id, 1) ) {
        return ( CIO_IOERR );
    }

    if ( class_id != FloatMatrixClass ) {
        return ( CIO_BADVERSION );
    }

    // read size
    if ( !stream->read(& nRows, 1) ) {
        return ( CIO_IOERR );
    }

    if ( !stream->read(& nColumns, 1) ) {
        return ( CIO_IOERR );
    }

    if ( values != NULL ) {
        freeDouble(values);
    }

    if ( nRows * nColumns ) {
        values = allocDouble(nRows * nColumns);
        allocatedSize = nRows * nColumns;
    } else {
        values = NULL;
        allocatedSize = 0;
    }

    // write raw data
    if ( !stream->read(values, nRows * nColumns) ) {
        return ( CIO_IOERR );
    }

    // return result back
    return CIO_OK;
}

/*
 * void
 * FloatMatrix :: assemble (FloatMatrix* src, IntArray* loc)
 * {
 * int i,j,ii,jj,size;
 *
 * if ((size=src->giveNumberOfRows()) != loc->giveSize()) {
 * printf ("dimensions of 'src' and 'loc' mismatch \n") ;
 * exit(0) ; }
 *
 * if (!src->isSquare()) {
 * printf ("'src' is not sqaure matrix\n") ;
 * exit(0) ; }
 *
 * for (i=1; i<= size; i++) {
 * if (ii = loc->at(i)) {
 * for (j=1; j<=size; j++) {
 *  if (jj = loc->at(j)) this->at(ii,jj)+= src->at(i,j);
 * }
 * }
 * }
 * }
 */

/*
 * FloatMatrix*
 * FloatMatrix :: GiveSubMatrix ( int topRow, int bottomRow , int topCol, int bottomCol)
 * // returns newly allocated submatrix of the receiver
 * // size of returned submatrix is determined from
 * // input parametrs
 * {
 #ifdef DEBUG
 * // check input params
 * if ((topRow < 1) || (bottomRow < 1) || (topCol < 1) || (bottomCol < 1)) {
 *  printf ("GiveSubMatrix : subindexes size mismatch\nfile %s, line %d\n",__FILE__,__LINE__);
 * exit(1);
 * }
 *
 * if ((nRows < bottomRow) || (nColumns < bottomCol) || ((bottomRow-topRow) > nRows) ||
 *    ((bottomCol-topCol) > nColumns)) {
 *  printf ("GiveSubMatrix : subindexes size mismatch\nfile %s, line %d\n",__FILE__,__LINE__);
 * exit(1);
 * }
 #endif
 *
 *
 * int i,j,topRm1,topCm1;
 * topRm1=topRow-1;
 * topCm1=topCol-1;
 *
 * // allocate return value
 * FloatMatrix * subM = new FloatMatrix (bottomRow-topRm1,bottomCol-topCm1);
 * for (i = topRow; i<= bottomRow; i++)
 *  for (j = topCol; j<= bottomCol; j++)
 *    subM->at(i-topRm1,j-topCm1) = this->at(i,j);
 *
 * return subM;
 *
 * }
 */
/*
 * FloatMatrix*
 * FloatMatrix :: GiveSubMatrix (IntArray *indx)
 *
 * // returns newly allocated (sub)matrix of the receiver.
 * // (sub)matrix has size of max value in indx
 * // and on its position (indx->at(i),indx->at(j)) are values from receiver at (i, j).
 * // if indx->at(i) or indx->at(j) are <= 0 then at position i,j is zero.
 * //
 * // Warning:
 * //   This method should produce also bigger matrix than receiver
 * //   Works only for square matrices.
 *
 * {
 * int size,n,i,j,ii,jj;
 * FloatMatrix *answer;
 *
 * if (indx == NULL) return NULL;
 * if ((n = indx->giveSize()) == 0) return NULL;
 #  ifdef DEBUG
 * if (! this->isSquare()) {
 * printf ("error : cannot return submatrix constructed from IntArray from %d by %d matrix ! \nfile %s, line %d\n",
 *    nRows,nColumns, __FILE__,__LINE__);
 * exit(0);
 * }
 # endif
 *
 * if (n!=nRows) {
 * printf("error: GiveSubMatrix size mismatch\n,file %s, line %d\n",__FILE__,__LINE__);
 * exit (0);
 * }
 *
 * for (size = 0, i=1; i<=n;i++) if (indx->at(i) > size) size = indx->at(i);
 * answer = new FloatMatrix (size,size);
 *
 * for (i=1; i<= n; i++)
 * for (j=1; j<= n; j++) {
 * if (((ii = indx->at(i))!=0) && ((jj = indx->at(j))!=0)) {
 *  answer->at(ii,jj) = this -> at(i,j);
 * }
 * }
 * return answer;
 * }
 *
 *
 *
 * FloatMatrix*
 * FloatMatrix :: GiveSubMatrixOfSize (IntArray *indx, int size)
 *
 * // returns newly allocated (sub)matrix of the receiver.
 * // (sub)matrix has size size
 * // and on its position (indx->at(i),indx->at(j)) are values from receiver at (i, j).
 * // if indx->at(i) or indx->at(j) are <= 0 then at position i,j is zero.
 * //
 * // Warning:
 * //   This method should produce also bigger matrix than receiver
 * //   Works only for square matrices.
 *
 * {
 * int tsize,n,i,j,ii,jj;
 * FloatMatrix *answer;
 *
 * if (indx == NULL) return NULL;
 * if ((n = indx->giveSize()) == 0) return NULL;
 #  ifdef DEBUG
 * if (! this->isSquare()) {
 * printf ("error : cannot return submatrix constructed from IntArray from %d by %d matrix ! \nfile %s, line %d\n",
 *    nRows,nColumns,__FILE__,__LINE__);
 * exit(0);
 * }
 # endif
 *
 * if (n!=nRows) {
 * printf("error: GiveSubMatrix size mismatch\nfile %s, line %d\n",__FILE__,__LINE__);
 * exit (0);
 * }
 *
 * for (tsize = 0, i=1; i<=n;i++) if (indx->at(i) > tsize) tsize = indx->at(i);
 * if (tsize > size) {
 * printf("error: GiveSubMatrixOfSize index in mask exceed size\nfile %s, line %d\n",
 *    __FILE__,__LINE__);
 * exit (0);
 * }
 *
 * answer = new FloatMatrix (size,size);
 *
 * for (i=1; i<= n; i++)
 * for (j=1; j<= n; j++) {
 * if (((ii = indx->at(i))!=0) && ((jj = indx->at(j))!=0)) {
 *  answer->at(ii,jj) = this -> at(i,j);
 * }
 * }
 * return answer;
 * }
 */


int
FloatMatrix :: jaco_(FloatArray &eval, FloatMatrix &v, int nf)
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
        OOFEM_ERROR("FloatMatrix::jaco_: Not square matrix");
    }

    // check for symetry
    for ( i = 1; i <= neq; i++ ) {
        for ( j = i + 1; j <= neq; j++ ) {
            if ( this->at(i, j) != this->at(j, i) ) {
                OOFEM_ERROR("FloatMatrix::jaco_: Not Symmetric matrix");
            }
        }
    }

#endif

    eval.resize(neq);
    v.resize(neq, neq);

    for ( i = 1; i <= neq; i++ ) {
        eval.at(i) = this->at(i, i);
    }

    tol = __OOFEM_POW(c_b2, nf);
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
                        this->at(k, i) = co * tt + si * this->at(k, j);
                        this->at(k, j) = -si * tt + co * this->at(k, j);
                        tt = v.at(k, i);
                        v.at(k, i) = co * tt + si *v.at(k, j);
                        v.at(k, j) = -si * tt + co *v.at(k, j);
                    }

                    // diagonal term (i,i)
                    tt = eval.at(i);
                    eval.at(i) = co * tt + si * this->at(i, j);
                    aij = -si * tt + co * this->at(i, j);
                    tt = v.at(i, i);
                    v.at(i, i) = co * tt + si *v.at(i, j);
                    v.at(i, j) = -si * tt + co *v.at(i, j);

                    for ( k = i + 1; k < j; ++k ) {
                        tt = this->at(i, k);
                        this->at(i, k) = co * tt + si * this->at(k, j);
                        this->at(k, j) = -si * tt + co * this->at(k, j);
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
                        this->at(i, k) = co * tt + si * this->at(j, k);
                        this->at(j, k) = -si * tt + co * this->at(j, k);
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
            OOFEM_ERROR("FloatMatrix::jaco_: too many iterations");
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



#ifdef __PARALLEL_MODE
int
FloatMatrix :: packToCommBuffer(CommunicationBuffer &buff) const
{
    int result = 1;
    // pack size
    result &= buff.packInt(nRows);
    result &= buff.packInt(nColumns);
    // pack data
    result &= buff.packArray(this->values, nRows * nColumns);

    return result;
}

int
FloatMatrix :: unpackFromCommBuffer(CommunicationBuffer &buff)
{
    int newNRows, newNColumns, result = 1;
    // unpack size
    result &= buff.unpackInt(newNRows);
    result &= buff.unpackInt(newNColumns);
    // resize yourself
    this->resize(newNRows, newNColumns);
    result &= buff.unpackArray(this->values, newNRows * newNColumns);

    return result;
}

int
FloatMatrix :: givePackSize(CommunicationBuffer &buff)
{
    return buff.givePackSize(MPI_INT, 1) + buff.givePackSize(MPI_INT, 1) +
           buff.givePackSize(MPI_DOUBLE, nRows * nColumns);
}


#endif
} // end namespace oofem
