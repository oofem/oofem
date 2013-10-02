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

#include "DenseMatrixArithmeticsNN.h"

DSS_NAMESPASE_BEGIN

long DenseMatrixArithmetics :: zero_pivots = 0;

DenseMatrixArithmetics *DenseMatrixArithmetics :: NewArithmetics(long block_size)
{
    if ( block_size == 1 ) {
        return new DenseMatrixArithmetics1x1(1);
    }

    if ( block_size == 2 ) {
        return new DenseMatrixArithmetics2x2(2);
    }

    if ( block_size == 3 ) {
        return new DenseMatrixArithmetics3x3(3);
    }

    if ( block_size == 4 ) {
        return new DenseMatrixArithmetics4x4(4);
    }

    if ( block_size == 5 ) {
        return new DenseMatrixArithmetics5x5(5);
    }

    if ( block_size == 6 ) {
        return new DenseMatrixArithmetics6x6(6);
    } else  {
        return new DenseMatrixArithmetics(block_size);
    }
}

/////////////////////////////////
// DenseMatrixArithmetics
DenseMatrixArithmetics :: DenseMatrixArithmetics(long bn)
{
    this->bn = bn;
    bnbn = bn * bn;
    bn1 = bn + 1;
    p = new double [ bn ];
    DenseMatrixArithmetics :: zero_pivots = 0;
    prefered_decomposition = eLDL_decomposition;
    eMT = & MT;
}

DenseMatrixArithmetics :: ~DenseMatrixArithmetics()
{
    if ( p ) {
        delete [] p;
    }

    p = NULL;
}

// b += B*x
void DenseMatrixArithmetics :: AddMultBlockByVectorSym(double *B, double *x, double *b, long bi_bn, long bj_bn)
{
    for ( long i = 0; i < bn; i++ ) {
        for ( long j = 0; j < bn; j++ ) {
            double DataIJ = B [ i + bn * j ];
            b [ bi_bn + i ] += DataIJ * x [ bj_bn + j ];
            b [ bj_bn + j ] += DataIJ * x [ bi_bn + i ];
        }
    }
}

// b -= A*x
void DenseMatrixArithmetics :: SubMultBlockByVector(double *B, double *x, double *b)
{
    for ( long i = 0; i < bn; i++ ) {
        for ( long j = 0; j < bn; j++ ) {
            b [ i ] -= B [ i + bn * j ] * x [ j ];
        }
    }
}

// b -= A^T * x
void DenseMatrixArithmetics :: SubMultTBlockByVector(double *B, double *x, double *b)
{
    for ( long i = 0; i < bn; i++ ) {
        for ( long j = 0; j < bn; j++ ) {
            b [ i ] -= B [ j + bn * i ] * x [ j ];
        }
    }
}

// b += B * x
void DenseMatrixArithmetics :: MultDiagonalBlockByVector(double *B, double *x, double *b)
{
    for ( long i = 0; i < bn; i++ ) {
        b [ i ] += B [ i + bn * i ] * x [ i ];
        for ( long j = i + 1; j < bn; j++ ) {
            double DataIJ = B [ i + bn * j ];
            b [ i ] += DataIJ * x [ j ];
            b [ j ] += DataIJ * x [ i ];
        }
    }
}

// C -= A^T * B
// This is the most SPARSEDIRECT SOLVER innerloop operation
void DenseMatrixArithmetics :: SubATBproduct(double *pC, double *pA, double *pB)
{
    long i, j, k;
    for ( j = 0; j < bn; j++, pA -= bnbn, pB += bn ) {
        for ( i = 0; i < bn; i++, pC++, pB -= bn ) {
            tmp = 0;
            for ( k = 0; k < bn; k++, pA++, pB++ ) {
                tmp += * pA * ( * pB );
            }

            * pC -= tmp;
        }
    }
}

/// <summary>
/// Solves LL' decomposition of symmetric pC column based matrix
/// The result is stored in lower half of the matrix pC
/// </summary>
/// <param name="pC"> matrix to be factorized</param>
/// <param name="p"> here we store the diagonal elements of L</param>
void DenseMatrixArithmetics :: Cholesky_Decomposition(double *pC, double * &p)
{
    double *pi, *pj, *a = pC;
    long i, j, k, n = bn;
    double sum;
    if ( p == NULL ) {
        p = new double [ bn ];
    }

    for ( i = 0; i < n; i++ ) {
        for ( j = i; j < n; j++ ) {
            //for (sum=a[i+n*j],k=i-1;k>=0;k--)
            //sum -= a[i+n*k]*a[j+n*k];
            for ( k = i - 1, pi = a + i + k * n, pj = a + j + k * n, sum = a [ i + n * j ]; k >= 0; k--, pi -= n, pj -= n ) {
                sum -= * pi * * pj;              //a[i+n*k]*a[j+n*k];
            }

            if ( i == j ) {
                if ( sum > 0.0 ) {
                    p [ i ] = sqrt(sum);
                } else {
                    p [ i ] = 1.0;
                    DenseMatrixArithmetics :: zero_pivots++;
                }
            } else   {
                a [ j + n * i ] = sum / p [ i ];
            }
        }
    }
}

/// <summary> Solve LL^T x = b </summary>
/// <param name="pC">factorized L matrix (lower triangle)</param>
/// <param name="p">diagonal elements of L</param>
/// <param name="b">right side</param>
/// <param name="x">unknowns</param>
void DenseMatrixArithmetics :: Cholesky_Solve(double *pC, double *p, double *b, double *x)
{
    double *a = pC;
    long i, k, n = bn;
    double sum;

    for ( i = 0; i < n; i++ ) { //Solve L y = b, storing y in x.
        for ( sum = b [ i ], k = i - 1; k >= 0; k-- ) {
            sum -= a [ i + n * k ] * x [ k ];
        }

        x [ i ] = sum / p [ i ];
    }

    for ( i = n - 1; i >= 0; i-- ) { //Solve LT  x = y.
        for ( sum = x [ i ], k = i + 1; k < n; k++ ) {
            sum -= a [ k + n * i ] * x [ k ];
        }

        x [ i ] = sum / p [ i ];
    }
}

void DenseMatrixArithmetics :: Cholesky_Linv(double *pC, double *p)
{
    double *a = pC;
    long i, j, k, n = bn;
    double sum;
    for ( i = 0; i < n; i++ ) {
        a [ i + n * i ] = 1.0 / p [ i ];
        for ( j = i + 1; j < n; j++ ) {
            sum = 0.0;
            for ( k = i; k < j; k++ ) {
                sum -= a [ j + n * k ] * a [ k + n * i ];
            }

            a [ j + n * i ] = sum / p [ j ];
        }
    }
}

void DenseMatrixArithmetics :: ComputeInversionByCholesky(double *C, double *Inv)
{
    Cholesky_Decomposition(C, this->p);
    Cholesky_Linv(C, this->p);

    double *a = C;
    long i, j, k, n = bn;
    double sum;

    for ( i = 0; i < n; i++ ) {
        for ( j = i; j < n; j++ ) {
            sum = 0;
            for ( k = j; k < n; k++ ) {
                sum += a [ k + n * i ] * a [ k + n * j ];
            }

            Inv [ j + n * i ] = Inv [ i + n * j ] = sum;
        }
    }
}


/// <summary>
/// Solves LL' decomposition of symmetric pC column based matrix
/// The result is stored in lower half of the matrix pC
/// </summary>
/// <param name="pC"> matrix to be factorized</param>
void DenseMatrixArithmetics :: LL_Decomposition(double *pC)
{
    double *pi, *pj, *a = pC;
    long i, j, k, n = bn;
    double sum;

    for ( i = 0; i < n; i++ ) {
        for ( j = i; j < n; j++ ) {
            //for (sum=a[i+n*j],k=i-1;k>=0;k--)
            // sum -= a[i+n*k]*a[j+n*k];
            for ( k = i - 1, pi = a + i + k * n, pj = a + j + k * n, sum = a [ i + n * j ]; k >= 0; k--, pi -= n, pj -= n ) {
                sum -= * pi * * pj;              //a[i+n*k]*a[j+n*k];
            }

            if ( i == j ) {
                if ( sum > 0.0 ) {
                    a [ bn * i ] = sqrt(sum);
                } else {
                    this->MT.Write("Matrix is not positive definite.");
                    a [ bn * i ] = 1.0;
                    DenseMatrixArithmetics :: zero_pivots++;
                }
            } else   {
                a [ j + n * i ] = sum / a [ bn * i ];
            }
        }
    }
}

/// <summary> Solve LL^T x = b </summary>
/// <param name="pC">factorized L matrix (lower triangle)</param>
/// <param name="b">right side</param>
/// <param name="x">unknowns</param>
void DenseMatrixArithmetics :: LL_Solve(double *pC, double *b, double *x)
{
    double *a = pC;
    long i, k, n = bn;
    double sum;

    for ( i = 0; i < n; i++ ) { //Solve L y = b, storing y in x.
        for ( sum = b [ i ], k = i - 1; k >= 0; k-- ) {
            sum -= a [ i + n * k ] * x [ k ];
        }

        x [ i ] = sum / a [ bn * i ];
    }

    for ( i = n - 1; i >= 0; i-- ) { //Solve LT  x = y.
        for ( sum = x [ i ], k = i + 1; k < n; k++ ) {
            sum -= a [ k + n * i ] * x [ k ];
        }

        x [ i ] = sum / a [ bn * i ];
    }
}

/// <summary> Solve LL^T x = b </summary>
/// <param name="pC">factorized L matrix (lower triangle)</param>
/// <param name="b">right side</param>
/// <param name="x">unknowns</param>
void DenseMatrixArithmetics :: L_BlockSolve(double *C, double *B)
{
    double *a = C;
    int i, j, n = bn, k;
    double sum;
    double *b = B;

    for ( j = 0; j < n; j++ ) {
        for ( i = 0; i < n; i++ ) { //Solve L y = b, storing y in x.
            for ( sum = b [ i ], k = i - 1; k >= 0; k-- ) {
                sum -= a [ i + n * k ] * b [ k ];
            }

            b [ i ] = sum / a [ bn * i ];
        }

        b += n;
    }
}

void DenseMatrixArithmetics :: SubstSolveL(double *a, double *x)
{
    long i, k, n = bn;
    double sum;

    for ( i = 0; i < n; i++ ) { //Solve L y = b, storing y in x.
        for ( sum = x [ i ], k = i - 1; k >= 0; k-- ) {
            sum -= a [ i + n * k ] * x [ k ];
        }

        x [ i ] = sum / a [ bn * i ];
    }
}


void DenseMatrixArithmetics :: SubstSolveLT(double *a, double *x)
{
    long i, k, n = bn;
    double sum;

    for ( i = n - 1; i >= 0; i-- ) { //Solve LT  x = y.
        for ( sum = x [ i ], k = i + 1; k < n; k++ ) {
            sum -= a [ k + n * i ] * x [ k ];
        }

        x [ i ] = sum / a [ bn * i ];
    }
}


void DenseMatrixArithmetics :: LU_Decomposition(double *C)
{
    double *a = C;
    //double[,] a = C.data;

    long n = bn;
    double TINY = 1e-10;
    long i, j, k;
    double dum, sum;

    for ( j = 0; j < n; j++ ) { //This is the loop over columns of Crout's method.
        for ( i = 0; i < j; i++ ) { //This is equation (2.3.12) except for i = j.
            sum = a [ i + n * j ];
            for ( k = 0; k < i; k++ ) {
                sum -= a [ i + n * k ] * a [ k + n * j ];
            }

            a [ i + n * j ] = sum;
        }

        for ( i = j; i < n; i++ ) { //This is i = j of equation (2.3.12) and i = j+1 : ::N of equation (2.3.13).
            sum = a [ i + n * j ];
            for ( k = 0; k < j; k++ ) {
                sum -= a [ i + n * k ] * a [ k + n * j ];
            }

            a [ i + n * j ] = sum;
        }

        if ( a [ j + n * j ] == 0.0 ) {
            a [ j + n * j ] = TINY;
            DenseMatrixArithmetics :: zero_pivots++;
        } else   {
            a [ j + n * j ] = 1.0 / a [ j + n * j ]; //invert block
        }

        if ( j != n ) {
            dum = a [ j + n * j ];
            for ( i = j + 1; i < n; i++ ) {
                a [ i + n * j ] *= dum;
            }
        }
    }
}

void DenseMatrixArithmetics :: LU_Solve(double *C, double *b)
{
    double *a = C;

    long i, ii = -1, j, n = bn;
    double sum;
    for ( i = 0; i < n; i++ ) { //When ii is set to a positive value, it will become the
        //ip=indx[i];
        sum = b [ i ];
        if ( ii != -1 ) {
            for ( j = ii; j <= i - 1; j++ ) {
                sum -= a [ i + n * j ] * b [ j ];
            }
        } else if ( sum != 0.0 ) {
            ii = i;
        }

        //A nonzero element was encountered, so from now on we will have to do the sums in the loop above.
        b [ i ] = sum;
    }

    for ( i = n - 1; i >= 0; i-- ) { //Now we do the backsubstitution, equation (2.3.7).
        sum = b [ i ];
        for ( j = i + 1; j < n; j++ ) {
            sum -= a [ i + n * j ] * b [ j ];
        }

        b [ i ] = sum * a [ i + n * i ]; //Store a component of the solution vector X.
    }     //All done!

}

// Solves system X*LU=B
// B is stored by rows
void DenseMatrixArithmetics :: ULT_BlockSolve(double *C, double *B)
{
    //double[,] a = C.data;
    double *a = C;
    int i, j, n = bn, k;
    double sum;

    for ( k = 0; k < n; k++ ) {
        for ( j = 0; j < n; j++ ) {
            sum = B [ k * n + j ];
            for ( i = 0; i < j; i++ ) {
                sum -= a [ n * j + i ] * B [ k * n + i ];
            }

            B [ k * n + j ] = sum * a [ j + n * j ];
        }

        for ( j = n - 2; j >= 0; j-- ) {
            sum = 0;
            for ( i = j + 1; i < n; i++ ) {
                sum += a [ n * j + i ] * B [ k * n + i ];
            }

            B [ k * n + j ] -= sum;
        }         //All done!

    }
}


void DenseMatrixArithmetics :: LDL_Decomposition(double *A)
{
    long i, j, k;
    double *Ai, *Aj, *Aii, *Ajj = A;
    double *A0j = A;
    long n = bn;

    for ( j = 0; j < n; j++, Ajj += bn1, A0j += n ) {
        // eliminate column
        for ( i = 1; i < j; i++ ) {
            double sum = 0;
            for ( k = 0, Ai = A + n * i, Aj = A0j; k < i; k++ ) {
                sum += * ( Ai++ ) * ( * ( Aj++ ) );
            }

            * Aj -= sum;
        }

        // divide by pivot
        for ( Aj = A0j, Aii = A; Aii < Ajj; Aj++, Aii += bn1 ) {
            double A_ij = * Aj;
            * Aj *= * Aii;
            * Ajj -= ( * Aj ) * A_ij;
        }

        if ( fabs(* Ajj) <= eMT->min_pivot ) {
            if ( eMT->stabil_pivot == 0.0 ) {
                eMT->act_row = eMT->act_block + j;
                if ( !eMT->CallUnstableDialog() ) {
                    return;
                }
            }

            DenseMatrixArithmetics :: zero_pivots++;
            * Ajj = eMT->stabil_pivot;
        }

        * Ajj = 1.0 / ( * Ajj );
    }
}

void DenseMatrixArithmetics :: LDL_Solve(double *A, double *b, const long startIndex)
{
    long i, j, n = bn;
    double tmp;

    for ( j = startIndex; j < n - 1; j++ ) {
        tmp = b [ j ];
        for ( i = j + 1; i < n; i++ ) {
            b [ i ] -= A [ j + i * n ] * tmp;   //L[i,j] is stored in A[j,i]
        }
    }

    for ( j = startIndex; j < n; j++ ) {
        b [ j ] *= A [ j + j * n ];
    }

    for ( j = n - 1; j > 0; j-- ) {
        tmp = b [ j ];
        for ( i = 0; i < j; i++ ) {
            b [ i ] -= A [ i + j * n ] * tmp;
        }
    }
}

// A - decomposed matrix
// B - block of right sides
void DenseMatrixArithmetics :: SubstSolveBlock(double *A, double *B)
{
    double *b = B;
    for ( int j = 0; j < bn; j++, b += bn ) {
        SubstSolve(A, b);
    }
}

void DenseMatrixArithmetics :: FactorizeBlock(double *A)
{
    switch ( prefered_decomposition ) {
    case eLDL_decomposition:        LDL_Decomposition(A);
        break;
    case eLL_decomposition:         LL_Decomposition(A);
        break;
    }
}

// Solves A b = b
void DenseMatrixArithmetics :: SubstSolve(double *A, double *b)
{
    switch ( prefered_decomposition ) {
    case eLDL_decomposition:        LDL_Solve(A, b);
        break;
    case eLL_decomposition:         LL_Solve(A, b, b);
        break;
    }
}

void DenseMatrixArithmetics :: ComputeInversionByLDL(double *C, double *Inv)
{
    long j;
    LDL_Decomposition(C);    //Decompose the matrix just once.
    memset( Inv, 0, bn * bn * sizeof( double ) );
    for ( j = 0; j < bn; j++, Inv += bn ) { //Find inverse by columns.
        Inv [ j ] = 1.0;
        LDL_Solve(C, Inv, 0);
    }
}

void DenseMatrixArithmetics :: GetInversion(double *C, double *blockA)
{
    ComputeInversionByLDL(C, blockA);
}

bool DenseMatrixArithmetics :: GaussElimination(double *C, double *blockA)
{
    * blockA = 1 / ( * C );
    return true;
    // C;blockA;
}

/////////////////////////////////
// DenseMatrixArithmetics1x1
void DenseMatrixArithmetics1x1 :: SubATBproduct(double *pC, double *pA, double *pB)
{
    * pC -= ( * pA ) * ( * pB );
}

void DenseMatrixArithmetics1x1 :: GetInversion(double *pC, double *pA)
{
    * pA = 1.0 / ( * pC );
}

void DenseMatrixArithmetics1x1 :: FactorizeBlock(double *A)
{
    if ( * A == 0.0 ) {
        DenseMatrixArithmetics :: zero_pivots++;
        * A = 1.0;
    } else {
        * A = 1.0 / ( * A );
    }
}
void DenseMatrixArithmetics1x1 :: SubstSolveBlock(double *A, double *B)
{
    * B *= * A;
}
void DenseMatrixArithmetics1x1 :: SubstSolve(double *A, double *b)
{
    * b *= * A;
}

/////////////////////////////////
// DenseMatrixArithmetics2x2
void DenseMatrixArithmetics2x2 :: SubATBproduct(double *pC, double *pa, double *pb)
{
    pC [ 0 ] -= pa [ 0 ] * pb [ 0 ] + pa [ 1 ] * pb [ 1 ];
    pC [ 1 ] -= pa [ 2 ] * pb [ 0 ] + pa [ 3 ] * pb [ 1 ];
    pC [ 2 ] -= pa [ 0 ] * pb [ 2 ] + pa [ 1 ] * pb [ 3 ];
    pC [ 3 ] -= pa [ 2 ] * pb [ 2 ] + pa [ 3 ] * pb [ 3 ];
}

void DenseMatrixArithmetics2x2 :: GetInversion(double *C, double *pA)
{
    double D = C [ 0 ] * C [ 3 ] - C [ 1 ] * C [ 1 ];
    pA [ 0 ] = C [ 3 ] / D;
    pA [ 1 ] = pA [ 2 ] = -C [ 1 ] / D;
    pA [ 3 ] = C [ 0 ] / D;
}

/////////////////////////////////
// DenseMatrixArithmetics3x3
void DenseMatrixArithmetics3x3 :: SubATBproduct(double *pC, double *pa, double *pb)
{
    pC [ 0 ] -= pa [ 0 ] * pb [ 0 ] + pa [ 1 ] * pb [ 1 ] + pa [ 2 ] * pb [ 2 ];
    pC [ 1 ] -= pa [ 3 ] * pb [ 0 ] + pa [ 4 ] * pb [ 1 ] + pa [ 5 ] * pb [ 2 ];
    pC [ 2 ] -= pa [ 6 ] * pb [ 0 ] + pa [ 7 ] * pb [ 1 ] + pa [ 8 ] * pb [ 2 ];

    pC [ 3 ] -= pa [ 0 ] * pb [ 3 ] + pa [ 1 ] * pb [ 4 ] + pa [ 2 ] * pb [ 5 ];
    pC [ 4 ] -= pa [ 3 ] * pb [ 3 ] + pa [ 4 ] * pb [ 4 ] + pa [ 5 ] * pb [ 5 ];
    pC [ 5 ] -= pa [ 6 ] * pb [ 3 ] + pa [ 7 ] * pb [ 4 ] + pa [ 8 ] * pb [ 5 ];

    pC [ 6 ] -= pa [ 0 ] * pb [ 6 ] + pa [ 1 ] * pb [ 7 ] + pa [ 2 ] * pb [ 8 ];
    pC [ 7 ] -= pa [ 3 ] * pb [ 6 ] + pa [ 4 ] * pb [ 7 ] + pa [ 5 ] * pb [ 8 ];
    pC [ 8 ] -= pa [ 6 ] * pb [ 6 ] + pa [ 7 ] * pb [ 7 ] + pa [ 8 ] * pb [ 8 ];
}

/////////////////////////////////
// DenseMatrixArithmetics4x4
void DenseMatrixArithmetics4x4 :: SubATBproduct(double *pC, double *pa, double *pb)
{
    pC [ 0 ]  -= pa [ 0 ] * pb [ 0 ] + pa [ 1 ] * pb [ 1 ] + pa [ 2 ] * pb [ 2 ] + pa [ 3 ] * pb [ 3 ];
    pC [ 1 ]  -= pa [ 4 ] * pb [ 0 ] + pa [ 5 ] * pb [ 1 ] + pa [ 6 ] * pb [ 2 ] + pa [ 7 ] * pb [ 3 ];
    pC [ 2 ]  -= pa [ 8 ] * pb [ 0 ] + pa [ 9 ] * pb [ 1 ] + pa [ 10 ] * pb [ 2 ] + pa [ 11 ] * pb [ 3 ];
    pC [ 3 ]  -= pa [ 12 ] * pb [ 0 ] + pa [ 13 ] * pb [ 1 ] + pa [ 14 ] * pb [ 2 ] + pa [ 15 ] * pb [ 3 ];

    pC [ 4 ]  -= pa [ 0 ] * pb [ 4 ] + pa [ 1 ] * pb [ 5 ] + pa [ 2 ] * pb [ 6 ] + pa [ 3 ] * pb [ 7 ];
    pC [ 5 ]  -= pa [ 4 ] * pb [ 4 ] + pa [ 5 ] * pb [ 5 ] + pa [ 6 ] * pb [ 6 ] + pa [ 7 ] * pb [ 7 ];
    pC [ 6 ]  -= pa [ 8 ] * pb [ 4 ] + pa [ 9 ] * pb [ 5 ] + pa [ 10 ] * pb [ 6 ] + pa [ 11 ] * pb [ 7 ];
    pC [ 7 ]  -= pa [ 12 ] * pb [ 4 ] + pa [ 13 ] * pb [ 5 ] + pa [ 14 ] * pb [ 6 ] + pa [ 15 ] * pb [ 7 ];

    pC [ 8 ]  -= pa [ 0 ] * pb [ 8 ] + pa [ 1 ] * pb [ 9 ] + pa [ 2 ] * pb [ 10 ] + pa [ 3 ] * pb [ 11 ];
    pC [ 9 ]  -= pa [ 4 ] * pb [ 8 ] + pa [ 5 ] * pb [ 9 ] + pa [ 6 ] * pb [ 10 ] + pa [ 7 ] * pb [ 11 ];
    pC [ 10 ] -= pa [ 8 ] * pb [ 8 ] + pa [ 9 ] * pb [ 9 ] + pa [ 10 ] * pb [ 10 ] + pa [ 11 ] * pb [ 11 ];
    pC [ 11 ] -= pa [ 12 ] * pb [ 8 ] + pa [ 13 ] * pb [ 9 ] + pa [ 14 ] * pb [ 10 ] + pa [ 15 ] * pb [ 11 ];

    pC [ 12 ] -= pa [ 0 ] * pb [ 12 ] + pa [ 1 ] * pb [ 13 ] + pa [ 2 ] * pb [ 14 ] + pa [ 3 ] * pb [ 15 ];
    pC [ 13 ] -= pa [ 4 ] * pb [ 12 ] + pa [ 5 ] * pb [ 13 ] + pa [ 6 ] * pb [ 14 ] + pa [ 7 ] * pb [ 15 ];
    pC [ 14 ] -= pa [ 8 ] * pb [ 12 ] + pa [ 9 ] * pb [ 13 ] + pa [ 10 ] * pb [ 14 ] + pa [ 11 ] * pb [ 15 ];
    pC [ 15 ] -= pa [ 12 ] * pb [ 12 ] + pa [ 13 ] * pb [ 13 ] + pa [ 14 ] * pb [ 14 ] + pa [ 15 ] * pb [ 15 ];
}

/////////////////////////////////
// DenseMatrixArithmetics5x5
void DenseMatrixArithmetics5x5 :: SubATBproduct(double *pC, double *pA, double *pB)
{
    long i, j;
    for ( j = 0; j < bn; j++, pB += bn, pA -= bnbn ) {
        for ( i = 0; i < bn; i++, pC++, pA += bn ) {
            * pC -= pA [ 0 ] * pB [ 0 ] + pA [ 1 ] * pB [ 1 ] + pA [ 2 ] * pB [ 2 ] + pA [ 3 ] * pB [ 3 ] + pA [ 4 ] * pB [ 4 ];
        }
    }
}

/////////////////////////////////
// DenseMatrixArithmetics6x6
void DenseMatrixArithmetics6x6 :: SubATBproduct(double *pC, double *pA, double *pB)
{
    long i, j;
    for ( j = 0; j < bn; j++, pB += bn, pA -= bnbn ) {
        for ( i = 0; i < bn; i++, pC++, pA += bn ) {
            * pC -= pA [ 0 ] * pB [ 0 ] + pA [ 1 ] * pB [ 1 ] + pA [ 2 ] * pB [ 2 ] + pA [ 3 ] * pB [ 3 ] + pA [ 4 ] * pB [ 4 ] + pA [ 5 ] * pB [ 5 ];
        }
    }
}

/////////////////////////////////
// DenseMatrixArithmetics_Fake
void DenseMatrixArithmetics_Fake :: SubATBproduct(double *C, double *blockA, double *blockB)
{
    * C = * blockA * * blockB;
    //C;blockA;blockB;
}

DSS_NAMESPASE_END

