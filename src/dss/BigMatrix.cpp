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

#include "BigMatrix.h"

DSS_NAMESPASE_BEGIN


LargeVectorAttach :: LargeVectorAttach()
{
    data = NULL;
    n = 0;
}

LargeVectorAttach :: LargeVectorAttach(long N, double *data)
{
    this->data = data;
    n = N;
}

LargeVectorAttach :: ~LargeVectorAttach()
{}

LargeVectorAttach LargeVectorAttach :: GetPermutedVector(long *perm) const
{
    LargeVector tmp(n);
    for ( long i = 0; i < n; i++ ) {
        tmp [ i ] = data [ perm [ i ] ];
    }

    return tmp;
}

LargeVectorAttach LargeVectorAttach :: GetPermuted_1Vector(long *perm) const
{
    LargeVector tmp(n);
    for ( long i = 0; i < n; i++ ) {
        tmp [ perm [ i ] ] = data [ i ];
    }

    return tmp;
}

void LargeVectorAttach :: GetPermutedVector(LargeVectorAttach *outV, IntArrayList &perm, long bn) const
{
    if ( outV == NULL ) {
        outV = new LargeVector(perm.Count * bn);
    }

    for ( long i = outV->N() - 1; i >= 0; i-- ) {
        ( * outV ) [ i ] = data [ perm [ i / bn ] * bn + ( i % bn ) ];
    }
}

void LargeVectorAttach :: GetPermuted_1Vector(LargeVectorAttach *outV, IntArrayList &perm, long bn) const
{
    if ( outV == NULL ) {
        outV = new LargeVector(perm.Count * bn);
    }

    for ( long i = 0; i < n; i++ ) {
        ( * outV ) [ perm [ i / bn ] * bn + ( i % bn ) ] = data [ i ];
    }
}

void LargeVectorAttach :: GetPermutedVector(LargeVectorAttach *outV, long *perm) const
{
    if ( outV == NULL ) {
        outV = new LargeVector(n);
    }

    for ( long i = outV->N() - 1; i >= 0; i-- ) {
        ( * outV ) [ i ] = data [ perm [ i ] ];
    }
}

void LargeVectorAttach :: GetPermuted_1Vector(LargeVectorAttach *outV, long *perm) const
{
    if ( outV == NULL ) {
        outV = new LargeVector(n);
    }

    for ( long i = 0; i < n; i++ ) {
        ( * outV ) [ perm [ i ] ] = data [ i ];
    }
}

void LargeVectorAttach :: Init(long N)
{
    if ( data ) {
        delete [] data;
    }

    data = new double [ N ];
    n = N;
}

void LargeVectorAttach :: Grow(long addN)
{
    double *oldData = this->data;
    data = new double [ n + addN ];
    Array :: Copy(oldData, data, n);
    if ( oldData ) {
        delete oldData;
    }

    n += addN;
}

/*
 * static double operator | (LargeVectorAttach& A,LargeVectorAttach& B)
 * {
 * //System.Diagnostics.Debug.Assert(A.N==B.N,"Vector sizes differ");
 * double* pa = A.DataPtr(),*pb = B.DataPtr();
 * return DenseMatrix::InnerProduct(pa,pb,A.N());
 * }
 */

void LargeVectorAttach :: Add(const LargeVectorAttach &A)
{
    for ( long i = n - 1; i >= 0; i-- ) {
        data [ i ] += A [ i ];
    }
}

void LargeVectorAttach :: Add(const LargeVectorAttach &A, long start_index, long length)
{
    for ( long i = start_index; i < start_index + length; i++ ) {
        data [ i ] += A [ i ];
    }
}

void LargeVectorAttach :: AddSmaler(const LargeVectorAttach &A)
{
    for ( long i = A.N() - 1; i >= 0; i-- ) {
        data [ i ] += A [ i ];
    }
}

void LargeVectorAttach :: Mult(double &alfa)
{
    for ( long i = n - 1; i >= 0; i-- ) {
        data [ i ] *= alfa;
    }
}

void LargeVectorAttach :: AddMult(double alfa, const LargeVectorAttach &B)
{
    /*
     * //System.Diagnostics.Debug.Assert(B.N==this.N,"Vector sizes differ");
     * double* pa,*pb;
     * double* pA = this->DataPtr(),*pB = B.DataPtr();
     * {
     *      pa = pA;pb = pB;
     *      for (long i=n-1; i>=0; i-- )
     * *pa++ += (alfa) * (*pb++);
     * }
     */
    double *dp1 = B.DataPtr();
    double *d = data;

    int i, len4, len;
    len4 = n / 4;
    len  = n % 4;

    for ( i = 0; i < len4; i++ ) {
        d [ 4 * i ]   += alfa * dp1 [ 4 * i ];
        d [ 4 * i + 1 ] += alfa * dp1 [ 4 * i + 1 ];
        d [ 4 * i + 2 ] += alfa * dp1 [ 4 * i + 2 ];
        d [ 4 * i + 3 ] += alfa * dp1 [ 4 * i + 3 ];
    }

    dp1 += 4 * len4;
    d += 4 * len4;

    for ( i = 0; i < len; i++ ) {
        d [ i ] += alfa * ( * dp1++ );
    }
}

void LargeVectorAttach :: LinComb(double alfa, LargeVectorAttach &A, double beta, LargeVectorAttach &B)
{
    for ( long i = n - 1; i >= 0; i-- ) {
        data [ i ] = alfa * A [ i ] + beta * B [ i ];
    }
}

void LargeVectorAttach :: LinComb(double alfa, LargeVectorAttach &A, LargeVectorAttach &B)
{
    for ( long i = n - 1; i >= 0; i-- ) {
        data [ i ] = alfa * A [ i ] + B [ i ];
    }

    /*
     *      double* dp1 = A.DataPtr();
     *      double* dp2 = B.DataPtr();
     *      double* d = data;
     *
     *      double s0,s1,s2,s3;
     *
     *      int i,len4,len;
     *      len4 = n / 4;
     *      len  = n % 4;
     *
     *      for ( i = 0; i < len4; i++ )
     *      {
     *              s0 = alfa*dp1[4*i]+dp2[4*i];
     *              s1 = alfa*dp1[4*i+1]+dp2[4*i+1];
     *              s2 = alfa*dp1[4*i+2]+dp2[4*i+2];
     *              s3 = alfa*dp1[4*i+3]+dp2[4*i+3];
     *              d[4*i] = s0;
     *              d[4*i+1] = s1;
     *              d[4*i+2] = s2;
     *              d[4*i+3] = s3;
     *      }
     *      dp1 += 4*len4; dp2 += 4*len4; d +=4*len4;
     *
     *      for ( i = 0; i < len; i++ )
     *              d[i] = alfa*(*dp1++)+(*dp2++);
     */
}

void LargeVectorAttach :: LinComb(double alfa, LargeVectorAttach &A)
{
    for ( long i = n - 1; i >= 0; i-- ) {
        data [ i ] = alfa * A [ i ];
    }
}

void LargeVectorAttach :: LinComb(LargeVectorAttach &A)
{
    for ( long i = n - 1; i >= 0; i-- ) {
        data [ i ] = A [ i ];
    }
}

void LargeVectorAttach :: Initialize(const LargeVectorAttach &A)
{
    Array :: Copy(A.data, data, n);
}

void LargeVectorAttach :: Initialize(const LargeVectorAttach &A, long start_index, long length)
{
    Array :: Copy(A.data, start_index, data, start_index, length);
}

void LargeVectorAttach :: Initialize(long start_index, long length)
{
    Array :: Clear(data, start_index, length);
}

void LargeVectorAttach :: Initialize(double val)
{
    for ( long i = n - 1; i >= 0; i-- ) {
        data [ i ] = val;
    }
}

void LargeVectorAttach :: Initialize()
{
    Array :: Clear(data, 0, n);
}

double LargeVectorAttach :: Norm()
{
    return sqrt( NormSqr() );
}

double LargeVectorAttach :: NormSqr()
{
    double tmp = 0;
    for ( long i = n - 1; i >= 0; i-- ) {
        tmp += pow(data [ i ], 2);
    }

    return tmp;
}

void LargeVectorAttach :: DiagonalSolve(double *b, double *x)
{
    for ( long i = n - 1; i >= 0; i-- ) {
        x [ i ] = b [ i ] / data [ i ];
    }
}

double LargeVectorAttach :: InnerProduct(double *dp1, double *dp2, int len)
{
    int i, len4;
    double sum0, sum1, sum2, sum3;

    sum0 = sum1 = sum2 = sum3 = 0.0;

    len4 = len / 4;
    len  = len % 4;

    for ( i = 0; i < len4; i++ ) {
        sum0 += dp1 [ 4 * i ] * dp2 [ 4 * i ];
        sum1 += dp1 [ 4 * i + 1 ] * dp2 [ 4 * i + 1 ];
        sum2 += dp1 [ 4 * i + 2 ] * dp2 [ 4 * i + 2 ];
        sum3 += dp1 [ 4 * i + 3 ] * dp2 [ 4 * i + 3 ];
    }

    sum0 += sum1 + sum2 + sum3;
    dp1 += 4 * len4;
    dp2 += 4 * len4;

    for ( i = 0; i < len; i++ ) {
        sum0 += ( * dp1++ ) * ( * dp2++ );
    }

    return sum0;
}

// A += alfa*B
void LargeVectorAttach :: AddMul(double *A, const double alfa, const double *B, int n)
{
    for ( long i = n - 1; i >= 0; i-- ) {
        * A++ += ( alfa ) * ( * B++ );
    }
}

void LargeVectorAttach :: LoadBinary(FILE *stream)
{
    if ( data ) {
        delete [] data;
    }

    data = NULL;

    fread(& n, sizeof( n ), 1, stream);
    data = new double [ n ];
    fread(data, sizeof( double ), n, stream);
}

///////////////////////////////////////////////////////////////////////
///  TraceableMatrix

LargeVector :: LargeVector()
{
    data = NULL;
    n = 0;
}

LargeVector :: LargeVector(long N)
{
    data = new double [ N ];
    memset( data, 0, N * sizeof( double ) );
    n = N;
}

LargeVector :: LargeVector(const LargeVectorAttach &B)
{
    data = new double [ B.N() ];
    this->n = B.N();
    Array :: Copy(B.DataPtr(), data, n);
}

LargeVector :: LargeVector(const LargeVector &B)
{
    data = new double [ B.N() ];
    this->n = B.N();
    Array :: Copy(B.DataPtr(), data, n);
}

void LargeVector :: Detach()
{
    data = NULL;
}

LargeVector :: ~LargeVector()
{
    if ( data ) {
        delete [] data;
    }

    data = NULL;
}


///////////////////////////////////////////////////////////////////////
///  TraceableMatrix

TraceableMatrix :: TraceableMatrix()
{
    eMT = & MT;
}

void TraceableMatrix :: Writeln(const char *cmd)
{
    if ( eMT ) {
        eMT->Writeln(cmd);
    }
}
void TraceableMatrix :: Write(const char *cmd)
{
    if ( eMT ) {
        eMT->Write(cmd);
    }
}

void TraceableMatrix :: CS()
{
    //temporary_measure_start = CTime::GetCurrentTime();
    time(& temporary_measure_start);
    clock_start = clock();
}

char *TraceableMatrix :: MC_()
{
    clock_t end = clock();
    double duration = ( double ) ( end - clock_start ) / CLOCKS_PER_SEC;
    sprintf(m_string, "...%0.3f s", duration);
    return m_string;
}

DSS_NAMESPASE_END
