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

#ifndef _DENSEMATRIX_H__
#define _DENSEMATRIX_H__

#include "DenseMatrixArithmeticsNN.h"

DSS_NAMESPASE_BEGIN

/**
 * Dense matrix stored by columns
 *
 * @author Richard Vondracek
 */

class DenseMatrix
{
public:
    long n;
    double *data;               // Matrix stored by columns

private:
    DenseMatrixArithmetics *dma;

public:
    DenseMatrixArithmetics &DMA()
    {
        if ( dma == NULL ) {
            dma = DenseMatrixArithmetics :: NewArithmetics(n);
        }

        return * dma;
    }

    /*public double this[long i, long j]
     * {
     * get{return data[i+n*j];}
     * set{data[i+n*j] = value;}
     * }*/

    void Add(long i, long j, double val)
    {
        data [ i + n * j ] += val;
    }

    DenseMatrix(long n)
    {
        dma = NULL;
        this->n = n;
        data = new double [ n * n ];
        memset( data, 0, n * n * sizeof( double ) );
    }


    DenseMatrix(long n, double d)
    {
        dma = NULL;
        this->n = n;
        data = new double [ n * n ];
        for ( long i = 0; i < n; i++ ) { data [ i + i * n ] = d; }
    }

    DenseMatrix(double *dataFrom, long start_idx, long n)
    {
        dma = NULL;
        this->n = n;
        data = new double [ n * n ];
        Array :: Copy(dataFrom, start_idx, data, 0, n * n);
    }

    ~DenseMatrix()
    {
        if ( dma ) {
            delete dma;
        }

        if ( data ) {
            delete [] data;
        }
    }



    void CopyTo(DenseMatrix &blockA, long bn)
    {
        Array :: Copy(this->data, blockA.data, bn * bn);
    }

    void CopyTo(double *dataTo, long start_idx)
    {
        Array :: Copy(data, 0, dataTo, start_idx, n * n);
    }

    void Clear()
    {
        Array :: Clear(data, 0, n * n);
    }

    static double InnerProduct(double *dp1, double *dp2, long len)
    {
        long i, len4;
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
};

DSS_NAMESPASE_END

#endif //_DENSEMATRIX_H__
