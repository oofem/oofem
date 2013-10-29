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

#ifndef _MATRIX_H__
#define _MATRIX_H__

#include "DenseMatrix.h"
#include "IntArrayList.h"
#include "SparseMatrixF.h"

DSS_NAMESPASE_BEGIN

/**
 * @author Richard Vondracek
 */

// Simple mathematical methods
class Array2
{
public:

    inline static double ScalarProduct(double *a, double *b, long n)
    {
        double tmp = 0;
        for ( long i = 0; i < n; i++ ) {
            tmp += a [ i ] * b [ i ];
        }

        return tmp;
    }

    inline static void Add(const double *a, double *c, long n)
    {
        for ( long i = 0; i < n; i++ ) {
            c [ i ] += a [ i ];
        }
    }

    inline static void Add(const double *a, const double *b, double *c, long n)
    {
        for ( long i = 0; i < n; i++ ) {
            c [ i ] = a [ i ] + b [ i ];
        }
    }

    inline static void Add(const double *a, const double coeff_a, const double *b, const double coeff_b, double *c, long n)
    {
        for ( long i = 0; i < n; i++ ) {
            c [ i ] = coeff_a * a [ i ] + coeff_b * b [ i ];
        }
    }

    inline static void Add(const double *a, const double coeff_a, double *c, long n)
    {
        for ( long i = 0; i < n; i++ ) {
            c [ i ] += coeff_a * a [ i ];
        }
    }

    inline static void Mul(double *a, const double coeff_a, long n)
    {
        for ( long i = 0; i < n; i++ ) {
            a [ i ] *= coeff_a;
        }
    }
};

// interface for large vector Attached to pointer
class LargeVectorAttach
{
protected:
    long n;
    double *data;

public:
    long N() const { return n; }
    inline double *DataPtr() const { return data; }

    LargeVectorAttach();
    LargeVectorAttach(long N, double *data);
    ~LargeVectorAttach();

    LargeVectorAttach GetPermutedVector(long *perm) const;
    LargeVectorAttach GetPermuted_1Vector(long *perm) const;
    void GetPermutedVector(LargeVectorAttach *outV, IntArrayList &perm, long bn) const;
    void GetPermuted_1Vector(LargeVectorAttach *outV, IntArrayList &perm, long bn) const;
    void GetPermutedVector(LargeVectorAttach *outV, long *perm) const;

    void GetPermuted_1Vector(LargeVectorAttach *outV, long *perm) const;
    void Init(long N);
    void Grow(long addN);

    inline double operator[](long i) const
    {
        return data [ i ];
    }

    inline double &operator[](long i)
    {
        return data [ i ];
    }

    /*
     * static double operator | (LargeVector& A,LargeVector& B)
     * {
     * //System.Diagnostics.Debug.Assert(A.N==B.N,"Vector sizes differ");
     * double* pa = A.DataPtr(),*pb = B.DataPtr();
     * return DenseMatrix::InnerProduct(pa,pb,A.N());
     * }
     */

    void Add(const LargeVectorAttach &A);
    void Add(const LargeVectorAttach &A, long start_index, long length);
    void AddSmaler(const LargeVectorAttach &A);
    void Mult(double &alfa);
    void AddMult(double alfa, const LargeVectorAttach &B);
    void LinComb(double alfa, LargeVectorAttach &A, double beta, LargeVectorAttach &B);
    void LinComb(double alfa, LargeVectorAttach &A, LargeVectorAttach &B);
    void LinComb(double alfa, LargeVectorAttach &A);
    void LinComb(LargeVectorAttach &A);
    void Initialize(const LargeVectorAttach &A);
    void Initialize(const LargeVectorAttach &A, long start_index, long length);
    void Initialize(long start_index, long length);
    void Initialize(double val);
    void Initialize();

    double Norm();
    double NormSqr();

    void DiagonalSolve(double *b, double *x);

    static double InnerProduct(double *dp1, double *dp2, int len);

    // A += alfa*B
    static void AddMul(double *A, const double alfa, const double *B, int n);

    void LoadBinary(FILE *stream);
};

// interface for large vector
class LargeVector : public LargeVectorAttach
{
public:
    LargeVector();
    LargeVector(long N);
    LargeVector(const LargeVectorAttach &B);
    LargeVector(const LargeVector &B);

    void Detach();
    ~LargeVector();
};


////////////////// ////////////////////////////////////////////////////////////////////////////////////
// interface for any matrix, not very usefull
struct IMatrix
{
    virtual long N() const = 0;
    virtual long Nonzeros() const = 0;
    virtual ~IMatrix() {};
};

////////////////// ////////////////////////////////////////////////////////////////////////////////////
// interface for large matrix
struct ILargeMatrix : public IMatrix
{
    //double this [long i, long j]{get;set;}
    virtual double &ElementAt(int i, int j) = 0;

    virtual void WriteStatistics(long no_init_blocks, long no_nonzeros) = 0;
    virtual long No_Multiplications() = 0;

    virtual void SolveLV(const LargeVector &x, LargeVector &y) = 0;
    virtual void Solve(double *b, double *x) = 0;
    virtual void MultiplyByVector(const LargeVectorAttach &x, LargeVectorAttach &y) = 0;
    virtual void Factorize() = 0;
    virtual void LoadZeros() = 0;
    virtual void LoadMatrixNumbers(SparseMatrixF &sm) = 0;
};

////////////////// ////////////////////////////////////////////////////////////////////////////////////

class TraceableMatrix
{
private:
    char m_string [ 128 ];

public:
    MathTracer MT;
    MathTracer *eMT;

    TraceableMatrix();

    void Writeln(const char *cmd);
    void Write(const char *cmd);

    void CS();

    char *MC_();

protected:
    time_t temporary_measure_start;
    clock_t clock_start;
};

////////////////// ////////////////////////////////////////////////////////////////////////////////////
struct IConectMatrix : public IMatrix
{
    virtual IntArrayList *GetIndexesAboveDiagonalInColumn(long j) = 0;
    virtual IntArrayList *DetachIndexesAboveDiagonalInColumn(long j) = 0;
};


DSS_NAMESPASE_END

#endif //_MATRIX_H__
