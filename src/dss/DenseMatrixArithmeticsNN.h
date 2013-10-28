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

#ifndef _DENSEMATRIXARITH_H__
#define _DENSEMATRIXARITH_H__

#include "MathTracer.h"
#include "Array.h"

DSS_NAMESPASE_BEGIN

enum ePreferedDecomposition {
    eLDL_decomposition,
    eLL_decomposition
};

/**
 * @author Richard Vondracek
 */

class DenseMatrixArithmetics
{
public:
    long bn;
    long bnbn;
    long bn1;
    static long zero_pivots;

    ePreferedDecomposition prefered_decomposition;
    MathTracer MT;
    MathTracer *eMT;

private:
    double tmp;
    double *p;

public:
    DenseMatrixArithmetics(long bn);

    virtual ~DenseMatrixArithmetics();

    static DenseMatrixArithmetics *NewArithmetics(long block_size);

    // b += B*x
    virtual void AddMultBlockByVectorSym(double *B, double *x, double *b, long bi_bn, long bj_bn);

    // b -= A*x
    virtual void SubMultBlockByVector(double *B, double *x, double *b);

    // b -= A^T * x
    virtual void SubMultTBlockByVector(double *B, double *x, double *b);

    // b += B * x
    virtual void MultDiagonalBlockByVector(double *B, double *x, double *b);

    // C -= A^T * B
    // This is the most SPARSEDIRECT SOLVER innerloop operation
    virtual void SubATBproduct(double *pC, double *pA, double *pB);

    /// <summary>
    /// Solves LL' decomposition of symmetric pC column based matrix
    /// The result is stored in lower half of the matrix pC
    /// </summary>
    /// <param name="pC"> matrix to be factorized</param>
    /// <param name="p"> here we store the diagonal elements of L</param>
    virtual void Cholesky_Decomposition(double *pC, double * &p);

    /// <summary> Solve LL^T x = b </summary>
    /// <param name="pC">factorized L matrix (lower triangle)</param>
    /// <param name="p">diagonal elements of L</param>
    /// <param name="b">right side</param>
    /// <param name="x">unknowns</param>
    virtual void Cholesky_Solve(double *pC, double *p, double *b, double *x);

    virtual void Cholesky_Linv(double *pC, double *p);

    virtual void ComputeInversionByCholesky(double *C, double *Inv);

    /// <summary>
    /// Solves LL' decomposition of symmetric pC column based matrix
    /// The result is stored in lower half of the matrix pC
    /// </summary>
    /// <param name="pC"> matrix to be factorized</param>
    virtual void LL_Decomposition(double *pC);

    /// <summary> Solve LL^T x = b </summary>
    /// <param name="pC">factorized L matrix (lower triangle)</param>
    /// <param name="b">right side</param>
    /// <param name="x">unknowns</param>
    virtual void LL_Solve(double *pC, double *b, double *x);

    /// <summary> Solve LL^T x = b </summary>
    /// <param name="pC">factorized L matrix (lower triangle)</param>
    /// <param name="b">right side</param>
    /// <param name="x">unknowns</param>
    virtual void L_BlockSolve(double *C, double *B);

    virtual void SubstSolveL(double *pC, double *x);

    virtual void SubstSolveLT(double *pC, double *x);

    virtual void LU_Decomposition(double *C);

    void LU_Solve(double *C, double *b);

    // Solves system X*LU=B
    // B is stored by rows
    void ULT_BlockSolve(double *C, double *B);

    void LDL_Decomposition(double *A);

    void LDL_Solve(double *A, double *b, const long startIndex = 0);

    // A - decomposed matrix
    // B - block of right sides
    virtual void SubstSolveBlock(double *A, double *B);

    virtual void FactorizeBlock(double *A);

    // Solves A b = b
    virtual void SubstSolve(double *A, double *b);

    void ComputeInversionByLDL(double *C, double *Inv);

    virtual void GetInversion(double *C, double *blockA);

    virtual bool GaussElimination(double *C, double *blockA);
};

class DenseMatrixArithmetics1x1 : public DenseMatrixArithmetics
{
public:
    DenseMatrixArithmetics1x1(long n) : DenseMatrixArithmetics(n) {}

    void SubATBproduct(double *pC, double *pA, double *pB);
    void GetInversion(double *pC, double *pA);
    virtual void FactorizeBlock(double *A);
    virtual void SubstSolveBlock(double *A, double *B);
    virtual void SubstSolve(double *A, double *b);
};

class DenseMatrixArithmetics2x2 : public DenseMatrixArithmetics
{
public:
    DenseMatrixArithmetics2x2(long n) : DenseMatrixArithmetics(n) {}
    void SubATBproduct(double *pC, double *pa, double *pb);
    void GetInversion(double *C, double *pA);
};

class DenseMatrixArithmetics3x3 : public DenseMatrixArithmetics
{
public:
    DenseMatrixArithmetics3x3(long n) : DenseMatrixArithmetics(n) {}
    void SubATBproduct(double *pC, double *pa, double *pb);
};

class DenseMatrixArithmetics4x4 : public DenseMatrixArithmetics
{
public:
    DenseMatrixArithmetics4x4(long n) : DenseMatrixArithmetics(n)     {}
    void SubATBproduct(double *pC, double *pa, double *pb);
};

class DenseMatrixArithmetics5x5 : public DenseMatrixArithmetics
{
public:
    DenseMatrixArithmetics5x5(long n) : DenseMatrixArithmetics(n)     {}
    void SubATBproduct(double *pC, double *pA, double *pB);
};

class DenseMatrixArithmetics6x6 : public DenseMatrixArithmetics
{
public:
    DenseMatrixArithmetics6x6(long n) : DenseMatrixArithmetics(n) {}

    void SubATBproduct(double *pC, double *pA, double *pB);
};

class DenseMatrixArithmetics_Fake : public DenseMatrixArithmetics
{
public:
    DenseMatrixArithmetics_Fake(long n) : DenseMatrixArithmetics(n) {}
    void SubATBproduct(double *C, double *blockA, double *blockB);
};

DSS_NAMESPASE_END

#endif //_DENSEMATRIXARITH_H__

