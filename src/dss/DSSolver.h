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

#ifndef _DSSOLVER_H__
#define _DSSOLVER_H__

#include "DSSAfx.h"
#include "SparseConectivityMtx.h"
#include "IntArrayList.h"
//#include "SparseGridMtx.h"
#include "SparseGridMtxPD.h"

DSS_NAMESPASE_BEGIN

enum eDSMatrixType {
    eDSSparseMatrix = 0,
    eDSSkylineMatrix = 1
};

enum eDSSolverType {
    eDSSFactorizationLDLT = 0,
    eDSSFactorizationLLT = 1,
    eDSSFactorizationLU = 2,
    eDSSFactorizationLDLTIncomplete = 3,
    eDSSFactorizationLLTIncomplete = 4,
    eDSSDiagonalScaling = 5,
    eDSSFastCG = 6
};

/**
 * ISolver is the user interface for the DirectSparseSolver
 *
 * @author Richard Vondracek
 */

struct ISolver
{
    // Makes necessary intialization
    // 0 - produces sparse symmetrical LDL^T factorization (default)
    // 1 - produces sparse symmetrical LLT Cholesky factorization
    // 2 - produces sparse LU factorization on symmetrical pattern -> undirected graph
    // 3 - produces sparse incomplete LDL^T factorization
    virtual long Initialize(unsigned char run_code, eDSSolverType solverType = eDSSFactorizationLDLT, eDSMatrixType matrixType = eDSSparseMatrix) = 0;

    // Sets the type of the reordering of the sparse matrix to something else than default
    // returns true if possible
    virtual bool SetOrderingType(Ordering :: Type otype) = 0;

    // Registers the source to the DSSolver
    virtual bool LoadMatrix(unsigned long neq, unsigned char block_size, double *a, unsigned long *ci, unsigned long *adr) = 0;
    virtual bool LoadMatrix(SparseMatrixF *sm, unsigned char block_size) = 0;
    virtual bool SetMatrixPattern(SparseMatrixF *smt, unsigned char block_size) = 0;

    // Loads the "MatrixCodeNumbers" vector[n_blocks*block_size]
    // each entry in the vector[i] can be :
    // vector[i] >=  0  normal unknown (means the corresponding row in SparseMatrixF [zero-based])
    // vector[i] == -1  this entry is not used (no corresponding row in SparseMatrixF)
    // vector[i] <= -2  unknown to be left uncondensed ( -vector[i]-2 is the zero-based row index in SparseMatrixF)

    //  Ex:
    // 0,1,2,....... code numbers
    // -1   ........ eliminated row
    // -2, -3, -4 .. rows (0,1,2) ment to be left uncondensed
    //  n_blocks - number of blocks

    virtual bool LoadMCN(unsigned long n_blocks, unsigned char block_size, long *mcn) = 0;


    // Creates the connectivity matrix
    // Computes the MinimumDegree ordering
    // Allocates the memory for the SparseGrid matrix
    virtual bool PreFactorize() = 0;

    // Loads or reloads the numeric values of the sparse matrix 'sm' to the already allocated SparseGrid matrix
    virtual bool LoadNumbers(SparseMatrixF *sm) = 0;

    // Runs the LDL factorization on the already allocated and loaded matrix
    virtual bool ReFactorize() = 0;

    // This function unites the three previous function into one step if necessary
    // PreFactorize ( ) + LoadNubers( ) + ReFactorize ( )
    // This is usefull for single-pass solutions
    virtual bool Factorize() = 0;

    // When the matrix was factorized we can solve several A*r=f equations
    // if the pointers are equal (r==f) the result will overwrite the RHS vector (In that case the soultion is faster)
    // It is recomended that both vectors have size in a whole multiple of the 'block_size' or more
    virtual bool Solve(double *r, double *f) = 0;

    // When the matrix was factorized we solve r from A*r=f equation
    // where f is vector with units. And we return (f-A*r)|(f-A*r)/N value.
    virtual double GetFactorizationError() = 0;

    // Disposes the solver internal data from the memory for us to be able to load another matrix.
    virtual void Dispose() = 0;

    // Contains calls Dispose function.
    virtual long Close() = 0;

    virtual void SetMT(MathTracer *MT) = 0;

    virtual bool decomp() = 0;
    virtual void changedecomp() = 0;

    virtual bool IsFactorized() = 0;

    virtual ~ISolver() {};

    // |A11 A12| |r1| = |f1|
    // |A21 A22| |r2|   |f2|

    //  a - condensed matrix
    //  lhs - left hand side
    //  rhs - right hand side
    //  tc - type of computation
    //  tc=1 - return condensed matrix a=(A22-A21*A11inv*A12) and modify rhs part f2 = (f2-A21*A11inv*f1)
    //  tc=2 - return part of left hand side
    //  tc=3 - return (Kii lhs = rhs) solution
    virtual void condense(double *a, double *lhs, double *rhs, long tc) = 0;

    // A few words to different ordering systems used during sparse PD solution
    //
    // A-order [neq]
    //    This is the order of data stored in user specified SparseMatrixF* sm.
    //    Usually, we get RHS and return LHS in this order system.
    //
    // B-order [neq] -> [n_blocks*block_size]
    //    You obtain B-order after applying MCN vector on A-order.
    //    MCN re-introduces removed lines and DOFs from certain nodes.
    //
    //  ex:    for (i=0;i<neq;i++) lhsA[i] = tmpB[mcn->order[i]];
    //      for (i=0;i<neq;i++) tmpB[mcn->order[i]] = lhsA[i];
    //
    // C-order [neq]
    //    This is such a permutaion of A order where all noncondensed nodes are on the end of
    //    the order.
    //          dom_order

    virtual void GetA12block(double *pA12) = 0;

    // c = A * b
    virtual void MulMatrixByVector(double *b, double *c) = 0;
    
    // Copute x by Preconditioned Conjugate gradient method
    virtual int PreCG(double *b, double *x, double epsilon, int max_iter) = 0;

    // Copute x by Conjugate gradient method
    virtual int CG(double *b, double *x, double epsilon, int max_iter) = 0;

    enum eState {
        None = 0,
        Initialized = 1,
        Allocated = 2,
        Factorized = 3,
        ErrorInMCN = -1
    };
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////

class DSSolver : public ISolver
{
    MathTracer MT;
    MathTracer *eMT;

    SparseMatrixF sm;
    SparseGridMtx *matrix;
    SparseGridMtxPD *matrixPD;
    char str [ 512 ];
    int decompid;
    ISolver :: eState m_eState;

    eDSSolverType SolverType;
    eDSMatrixType MatrixType;
    unsigned char run_code;

    // PDD data
    double *tmpR;
    IntArrayList *dom_order;

    // blocked matrix
    unsigned char blockSize;
    long n_blocks;
    long neq;
    Ordering *mcn;

    // Diagonal precoditioner
    LargeVector matrixD;
    SparseGridMtx *orig_matrix;     // blockmatrix for CG multiplication

    // Dirichlets conditions - noncondensed DOFs
    IntArrayList *fixed;                // noncondensed blocks
    IntArrayList *lncn;                         // noncondensed DOFs

    Ordering :: Type OrderingType;      //MinimumDegree,ApproxMinimumDegree,ApproxMinimumDegreeIncomplete,...

public:
    DSSolver(MathTracer *pMT = NULL);

    virtual ~DSSolver();
    virtual long Initialize(unsigned char run_code, eDSSolverType solverType = eDSSFactorizationLDLT, eDSMatrixType matrixType  = eDSSparseMatrix);
    virtual void SetMT(MathTracer *pMT);
    virtual void Dispose();

    virtual bool SetOrderingType(Ordering :: Type otype);
    virtual bool LoadMatrix(unsigned long neq, unsigned char block_size, double *a, unsigned long *ci, unsigned long *adr);
    virtual bool LoadMatrix(SparseMatrixF *smt, unsigned char block_size);
    virtual bool SetMatrixPattern(SparseMatrixF *smt, unsigned char block_size);

    virtual bool decomp();
    virtual void changedecomp();
    virtual bool IsFactorized();

    //multiply by a scalar
    virtual void times(double x);
    
    virtual bool Factorize();
    virtual bool PreFactorize();

    virtual void LoadZeros();
    virtual bool LoadNumbers(SparseMatrixF *sm);
    virtual double &ElementAt(int i, int j);

    virtual bool ReFactorize();

    virtual bool Solve(double *r, double *f);
    virtual long Close();

    void StartSolverWriteInfo();
    void EndSolverWriteInfo();
    bool StartSolver();

    void SetSM(SparseMatrixF *sm);
    virtual double GetFactorizationError();

    // Loads the "MatrixCodeNumbers" vector[n_blocks*block_size]
    // each entry in the vector[i] can be :
    // vector[i] >=  0  normal unknown      (means the corresponding row in SparseMatrixF)
    // vector[i] == -1  this entry is not used  (no corresponding row in SparseMatrixF)
    virtual bool LoadMCN(unsigned long n_blocks, unsigned char block_size, long *mcn);

    virtual bool LoadMCN(IntArrayList &mcn);

    IntArrayList *GetFixedBlocks() { return this->fixed; }
    IntArrayList *GetFixedDOFs()   { return this->lncn; }

    virtual void condense(double *a, double *lhs, double *rhs, long tc);

    void GetA12block(double *pA12);

    // c = A * b
    void MulMatrixByVector(double *b, double *c);

    int PreCG(double *b, double *x, double epsilon, int max_iter);

    // Copute x by Conjugate gradient method
    int CG(double *b, double *x, double epsilon, int max_iter);

private:
    bool LoadMCN_int(IntArrayList *mcn_order);
    void ExpandMCN(IntArrayList &mcn);
    SparseGridMtx *CreateNewSparseGridMtx(IntArrayList *fixed = NULL);
    void WriteFactorizationInfo();
    bool CreateFixedArray(long no_noncondensed_DOFs);
    bool FactorizeSchur();
    void StoreFixedLastPermutation_dom_order();
};

DSS_NAMESPASE_END

#endif //_DSSOLVER_H__
