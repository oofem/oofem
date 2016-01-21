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

#ifndef _SPARSEGRIDMTX_H__
#define _SPARSEGRIDMTX_H__

#include "SparseGridColumn.h"
#include "SparseMatrixF.h"
#include "Ordering.h"
#include "BigMatrix.h"

DSS_NAMESPASE_BEGIN

/**
 * @author: Richard Vondracek
 */

class SparseGridMtx :
    public TraceableMatrix,
    public ILargeMatrix
{
public:
    // Allocates new space according to bskl and reads old matrix with respect
    // to permutation blockP
    SparseGridMtx(SparseMatrixF &sm, long block_size, Ordering *block_order, MathTracer *eMT);

    // Allocates new space according to bskl and reads old matrix with respect
    // to permutation blockP
    SparseGridMtx(SparseMatrixF &sm, long block_size, Ordering *block_order, Ordering *node_order, MathTracer *eMT);

    virtual ~SparseGridMtx();

public:
    long block_size;

    // this is block_size*block_size
    long block_storage;

    DenseMatrixArithmetics *BlockArith;
    SparseGridColumn **Columns;
    Ordering *node_order;

    Ordering *block_order;
    // permutation dof in block i is
    // in this matrix stored on block_permutation[i] block position

protected:
    long n;
    long n_blocks;

    long noDummyDOFs;

    long nonzeros;
    long blocks;

    // tells how many multipication have been done during the factorization
    long no_multiplications;

public:
    long N() const { return n; }
    long Nonzeros() const { return ( long ) columns_data_length; }

    // Number of nonempty blocks in whole matrix
    long Blocks()   { return blocks; }
    long BlockSize() { return block_size; }

    // Number of blocks in one row (n/blocksize)
    long N_blocks() { return n_blocks; }

    // How many block multipication have been done during the factorization
    long No_Multiplications()       { return no_multiplications; }

    // This data is used in the Sealed state
    long columns_data_length;


    void ComputeBlocks();
    double GetWaste();

    void AlocateMemoryByPattern(IConectMatrix *bskl);


    //ILargeMatrix
    void WriteStatistics(long no_init_blocks, long no_nonzeros);
    virtual void LoadMatrixNumbers(SparseMatrixF &sm) = 0;
    virtual void SolveLV(const LargeVector &b, LargeVector &x) = 0;
    virtual void Solve(double *b, double *x) = 0;
    virtual void MultiplyByVector(const LargeVectorAttach &x, LargeVectorAttach &y);
    virtual void times(double x);
    virtual void Factorize() = 0;

public:

    virtual void SchurComplementFactorization(int fixed_blocks) = 0;
    virtual void SolveA11(double *x, long fixed_blocks) = 0;
    virtual void Sub_A21_A11inv(double *x, long fixed_blocks) = 0;
    virtual void Sub_A11inv_A12(double *x, long fixed_blocks) = 0;
    virtual void WriteCondensedMatrixA22(double *a, Ordering *mcn, IntArrayList *lncn) = 0;
}; //class SparseGridMtx

DSS_NAMESPASE_END

#endif // _SPARSEGRIDMTX_H__
