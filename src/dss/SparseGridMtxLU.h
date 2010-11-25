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
 * Author: Richard Vondracek, <richard.vondracek@seznam.cz>
 */

// SparseGridMtxLU.h

#ifndef _SPARSEGRIDMTXLU_H__
#define _SPARSEGRIDMTXLU_H__

#include "SparseGridMtx.h"

DSS_NAMESPASE_BEGIN

/// <summary>
/// Summary description for SparseGridMtx.
/// </summary>
class SparseGridMtxLU :
    public SparseGridMtx
{
public:
    // Allocates new space according to bskl and reads old matrix with respect
    // to permutation blockP
    SparseGridMtxLU(SparseMatrixF &sm, BYTE block_size, Ordering *block_order, MathTracer *eMT, BOOL load_data = TRUE);

    // Allocates new space according to bskl and reads old matrix with respect
    // to permutation blockP
    SparseGridMtxLU(SparseMatrixF &sm, BYTE block_size, Ordering *block_order, Ordering *node_order, MathTracer *eMT, BOOL load_data = TRUE);

    virtual ~SparseGridMtxLU();

public:
    // This data is used in the Sealed state
    double *Columns_data;                               // data of the columns above the diagonal
    double *Rows_data;                                          // data of the rows below the diagonal
    double *Diagonal_data;                              // data of the diagonal blocks

    void AlocateMemoryByPattern(IConectMatrix *bskl);

    double GetValue(long bi, long bj, long si, long sj, long &aux_bi_idx, long &aux_bj_idx)
    {
        if ( bi == bj ) {
            return this->Diagonal_data [ bi * block_storage + si + sj * block_size ];
        } else
        if ( bj > bi ) {
            return Columns [ bj ]->GetValue(block_size, bi, si, sj, Columns_data, aux_bj_idx);
        } else {
            return Columns [ bi ]->GetValue(block_size, bj, sj, si, Rows_data, aux_bi_idx);
        }
    }

    void Solve(double *b, double *x);

    //ILargeMatrix
    virtual double &ElementAt(int i, int j);
    virtual void LoadZeros();
    virtual void LoadMatrixNumbers(SparseMatrixF &sm);
    virtual void SolveLV(const LargeVector &b, LargeVector &x);
    virtual void MultiplyByVector(const LargeVectorAttach &x, LargeVectorAttach &y);
    virtual void Factorize();

    void BackSubstU(double *x, long fixed_blocks);
    void ForwardSubstL(double *x, long fixed_blocks);

    // x = A^(-1) * b
    void SolveLU(double *x, long fixed_blocks = 0);

    void SubMultU12(double *px, double *py, long fixed_blocks);
    void SubMultL21(double *px, double *py, long fixed_blocks);

public:
    // Schur complement solution methods
    virtual void SchurComplementFactorization(int fixed_blocks);
    virtual void SolveA11(double *x, long fixed_blocks);
    virtual void Sub_A21_A11inv(double *x, long fixed_blocks);
    virtual void Sub_A11inv_A12(double *x, long fixed_blocks);
    virtual void WriteCondensedMatrixA22(double *a, Ordering *mcn, IntArrayList *lncn);
}; //class SparseGridMtx

DSS_NAMESPASE_END

#endif // _SPARSEGRIDMTXLU_H__
