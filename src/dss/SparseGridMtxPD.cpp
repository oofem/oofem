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

#include "SparseGridMtxPD.h"

DSS_NAMESPASE_BEGIN

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxPD :: SparseGridMtxPD(SparseMatrixF &sm, long block_size, Ordering *block_order, long fixed_blocks, MathTracer *MT)
{
    this->m_pMatrix = new SparseGridMtxLDL(sm, block_size, block_order, MT);
    this->m_lFixed_blocks = fixed_blocks;
}

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxPD :: SparseGridMtxPD(SparseMatrixF &sm, long block_size, Ordering *block_order, Ordering *node_order, long fixed_blocks, MathTracer *MT)
{
    this->m_pMatrix = new SparseGridMtxLDL(sm, block_size, block_order, node_order, MT);
    this->m_lFixed_blocks = fixed_blocks;
}

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SparseGridMtxPD :: SparseGridMtxPD(SparseGridMtx *pMatrix, long fixed_blocks)
{
    this->m_pMatrix = pMatrix;
    this->m_lFixed_blocks = fixed_blocks;
}

SparseGridMtxPD :: ~SparseGridMtxPD()
{
    delete m_pMatrix;
}

SparseGridMtx *SparseGridMtxPD :: Matrix()
{
    return m_pMatrix;
}

// This functions computes LDL' decomposition of first (nblocks-fixed_bn) columns
// The resulting submatrix will contain the schur complement A22 - A21 * A11^(-1) * A12
void SparseGridMtxPD :: SchurComplementFactorization()
{
    m_pMatrix->SchurComplementFactorization(m_lFixed_blocks);
}

void SparseGridMtxPD :: SolveA11(double *x)
{
    m_pMatrix->SolveA11(x, m_lFixed_blocks);
}

void SparseGridMtxPD :: Sub_A21_A11inv(double *x)
{
    m_pMatrix->Sub_A21_A11inv(x, m_lFixed_blocks);
}

void SparseGridMtxPD :: Sub_A11inv_A12(double *x)
{
    m_pMatrix->Sub_A11inv_A12(x, m_lFixed_blocks);
}

void SparseGridMtxPD :: WriteCondensedMatrixA22(double *a, Ordering *mcn, IntArrayList *lncn)
{
    m_pMatrix->WriteCondensedMatrixA22(a, mcn, lncn);
}

DSS_NAMESPASE_END

