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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef _SKYLINEMTXLDL_H__
#define _SKYLINEMTXLDL_H__

#include "SkyLineMtx.h"

DSS_NAMESPASE_BEGIN

/**
 * @author: Richard Vondracek
 */

class SkyLineMtxLDL :
    public SkyLineMtx
{
public:
    SkyLineMtxLDL(SparseMatrixF &sm, Ordering *order, MathTracer *eMT);
    virtual ~SkyLineMtxLDL();

    void LoadMatrixData(SparseMatrixF &sm);


public:
    virtual void Solve(double *b, double *x) override;

    //ILargeMatrix
    virtual double &ElementAt(int i, int j) override;
    virtual void LoadZeros() override;
    virtual void LoadMatrixNumbers(SparseMatrixF &sm) override;
    virtual void SolveLV(const LargeVector &b, LargeVector &x) override;
    virtual void Factorize() override;
    virtual void MultiplyByVector(const LargeVectorAttach &x, LargeVectorAttach &y) override;

public:
    virtual void SchurComplementFactorization(int fixed_blocks) override;
    virtual void SolveA11(double *x, long fixed_blocks) override;
    virtual void Sub_A21_A11inv(double *x, long fixed_blocks) override;
    virtual void Sub_A11inv_A12(double *x, long fixed_blocks) override;
    virtual void WriteCondensedMatrixA22(double *a, Ordering *mcn, IntArrayList *lncn) override;
}; //class SkyLineMtxLDL

DSS_NAMESPASE_END

#endif // _SKYLINEMTXLDL_H__
