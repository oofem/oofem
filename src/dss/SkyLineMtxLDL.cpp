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

#include "SkyLineMtxLDL.h"

DSS_NAMESPASE_BEGIN

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SkyLineMtxLDL :: SkyLineMtxLDL(SparseMatrixF &sm, Ordering *order, MathTracer *eMT) : SkyLineMtx(sm, order, eMT)
{
    nonzeros = 0;
    this->eMT = eMT;
    this->order = order;
}

void SkyLineMtxLDL :: LoadZeros()
{}

void SkyLineMtxLDL :: Solve(double *b, double *x)
{}

double dummy = 0;
double &SkyLineMtxLDL :: ElementAt(int i, int j)
{
    return dummy;
}

void SkyLineMtxLDL :: LoadMatrixData(SparseMatrixF &sm)
{
    long *nodeP = this->order ? order->perm->Items : NULL;

    // diagonal
    for ( unsigned long j = 0; j < sm.neq; j++ ) {
        long nj = ( nodeP == NULL ) ? j : nodeP [ j ];
        int col_start = column_starts [ nj ] + nj - 1;

        for ( unsigned long ad = sm.Adr(j); ad < sm.Adr(j + 1); ad++ ) {
            long i = ( long ) sm.Ci(ad);
            long ni = ( nodeP == NULL ) ? i : nodeP [ i ];
            double val = sm.a [ ad ];

            if ( ni < nj ) {
                columndata [ col_start - ni ] = val;
            } else if ( ni == nj ) {
                D [ ni ] = val;
            }
        }
    }
}


SkyLineMtxLDL :: ~SkyLineMtxLDL()
{}

void SkyLineMtxLDL :: LoadMatrixNumbers(SparseMatrixF &sm)
{}
void SkyLineMtxLDL :: SolveLV(const LargeVector &b, LargeVector &x)
{}
void SkyLineMtxLDL :: Factorize()
{}
void SkyLineMtxLDL :: MultiplyByVector(const LargeVectorAttach &x, LargeVectorAttach &y)
{}

void SkyLineMtxLDL :: SchurComplementFactorization(int fixed_blocks)
{}

void SkyLineMtxLDL :: SolveA11(double *x, long fixed_blocks)
{}

void SkyLineMtxLDL :: Sub_A21_A11inv(double *x, long fixed_blocks)
{}

void SkyLineMtxLDL :: Sub_A11inv_A12(double *x, long fixed_blocks)
{}
void SkyLineMtxLDL :: WriteCondensedMatrixA22(double *a, Ordering *mcn, IntArrayList *lncn)
{}

DSS_NAMESPASE_END
