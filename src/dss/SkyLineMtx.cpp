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

#include "SkyLineMtx.h"

DSS_NAMESPASE_BEGIN

// Allocates new space according to bskl and reads old matrix with respect
// to permutation blockP
SkyLineMtx :: SkyLineMtx(SparseMatrixF &sm, Ordering *order, MathTracer *eMT)
{
    this->n = sm.neq;
    this->nonzeros = 0;
    this->eMT = eMT;
    this->order = order;
    this->column_starts = NULL;
    this->columndata = NULL;
    this->D = NULL;
}

SkyLineMtx :: ~SkyLineMtx()
{
    if ( column_starts ) {
        delete[] column_starts;
    }

    if ( columndata ) {
        delete[] columndata;
    }

    if ( D ) {
        delete[] D;
    }
}

void SkyLineMtx :: GrowSkyline(int i, int j, long *column_ns)
{
    if ( i > j ) {
        int l = i - j;
        if ( column_ns [ i ] < l ) {
            column_ns [ i ] = l;
        }
    } else {
        int l = j - i;
        if ( column_ns [ j ] < l ) {
            column_ns [ j ] = l;
        }
    }
}


void SkyLineMtx :: AllocateMemory(IConectMatrix *spm, int neq)
{
    long *perm = order->perm->Items;
    column_starts = new int [ n + 1 ];
    long *column_ns = new long [ n ];
    D = new double [ n ];

    for ( int j = 0; j < n; j++ ) {
        IntArrayList *ColumnIndexes = spm->GetIndexesAboveDiagonalInColumn(j);
        int nj = perm != NULL ? perm [ j ] : j;
        if ( nj >= neq ) {
            continue;
        }

        int cnt = ColumnIndexes->Count;
        for ( int idx = 0; idx < cnt; idx++ ) {
            int ni = perm != NULL ? perm [ ColumnIndexes->Items [ idx ] ] : ColumnIndexes->Items [ idx ];
            if ( ni >= neq ) {
                continue;
            }

            GrowSkyline(ni, nj, column_ns);
        }
    }

    long column_field_lenght = 0;
    for ( int j = 0; j < n; j++ ) {
        column_starts [ j ] = column_field_lenght;
        column_field_lenght += column_ns [ j ];
    }

    columns_data_length = column_starts [ n ] = column_field_lenght;

    if ( column_field_lenght == 0 ) {   // This array seems to need nonzero size, needed for the fixed statement
        column_field_lenght = 1;
    }

    columndata = new double [ column_field_lenght ];

    delete[] column_ns;
}

void SkyLineMtx :: WriteStatistics(long no_init_blocks, long no_nonzeros)
{
    char str [ 512 ];

    sprintf( str, " number of nonzeros  : %ld", Nonzeros() );
    eMT->Writeln(str);
}

DSS_NAMESPASE_END
