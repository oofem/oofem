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

#ifndef _GRIDCOLUMN_H__
#define _GRIDCOLUMN_H__

#include "IntArrayList.h"
#include "Array.h"

DSS_NAMESPASE_BEGIN

/**
 * @author: Richard Vondracek
 */

class SparseGridColumn
{
    // used during sealed state
public:
    long column_start_idx;      // pointer into Columns_data array
    IntArrayList *IndexesUfa;
    long Entries;

    SparseGridColumn(IntArrayList *indexes)
    {
        Entries = 0;
        column_start_idx = -1;
        IndexesUfa = indexes;
        Entries = indexes->Count;
    }

    ~SparseGridColumn()
    {
        if ( IndexesUfa ) { delete IndexesUfa;
                            IndexesUfa = NULL; }
    }

    // Sealed procedures
    //---------------------------

    // Column must not be sealed
    long FindExistingBlockIndex(long bi)        // returns the pointer in the data
    {
        long myIndex = IndexesUfa->BinarySearch(bi);
        //System.Diagnostics.Debug.Assert(myIndex>=0,"Index not found !");
        //if (myIndex<0)
        //return 0;
        return myIndex;
    }

    void SetValue(long bn, long bi, long si, long sj, double val, double *Columns_data, long &aux_idx)
    {
        if ( ( aux_idx >= IndexesUfa->Count ) || ( IndexesUfa->Items [ aux_idx ] != bi ) ) {
            aux_idx = FindExistingBlockIndex(bi);
        }

        if ( aux_idx >= 0 ) {
            Columns_data [ column_start_idx + bn * ( bn * aux_idx + sj ) + si ] = val;
        }
    }

    double GetValue(long bn, long bi, long si, long sj, double *Columns_data, long &aux_idx)
    {
        if ( ( aux_idx >= IndexesUfa->Count ) || aux_idx < 0 || ( IndexesUfa->Items [ aux_idx ] != bi ) ) {
            aux_idx = FindExistingBlockIndex(bi);
        }

        if ( aux_idx < 0 ) {
            return 0;
        }

        return Columns_data [ column_start_idx + bn * ( bn * aux_idx + sj ) + si ];
    }

    void AddValue(long bn, long bi, long si, long sj, double val, double *Columns_data)
    {
        long newBlock = FindExistingBlockIndex(bi);
        Columns_data [ column_start_idx + bn * bn * newBlock + si + sj * bn ] += val;
    }
};

DSS_NAMESPASE_END

#endif //GRIDCOLUMN_H__
