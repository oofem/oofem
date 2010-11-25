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

using System;
using System.Collections;
using System.Diagnostics;
using System.Runtime.InteropServices;

using Golem;

namespace MathKer
{
/// <summary>
/// This is a fixed sparse matrix format
/// There is full set of rows but in each row are stored only nonzeros
/// </summary>
public class SparseConectivityMtxFixed : IConectMatrix
{
public:
    IntArrayList **ColumnsIndexes;
    int n;
    int nonzeros;

public:
    int N()                 { return n; }
    int Nonzeros()  { return nonzeros; }


    public void Init(int N)
    {
        n = N;
        ColumnsIndexes = new IntArrayList * [ n ];
    }

    ~SparseConectivityMtxFixed()
    {
        for ( int i = 0; i < n; i++ ) {
            if ( ColumnsIndexes [ i ] ) { delete ColumnsIndexes [ i ];
                                          ColumnsIndexes [ i ] = NULL; } }

        delete [] ColumnsIndexes;
    }

    IntArrayList *GetIndexesAboveDiagonalInColumn(int j)
    {
        return ColumnsIndexes [ j ];
    }
}
}
