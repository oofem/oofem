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

#ifndef _COLHASH_H__
#define _COLHASH_H__

#include "Array.h"
#include "IntArrayList.h"

DSS_NAMESPASE_BEGIN

/**
 * @author Richard Vondracek
 */

class ColHash
{
public:
    IntArrayList **buckets;
    IntArrayList occupied_buckets;
    int n;

    ColHash()
    {
        buckets = NULL;
    }

    ColHash(long n)
    {
        buckets = new IntArrayList * [ n ];
        memset( buckets, 0, n * sizeof( IntArrayList * ) );
        this->n = n;
        //for (long i=0; i<n; i++)
        //buckets[i] = new IntArrayList();
    }

    void Init(long n)
    {
        if ( buckets ) { delete [] buckets; }

        buckets = new IntArrayList * [ n ];
        memset( buckets, 0, n * sizeof( IntArrayList * ) );
        this->n = n;
    }

    ~ColHash()
    {
        for ( int i = 0; i < n; i++ ) {
            if ( buckets [ i ] ) { delete buckets [ i ];
                                   buckets [ i ] = NULL; } }

        if ( buckets ) { delete [] buckets; }

        buckets = NULL;
    }

    void AddValue(long bucket, long obj)
    {
        if ( buckets [ bucket ] == NULL ) {
            buckets [ bucket ] = new IntArrayList();
        }

        //buckets[bucket].Add(obj);

        long i = buckets [ bucket ]->Add(obj);
        if ( i == 0 ) {
            occupied_buckets.Add(bucket);
        }
    }

    void Clear()
    {
        long *pI = occupied_buckets.Items;
        for ( long *p = pI + occupied_buckets.Count - 1; p >= pI; p-- ) {
            buckets [ * p ]->Clear();
        }

        occupied_buckets.Clear();
    }
};

DSS_NAMESPASE_END

#endif //_COLHASH_H__
