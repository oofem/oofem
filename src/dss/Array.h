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

#ifndef _ARRAY_H__
#define _ARRAY_H__

#include "DSSAfx.h"

DSS_NAMESPASE_BEGIN

/**
 * @author Richard Vondracek
 */

class Array
{
public:
    /*
     * template <class ItemType>
     * static void Copy(ItemType* source,ItemType* dest,long count)
     * {
     * memcpy(dest,source,count*sizeof(ItemType));
     * }
     */

    static void Copy(long *source, long *dest, long count)
    {
        if ( count > 0 ) {
            memcpy( dest, source, count * sizeof( long ) );
        }
    }

    static void Copy(int *source, int *dest, int count)
    {
        if ( count > 0 ) {
            memcpy( dest, source, count * sizeof( int ) );
        }
    }

    static void Copy(unsigned long *source, unsigned long *dest, unsigned long count)
    {
        if ( count > 0 ) {
            memcpy( dest, source, count * sizeof( unsigned long ) );
        }
    }

    static void Copy(double *source, double *dest, long count)
    {
        if ( count > 0 ) {
            memcpy( dest, source, count * sizeof( double ) );
        }
    }

    /*
     * template <class ItemType>
     * static void Copy(ItemType* source,long sindex,ItemType* dest,long dindex,long count)
     * {
     * memcpy(dest+dindex,source+sindex,count*sizeof(ItemType));
     * }
     */

    static void Copy(long *source, long sindex, long *dest, long dindex, long count)
    {
        if ( count > 0 ) {
            memcpy( dest + dindex, source + sindex, count * sizeof( long ) );
        }
    }

    static void Copy(double *source, long sindex, double *dest, long dindex, long count)
    {
        if ( count > 0 ) {
            memcpy( dest + dindex, source + sindex, count * sizeof( double ) );
        }
    }

    /*
     * template <class ItemType>
     * static void Clear(ItemType* dest,long dindex,long count)
     * {
     * memset(dest+dindex,0,count*sizeof(ItemType));
     * }
     */

    static void Clear(long *dest, long dindex, long count)
    {
        if ( count > 0 ) {
            memset( dest + dindex, 0, count * sizeof( long ) );
        }
    }

    static void Clear(double *dest, long dindex, long count)
    {
        if ( count > 0 ) {
            memset( dest + dindex, 0, count * sizeof( double ) );
        }
    }

    /*
     * template <class ItemType>
     * static void Reverse(ItemType* dest,long count)
     * {
     * ItemType tmp;
     * for (long i=count/2; i>=0; i--)
     *    {
     *      tmp = dest[i];
     *      dest[i] = dest[count-i];
     *      dest[count-i] = tmp;
     *    }
     * }
     */

    static void Reverse(long *dest, long count)
    {
        long tmp;
        for ( long i = count / 2 - 1; i >= 0; i-- ) {
            tmp = dest [ i ];
            dest [ i ] = dest [ count - i - 1 ];
            dest [ count - i - 1 ] = tmp;
        }
    }
};

DSS_NAMESPASE_END

#endif
