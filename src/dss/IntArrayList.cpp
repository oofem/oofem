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

#include "IntArrayList.h"

DSS_NAMESPASE_BEGIN

IntArrayList :: IntArrayList()
{
    Items = new long [ Capacity = 16 ];
    Count = 0;
}


IntArrayList :: IntArrayList(long capacity)
{
    Items = new long [ Capacity = capacity ];
    Count = 0;
}

IntArrayList :: IntArrayList(long n, long *data)
{
    Items = new long [ n ];
    Array :: Copy(data, Items, n);
    Capacity = Count = n;
}

//IntArrayList(ArrayList al)
//{
//Items = (long[])al.ToArray(typeof(long));
//Capacity = Count = al.Count;
//}

IntArrayList :: IntArrayList(const IntArrayList &al)
{
    Items = new long [ al.Count ];
    Array :: Copy(al.Items, 0, Items, 0, al.Count);
    Capacity = Count = al.Count;
}

IntArrayList :: ~IntArrayList()
{
    if ( Items ) {
        delete [] Items;
    }

    Items = NULL;
}

long IntArrayList :: Add(long a)
{
    long cn = Count;
    if ( Capacity <= cn ) {
        long *old_items = Items;
        Capacity *= 2;
        if ( Capacity < 4 ) {
            Capacity = 4;
        }

        Items = new long [ Capacity ];
        Array :: Copy(old_items, Items, Count);
        if ( old_items ) {
            delete [] old_items;
        }
    }

    Items [ Count++ ] = a;
    return cn;
}

void IntArrayList :: AddRange(long count, long val)
{
    long cn = Count + count;
    if ( Capacity < cn ) {
        long *old_items = Items;
        Capacity = Count + count;
        Items = new long [ Capacity ];
        Array :: Copy(old_items, Items, Count);
        if ( old_items ) {
            delete [] old_items;
        }
    }

    for ( long i = Count; i < Count + count; i++ ) {
        Items [ i ] = val;
    }

    Count += count;
}

void IntArrayList :: Fill(const IntArrayList &al)
{
    if ( Items ) {
        delete [] Items;
    }

    Items = new long [ al.Count ];
    Array :: Copy(al.Items, 0, Items, 0, al.Count);
    Capacity = Count = al.Count;
}

/// <summary>If the last item of the array is not Item, then add it.</summary>
/// <param name="item">Item</param>
/// <returns>index of item</returns>
long IntArrayList :: AddIfNotLast(long item)
{
    if ( Count == 0 || Items [ Count - 1 ] != item ) {
        return Add(item);
    }

    return -1;
}

void IntArrayList :: Clear()
{
    Count = 0;
}

void IntArrayList :: TrimToSize()
{
    if ( Capacity > Count ) {
        long *old_items = Items;
        Capacity = Count;
        Items = new long [ Capacity ];
        Array :: Copy(old_items, Items, Count);
        if ( old_items ) {
            delete [] old_items;
        }
    }
}

void IntArrayList :: Alloc()
{
    Count = Capacity;
    Initialize();
}

void IntArrayList :: InitIdentity()
{
    Count = Capacity;
    for ( int i = 0; i < Count; i++ ) {
        Items [ i ] = i;
    }
}

void IntArrayList :: Alloc(long size)
{
    if ( Capacity < size ) {
        Count = Capacity = size;
        if ( Items ) {
            delete [] Items;
        }

        Items = new long [ Capacity ];
        Initialize();
    } else {
        Count = size;
        Initialize();
    }
}

void IntArrayList :: Initialize()
{
    memset( Items, 0, Count * sizeof( long ) );
}

void IntArrayList :: Initialize(long val)
{
    for ( long *p = Items + Count - 1; p >= Items; p-- ) {
        * p = val;
    }
}

void IntArrayList :: SortInsert(long item)
{
    if ( Capacity <= Count ) {
        long *old_items = Items;
        Capacity *= 2;
        if ( Capacity < 4 ) {
            Capacity = 4;
        }

        Items = new long [ Capacity ];
        Array :: Copy(old_items, Items, Count);
        if ( old_items ) {
            delete [] old_items;
        }
    }

    long *pI = this->Items;
    {
        long *ptr = pI + Count;
        for ( ; ptr > pI && ptr [ -1 ] > item; ptr-- ) {
            * ptr = ptr [ -1 ];
        }

        * ptr = item;
    }
    Count++;
}

long IntArrayList :: SortInsertIfNot(long item)
{
    long myIndex = BinarySearch(item);
    if ( myIndex < 0 ) {
        Insert(~myIndex, item);
        return ~myIndex;
    }

    return -1;
}

void IntArrayList :: Insert(long pos, long item)
{
    if ( pos < 0 || pos > Count ) {
        return;
    }

    if ( Capacity <= Count ) {
        long *old_items = Items;
        Capacity *= 2;
        if ( Capacity < 4 ) {
            Capacity = 4;
        }

        Items = new long [ Capacity ];
        Array :: Copy(old_items, 0, Items, 0, pos);
        Array :: Copy(old_items, pos, Items, pos + 1, Count - pos);
        if ( old_items ) {
            delete [] old_items;
        }
    } else {
        //Array::Copy(Items,pos,Items,pos+1,Count-pos);
        for ( int i = Count; i > pos; i-- ) {
            Items [ i ] = Items [ i - 1 ];
        }
    }

    Count++;
    Items [ pos ] = item;
}

void IntArrayList :: RemoveAt(long pos)
{
    Count--;
    Array :: Copy(Items, pos + 1, Items, pos, Count - pos);
}

long IntArrayList :: BinarySearch(long item)
{
    return BinarySearch(Count, item);
}

//long BinarySearch(long startindex,long count, long item)
long IntArrayList :: BinarySearch(long count, long item)
{
    if ( count == 0 ) {
        return ~0;
    }

    long h = 0;
    long e = count - 1;

    if ( item < Items [ h ] ) {
        return ~0;
    }

    if ( item == Items [ h ] ) {
        return 0;
    }

    if ( item > Items [ e ] ) {
        return ~count;
    }

    if ( item == Items [ e ] ) {
        return e;
    }

    while ( ( h + 1 ) < e ) {
        long m = ( h + e ) / 2;

        if ( item == Items [ m ] ) {
            return m;
        } else if ( item < Items [ m ] ) {
            e = m;
        } else {
            h = m;
        }
    }

    return ~e;
    //startindex;
}

long *IntArrayList :: ToArray()
{
    long *arr = new long [ Count ];
    memcpy(arr, Items, Count);
    //Array::Copy(Items,arr,Count);
    return arr;
}

/// <summary>Sums all contained integers.</summary>
/// <returns>The sum.</returns>
long IntArrayList :: SumElements()
{
    long sum = 0;
    long *pI = this->Items;
    if ( pI ) {
        for ( long *p = pI + Count - 1; p >= pI; p-- ) {
            sum += * p;
        }
    }

    return sum;
}

/// <summary>Sets items of a ptr* array using all indexes stored in this array.</summary>
/// <param name="ptr">update items of this array</param>
/// <param name="val">set value</param>
void IntArrayList :: SetIndexesTo(long *ptr, long val)
{
    long *pI = this->Items;
    if ( pI ) {
        for ( long *p = pI + Count - 1; p >= pI; p-- ) {
            ptr [ * p ] = val;
        }
    }
}

/// <summary>Sets ~indexes of a ptr* array using all values stored in this array.</summary>
/// <param name="ptr">update items of this array</param>
void IntArrayList :: SetIndexesToMap(long *ptr)
{
    if ( Items ) {
        long *pI = Items;
        long idx = Count - 1;
        for ( long *p = pI + idx; p >= pI; p--, idx-- ) {
            ptr [ * p ] = ~idx;
        }
    }
}


bool IntArrayList :: TestSetIndexesTo(long *ptr, long valtest, long valset, bool &ind)
{
    long *pI = this->Items;
    if ( pI ) {
        for ( long *p = pI + Count - 1; p >= pI; p-- ) {
            if ( ptr [ * p ] != valtest ) {
                ind = false;
            }

            ptr [ * p ] = valset;
        }
    }

    return ind;
}

/// <summary>If ptr[element] is nonzero, then element will be removed from the array. </summary>
/// <param name="ptr">Mask of which elements must be removed.</param>
void IntArrayList :: RemoveByBattern(long *ptr)
{
    long *pI = this->Items;
    if ( pI ) {
        long *src = pI;
        long *des = pI;
        for ( long k = Count; k > 0; k-- ) {
            long e = * src++;
            if ( ptr [ e ] == 0 ) {
                * des++ = e;
            } else {
                Count--;
            }
        }
    }
}

/// <summary>If ptr[element] is nonzero, then element will be removed from the array.
/// All remaining elements update ptr[element] to 1 </summary>
/// <param name="ptr">Mask of which elements must be removed.</param>
void IntArrayList :: RemoveMarkByBattern(long *ptr)
{
    long *src = Items;
    long *des = Items;
    for ( long k = Count; k > 0; k-- ) {
        long e = * src++;
        if ( ptr [ e ] == 0 ) {
            * des++ = e;
            ptr [ e ] = 1;
        } else {
            Count--;
        }
    }
}

/// <summary>Computes how many elements are not set in the mask.</summary>
/// <param name="ptr">mask</param>
long IntArrayList :: CountWithoutMask(long *ptr)
{
    long cnt = 0;
    if ( Items ) {
        for ( long *p = Items + Count - 1; p >= Items; p-- ) {
            if ( ptr [ * p ] == 0 ) {
                cnt++;
            }
        }
    }

    return cnt;
}

DSS_NAMESPASE_END
