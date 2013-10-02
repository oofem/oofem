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

#ifndef _INTARRAYLIST_H__
#define _INTARRAYLIST_H__

#include "Array.h"


DSS_NAMESPASE_BEGIN

/**
 * @author Richard Vondracek
 */

class IntArrayList
{
public:
    long *Items;
    long Count;
    long Capacity;

    IntArrayList();

    IntArrayList(long capacity);

    IntArrayList(long n, long *data);

    //IntArrayList(ArrayList al)
    //{
    //  Items = (long[])al.ToArray(typeof(long));
    //  Capacity = Count = al.Count;
    //}

    IntArrayList(const IntArrayList &al);

    virtual ~IntArrayList();

    inline long operator[](const long i) const
    {
        return Items [ i ];
    }

    inline long &operator[](const long i)
    {
        return Items [ i ];
    }

    long Add(long a);

    void AddRange(long count, long val);

    void Fill(const IntArrayList &al);

    /// <summary>If the last item of the array is not Item, then add it.</summary>
    /// <param name="item">Item</param>
    /// <returns>index of item</returns>
    long AddIfNotLast(long item);

    void Clear();
    void TrimToSize();
    void Alloc();
    void Alloc(long size);
    void Initialize();
    void Initialize(long val);
    void InitIdentity();

    void SortInsert(long item);
    long SortInsertIfNot(long item);
    void Insert(long pos, long item);
    void RemoveAt(long pos);
    long BinarySearch(long item);
    //  long BinarySearch(long startindex,long count, long item)
    long BinarySearch(long count, long item);

    long *ToArray();

    /// <summary>Sums all contained integers.</summary>
    /// <returns>The sum.</returns>
    long SumElements();

    /// <summary>Sets items of a ptr* array using all indexes stored in this array.</summary>
    /// <param name="ptr">update items of this array</param>
    /// <param name="val">set value</param>
    void SetIndexesTo(long *ptr, long val);

    /// <summary>Sets ~indexes of a ptr* array using all values stored in this array.</summary>
    /// <param name="ptr">update items of this array</param>
    void SetIndexesToMap(long *ptr);

    bool TestSetIndexesTo(long *ptr, long valtest, long valset, bool &ind);

    /// <summary>If ptr[element] is nonzero, then element will be removed from the array. </summary>
    /// <param name="ptr">Mask of which elements must be removed.</param>
    void RemoveByBattern(long *ptr);

    /// <summary>If ptr[element] is nonzero, then element will be removed from the array.
    /// All remaining elements update ptr[element] to 1 </summary>
    /// <param name="ptr">Mask of which elements must be removed.</param>
    void RemoveMarkByBattern(long *ptr);

    /// <summary>Computes how many elements are not set in the mask.</summary>
    /// <param name="ptr">mask</param>
    long CountWithoutMask(long *ptr);
};

DSS_NAMESPASE_END

#endif //INTARRAYLIST_H__
