/* $Header: /home/cvs/bp/oofem/oofemlib/src/intarray.h,v 1.9.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

//   ***********************
//   *** CLASS INT ARRAY ***
//   ***********************


#ifndef intarray_h
#define intarray_h

#include "freestor.h"
#include "debug.h"
#include "cltypes.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#endif

class DataStream;
#ifdef __PARALLEL_MODE
class CommunicationBuffer;
#endif
/**
 * Class implementing an array of integers. Array can grow or shrink to desired dimension.
 * The lower value index of array is 1, upper depends on array size.
 */
class IntArray
{
    /*
     * This class implements an array of integers.
     * DESCRIPTION :
     * An IntArray stores its coefficients in an array 'value' of size 'size'.
     * TASKS :
     * - storing and reterning coefficients (method 'at') ;
     * - appending another IntArray to itself ;
     * - introduced allocatedSize variable to allow dynamic rescaling of array
     *   size possibly without memory realocation. At startup array occupies space
     *   given by allocatedSpace = size. Then there can be
     *   1) further request for resizeing array to smaller dimension
     *      then we only change size wariable, but allocatedSize
     *  variable remain untouched - expecting possible array grow and then re-using
     *  previously allocated space.
     *   2) if further request for growing then is necessary memory realocation.
     *   This process is controlled in resize member function.
     * REMARK :
     * see Remark 2 in file "floatarry.hxx".
     */


private:
    /// size of array
    int size;
    /// allocated size for array
    int allocatedSize;
    /// stored values
    int *values;



public:
    /// Constructor for zero sized array
    IntArray(int = 0);                                   // constructor
    /** Copy constructor. Creates the array from another array.
     */
    IntArray(const IntArray &);                      // copy constructor
    /// Destructor.
    ~IntArray()  { if ( values ) { freeInt(values); } } // destructor

    /// Assingnment operator
    IntArray &operator=(const IntArray &);              // assignment: cleanup and copy


#     ifdef DEBUG
    /** Coefficient access function. Returns l-value of coeffiicient at given
     * position of the receiver.
     * @param i position of coefficient in array
     */
    int &at(int i);
    /** Coefficient access function. Returns value of coeffiicient at given
     * position of the receiver.
     * @param i position of coefficient in array
     */
    int     at(int i) const;
#     else
    int &at(int i)                  { return values [ i - 1 ]; }
    int     at(int i) const { return values [ i - 1 ]; }
#     endif
    /**
     * Coefficient access function. Returns value of coeffiicient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i position of coefficient in array
     */
    int &operator()(int i)
    {
#       ifdef DEBUG
        assert(i < size);
#       endif
        return values [ i ];
    }
    /**
     * Coefficient access function. Returns value of coeffiicient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i position of coefficient in array
     */
    const int &operator()(int i) const
    {
#       ifdef DEBUG
        assert(i < size);
#       endif
        return values [ i ];
    }

    /** Checks size of receiver towards requested bounds.
     * Current implementation will call exit(1), if dimension
     * mismatch found.
     * @param i required size of receiver
     */
    void       checkBounds(int i) const;
    /** Checks size of receiver towards requested bounds.
     * If dimension mismatch, size is adjusted accordingly.
     * Warning: after this operation array values are in undefined state, programmer should
     * zero receiver
     * @param allocChunk if reallocation needed, an aditional space for allocChunk values will be allocated
     */
    void       resize(int n, int allocChunk = 0);
    /**
     * Appends array b at the end of receiver.
     * @param b array to be appended at the end of receiver
     * @param allocChunk if reallocation needed, an aditional space for allocChunk values will be allocated
     * will be allocated to prevent excessive realocation
     */
    void       followedBy(const IntArray &b, int allocChunk = 0);
    /**
     * Appends given Number at the end of receiver.
     * @param b value to be appended
     * @param allocChunk if reallocation needed, an aditional space for allocChunk values will be allocated
     * will be allocated to prevent excessive realocation
     */
    void       followedBy(const int b, int allocChunk = 0);
    /// Returns the size of receiver.
    int        giveSize() const { return size; }
    /// Checks if receiver is empty (i.e., zero sized).
    int        isEmpty()  const { return size == 0; }
    /// Returns true if receiver contains only zeroes
    int         containsOnlyZeroes() const;

    /** finds the first occurence of given value, assuming that the receiver is sorted.
     *  Returns its index in (1-based) indexing, if not present, 0 is returned.
     */
    int findSorted(int value) const;
    /**
     * returns true if receiver contains given value
     */
    bool containsSorted(int value) const { return ( findSorted(value) > 0 ); }

    /** Inserts given value into a receiver, which is assumed to be sorted.
     *  The size of receiver is changed accordingly.
     *  @param value value to insert
     *  @param allocChunk if reallocation needed, an aditional space for allocChunk values will be allocated
     *  @return index of inserted (or existing) value
     */
    int insertSorted(int value, int allocChunk = 0);
    /** Inserts given value into a receiver, which is assumed to be sorted.
     *  The value is inserted only if it does not exist.
     *  The size of receiver is changed accordingly.
     *  @param value value to insert
     *  @param allocChunk if reallocation needed, an aditional space for allocChunk values will be allocated
     *  @return index of inserted (or existing) value
     */
    int insertSortedOnce(int value, int allocChunk = 0);




    /**
     * Finds index of first occurence of given value in array. If such value is not presented,
     * returns zero value.
     * @param value scanned value
     * @return index of first value in array, otherwise zero
     */
    int        findFirstIndexOf(int value)  const;
    /**
     * returns true if receiver contains given value
     */
    bool contains(int value) const { return ( findFirstIndexOf(value) > 0 ); }
    /**
     * Insert once (does not make any assumption about receiver state or ordering, quite
     * inefficient). More efficient version InsertSortedOnce exist. Insert _p  value into
     * receiver if it does not exist.
     * @return index of sorted value;
     */
    int insertOnce(int _p);
    /**
     * Erase the element at given position (1-based index)
     * Receiver will shrink accordingly, the values at positions (_pos+1,...,size)
     * will be moved to positions (_pos,...,size-1)
     */
    void       erase(int _pos);

    /// add given subvector to receiver values starting at position si
    void addSubVector(const IntArray &src, int si);
    /// copy given subvector to receiver values starting at position si
    void copySubVector(const IntArray &src, int si);

    /// Prints receiver on stdin.
    void       printYourself() const;
    /** Stores array  to output stream.
     * @see FEMComponent class */
    contextIOResultType          storeYourself(DataStream *stream, ContextMode mode) const;
    /** Restores array from image on stream.
     * @see FEMComponent class */
    contextIOResultType          restoreYourself(DataStream *stream, ContextMode mode);
    /// Sets all component to zero.
    void       zero();
    int *givePointer()  const { return values; }
    /// adds given value to all values of receiver
    void       add(int val);

#ifdef __PARALLEL_MODE
    /**@name Methods for  packing/unpacking to/from communication buffer */
    //@{
    /**
     * Packs receiver into communication buffer.
     * @param buff buffer to pack itself into.
     * @return nonzero if succesfull
     */
    int packToCommBuffer(CommunicationBuffer &buff) const;
    /**
     * Unpacks receiver from communication buffer.
     * @param buff buffer from which unpack itself.
     * @return nonzero if succesfull
     */
    int unpackFromCommBuffer(CommunicationBuffer &buff);
    /**
     * Returns how much space is needed to pack receivers message.
     * @param buff buffer used for packing
     */
    int givePackSize(CommunicationBuffer &buff);
    //@}
#endif
};


template< class operation > int
quickSortPartition(IntArray &arry, int l, int r, operation op) {
    int i = l - 1, j = r;
    int v = arry.at(r);
    int swap;

    for ( ; ; ) {
        while ( ( op(arry.at(++i), v) ) < 0 ) {
            ;
        }

        while ( ( op( v, arry.at(--j) ) ) < 0 ) {
            if ( j == l ) {
                break;
            }
        }

        if ( i >= j ) {
            break;
        }

        swap = arry.at(i);
        arry.at(i) = arry.at(j);
        arry.at(j) = swap;
    }

    swap = arry.at(i);
    arry.at(i) = arry.at(r);
    arry.at(r) = swap;
    return i;
}



template< class operation > void quickSort(IntArray &arry, int l, int r, operation op) {
    if ( r <= l ) {
        return;
    }

    int i = quickSortPartition(arry, l, r, op);
    quickSort(arry, l, i - 1, op);
    quickSort(arry, i + 1, r, op);
}


/**
 * Sorts the receiver using quiksort algorithm.
 * @param op is Function object, required to have member function int class::operator() (int, int),
 * must return a negative value if first argument is less than the second,
 * zero if the arguments are equal, and a positive number otherwise.
 */
template< class operation > void sort(IntArray &arry, operation op) { quickSort(arry, 1, arry.giveSize(), op); }

#endif // intarray_h


