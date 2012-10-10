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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef intarray_h
#define intarray_h

#include "freestor.h"
#include "contextioresulttype.h"
#include "contextmode.h"

namespace oofem {
class DataStream;
#ifdef __PARALLEL_MODE
class CommunicationBuffer;
#endif
/**
 * Class implementing an array of integers. Array can grow or shrink to desired dimension.
 * The lower value index of array is 1, upper depends on array size.
 *
 * Tasks:
 * - Storing and returning coefficients (method 'at')
 * - Appending another IntArray to itself
 *
 * The allocatedSize variable to allow dynamic rescaling of array
 * size possibly without memory reallocation. At startup an array occupies space
 * given by allocatedSpace = size. Then there can be
 * further request for resizing array to smaller dimension
 * then we only change size variable, but allocatedSize
 * variable remain untouched - expecting possible array grow and then re-using
 * previously allocated space.
 * if further request for growing then is necessary memory reallocation.
 * This process is controlled in resize member function.
 */
class IntArray
{
private:
    /// Size of array.
    int size;
    /// Allocated size for array.
    int allocatedSize;
    /// Stored values.
    int *values;

public:
    /// Constructor for zero sized array
    IntArray(int = 0);
    /// Copy constructor. Creates the array from another array.
    IntArray(const IntArray &);
    /// Destructor.
    ~IntArray() { if ( values ) { freeInt(values); } }

    /// Assignment operator
    IntArray & operator=(const IntArray &);

    /**
     * Coefficient access function. Returns l-value of coefficient at given
     * position of the receiver.
     * @param i Position of coefficient in array.
     * @return Value at position.
     */
#ifdef DEBUG
    int &at(int i);
#else
    int &at(int i) { return values [ i - 1 ]; }
#endif
    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver.
     * @param i position of coefficient in array.
     * @return Value at position.
     */
#ifdef DEBUG
    int at(int i) const;
#else
    int at(int i) const { return values [ i - 1 ]; }
#endif
    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i Position of coefficient in array.
     * @return Value at position.
     */
    int &operator()(int i)
    {
#ifdef DEBUG
        checkBounds(i);
#endif
        return values [ i ];
    }
    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i position of coefficient in array.
     * @return Value at position.
     */
    const int &operator()(int i) const
    {
#ifdef DEBUG
        checkBounds(i);
#endif
        return values [ i ];
    }

#ifdef DEBUG
    /**
     * Checks size of receiver towards requested bounds.
     * Current implementation will call exit(1) if dimension
     * mismatch found.
     * @param i Required size of receiver
     */
    void checkBounds(int i) const;
#endif
    /**
     * Checks size of receiver towards requested bounds.
     * If dimension mismatch, size is adjusted accordingly.
     * @note{After this operation array values are in undefined state, programmer should zero receiver.}
     * @param n New size of array.
     * @param allocChunk If reallocation needed, an additional space for allocChunk values will be allocated
     * to prevent excessive reallocation.
     */
    void resize(int n, int allocChunk = 0);
    /**
     * Preallocates receiver to given futureSize if larger then allocatedSize.
     * @note{After this operation array values are in undefined state, programmer should zero receiver.}
     * @param futureSize Size to be allocated.
     */
    void preallocate(int futureSize);
    /**
     * Appends array b at the end of receiver.
     * @param b Array to be appended at the end of receiver
     * @param allocChunk If reallocation needed, an additional space for allocChunk values will be allocated
     * to prevent excessive reallocation.
     */
    void followedBy(const IntArray &b, int allocChunk = 0);
    /**
     * Appends given Number at the end of receiver.
     * @param b value to be appended.
     * @param allocChunk If reallocation needed, an additional space for allocChunk values will be allocated
     * to prevent excessive reallocation.
     */
    void followedBy(const int b, int allocChunk = 0);
    /// @return Size of receiver.
    int giveSize() const { return size; }
    /**
     * Checks if receiver is empty (i.e., zero sized).
     * @return True is size is zero.
     */
    bool isEmpty() const { return size == 0; }
    /**
     * Checks if receiver is all zero.
     * @return True is receiver contains only zeroes.
     */
    bool containsOnlyZeroes() const;
    /**
     * Finds the first occurrence of given value, assuming that the receiver is sorted.
     * @return Index (1-based), if not present, 0 is returned.
     */
    int findSorted(int value) const;
    /**
     * Finds the minimum component in the array.
     * @return Minimum of array or prints error is array is empty.
     */
    int minimum() const;
    /**
     * Finds the maximum component in the array.
     * @return Maximum of array or prints error is array is empty.
     */
    int maximum() const;
    /**
     * Checks if sorted receiver contains a given value.
     * @return True if receiver contains given value.
     */
    bool containsSorted(int value) const { return ( findSorted(value) > 0 ); }

    /**
     * Inserts given value into a receiver, which is assumed to be sorted.
     * The size of receiver is changed accordingly.
     * @param value value to insert.
     * @param allocChunk If reallocation needed, an additional space for allocChunk values will be allocated.
     * @return Index of inserted (or existing) value.
     */
    int insertSorted(int value, int allocChunk = 0);
    /**
     * Inserts given value into a receiver, which is assumed to be sorted.
     * The value is inserted only if it does not exist.
     * The size of receiver is changed accordingly.
     * @param value Value to insert.
     * @param allocChunk If reallocation needed, an additional space for allocChunk values will be allocated.
     * @return Index of inserted (or existing) value.
     */
    int insertSortedOnce(int value, int allocChunk = 0);
    /**
     * Erase the element of given value.
     * If the value is found receiver will shrink accordingly,
     * while preserving a sorted state.
     * @param value Value to erase.
     */
    void eraseSorted(int value);
    /**
     * Extracts common values in receiver and iarray.
     * Assumes that receiver as well as iarray are sorted.
     * The size of array common is changed accordingly.
     * @param iarray Array to search for values common with receiver.
     * @param common Array of common values.
     * @param allocChunk If reallocation needed, an additional space for allocChunk values will be allocated.
     * @return Number of common values.
     */
    int findCommonValuesSorted(const IntArray &iarray, IntArray &common, int allocChunk = 0) const;

    /**
     * Finds index of first occurrence of given value in array.
     * @param value Scanned value.
     * @return Index of first value in array, otherwise zero.
     */
    int findFirstIndexOf(int value) const;
    /**
     * @return True if receiver contains given value.
     */
    bool contains(int value) const { return ( findFirstIndexOf(value) > 0 ); }
    /**
     * Insert once (does not make any assumption about receiver state or ordering, quite
     * inefficient). More efficient version insertSortedOnce exist.
     * @param p Value to insert.
     * @return Index of sorted value.
     */
    int insertOnce(int p);
    /**
     * Erase the element at given position (1-based index)
     * Receiver will shrink accordingly, the values at positions (pos+1,...,size)
     * will be moved to positions (pos,...,size-1)
     * @param pos Position to erase.
     */
    void erase(int pos);

    /**
     * Add given array to receiver values starting at position si.
     * @param src Array to add from.
     * @param si Index to start adding from.
     */
    void addSubVector(const IntArray &src, int si);
    /**
     * Copy given array to receiver values starting at position si.
     * @param src Array to copy from.
     * @param si Index to start copying from.
     */
    void copySubVector(const IntArray &src, int si);
    /**
     * Sets the array to be a sequence of values.
     * @param x Number of integers in the sequence to follow.
     */
    void setValues(int n, ...);
    /**
     * Adds given scalar to all values of receiver
     * @param val Value to add.
     */
    void add(int val);

    /// Sets all component to zero.
    void zero();

    /// Prints receiver on stdout.
    void printYourself() const;

    /**
     * Breaks encapsulation. Avoid using this unless absolutely necessary.
     * @return Internal pointer to stored values.
     */
    int *givePointer() const { return values; }

    /**
     * Stores array to output stream.
     * @see FEMComponent
     */
    contextIOResultType storeYourself(DataStream *stream, ContextMode mode) const;
    /**
     * Restores array from image on stream.
     * @see FEMComponent
     */
    contextIOResultType restoreYourself(DataStream *stream, ContextMode mode);

#ifdef __PARALLEL_MODE
    /**@name Methods for packing/unpacking to/from communication buffer */
    //@{
    /**
     * Packs receiver into communication buffer.
     * @param buff Buffer to pack itself into.
     * @return Nonzero if successful.
     */
    int packToCommBuffer(CommunicationBuffer &buff) const;
    /**
     * Unpacks receiver from communication buffer.
     * @param buff Buffer from which unpack itself.
     * @return Nonzero if successful.
     */
    int unpackFromCommBuffer(CommunicationBuffer &buff);
    /**
     * Returns how much space is needed to pack receivers message.
     * @param buff Buffer used for packing.
     */
    int givePackSize(CommunicationBuffer &buff);
    //@}
#endif

#ifdef BOOST_PYTHON
    void __setitem__(int i, int val) { this->at(i+1) = val; }
    int __getitem__(int i) { return this->at(i+1); }
    IntArray copy() { IntArray result = *this; return result; }
    void beCopyOf(IntArray &src) { this->operator=(src); }
#endif
};


template< class operation >int
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



template< class operation >void quickSort(IntArray &arry, int l, int r, operation op) {
    if ( r <= l ) {
        return;
    }

    int i = quickSortPartition(arry, l, r, op);
    quickSort(arry, l, i - 1, op);
    quickSort(arry, i + 1, r, op);
}


/**
 * Sorts the receiver using quicksort algorithm.
 * @param op Function object, required to have member function int class::operator() (int, int),
 * must return a negative value if first argument is less than the second,
 * zero if the arguments are equal, and a positive number otherwise.
 * @param arry Array to sort.
 */
template< class operation >void sort(IntArray &arry, operation op) { quickSort(arry, 1, arry.giveSize(), op); }
} // end namespace oofem
#endif // intarray_h


