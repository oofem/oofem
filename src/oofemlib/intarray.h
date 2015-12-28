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

#ifndef intarray_h
#define intarray_h

#include "oofemcfg.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#include <cstdio>
#include <vector>
#include <iosfwd>
#include <initializer_list>

namespace oofem {
class DataStream;

/**
 * Class implementing an array of integers.
 * Acts as a wrapper for std::vector of integers, but adds many convenience functions that are commonly used in OOFEM.
 * Array can grow or shrink to desired dimension.
 * The lower value index of array is 1, upper depends on array size.
 * 
 * @author Mikael Öhman
 * @author Jim Brouzoulis 
 * @author many others (please add yourselves)
 */
class OOFEM_EXPORT IntArray
{
private:
    /// Stored values.
    std::vector< int > values;

public:
    /// @name Iterator for for-each loops:
    //@{
    std::vector< int > :: iterator begin() { return this->values.begin(); }
    std::vector< int > :: iterator end() { return this->values.end(); }
    std::vector< int > :: const_iterator begin() const { return this->values.begin(); }
    std::vector< int > :: const_iterator end() const { return this->values.end(); }
    //@}

    /// Constructor for sized array. Data is zeroed.
    IntArray(int n = 0) : values(n) { }
    /// Copy constructor. Creates the array from another array.
    IntArray(const IntArray &src) : values(src.values) { }
    /// Move constructor. Creates the array from another array.
    IntArray(IntArray &&src) : values(std::move(src.values)) { }
    /// Initializer list constructor.
    inline IntArray(std :: initializer_list< int >list) : values(list) { }
    /// Destructor.
    ~IntArray() {};

    /// Assignment operator
    IntArray &operator = (const IntArray &src) { values = src.values; return *this; }
    /// Move operator
    IntArray &operator = (IntArray &&src) { values = std::move(src.values); return *this; }
    /// Assignment operator.
    inline IntArray &operator = (std :: initializer_list< int >list) { values = list; return *this; }

    /**
     * Coefficient access function. Returns l-value of coefficient at given
     * position of the receiver.
     * @param i Position of coefficient in array.
     * @return Value at position.
     */
#ifdef DEBUG
    int &at(int i);
#else
    inline int &at(int i) { return values [ i - 1 ]; }
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
    inline int at(int i) const { return values [ i - 1 ]; }
#endif
    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i Position of coefficient in array.
     * @return Value at position.
     */
#ifdef DEBUG
    int &operator() (int i);
#else
    inline int &operator() (int i) { return values [ i ]; }
#endif
    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i position of coefficient in array.
     * @return Value at position.
     */
#ifdef DEBUG
    const int &operator() (int i) const;
#else
    inline const int &operator() (int i) const { return values [ i ]; }
#endif

#ifdef DEBUG
    int &operator[] ( int i );
#else
    inline int &operator[] ( int i ) { return values [ i ]; }
#endif

#ifdef DEBUG
    const int &operator[] ( int i ) const;
#else
    inline const int &operator[] ( int i ) const { return values [ i ]; }
#endif

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
     * If dimension mismatch, size is adjusted accordingly and memory is copied over.
     * @param n New size of array.
     * @param allocChunk If reallocation needed, an additional space for allocChunk values will be allocated
     * to prevent excessive reallocation.
     */
    void resizeWithValues(int n, int allocChunk = 0);
    /**
     * Checks size of receiver towards requested bounds.
     * Data is always zeroed.
     * @param n New size of array.
     */
    void resize(int n);
    /**
     * Clears the array (zero size).
     */
    void clear() { this->values.clear(); }
    /**
     * Preallocates receiver to given futureSize if larger then allocatedSize.
     * @param futureSize Size to be allocated.
     */
    void preallocate(int futureSize);
    /**
     * Resizes receiver and enumerates from 1 to the maximum value given.
     * @param maxVal to enumerate to.
     */
    void enumerate(int maxVal);
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
    void followedBy(int b, int allocChunk = 0);
    /// @return Size of receiver.
    int giveSize() const { return (int)values.size(); }
    /**
     * Checks if receiver is empty (i.e., zero sized).
     * @return True is size is zero.
     */
    bool isEmpty() const { return values.size() == 0; }
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
     * Finds all indices where the input array is nonzero
     * @param logical Array of logical values (0 = false, nonzero = true) to have index extracted.
     */
    void findNonzeros(const IntArray &logical);
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
     */
    void insertSorted(int value, int allocChunk = 0);
    /**
     * Inserts given value into a receiver, which is assumed to be sorted.
     * The value is inserted only if it does not exist.
     * The size of receiver is changed accordingly.
     * @param value Value to insert.
     * @param allocChunk If reallocation needed, an additional space for allocChunk values will be allocated.
     */
    void insertSortedOnce(int value, int allocChunk = 0);
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
     */
    void insertOnce(int p);
    /**
     * Sorts array.
     */
    void sort();
    /**
     * Erase the element at given position (1-based index)
     * Receiver will shrink accordingly, the values at positions (pos+1,...,size)
     * will be moved to positions (pos,...,size-1)
     * @param pos Position to erase.
     */
    void erase(int pos);
    /**
     * Adds given scalar to all values of receiver
     * @param val Value to add.
     */
    void add(int val);

    /// Sets all component to zero.
    void zero();

    /// Prints receiver on stdout.
    void printYourself() const;

    /// Abbreviation for printYourself().
    void pY() const;

    /**
     * Prints receiver on stdout with custom name.  
     * @param name Display name of reciever.
     */
    void printYourself(const std::string name) const;    
    /**
     * Breaks encapsulation. Avoid using this unless absolutely necessary.
     * @return Internal pointer to stored values.
     */
    inline const int *givePointer() const { return values.data(); }
    inline int *givePointer() { return values.data(); }

    /**
     * Stores array to output stream.
     */
    contextIOResultType storeYourself(DataStream &stream) const;
    /**
     * Restores array from image on stream.
     */
    contextIOResultType restoreYourself(DataStream &stream);
    /**
     * Returns how much space is needed to pack receivers message.
     * @param buff Buffer used for packing.
     */
    int givePackSize(DataStream &buff) const;


    friend std :: ostream &operator << ( std :: ostream & out, const IntArray & x );

#ifdef BOOST_PYTHON
    void __setitem__(int i, int val) { this->at(i + 1) = val; }
    int __getitem__(int i) { return this->at(i + 1); }
    void beCopyOf(const IntArray &src) { this->operator = ( src ); }
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
 * Sorts the receiver using quicksort algorithm.
 * @param op Function object, required to have member function int class::operator() (int, int),
 * must return a negative value if first argument is less than the second,
 * zero if the arguments are equal, and a positive number otherwise.
 * @param arry Array to sort.
 */
template< class operation > void sort(IntArray &arry, operation op) { quickSort(arry, 1, arry.giveSize(), op); }
} // end namespace oofem
#endif // intarray_h
