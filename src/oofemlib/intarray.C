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

#include "intarray.h"
#include "freestor.h"
#include "error.h"
#include "datastream.h"
#include "classtype.h"

#ifndef __MAKEDEPEND
 #include <cstdarg>
#endif

#ifdef __PARALLEL_MODE
 #include "combuff.h"
#endif

namespace oofem {
IntArray :: IntArray(int n)
// Constructor : creates an array of size n (filled with garbage).
{
    allocatedSize = size = n;
    if ( n ) {
        values = allocInt(size);
    } else {
        values = NULL;
    }
}


IntArray :: IntArray(const IntArray &src)
{
    // copy constructor
    allocatedSize = size = src.size;
    if ( size ) {
        values = allocInt(size);
    } else {
        values = NULL;
    }

    int *p2        = src.values;
    int *p3        = values;

    int i = size;
    while ( i-- ) {
        * p3++ = * p2++;
    }
}


IntArray :: ~IntArray()
{
    if ( values ) {
        freeInt(values);
    }
}


IntArray &IntArray :: operator = ( const IntArray & src )
{
    // assignment: cleanup and copy
    if ( values ) {
        freeInt(values);
    }


    allocatedSize = size = src.size;
    if ( size ) {
        values = allocInt(size);

        int *p2        = src.values;
        int *p3        = values;

        int i = size;
        while ( i-- ) {
            * p3++ = * p2++;
        }
    } else {
        values = NULL;
    }

    return * this;
}


void IntArray :: zero()
{
    int *p1 = values;

    int i = size;
    while ( i-- ) {
        * p1++ = 0;
    }
}


void IntArray :: add(int value)
{
    int *p1 = values;

    int i = size;
    while ( i-- ) {
        * p1 += value;
        p1++;
    }
}


#ifdef DEBUG
int &IntArray :: at(int i)
// Returns the i-th coefficient of the receiver. Slow but safe.
{
    this->checkBounds(i);
    return values [ i - 1 ];
}

int IntArray :: at(int i) const
// Returns the i-th coefficient of the receiver. Slow but safe.
{
    this->checkBounds(i);
    return values [ i - 1 ];
}
#endif


#ifdef DEBUG
void IntArray :: checkBounds(int i) const
// Checks that the receiver includes an index i.
{
    if ( i < 0 ) {
        OOFEM_ERROR2("IntArray::checkBounds : array error on index : %d < 0", i);
    }

    if ( i > size ) {
        OOFEM_ERROR3("IntArray::checkBounds : array error on index : %d > %d", i, size);
    }
}
#endif


void IntArray :: resize(int n, int allocChunk)
{
    int *p1, *p2, *newValues, i;

    if ( n <= allocatedSize ) {
        size = n;
        return;
    }

    if ( allocChunk < 0 ) {
        allocChunk = 0;
    }

    newValues = allocInt(n + allocChunk);

    p1 = values;
    p2 = newValues;
    i  = size;
    while ( i-- ) {
        * p2++ = * p1++;
    }

    if ( values ) {
        freeInt(values);
    }

    values = newValues;
    allocatedSize = n + allocChunk;
    size   = n;
}


void IntArray :: preallocate(int futureSize)
{
    int *p1, *p2, *newValues, i;

    if ( allocatedSize >= futureSize ) {
        return;
    }

    newValues = allocInt(futureSize);

    p1 = values;
    p2 = newValues;
    i  = size;
    while ( i-- ) {
        * p2++ = * p1++;
    }

    if ( values ) {
        freeInt(values);
    }

    values = newValues;
    allocatedSize = futureSize;
}


void IntArray :: followedBy(const IntArray &b, int allocChunk)
// Appends the array 'b' the receiver. Returns the receiver.
{
    register int i;
    int newSize;
    int *newValues, *p1, *p2, *p3;

    newSize = size + b.size;
    if ( newSize == size ) {
        return;
    }

    if ( allocChunk < 0 ) {
        allocChunk = 0;
    }

    if ( newSize > allocatedSize ) {
        newValues = allocInt(newSize + allocChunk);
        p3        = newValues;

        p1        = values;
        p2        = b.values;
        p3        = newValues;

        i = size;
        while ( i-- ) {
            * p3++ = * p1++;
        }

        i = b.size;
        while ( i-- ) {
            * p3++ = * p2++;
        }

        if ( values ) {
            freeInt(values);
        }

        values = newValues;
        allocatedSize = newSize + allocChunk;
        size   = newSize;
    } else {
        p1        = values + size;
        p2        = b.values;

        i = b.size;
        while ( i-- ) {
            * p1++ = * p2++;
        }

        size   = newSize;
    }
}


void IntArray :: followedBy(const int b, int allocChunk)
// Appends the array 'b' the receiver. Returns the receiver.
{
    register int i;
    int newSize;
    int *newValues, *p1, *p3;

    newSize = size + 1;

    if ( newSize > allocatedSize ) {
        newValues = allocInt(newSize + allocChunk);

        p1        = values;
        p3        = newValues;

        i = size;
        while ( i-- ) {
            * p3++ = * p1++;
        }

        * p3 = b;

        if ( values ) {
            freeInt(values);
        }

        values = newValues;
        allocatedSize = newSize + allocChunk;
        size   = newSize;
    } else {
        * ( values + size ) = b;
        size   = newSize;
    }
}


void IntArray :: erase(int _pos)
{
    // this will erase the element at given position (1-based index)
    // receiver size will shrink accordingly
    // the (_pos+1, size) elements will become (_pos, size-1) elements;

#ifdef DEBUG
    this->checkBounds(_pos);
#endif
    // adjust size; keep allocated size untouched
    size--;
    for ( int _i = _pos - 1; _i < size; _i++ ) {
        values [ _i ] = values [ _i + 1 ];
    }
}


bool IntArray :: containsOnlyZeroes() const
{
    for ( int i = 0; i < size; i++ ) {
        if ( values [ i ] ) {
            return false;
        }
    }

    return true;
}


int IntArray :: minimum() const
{
    if (size == 0) {
        OOFEM_ERROR("IntArray :: minimum - Empty array.");
    }
    int x = values[0];
    for ( int i = 1; i < size; ++i ) {
        if (values[i] < x) {
            x = values[i];
        }
    }
    return x;
}


int IntArray :: maximum() const
{
    if (size == 0) {
        OOFEM_ERROR("IntArray :: maximum - Empty array.");
    }
    int x = values[0];
    for ( int i = 1; i < size; ++i ) {
        if (values[i] > x) {
            x = values[i];
        }
    }
    return x;
}


void IntArray :: printYourself() const
// Prints the receiver on screen.
{
    printf("IntArray of size : %d\n", size);
    for ( int i = 1; i <= size; ++i ) {
        if ( i > 35 ) {
            printf("   (other components not printed)");
            break;
        } else {
            printf( "%d  ", this->at(i) );
        }
    }

    printf("\n");
}


contextIOResultType IntArray :: storeYourself(DataStream *stream, ContextMode mode) const
// writes receiver's binary image into stream
// use id to distinguish some instances
// return value >0 succes
//              =0 file i/o error
{
    int type_id = IntArrayClass;
    // write class header
    if ( !stream->write(& type_id, 1) ) {
        return ( CIO_IOERR );
    }

    // write size
    if ( !stream->write(& size, 1) ) {
        return ( CIO_IOERR );
    }

    // write raw data
    if ( !stream->write(values, size) ) {
        return ( CIO_IOERR );
    }

    // return result back
    return CIO_OK;
}

contextIOResultType IntArray :: restoreYourself(DataStream *stream, ContextMode mode)
// reads receiver from stream
// warning - overwrites existing data!
// returns 0 if file i/o error
//        -1 if id od class id is not correct
{
    int class_id;
    // read class header
    if ( !stream->read(& class_id, 1) ) {
        return ( CIO_IOERR );
    }

    if ( class_id != IntArrayClass ) {
        return CIO_BADVERSION;
    }

    // read size
    if ( !stream->read(& size, 1) ) {
        return ( CIO_IOERR );
    }

    if ( values != NULL ) {
        freeInt(values);
    }

    if ( size ) {
        values = allocInt(size);
        allocatedSize = size;
    } else {
        values = NULL;
        allocatedSize = 0;
    }

    // write raw data
    if ( !stream->read(values, size) ) {
        return ( CIO_IOERR );
    }

    // return result back
    return CIO_OK;
}


int IntArray :: findFirstIndexOf(int value)   const
{
    // finds index of value in receiver
    // if such value  does not exists, returns zero index
    int i;
    for ( i = 0; i < size; i++ ) {
        if ( values [ i ] == value ) {
            return i + 1;
        }
    }

    // nothing found
    return 0;
}


void IntArray :: addSubVector(const IntArray &src, int si)
{
    int i, reqSize, n = src.giveSize();

    si--;
    reqSize = si + n;
    if ( this->giveSize() < reqSize ) {
        this->resize(reqSize);
    }

    for ( i = 1; i <= n; i++ ) {
        this->at(si + i) += src.at(i);
    }
}


void IntArray :: copySubVector(const IntArray &src, int si)
{
    int i, reqSize, n = src.giveSize();

    si--;
    reqSize = si + n;
    if ( this->giveSize() < reqSize ) {
        this->resize(reqSize);
    }

    for ( i = 1; i <= n; i++ ) {
        this->at(si + i) = src.at(i);
    }
}


void IntArray :: setValues(int n, ...)
{
    va_list vl;
    va_start(vl,n);
    this->resize(n);
    for (int i = 0; i < n; i++ ) {
        this->values [ i ] = va_arg(vl, int);
    }
    va_end(vl);
}


int IntArray :: findSorted(int _val)   const
{
    int first = 0;
    int last = size - 1;
    int mid;

    while ( first <= last ) { //while we haven't reached the end of
        mid = ( first + last ) / 2; //rule out half of the data by spliting the array
        if ( values [ mid ] == _val ) { //if we have found the target
            return mid + 1;
        } else if ( values [ mid ] > _val ) {
            last = mid - 1; //the desired value is in the lower part of the array
        } else {
            first = mid + 1; //the desired value is in the upper part of the array
        }
    }

    return 0; //there is no such element in the array
}


int IntArray :: insertSorted(int _val, int allocChunk)
{
    int pos, i = size;
    int newSize = size + 1;
    int *newValues = NULL, *p1, *p2;

    if ( newSize > allocatedSize ) { // realocate if needed
        newValues = allocInt(newSize + allocChunk);

        p1        = values;
        p2        = newValues;
    } else {
        p1 = p2 = values;
    }

    while ( ( i > 0 ) && values [ i - 1 ] > _val ) { // copy values larger than _val into destination
        p2 [ i ] = p1 [ i - 1 ];
        i--;
    }

    pos = i;
    p2 [ pos ] = _val; // insert _val
    i--;

    // if allocated copy smaller values to destination
    if ( newSize > allocatedSize ) {
        // if allocated copy smaller values to destination
        while ( i >= 0 ) {
            p2 [ i ] = p1 [ i ];
            i--;
        }

        // dealocate original space
        if ( values ) {
            freeInt(values);
        }

        values = newValues;
        allocatedSize = newSize + allocChunk;
    }

    // update size
    size   = newSize;
    return pos + 1; // return 1-based index
}


int IntArray :: insertSortedOnce(int _val, int allocChunk)
{
    int res;
    if ( ( res = findSorted(_val) ) ) {
        return res;                          // value already present
    }

    return insertSorted(_val, allocChunk);
}


void IntArray :: eraseSorted(int value)
{
    int pos;

    if ( ( pos = findSorted(value) ) ) {
        erase(pos);
    }
}


int IntArray :: findCommonValuesSorted(const IntArray &iarray, IntArray &common, int allocChunk) const
{
    int i = 0, j, val;

    for ( j = 1; j <= iarray.giveSize(); j++ ) {
        val = iarray.at(j);

        while ( i < size ) {
            if ( values [ i ] == val ) {
                common.followedBy(val, allocChunk);
                i++;
                break;
            }

            if ( values [ i ] > val ) {
                break;
            }

            i++;
        }

        if ( i == size ) {
            break;
        }
    }

    return ( common.giveSize() );
}


int IntArray :: insertOnce(int _p)
{
    if ( !this->findFirstIndexOf(_p) ) {
        this->followedBy(_p, 2);
    }

    return size;
}


#ifdef __PARALLEL_MODE
int IntArray :: packToCommBuffer(CommunicationBuffer &buff) const
{
    int result = 1;
    // pack size
    result &= buff.packInt(size);
    // pack data
    result &= buff.packArray(this->values, size);

    return result;
}


int IntArray :: unpackFromCommBuffer(CommunicationBuffer &buff)
{
    int newSize, result = 1;
    // unpack size
    result &= buff.unpackInt(newSize);
    // resize yourself
    this->resize(newSize);
    result &= buff.unpackArray(this->values, newSize);

    return result;
}


int IntArray :: givePackSize(CommunicationBuffer &buff)
{
    return buff.givePackSize(MPI_INT, 1) + buff.givePackSize(MPI_INT, this->size);
}
#endif
} // end namespace oofem
