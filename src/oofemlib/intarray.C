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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#include "intarray.h"
#include "error.h"
#include "datastream.h"

#include <cstdarg>
#include <cstdlib>
#include <cstring>
#if __cplusplus > 199711L
 #include <memory>
#endif

#ifdef __PARALLEL_MODE
 #include "combuff.h"
#endif

#define ALLOC(size) new int [ size ];

#define RESIZE(n) \
    { \
        size = n; \
        if ( n > allocatedSize ) { \
            allocatedSize = n; \
            if ( values ) { delete[] values; } \
            values = ALLOC(size); \
        } \
    }

namespace oofem {
IntArray :: iterator :: iterator(const IntArray *vec, int pos) : pos(pos), vec(vec) { }

bool IntArray :: iterator :: operator != ( const IntArray :: iterator & other ) const {
    return pos != other.pos;
}

int IntArray :: iterator :: operator *( ) const {
    return ( * vec )( pos );
}

const IntArray :: iterator &IntArray :: iterator :: operator++ ( ) {
    ++pos;
    return * this;
}

IntArray :: iterator IntArray :: begin() { return iterator(this, 0); }
IntArray :: iterator IntArray :: end() { return iterator(this, this->size); }


IntArray :: IntArray() :
    size(0),
    allocatedSize(0),
    values(NULL)
{ }


IntArray :: IntArray(int n) :
    size(n),
    allocatedSize(n)
    // Constructor : creates an array of size n (filled with garbage).
{
    if ( size ) {
        values = ALLOC(size);
        memset( values, 0, size * sizeof( int ) );
    } else {
        values = NULL;
    }
}


IntArray :: IntArray(const IntArray &src) :
    size(src.size),
    allocatedSize(src.size)
    // copy constructor
{
    if ( size ) {
        values = ALLOC(this->size);
        memcpy( values, src.values, size * sizeof( int ) );
    } else {
        values = NULL;
    }
}


#if __cplusplus > 199711L
IntArray :: IntArray(std :: initializer_list< int >list)
{
    this->size = this->allocatedSize = list.size();
    if ( this->size ) {
        this->values = ALLOC(this->size);
        std :: uninitialized_copy(list.begin(), list.end(), this->values);
    } else {
        this->values = NULL;
    }
}


IntArray &IntArray :: operator = ( std :: initializer_list< int >list )
{
    RESIZE( ( int ) list.size() );
    std :: uninitialized_copy(list.begin(), list.end(), this->values);
    return * this;
}

#endif


IntArray :: ~IntArray()
{
    if ( values ) {
        delete[] values;
    }
}


IntArray &IntArray :: operator = ( const IntArray & src )
{
    // assignment: cleanup and copy
    if ( values ) {
        delete[] values;
    }

    allocatedSize = size = src.size;
    if ( size ) {
        values = ALLOC(this->size);
        memcpy( values, src.values, size * sizeof( int ) );
    } else {
        values = NULL;
    }

    return * this;
}


void IntArray :: zero()
{
    memset( values, 0, size * sizeof( int ) );
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

int &IntArray :: operator() (int i)
{
    this->checkBounds(i);
    return values [ i ];
}

const int &IntArray :: operator() (int i) const
{
    this->checkBounds(i);
    return values [ i ];
}

int &IntArray :: operator[] ( int i )
{
    this->checkBounds(i);
    return values [ i ];
}

const int &IntArray :: operator[] ( int i ) const
{
    this->checkBounds(i);
    return values [ i ];
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


void IntArray :: resizeWithValues(int n, int allocChunk)
{
#ifdef DEBUG
    if ( allocChunk < 0 ) {
        OOFEM_FATAL2("IntArray :: resize - allocChunk must be non-negative; %d", allocChunk);
    }

#endif

    if ( n <= allocatedSize ) {
        size = n;
        return;
    }
    allocatedSize = n + allocChunk;

    // For the typical small sizes we have, realloc doesn't seem to be worth it, better with just malloc.
    int *newValues = ALLOC(allocatedSize);
#ifdef DEBUG
    if ( !newValues ) {
        OOFEM_FATAL2("IntArray :: resizeWithValues - Failed in allocating %d doubles", allocatedSize);
    }
#endif
    memcpy( newValues, values, size * sizeof( int ) );
    memset( & newValues [ size ], 0, ( allocatedSize - size ) * sizeof( int ) );

    if ( values ) {
        delete[] values;
    }
    values = newValues;
    size = n;
}

void IntArray :: resize(int n)
{
    size = n;
    if ( n <= allocatedSize ) {
        memset( values, 0, size * sizeof( int ) );
        return;
    }
    allocatedSize = n;

    if ( values ) {
        delete[] values;
    }
    values = ALLOC(allocatedSize);
    memset( values, 0, allocatedSize * sizeof( int ) );
#ifdef DEBUG
    if ( !values ) {
        OOFEM_FATAL2("IntArray :: resize - Failed in allocating %d doubles", allocatedSize);
    }
#endif
}


void IntArray :: preallocate(int futureSize)
{
    if ( allocatedSize >= futureSize ) {
        return;
    }
    allocatedSize = futureSize;

    int *newValues = ALLOC(allocatedSize);
#ifdef DEBUG
    if ( !newValues ) {
        OOFEM_FATAL2("IntArray :: preallocate - Failed in allocating %d doubles", allocatedSize);
    }
#endif
    memcpy( newValues, values, size * sizeof( int ) );
    memset( & newValues [ size ], 0, ( allocatedSize - size ) * sizeof( int ) );

    if ( values ) {
        delete[] values;
    }
    values = newValues;
}


void IntArray :: enumerate(int maxValue)
{
    this->resize(maxValue);
    for ( int i = 1; i <= maxValue; ++i ) {
        this->at(i) = i;
    }
}


void IntArray :: followedBy(const IntArray &b, int allocChunk)
// Appends the array 'b' the receiver. Returns the receiver.
{
    int newSize = size + b.size;
    if ( newSize == size ) {
        return;
    }

    if ( allocChunk < 0 ) {
        allocChunk = 0;
    }

    if ( newSize > allocatedSize ) {
        int *newValues = ALLOC(newSize + allocChunk);

        memcpy( newValues, values, size * sizeof( int ) );
        memcpy( newValues + size, b.values, b.size * sizeof( int ) );
        ///@todo Do we zero or just leave the last part uninitialized?
        memset( newValues + newSize, 0, allocChunk * sizeof( int ) );

        if ( values ) {
            delete[] values;
        }
        values = newValues;
        allocatedSize = newSize + allocChunk;
        size = newSize;
    } else {
        memcpy( values + size, b.values, b.size * sizeof( int ) );
        size = newSize;
    }
}


void IntArray :: followedBy(int b, int allocChunk)
// Appends the array 'b' the receiver. Returns the receiver.
{
    int newSize = size + 1;

    if ( newSize > allocatedSize ) {
        int *newValues = ALLOC(newSize + allocChunk);
        memset( newValues, 0, ( newSize + allocChunk ) * sizeof( int ) );

        memcpy( newValues, values, size * sizeof( int ) );
        newValues [ size ] = b;

        if ( values ) {
            delete[] values;
        }

        values = newValues;
        allocatedSize = newSize + allocChunk;
        size = newSize;
    } else {
        * ( values + size ) = b;
        size = newSize;
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
#ifdef DEBUG
    if ( size == 0 ) {
        OOFEM_ERROR("IntArray :: minimum - Empty array.");
    }
#endif
    int x = values [ 0 ];
    for ( int i = 1; i < size; ++i ) {
        if ( values [ i ] < x ) {
            x = values [ i ];
        }
    }
    return x;
}


int IntArray :: maximum() const
{
#ifdef DEBUG
    if ( size == 0 ) {
        OOFEM_ERROR("IntArray :: maximum - Empty array.");
    }
#endif
    int x = values [ 0 ];
    for ( int i = 1; i < size; ++i ) {
        if ( values [ i ] > x ) {
            x = values [ i ];
        }
    }
    return x;
}


void IntArray :: findNonzeros(const IntArray &logical)
{
    int newsize = 0;
    for ( int i = 1; i <= logical.size; ++i ) {
        if ( logical.at(i) ) {
            ++newsize;
        }
    }
    this->resize(newsize);

    int pos = 1;
    for ( int i = 1; i <= logical.size; ++i ) {
        if ( logical.at(i) ) {
            this->at(pos++) = i;
        }
    }
}


void IntArray :: printYourself() const
// Prints the receiver on screen.
{
    printf("IntArray of size : %d\n", size);
    for ( int i = 1; i <= size; ++i ) {
        if ( i > 42 ) {
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
    // read size
    if ( !stream->read(& size, 1) ) {
        return ( CIO_IOERR );
    }

    if ( values != NULL ) {
        delete[] values;
    }

    if ( size ) {
        values = ALLOC(size);
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
    for ( int i = 0; i < size; i++ ) {
        if ( values [ i ] == value ) {
            return i + 1;
        }
    }

    // nothing found
    return 0;
}


void IntArray :: setValues(int n, ...)
{
    va_list vl;
    va_start(vl, n);
    this->resize(n);
    for ( int i = 0; i < n; i++ ) {
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
        newValues = ALLOC(newSize + allocChunk);
        memset( values, 0, ( newSize + allocChunk ) * sizeof( int ) );

        p1 = values;
        p2 = newValues;
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
            delete[] values;
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
    int first = 0;
    int last = size - 1;

    if ( size == 0 || values [ first ] > _val ) {
        // first = 0 is correct already
    } else if ( values [ last ] < _val ) {
        first = size;
    } else {
        while ( first <= last ) { //while we haven't reached the end of
            int mid = ( first + last ) / 2; //rule out half of the data by spliting the array
            if ( values [ mid ] == _val ) { //if we have found the target
                return mid + 1;
            } else if ( values [ mid ] > _val ) {
                last = mid - 1; //the desired value is in the lower part of the array
            } else {
                first = mid + 1; //the desired value is in the upper part of the array
            }
        }
    }
    // After this loop, "first" is the new position to fill in

    int newSize = size + 1;
    if ( newSize > allocatedSize ) {
        int *newValues = ALLOC(newSize + allocChunk);
        if ( values ) {
            memcpy(newValues, values, sizeof( int ) * first); // Copy over the values before _val
            memcpy( & newValues [ first + 1 ], & values [ first ], sizeof( int ) * ( size - first ) ); // Copy over the values after _val, leaving a gap
            // Replace the old data:
            delete[] values;
        }
        values = newValues;
        allocatedSize = newSize + allocChunk;
        size = newSize;
    } else {
        memmove( & values [ first + 1 ], & values [ first ], sizeof( int ) * ( size - first ) );
        size = newSize;
    }
    // Finally we can insert the new value:
    values [ first ] = _val;

    return first;
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
    int i = 0, val;

    for ( int j = 1; j <= iarray.giveSize(); j++ ) {
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

std :: ostream &operator << ( std :: ostream & out, const IntArray & x )
{
    out << x.size;
    for ( int i = 0; i < x.size; ++i ) {
        out << " " << x(i);
    }
    return out;
}
} // end namespace oofem
