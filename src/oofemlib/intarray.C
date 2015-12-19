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

#include "intarray.h"
#include "error.h"
#include "datastream.h"

#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <memory>

namespace oofem {

void IntArray :: zero()
{
    std::fill(values.begin(), values.end(), 0);
}


void IntArray :: add(int value)
{
    for (int &x: values) x += value;
}


#ifdef DEBUG
int &IntArray :: at(int i)
{
    this->checkBounds(i);
    return values [ i - 1 ];
}

int IntArray :: at(int i) const
{
    this->checkBounds(i);
    return values [ i - 1 ];
}

int &IntArray :: operator()(int i)
{
    this->checkBounds(i);
    return values [ i ];
}

const int &IntArray :: operator()(int i) const
{
    this->checkBounds(i);
    return values [ i ];
}

int &IntArray :: operator[](int i)
{
    this->checkBounds(i);
    return values [ i ];
}

const int &IntArray :: operator[](int i) const
{
    this->checkBounds(i);
    return values [ i ];
}

void IntArray :: checkBounds(int i) const
// Checks that the receiver includes an index i.
{
    if ( i < 0 ) {
        OOFEM_ERROR("array error on index : %d < 0", i);
    }

    if ( i > this->giveSize() ) {
        OOFEM_ERROR("array error on index : %d > %d", i, this->giveSize());
    }
}
#endif


void IntArray :: resizeWithValues(int n, int allocChunk)
{
    if ( allocChunk > 0 && (int)this->values.capacity() < n ) {
        this->values.reserve(n + allocChunk);
    }
    this->values.resize(n);
}


void IntArray :: resize(int n)
{
    this->values.assign(n, 0);
}


void IntArray :: preallocate(int futureSize)
{
    values.reserve(futureSize);
}


void IntArray :: enumerate(int maxValue)
{
    this->values.resize(maxValue);
    for ( int i = 1; i <= maxValue; ++i ) {
        this->at(i) = i;
    }
}


void IntArray :: followedBy(const IntArray &b, int allocChunk)
{
    if ( allocChunk && (int)values.capacity() < this->giveSize() + b.giveSize() ) {
        values.reserve(values.capacity() + allocChunk + b.giveSize());
    }
    values.insert(values.end(), b.values.begin(), b.values.end());
}


void IntArray :: followedBy(int b, int allocChunk)
{
    if ( allocChunk && (int)values.capacity() < this->giveSize() + 1 ) {
        values.reserve(values.capacity() + allocChunk + 1);
    }
    values.push_back(b);
}


void IntArray :: erase(int _pos)
{
#ifdef DEBUG
    this->checkBounds(_pos);
#endif
    values.erase(values.begin() + _pos - 1);
}


bool IntArray :: containsOnlyZeroes() const
{
    for ( auto x: values ) {
        if ( x ) {
            return false;
        }
    }

    return true;
}


int IntArray :: minimum() const
{
#ifdef DEBUG
    if ( this->isEmpty() ) {
        OOFEM_ERROR("Empty array.");
    }
#endif
    return *std::min_element(values.begin(), values.end());
}


int IntArray :: maximum() const
{
#ifdef DEBUG
    if ( this->isEmpty() ) {
        OOFEM_ERROR("Empty array.");
    }
#endif
    return *std::max_element(values.begin(), values.end());
}


void IntArray :: findNonzeros(const IntArray &logical)
{
    int newsize = 0;
    for ( const int &x: logical.values) {
        if ( x ) {
            ++newsize;
        }
    }
    this->values.resize(newsize);

    int pos = 1;
    for ( int i = 1; i <= logical.giveSize(); ++i ) {
        if ( logical.at(i) ) {
            this->at(pos++) = i;
        }
    }
}


void IntArray :: printYourself() const
// Prints the receiver on screen.
{
    printf("IntArray of size : %d\n", this->giveSize());
    for ( int i = 1; i <= this->giveSize(); ++i ) {
        if ( i > 42 ) {
            printf("   (other components not printed)");
            break;
        } else {
            printf( "%d  ", this->at(i) );
        }
    }

    printf("\n");
}


void IntArray :: printYourself(const std::string name) const
// Prints the receiver on screen.
{
    printf("%s (%d): ", name.c_str(), this->giveSize());
    for ( int i = 1; i <= this->giveSize(); ++i ) {
        if ( i > 42 ) {
            printf("   (other components not printed)");
            break;
        } else {
            printf( "%d  ", this->at(i) );
        }
    }

    printf("\n");
}

void IntArray :: pY() const {
    printYourself();
}


contextIOResultType IntArray :: storeYourself(DataStream &stream) const
{
    // write size
    if ( !stream.write(this->giveSize()) ) {
        return ( CIO_IOERR );
    }

    // write raw data
    if ( !stream.write(values.data(), this->giveSize()) ) {
        return ( CIO_IOERR );
    }

    // return result back
    return CIO_OK;
}

contextIOResultType IntArray :: restoreYourself(DataStream &stream)
{
    // read size
    int size;
    if ( !stream.read(size) ) {
        return ( CIO_IOERR );
    }

    values.resize(size);

    // read raw data
    if ( !stream.read(values.data(), size) ) {
        return ( CIO_IOERR );
    }

    // return result back
    return CIO_OK;
}


int IntArray :: givePackSize(DataStream &buff) const
{
    return buff.givePackSizeOfInt(1) + buff.givePackSizeOfInt(this->giveSize());
}


int IntArray :: findFirstIndexOf(int value) const
{
    // finds index of value in receiver
    auto it = std::find(values.begin(), values.end(), value);
    // if such value  does not exists, returns zero index
    if ( it == values.end() ) {
        return 0;
    } else {
        return (int)(it - values.begin() + 1);
    }
}


int IntArray :: findSorted(int _val)   const
{
    return std::binary_search (values.begin(), values.end(), _val);
}


void IntArray :: insertSorted(int val, int allocChunk)
{
    if ( allocChunk > 0 && values.size() + 1 >= values.capacity() ) {
        values.reserve(allocChunk + values.capacity());
    }
    auto low = std::lower_bound(values.begin(), values.end(), val);
    values.insert(low, val);
}


void IntArray :: insertSortedOnce(int val, int allocChunk)
{
    if ( allocChunk > 0 && values.size() + 1 >= values.capacity() ) {
        values.reserve(allocChunk + values.capacity());
    }
    auto low = std::lower_bound(values.begin(), values.end(), val);
    if ( low == values.end() || *low != val ) {
        values.insert(low, val);
    }
}


void IntArray :: eraseSorted(int value)
{
    auto low = std::lower_bound(values.begin(), values.end(), value);
    if ( *low == value ) {
        values.erase(low);
    }
}


int IntArray :: findCommonValuesSorted(const IntArray &iarray, IntArray &common, int allocChunk) const
{
    int i = 0;

    for ( int val: iarray ) {

        while ( i < this->giveSize() ) {
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

        if ( i == this->giveSize() ) {
            break;
        }
    }

    return ( common.giveSize() );
}


void IntArray :: insertOnce(int _p)
{
    if ( !this->findFirstIndexOf(_p) ) {
        this->followedBy(_p, 2);
    }
}


void IntArray :: sort()
{
    std::sort(this->begin(), this->end());
}


std :: ostream &operator<<(std :: ostream &out, const IntArray &x)
{
    out << x.giveSize();
    for ( const int &val: x ) {
        out << " " << val;
    }
    return out;
}
} // end namespace oofem
