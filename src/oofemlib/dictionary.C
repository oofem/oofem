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

#include "dictionary.h"
#include "logger.h"
#include "datastream.h"
#include "contextioerr.h"
#include "contextmode.h"

#include <cstdlib>
#include <ostream>

namespace oofem {
Dictionary :: ~Dictionary()
// Destructor.
{
    this->clear();
}

void
Dictionary :: clear()
{
    Pair *Next;

    while ( first ) {
        Next = first->giveNext();
        delete first;
        first = Next;
    }

    first = NULL;
    last = NULL;
}

int
Dictionary :: giveSize()
{
    Pair *next;
    int size = 0;

    next = first;
    while ( next ) {
        size++;
        next = next->giveNext();
    }

    return size;
}

Pair *Dictionary :: add(int k, double v)
// Adds the pair (k,v) to the receiver. Returns this new pair.
{
    Pair *newPair;

#  ifdef DEBUG
    if ( this->includes(k) ) {
        OOFEM_ERROR("key (%d) already exists", k);
    }

#  endif

    newPair = new Pair(k, v);
    if ( last ) {
        last->append(newPair);
    } else {                              // empty dictionary
        first = newPair;
    }

    last = newPair;

    return newPair;
}


double &Dictionary :: at(int aKey)
// Returns the value of the pair which key is aKey. If such pair does
// not exist, creates it and assign value 0.
{
    Pair *next, *newPair;

    next = first;
    while ( next ) {
        if ( next->giveKey() == aKey ) {
            return next->giveValue();
        }

        next = next->giveNext();
    }

    newPair = this->add(aKey, 0);         // pair does not exist yet
    return newPair->giveValue();
}


bool Dictionary :: includes(int aKey)
// Returns True if the receiver contains a pair which key is aKey, else
// returns False.
{
    Pair *next;

    next = first;
    while ( next ) {
        if ( next->giveKey() == aKey ) {
            return true;
        }

        next = next->giveNext();
    }

    return false;
}


void Dictionary :: printYourself()
// Prints the receiver on screen.
{
    Pair *next;

    printf("Dictionary : \n");

    next = first;
    while ( next ) {
        next->printYourself();
        next = next->giveNext();
    }
}


void
Dictionary :: formatAsString(std :: string &str)
{
    Pair *next;
    char buffer [ 64 ];

    next = first;
    while ( next ) {
        sprintf( buffer, " %c %e", next->giveKey(), next->giveValue() );
        str += buffer;
        next = next->giveNext();
    }
}


contextIOResultType Dictionary :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
    int nitems = 0;
    int key;
    double value;
    Pair *next;

    next = first;
    while ( next ) {
        nitems++;
        next = next->giveNext();
    }

    // write size
    if ( !stream.write(nitems) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write raw data
    next = first;
    while ( next ) {
        key = next->giveKey();
        value = next->giveValue();
        if ( !stream.write(key) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.write(value) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        next = next->giveNext();
    }

    // return result back
    return CIO_OK;
}


contextIOResultType Dictionary :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
    int size;
    int key;
    double value;

    // delete currently occupied space
    this->clear();

    // read size
    if ( !stream.read(size) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read particular pairs
    for ( int i = 1; i <= size; i++ ) {
        if ( !stream.read(key) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(value) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        this->at(key) = value;
    }

    return CIO_OK;
}


std :: ostream &operator << ( std :: ostream & out, const Dictionary & r )
{
    int count = 0;
    Pair *next = r.first;
    while ( next ) {
        count++;
        next = next->giveNext();
    }

    out << count;
    next = r.first;
    while ( next ) {
        out << " " << next->giveKey() << " " << next->giveValue();
        next = next->giveNext();
    }
    return out;
}
} // end namespace oofem
