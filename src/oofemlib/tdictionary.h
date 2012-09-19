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

#ifndef tdictionary_h
#define tdictionary_h

#include "compiler.h"
#include "error.h"
#ifndef __MAKEDEPEND
 #include <cstdlib> // for NULL
#endif

namespace oofem {
template< class Key, class T >class TDictionaryIterator;
template< class Key, class T >class TDictionary;

/**
 * This class implements key/value associations,
 * A pair is used as an entry in a dictionary.

 * A pair has three components : its name (a character), its value (a class FEMComponent),
 * a pointer to the next pair in the dictionary.
 * Tasks:
 * - Returning its key, or its value, or the next pair ;
 * - Appending another pair to itself.
 * @note{ FEMComponent is assumed to have overloaded operator << }
 * @note{ FEMComponent is assumed to have printYourself() function }
 */
template< class Key, class T >class TPair
{
protected:
    /// Key, from client.
    Key key;
    /// Data, from client.
    T *data;
    /// Link to next TPair.
    TPair< Key, T > *next;
public:
    TPair(Key k, T *d) {
        key = k;
        data = d;
        next = NULL;
    }
    ~TPair() { delete data; }

    void append(TPair< Key, T > *p) { next = p; }
    int  giveKey()               { return key; }
    TPair< Key, T > *giveNext()      { return next; }
    T *giveValue()         { return data; }
};

/**
 * This class implements a linked list whose entries are TPairs (see below).
 * Dictionaries are typically used by degrees of freedom for storing their
 * unknowns.
 *
 * A dictionary stores its pairs in a linked list form. It knows the first
 * pair (attribute 'first') of the list. It also knows the last one (attribute
 * last) in order to append an additional pair quickly.
 * Tasks:
 * - Storing pairs (method add) and returning the value of a pair (method at).
 */
template< class Key, class T >class TDictionary
{
protected:
    TPair< Key, T > *first;
    TPair< Key, T > *last;

public:
    TDictionary() {
        first = NULL;
        last = NULL;
    }
    ~TDictionary();

    TPair< Key, T > *add(Key, T *);
    T *at(Key);
    bool includes(Key);
    void clear();

protected:
    friend class TDictionaryIterator< Key, T >;
};

template< class Key, class T >class TDictionaryIterator
{
protected:
    TPair< Key, T > *curr;
    TDictionary< Key, T > *list;
public:
    TDictionaryIterator(TDictionary< Key, T > *l);
    void initialize() { curr = list->first; }
    void initialize(TDictionary< Key, T > *l) {
        list = l;
        curr = l->first;
    }
    T *next();
};


template< class Key, class T >TDictionary< Key, T > :: ~TDictionary()
// Destructor.
{
    TPair< Key, T > *Next;

    while ( first ) {
        Next = first->giveNext();
        delete first;
        first = Next;
    }
}

template< class Key, class T >TPair< Key, T > *TDictionary< Key, T > :: add(Key k, T *v)
// Adds the pair (k,v) to the receiver. Returns this new pair.
{
    TPair< Key, T > *newPair;

#  ifdef DEBUG
    if ( this->includes(k) ) {
        OOFEM_ERROR("TDictionary::add : key already exists");
    }

#  endif

    newPair = new TPair< Key, T >(k, v);
    if ( last ) {
        last->append(newPair);
    } else {                              // empty dictionary
        first = newPair;
    }

    last = newPair;

    return newPair;
}

template< class Key, class T >T *TDictionary< Key, T > :: at(Key aKey)
// Returns the value of the pair which key is aKey. If such pair does
// not exist, creates it and assign value 0.
{
    TPair< Key, T > *next;

    next = first;
    while ( next ) {
        if ( next->giveKey() == aKey ) {
            return next->giveValue();
        }

        next = next->giveNext();
    }

    /*   newPair = this->add(aKey,0) ;          // pair does not exist yet
     *   return newPair->giveValue() ; */
    return NULL;
}

template< class Key, class T >bool TDictionary< Key, T > :: includes(Key aKey)
// Returns True if the receiver contains a pair which key is aKey, else
// returns False.
{
    TPair< Key, T > *next;

    next = first;
    while ( next ) {
        if ( next->giveKey() == aKey ) {
            return true;
        }

        next = next->giveNext();
    }

    return false;
}


template< class Key, class T >void
TDictionary< Key, T > :: clear()
{
    TPair< Key, T > *Next;

    while ( first ) {
        Next = first->giveNext();
        delete first;
        first = Next;
    }
}

template< class Key, class T >TDictionaryIterator< Key, T > :: TDictionaryIterator(TDictionary< Key, T > *l)
{
    list = l;
    curr = l->first;
}

template< class Key, class T >T *TDictionaryIterator< Key, T > :: next() {
  if (curr) {
    TPair< Key, T > *ret = curr;
    curr = curr->giveNext();
    if ( curr == list->last ) {
      curr = 0;
    }
    
    return ret->giveValue();
  } else {
    return 0;
  }
}
} // end namespace oofem
#endif // tdictionary_h
