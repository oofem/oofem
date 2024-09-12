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

#ifndef dictionr_h
#define dictionr_h

#include "oofemenv.h"
#include "pair.h"
#include "error.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#include <string>
#include <iosfwd>

namespace oofem {
class DataStream;

/**
 * This class implements a linked list whose entries are Pairs (see below).
 *
 * Dictionaries are typically used by degrees of freedom for storing their unknowns.
 * A dictionary stores its pairs in a linked list form. It knows the first
 * pair (attribute 'first') of the list. It also knows the last one (attribute
 * 'last') in order to append an additional pair fast.
 */
class OOFEM_EXPORT Dictionary
{
protected:
    /// First pair
    Pair *first;
    /// Last pair
    Pair *last;

public:
    /// Constructor, creates empty dictionary
    Dictionary() : first(NULL), last(NULL) { }
    /// Destructor
    ~Dictionary();

    /// copy constructor
    Dictionary(const Dictionary& other) {
        Pair *next = other.first;
        while ( next ) {
            this->add(next->giveKey(), next->giveValue());
            next = next->giveNext();
        }    
    }

    /// copy assignment constructor
    Dictionary& operator=(const Dictionary& other) { 
        this->clear();
        Pair *next = other.first;
        while ( next ) {
            this->add(next->giveKey(), next->giveValue());
            next = next->giveNext();
        }

        return *this;
    }

    /// Clears the receiver.
    void clear();
    /**
     * Adds a new Pair with given keyword and value into receiver.
     * @param aKey key of new pair
     * @param value value of new pair
     * @return New Pair with given keyword and value
     */
    Pair *add(int aKey, double value);
    /**
     * Returns the value of the pair which key is aKey.
     * If requested key doesn't exist, it is created with assigned value 0.
     * @param aKey Key for pair.
     * @return Reference to value of pair with given key
     */
    double &at(int aKey);
    double at(int aKey) const;
    /**
     * Checks if dictionary includes given key
     * @param aKey Dictionary key.
     * @return True if receiver contains pair with given key, otherwise false.
     */
    bool includes(int aKey) const;
    /// Prints the receiver on screen.
    void printYourself();
    /// Formats itself as string.
    void formatAsString(std :: string &str);
    /// Returns number of pairs of receiver.
    int giveSize();

    /**
     * Saves the receiver contends (state) to given stream.
     * @exception throws an ContextIOERR exception if error encountered
     */
    void saveContext(DataStream &stream);
    /**
     * Restores the receiver contents (state) from given stream.
     * @exception throws an ContextIOERR exception if error encountered
     */
    void restoreContext(DataStream &stream);

    friend std :: ostream &operator << ( std :: ostream & out, const Dictionary & r );
};
} // end namespace oofem
#endif // dictionr_h
