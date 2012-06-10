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

#ifndef dictionr_h
#define dictionr_h

#include "compiler.h"
#include "pair.h"

#include "error.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#ifndef __MAKEDEPEND
 #include <string>
#endif

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
class Dictionary
{
protected:
    /// First pair
    Pair *first;
    /// Last pair
    Pair *last;

public:
    /// Constructor, creates empty dictionary
    Dictionary()  {
        first = NULL;
        last = NULL;
    }
    /// Destructor
    ~Dictionary();

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
    /**
     * Checks if dictionary includes given key
     * @param aKey Dictionary key.
     * @return True if receiver contains pair with given key, otherwise false.
     */
    bool includes(int aKey);
    /// Prints the receiver on screen.
    void printYourself();
    /// Formats itself as string.
    void formatAsString(std :: string &str);
    /// Returns number of pairs of receiver.
    int giveSize();

    /**
     * Saves the receiver contends (state) to given stream.
     * @return contextIOResultType value.
     * @exception throws an ContextIOERR exception if error encountered
     */
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the receiver contents (state) from given stream.
     * @return contextIOResultType value.
     * @exception throws an ContextIOERR exception if error encountered
     */
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};
} // end namespace oofem
#endif // dictionr_h
