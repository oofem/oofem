/* $Header: /home/cvs/bp/oofem/oofemlib/src/dictionr.h,v 1.9.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   ************************
//   *** CLASS DICTIONARY ***
//   ************************


#ifndef dictionr_h
#define dictionr_h

#include "compiler.h"
#include "pair.h"

#include "error.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <string>
#endif

namespace oofem {
class DataStream;

/**
 * This class implements a linked list whose entries are Pairs (see below).
 * A dictionary stores its pairs in a linked list form. It knows the first
 * pair (attribute 'first') of the list. It also knows the last one (attri-
 * bute 'last') in order to append fastly an additional pair.
 */
class Dictionary
{
    /*
     * This class implements a linked list whose entries are Pairs (see below).
     * Dictionaries are typically used by degrees of freedom for storing their
     * unknowns.
     * DESCRIPTION :
     * A dictionary stores its pairs in a linked list form. It knows the first
     * pair (attribute 'first') of the list. It also knows the last one (attri-
     * bute 'last') in order to append fastly an additional pair.
     * TASK :
     * Storing pairs (method 'add') and returning the value of a pair (method
     * 'at').
     */

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

    /**
     * Clears the receiver
     */
    void clear();
    /**
     * Adds a new Pair with given keyword and value into receiver.
     * @return new Pair with given keyword and value
     */
    Pair *add(int, double);
    /**
     * Returns the value of the pair which key is aKey. If such pair does
     * not exist, creates it and assign value 0.
     */
    double &at(int);
    /**
     * Returns True if the receiver contains a pair which key is aKey, else
     * returns False.
     */
    int          includes(int);
    /**
     * Prints the receiver on screen.
     */
    void         printYourself();
    /**
     * Formats itself as string.
     */
    void         formatAsString(std :: string &str);
    /// Returns number of pairs of receiver.
    int giveSize();

    /**
     * Saves the receiver contens (state) to given stream.
     * @return contextIOResultType value.
     * @exception throws an ContextIOERR exception if error encountered
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the receiver contens (state) from given stream.
     * @return contextIOResultType value.
     * @exception throws an ContextIOERR exception if error encountered
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};
} // end namespace oofem
#endif // dictionr_h








