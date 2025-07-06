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

#ifndef pair_h
#define pair_h

#include <cstdio>

namespace oofem {
/**
 * This class implements key/value associations - the key and its associated value.
 * An instance of Pair is used as an entry in a dictionary.
 * Pair has three components - its key, its value and pointer to the next Pair in the dictionary.
 *
 * Tasks:
 * - Returning its key, or its value, or the next pair ;
 * - Appending another pair to itself.
 *
 */
class OOFEM_EXPORT Pair
{
private:
    /// Key.
    int key;
    /// Associate value.
    double value;
    /// Pointer to the next Pair.
    Pair *next;

public:
    /// Constructor - creates the new Pair with given key k and value v.
    Pair(int k, double v) : key(k), value(v), next(NULL) { }
    /// Destructor
    ~Pair() { }

    /// Appends a given pair to itself (sets the pointer to next pair to given Pair).
    void append(Pair *p) { next = p; }
    /// Returns the receiver key.
    int giveKey() { return key; }
    /// Returns pointer to the next pair.
    Pair *giveNext() { return next; }
    /// Returns associated value.
    double &giveValue() { return value; }
    /// Prints receiver to screen.
    void printYourself() { printf("   Pair (%d,%f)\n", key, value); }
};
} // end namespace oofem
#endif // pair_h
