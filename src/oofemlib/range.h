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

#ifndef range_h
#define range_h

//#include <iosfwd>
#include <ostream>

#include "oofemcfg.h"

namespace oofem {
/**
 * Class Range is an abstraction for interval of integer numbers. It is described using its start and end values of interval
 * it represents. The interval is defined to represent all values between start and end values, including start and end values.
 * Function for testing if number is in interval is provided. Used by OutputManager
 * to efficiently maintain intervals.
 */
class OOFEM_EXPORT Range
{
protected:
    /// Interval start value.
    int startIndx;
    /// Interval end value.
    int endIndx;

public:
    /// Constructor. Creates Range containing only given single number
    Range(int indx) {
        startIndx = endIndx = indx;
    }
    /// Constructor. Creates range <li, hi>
    Range(int li, int hi) {
        startIndx = li;
        endIndx = hi;
    }
    /// Empty range constructor.
    Range() {
        startIndx = 0;
        endIndx = -1;
    }

    /// Returns the start index (inclusive).
    int giveStart() { return startIndx; }
    /// Returns the end index (inclusive).
    int giveEnd() { return endIndx; }

    /// Tests if number is in range.
    bool test(int i) { return ( i >= startIndx ) && ( i <= endIndx ); }

    friend std :: ostream &operator << ( std :: ostream & out, const Range & r ) {
        return out << r.startIndx << " " << r.endIndx;
    }
};
} // end namespace oofem
#endif // range_h
