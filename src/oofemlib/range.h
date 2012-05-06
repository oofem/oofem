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

#ifndef range_h
#define range_h

namespace oofem {
/**
 * Class Range is an abstraction for interval of integer numbers. It is described using its start and end values of interval
 * it represents. The interval is defined to represent all values between start and end values, including start and end values.
 * Function for testing if number is in interval is provided. Used by OutputManager
 * to efficiently maintain intervals.
 */
class Range
{
protected:
    /// Interval start value.
    int startIndx;
    /// Interval end value.
    int endIndx;

public:
    /// Constructor. Creates Range containing only given single number
    Range(int indx) { startIndx = endIndx = indx; }
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

    /// Tests if number is in range.
    bool test(int i) { return ( i >= startIndx ) && ( i <= endIndx ); }
};
} // end namespace oofem
#endif // range_h
