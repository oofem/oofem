/* $Header: /home/cvs/bp/oofem/oofemlib/src/column.h,v 1.5 2003/04/06 14:08:23 bp Exp $ */
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

//   ********************
//   *** CLASS COLUMN ***
//   ********************


#ifndef column_h
#define column_h

#include "flotarry.h"
class IntArray;
class FloatMatrix;
class Skyline;

namespace oofem {

class Column : public FloatArray
{
    /*
     * This class implements a column in a matrix stored in segmented form
     * (symmetric skyline). A column is a particular kind of FloatArray.
     * DESCRIPTION :
     * A column n stores in 'values' its 'size' coefficients, upwards :
     * .values[0]      = diagonal coefficient A(n,n)
     * .values[1]      = off-diagonal coefficient A(n-1,n)
     * .values[size-1] = highest non-0 coefficient of the n-th column.
     * TASKS :
     * Those inherited from FloatArray.
     */

private:
    int number;
    Skyline *matrix;

public:
    Column(int n, int size, Skyline *m) : FloatArray(size)
    { number = n;
      matrix = m; }
    ~Column() { }
#ifdef DEBUG
    double &at(int i);
    double   at(int i) const;
#else
    double &at(int i)                      { return values [ number - i ]; }
    double   at(int i)  const { return values [ number - i ]; }
#endif

    void     checkSizeTowards(IntArray *, int);
    double   dot(Column *, int, int);
    Column *GiveCopy();
    int      giveHighestRow()               { return number - size + 1; }

protected:
    void     checkSizeTowards(const IntArray &i) { FloatArray :: checkSizeTowards(i); }
    void     resize(int i) { FloatArray :: resize(i); }
};

} // end namespace oofem
#endif // column_h
