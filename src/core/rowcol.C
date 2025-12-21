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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "rowcol.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"

#include <cstdlib>
#include <algorithm>

namespace oofem {

RowColumn :: RowColumn(int n):
    number(n),
    start(n),
    row(0),
    column(0),
    diag(0.)
{
}


RowColumn :: RowColumn(int n, int start):
    number(n),
    start(start),
    row(n-start),
    column(n-start),
    diag(0.)
{
}


#ifdef DEBUG
double &RowColumn :: atU(int i)
{
    this->checkBounds(i);
    return column [ i - start ];
}

double RowColumn :: atU(int i) const
{
    this->checkBounds(i);
    return column [ i - start ];
}

double &RowColumn :: atL(int i)
{
    this->checkBounds(i);
    return row [ i - start ];
}

double RowColumn :: atL(int i) const
{
    this->checkBounds(i);
    return row [ i - start ];
}
#endif


void RowColumn :: checkBounds(int i) const
{
    if ( i < start || i >= number ) {
        OOFEM_ERROR("error in bounds of RowColumn %d (start=%d) : no entry %d", number, start, i);
    }
}


void RowColumn :: checkSizeTowards(const IntArray &loc)
{
    int i, c, first = start;

    int n = loc.giveSize();
    if ( !n ) {
        OOFEM_ERROR("0-size location array");
    }

    // get first non-zero coefficient in loc
    for ( i = 1; i <= n; i++ ) {
        if ( ( first = loc.at(i) ) ) {
            break;
        }
    }

    // get min coefficient in loc
    for ( int j = i + 1; j <= n; j++ ) {
        if ( ( c = loc.at(j) ) ) {
            first = min(c, first);
        }
    }

    // enlarge the receiver if necessary
    if ( first < start && first ) {
        this->growTo(first);
    }
}


double RowColumn :: dot(const FloatArray &b, char c, int first, int last) const
// Returns the dot product of the row (c=='R') or the column (c=='C') of
// the receiver with array 'b', from indices 'first' to 'last'.
{
    double answer = 0.;
    int i = last - first + 1;

    const double *p1;
    if ( c == 'R' ) {
        p1 = row.data() + first - start;
    } else {
        p1 = column.data() + first - start;
    }

    // I do the bad practice, because it speed up debug runs significantly / Mikael
    const double *p2 = b.givePointer() + first-1 ;            // bad practice !
    //int pos = first;
    while ( i-- ) {
        answer += *p1++ * *p2++ ;
        //answer += * p1.at(pos) * b.at(pos++);
    }

    return answer;
}


void
RowColumn :: growTo(int newStart)
{
#  ifdef DEBUG
    if ( newStart <= 0 || newStart > start ) {
        OOFEM_ERROR("cannot enlarge RowCol %d (start=%d) to %d", number, start, newStart);
    }
#  endif
    int increase = start - newStart;
    row.resize(this->number-newStart);
    column.resize(this->number-newStart);
    if ( this->number - start > 0 ) {
        ///@todo Might be better to just store the data in the opposite order instead, so that a resize "does the right thing" directly.
        std::rotate(row.rbegin(), row.rbegin()+increase, row.rend());
        std::rotate(column.rbegin(), column.rbegin()+increase, column.rend());
    }
    start = newStart;
}


void RowColumn :: printYourself() const
{
    printf("Row-column %d : start = %d, diag = %.5f\n   col : ",
           number, start, diag);
    for ( int i = 0; i < number - start; i++ ) {
        printf(" % .5f", column [ i ]);
    }

    printf("\n   row : ");
    for ( int i = 0; i < number - start; i++ ) {
        printf(" % .5f", row [ i ]);
    }

    printf("\n");
}


void
RowColumn :: zero()
{
    std::fill(row.begin(), row.end(), 0.);
    std::fill(column.begin(), column.end(), 0.);
    diag = 0.;
}

} // end namespace oofem
