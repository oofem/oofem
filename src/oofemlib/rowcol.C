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

#include "rowcol.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"
#include "freestor.h"
#include "error.h"

namespace oofem {
RowColumn :: RowColumn(int n, int st)
// Constructor. Creates a row-column with number n, starting at index st.
{
    int size;

    number = n;
    start  = st;
    diag   = 0.;
    size   = number - start;
    if ( size ) {
        row    = allocDouble(size);
        column = allocDouble(size);
    } else {
        row    = NULL;
        column = NULL;
    }
}


RowColumn :: ~RowColumn()
// Destructor.
{
    int size;

    size = number - start;
    if ( size ) {
        freeDouble(row);
        freeDouble(column);
    }
}


#ifdef DEBUG
double &RowColumn :: atU(int i)
// Returns the i-th coefficient of the column of the receiver. Slow but
// safe.
{
    this->checkBounds(i);
    return column [ i - start ];
}
#endif


#ifdef DEBUG
double &RowColumn :: atL(int i)
// Returns the i-th coefficient of the row of the receiver. Slow but safe.
{
    this->checkBounds(i);
    return row [ i - start ];
}
#endif


void RowColumn :: checkBounds(int i)
// Checks that the receiver (k) possesses coefficients (i,k) and (k,i).
{
    if ( i < start || i >= number ) {
        OOFEM_ERROR4("RowColumn::checkBounds : error in bounds of RowColumn %d (start=%d) : no entry %d", number, start, i);
    }
}


void RowColumn :: checkSizeTowards(const IntArray &loc)
// If loc points to a coefficient which out of the receiver's range,
// enlarges the receiver.
{
    int n, i, j, c, first = start;

    n = loc.giveSize();
    if ( !n ) {
        OOFEM_ERROR("RowColumn::checkSizeTowards : error in RowCol::checkBounds : 0-size location array");
    }

    // get first non-zero coefficient in loc
    for ( i = 1; i <= n; i++ ) {
        if ( ( first = loc.at(i) ) ) {
            break;
        }
    }

    // get min coefficient in loc
    for ( j = i + 1; j <= n; j++ ) {
        if ( ( c = loc.at(j) ) ) {
            first = min(c, first);
        }
    }

    // enlarge the receiver if necessary
    if ( first < start && first ) {
        this->growTo(first);
    }
}


double RowColumn :: dot(const FloatArray &b, char c, int first, int last)
// Returns the dot product of the row (c=='R') or the column (c=='C') of
// the receiver with array 'b', from indices 'first' to 'last'.
{
    double answer, *p1;
    // double *p2 ;
    int i;

    answer = 0.;
    i      = last - first + 1;
    if ( i > 0 ) {
        if ( c == 'R' ) {
            p1 = row + first - start;
        } else {
            p1 = column + first - start;
        }

        //p2 = b->givePointer() + first-1 ;            // bad practice !
        int pos = first;
        while ( i-- ) {
            //answer += *p1++ * *p2++ ;
            answer += * p1++ * b.at(pos++);
        }
    }

    return answer;
}


void
RowColumn :: growTo(int newStart)
// Enlarges the receiver (the column upwards, the row leftwards) to
// 'start' = newStart.
{
    int newSize, size, i;
    double *newRow, *newColumn, *p1, *p2;

    if ( newStart == start ) {
        return;
    }

#  ifdef DEBUG
    if ( newStart <= 0 || newStart > start ) {
        OOFEM_ERROR4("RowColumn::growTo : cannot enlarge RowCol %d (start=%d) to %d", number, start, newStart);
    }

#  endif

    newSize   = number - newStart;
    newRow    = allocDouble(newSize);
    newColumn = allocDouble(newSize);
    size      = number - start;

    if ( size ) {
        p1  = row;
        p2  = newRow + ( start - newStart );
        i   = size;
        while ( i-- ) {
            * p2++ = * p1++;
        }

        freeDouble(row);

        p1  = column;
        p2  = newColumn + ( start - newStart );
        i   = size;
        while ( i-- ) {
            * p2++ = * p1++;
        }

        freeDouble(column);
    }

    row    = newRow;
    column = newColumn;
    start  = newStart;
}


void RowColumn :: printYourself()
// Prints the receiver on screen.
{
    int i;

    printf("Row-column %d : start = %d, diag = %.5f\n   col : ",
           number, start, diag);
    for ( i = 0; i < number - start; i++ ) {
        printf(" % .5f", column [ i ]);
    }

    printf("\n   row : ");
    for ( i = 0; i < number - start; i++ ) {
        printf(" % .5f", row [ i ]);
    }

    printf("\n");
}


void
RowColumn :: zero()
// Sets all coefficients of the receiver to 0.
{
    int i, size = number - start;

    diag = 0.;
    if ( size ) {
        for ( i = 0; i < size; i++ ) {
            column [ i ] = 0.0;
            row [ i ]    = 0.0;
        }
    }
}

RowColumn *
RowColumn :: GiveCopy()
{
    int i;
    RowColumn *answer;
    double *newRow, *newColumn, *p1, *p2;
    int size = number - start;

    if ( size ) {
        newRow    = allocDouble(size);
        newColumn = allocDouble(size);


        p1  = row;
        p2  = newRow;
        i   = size;
        while ( i-- ) {
            * p2++ = * p1++;
        }

        p1  = column;
        p2  = newColumn;
        i   = size;
        while ( i-- ) {
            * p2++ = * p1++;
        }
    } else {
        newRow = newColumn = NULL;
    }

    answer = new RowColumn(this->number, this->start, newRow, newColumn, this->diag);
    return answer;
}


RowColumn :: RowColumn(int inumber, int istart, double *irow, double *icol, double idiag)
{
    this->number = inumber;
    this->start  = istart;
    this->row    = irow;
    this->column = icol;
    this->diag   = idiag;
}
} // end namespace oofem
