/* $Header: /home/cvs/bp/oofem/oofemlib/src/column.C,v 1.4.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   file COLUMN.CC

#include "column.h"
#include "skyline.h"
#include "intarray.h"
#include "flotmtrx.h"
#include "mathfem.h"
#include "freestor.h"
#include "debug.h"
#include "error.h"

namespace oofem {
#ifdef DEBUG
double &Column :: at(int i)
//  Returns the i-th coefficient of the receiver. Slow but safe.
{
    int j = number - i;
    if ( j < 0 ) {
        OOFEM_ERROR2("Column::at : column error on index : %d < 0", j);
    }

    if ( j >= this->giveSize() ) {
        OOFEM_ERROR3( "Column::at : column error on index : %d (size=%d)", j, this->giveSize() );
    }

    return values [ j ];
}
#endif



void Column :: checkSizeTowards(IntArray *loc, int matrixIsDiagonal)
// Enlarges the receiver if 'loc' points to coefficients outside of the
// receiver. This test if trivial if the matched matrix is diagonal, since
// a column should always contin at least the diagonal coefficient.
{
    int i, coeff, highRow;

    if ( matrixIsDiagonal ) {
        return;
    }

    highRow = 32000;
    for ( i = 1; i <= loc->giveSize(); i++ ) { // find the most off-diagonal
        coeff = loc->at(i);                 //       coefficient in 'loc'
        if ( coeff ) {
            highRow = min(highRow, coeff);
        }
    }

    if ( highRow < this->giveHighestRow() ) {
        this->resize(number - highRow + 1);
    }
}


double Column :: dot(Column *b, int start, int stop)
// Returns the scalar product of the receiver by 'b', from index 'start'
// (included) to index 'stop' (included).
{
    double *P1, *P2;
    int i;

#  ifdef DEBUG
    if ( start <= number - size || start <= b->number - b->size || stop > number || stop > b->number ) {
        fprintf(stderr, "Column::dot : cannot make dot product from %d to %d\nof Column %d (size %d) and Column %d (size %d)\n",
                start, stop, number, size, b->number, b->size);
    }

#  endif

    P1 = values + number - stop;             // P1 points to a(stop)
    P2 = b->values + b->number - stop;       // P2 points to b(stop)
    i  = stop - start + 1;

    return dotProduct(P1, P2, i);
}


Column *Column :: GiveCopy()
// Returns a copy of the receiver.
{
    Column *answer;
    double *p1, *p2;
    register int i;

    answer = new Column(number, size, matrix);
    p1 = answer->values;
    p2 = values;
    i  = size;
    while ( i-- ) {
        * p1++ = * p2++;
    }

    return answer;
}
} // end namespace oofem
