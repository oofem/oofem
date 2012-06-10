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

#include "freestor.h"
#include "error.h"

#ifndef __MAKEDEPEND
 #include <cstdlib>
#endif

namespace oofem {
double *allocDouble(int n)
// Allocates and returns an array of n doubles. Issues an error message if
// the free store is exhausted.
{
    double *answer;

    if ( !n ) {
        OOFEM_FATAL("allocDouble : cannot allocate 0 coefficients");
    }

    answer = ( double * ) calloc( n, sizeof( double ) );
    if ( answer ) {
        return answer;
    } else {
        OOFEM_FATAL2("allocDouble : free store exhausted (tried to allocate %d doubles)",n);
    }

    return NULL;
}


int *allocInt(int n)
// Allocates and returns an array of n integers. Issues an error message
// if the free store is exhausted.
{
    int *answer;

    if ( !n ) {
        OOFEM_FATAL("allocInt: cannot allocate 0 coefficients");
    }

    answer = ( int * ) calloc( n, sizeof( int ) );
    if ( answer ) {
        return answer;
    } else {
        OOFEM_FATAL2("allocInt : free store exhausted (tried to allocate %d ints)",n);
    }

    return NULL;
}


void freeStoreError()
// This function is called whenever operator "new" is unable to allocate
// memory.
{
    OOFEM_FATAL("freeStoreError : free store exhausted");
}


void freeDouble(double *a)
// Deallocates the list of decimals 'a'.
{
    //#ifdef SUN_STATION
    //   free ((char*)a) ;        /* the Sun compiler is really lousy */
    //#else
    free(a);
    //#endif
}


void freeInt(int *a)
// Deallocates the list of integers 'a'.
{
    //#ifdef SUN_STATION
    //   free ((char*)a) ;        /* see above */
    //#else
    free(a);
    //#endif
}
} // end namespace oofem
