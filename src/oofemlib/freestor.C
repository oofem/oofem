/* $Header: /home/cvs/bp/oofem/oofemlib/src/freestor.C,v 1.5.4.1 2004/04/05 15:19:43 bp Exp $ */
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

/*   file FREESTOR.C */

#include "freestor.h"
#include "compiler.h"
#include "debug.h"
#include "error.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <stdlib.h>
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
        freeStoreError();
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
        freeStoreError();
    }

    return NULL;
}


void  freeStoreError()
// This function is called whenever operator "new" is unable to allocate
// memory.
{
    OOFEM_FATAL("freeStoreError : free store exhausted");
}


void  freeDouble(double *a)
// Deallocates the list of decimals 'a'.
{
    //#ifdef SUN_STATION
    //   free ((char*)a) ;        /* the Sun compiler is really lousy */
    //#else
    free(a);
    //#endif
}


void  freeInt(int *a)
// Deallocates the list of integers 'a'.
{
    //#ifdef SUN_STATION
    //   free ((char*)a) ;        /* see above */
    //#else
    free(a);
    //#endif
}

} // end namespace oofem
