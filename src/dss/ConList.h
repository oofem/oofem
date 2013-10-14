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

#ifndef CONLIST_H__
#define CONLIST_H__

#include "Array.h"

DSS_NAMESPASE_BEGIN

/**
 * @author Richard Vondracek
 */

class IntLinkArray
{
public:
    long *last;
    long *next;

    long first;
    long N;

    IntLinkArray()
    {
        last = next = NULL;
        first = 0;
        N = 0;
    }

    IntLinkArray(long *l, long *n, long N)
    {
        last = l;
        next = n;
        first = 0;
        this->N = N;
    }

    void Init(long *l, long *n, long N)
    {
        last = l;
        next = n;
        first = 0;
        this->N = N;
    }

    void Remove(long i)
    {
        //if(!(last[i]!=-2)) {i =i;}

        if ( last [ i ] >= 0 ) { next [ last [ i ] ] = next [ i ]; }

        if ( next [ i ] >= 0 ) { last [ next [ i ] ] = last [ i ]; }

        if ( i == first ) {
            first = next [ i ];
        }

        N--;

        //if(!((last[i]=-2)==-2)) {i =i;}
    }
};

DSS_NAMESPASE_END

#endif //CONLIST_H__
