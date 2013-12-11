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

#ifndef _ORDERING_H__
#define _ORDERING_H__

#include "IntArrayList.h"
#include "BigMatrix.h"

DSS_NAMESPASE_BEGIN

struct IConectMatrix;

/**
 * @author Richard Vondracek
 */

class Ordering
{
public:
    enum Type {
        None = 0,
        ReverseCuthillMcKee = 1,
        CuthillMcKee = 2,
        MinimumDegree = 3,
        ApproxMinimumDegree = 4,
        ApproxMinimumDegreeIncomplete = 5,
        NestedGraphBisection = 6,
        MetisND = 7,
        ColAMD = 8,
        ApproxMinimumDegreeAA = 9,
    };


public:
    IntArrayList *perm;
    IntArrayList *order;
    Type type;
    IConectMatrix *cm;

    Ordering(IntArrayList *perm, IntArrayList *order)
    {
        cm = NULL;
        type = None;
        this->perm  = perm;
        this->order = order;
    }

    Ordering(IntArrayList *order)
    {
        cm = NULL;
        this->order = order;
        perm = new IntArrayList(order->Count);
        perm->Alloc();
        for ( long i = 0; i < perm->Count; i++ ) {
            perm->Items [ order->Items [ i ] ] = i;
        }
    }

    virtual ~Ordering()
    {
        if ( perm ) { delete perm;
                      perm = NULL; }

        if ( order ) { delete order;
                       order = NULL; }

        if ( cm ) { delete cm;
                    cm = NULL; }
    }
};

DSS_NAMESPASE_END

#endif //_ORDERING_H__
