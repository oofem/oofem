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

#ifndef _BISECTION_H__
#define _BISECTION_H__

#include "SparseConectivityMtx.h"

DSS_NAMESPASE_BEGIN

/**
 * @author Richard Vondracek
 */

class CMcKee
{
private:
    long n;

    SparseConectivityMtxII *mtx;
    long *p_node_level;
    long *p_order;     //out
    // 0 - unsorted available
    //


public:
    long *nodes;      //in
    long size;
    long domA, domB;

    CMcKee();
    ~CMcKee();

    void Init(SparseConectivityMtxII *mtx);

    bool IsAvailable(int v);

    void PrepareValid();

    long FindFirstNode();
    void ComputeLevels();
    void DivideByMidLevel();
};


class CBiSection
{
private:
    SparseConectivityMtxII *mtx;
    CMcKee mck;

public:
    CBiSection(SparseConectivityMtxII *mtx);
    void RecurBiSectOrder(IntArrayList *order);

private:
    void RecurBiSect(long *nodes, long size);
    void BiSect(long *nodes, long size, long &domA, long &domB);
};

DSS_NAMESPASE_END

#endif //_BISECTION_H__
