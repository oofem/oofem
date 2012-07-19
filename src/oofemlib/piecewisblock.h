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

#ifndef piecewisblock_h
#define piecewisblock_h

#include "flotarry.h"
#include "loadtime.h"
#include "piecewis.h"

namespace oofem {

/**
 * This class implements a piecewise linear function. The format is in a block with two columns.
 * The block is defined by 'numberOfPoints' rows. The first column corresponds to (t), the second 
 * column contains f(t) values.
 */
class PiecewiseLinFunctionBlock : public PiecewiseLinFunction
{
public:
    PiecewiseLinFunctionBlock(int i, Domain *d) : PiecewiseLinFunction(i, d)
    { }
    virtual ~PiecewiseLinFunctionBlock() { }

    virtual IRResultType initializeFrom(InputRecord *ir) { OOFEM_ERROR("Not implemented"); return IRRT_NOTFOUND; };
    virtual IRResultType initializeFrom(InputRecord *ir, DataReader *dr);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual classType giveClassID() const { return PiecewiseLinFunctionBlockClass; }
    virtual const char *giveClassName() const { return "PiecewiseLinFunctionBlockClass"; }
    virtual const char *giveInputRecordName() const { return "PiecewiseLinFunctionBlock"; }
};
} // end namespace oofem
#endif // piecewisblock_h
