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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef piecewis_h
#define piecewis_h

#include "flotarry.h"
#include "loadtime.h"

namespace oofem {

/**
 * This class implements a piecewise linear function.
 * The function is defined by 'numberOfPoints' points. 'dates' and 'values'
 * store respectively the abscissas (t) and the values (f(t)) of the points
 */
class PiecewiseLinFunction : public LoadTimeFunction
{
protected:
    int numberOfPoints;
    FloatArray dates;
    FloatArray values;

public:
    PiecewiseLinFunction(int i, Domain *d) : LoadTimeFunction(i, d), dates(), values()
    { numberOfPoints = 0; }
    ~PiecewiseLinFunction() { }

    IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    classType giveClassID() const { return PiecewiceClass; }
    const char *giveClassName() const { return "PiecewiceClass"; }
    const char *giveInputRecordName() const { return "PiecewiseLinFunction"; }

    virtual double __at(double);
    virtual double __derAt(double);
    virtual double __accelAt(double) { return 0.; }
};
} // end namespace oofem
#endif // piecewis_h
