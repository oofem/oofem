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

#ifndef piecewiselinfunction_h
#define piecewiselinfunction_h

#include "floatarray.h"
#include "loadtimefunction.h"

///@name Input fields for PiecewiseLinFunction
//@{
#define _IFT_PiecewiseLinFunction_npoints "npoints"
#define _IFT_PiecewiseLinFunction_t "t"
#define _IFT_PiecewiseLinFunction_ft "f(t)"
#define _IFT_PiecewiseLinFunction_dataFile "datafile"
//@}

namespace oofem {

/**
 * This class implements a piecewise linear function.
 * The function is defined by 'numberOfPoints' points. 'dates' and 'values'
 * store respectively the abscissas (t) and the values (f(t)) of the points
 */
class PiecewiseLinFunction : public LoadTimeFunction
{
protected:
    FloatArray dates;
    FloatArray values;

public:
    PiecewiseLinFunction(int i, Domain *d);
    virtual ~PiecewiseLinFunction() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual classType giveClassID() const { return PiecewiceClass; }
    virtual const char *giveClassName() const { return "PiecewiceClass"; }
    virtual const char *giveInputRecordName() const { return "PiecewiseLinFunction"; }

    virtual double __at(double);
    virtual double __derAt(double);
    virtual double __accelAt(double) { return 0.; }
};
} // end namespace oofem
#endif // piecewiselinfunction_h
