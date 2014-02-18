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

#ifndef periodicpiecewiselinfunction_h
#define periodicpiecewiselinfunction_h

#include "floatarray.h"
#include "piecewiselinfunction.h"

///@name Input fields for PeriodicPiecewiseLinFunction
//@{
#define _IFT_PeriodicPiecewiseLinFunction_Name "periodicpiecewiselinfunction"
#define _IFT_PeriodicPiecewiseLinFunction_period "period"
#define _IFT_PeriodicPiecewiseLinFunction_addtf "addtf"
//@}

namespace oofem {
/**
 * This class implements an enhanced piecewise linear function with periodicity.
 * and possibility to add another arbitrary time function.
 *
 * The function is defined by 'numberOfPoints' points. 'dates' and 'values'
 * store respectively the abscissas (t) and the values (f(t)) of the points
 * The values are repeated after 'period'. 'AddTF' parameter specifies number
 * of function to add.
 */
class OOFEM_EXPORT PeriodicPiecewiseLinFunction : public PiecewiseLinFunction
{
private:
    /// If nonzero, the value of time function specified by addTF is added to computed value.
    int addTF;
    /**
     * If less than zero no periodicity, if >=0 date time is computed as
     * given time%period.
     * If points span more than period, span of LAST period is repeated
     */
    double period;

public:
    PeriodicPiecewiseLinFunction(int i, Domain * d) : PiecewiseLinFunction(i, d)
    {
        period = -1.0;
        addTF = 0;
    }
    virtual ~PeriodicPiecewiseLinFunction() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveClassName() const { return "PeriodicPiecewiseClass"; }
    virtual const char *giveInputRecordName() const { return _IFT_PeriodicPiecewiseLinFunction_Name; }

    virtual double evaluateAtTime(double);
    virtual double evaluateVelocityAtTime(double);
    virtual double evaluateAccelerationAtTime(double) { return 0.; }
};
} // end namespace oofem
#endif // periodicpiecewiselinfunction_h
