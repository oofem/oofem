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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef pulsefunction_h
#define pulsefunction_h

#include "function.h"

///@name Input fields for PulseFunction
//@{
#define _IFT_PulseFunction_Name "pulsefunction"
#define _IFT_PulseFunction_tsteptime "tstime"
#define _IFT_PulseFunction_tmin "tmin"
#define _IFT_PulseFunction_tmax "tmax"
#define _IFT_PulseFunction_value "value"
//@}

namespace oofem {
/**
 * This class implements a function that is 0 outside the user specified interval [t_min, t_max], where it has a specified
 * peak value. A variant allows to specify a single time and function is nonzero within a timestep containing given time.
 */
class OOFEM_EXPORT PulseFunction : public Function
{
private:
    /// mode of operation: 0 = interval [tmin,tmax], 1 = single time tstep
    int mode;
    /// Specific time when function is nonzero.
    double time;
    /// Start time of interval when function is nonzero.
    double tmin;
    /// End time of interval when function is nonzero.
    double tmax;
    /// Value of function within the interval.
    double value;

public:
    PulseFunction(int i, Domain * d) : Function(i, d)
    {
        mode = 0;
        time = 0.0;
        tmin = 0.0;
        tmax = 0.0;
        value = 0.0;
    }
    
    virtual ~PulseFunction() { }

    void initializeFrom(InputRecord &ir) override;
    const char *giveClassName() const override { return "PulseFunction"; }
    const char *giveInputRecordName() const override { return _IFT_PulseFunction_Name; }

    double evaluate(TimeStep *tStep, ValueModeType mode) override;
    double evaluateAtTime(double) override;
    double evaluateVelocityAtTime(double t) override { return 0.; }
    double evaluateAccelerationAtTime(double t) override { return 0.; }
};
} // end namespace oofem
#endif // pulsefunction_h
