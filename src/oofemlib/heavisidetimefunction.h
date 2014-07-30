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

#ifndef heavisidetimefunction_h
#define heavisidetimefunction_h

#include "floatarray.h"
#include "function.h"

///@name Input fields for HeavisideTimeFunction
//@{
#define _IFT_HeavisideTimeFunction_Name "heavisideltf"
#define _IFT_HeavisideTimeFunction_origin "origin"
#define _IFT_HeavisideTimeFunction_value "value"
//@}

namespace oofem {
/**
 * This class implements a Heaviside step load time function.
 *
 * The function is defined by the origin of step and value.
 * The result is value*H(t-origin),
 * where
 * @f[
 * H(t) = \begin{cases} 0,& t\leq 0  \\ 1, & t>0 \end{cases}
 * @f]
 */
class OOFEM_EXPORT HeavisideTimeFunction : public Function
{
private:
    double origin, value;

public:
    HeavisideTimeFunction(int i, Domain * d) : Function(i, d)
    {
        origin = value = 0.;
    }
    virtual ~HeavisideTimeFunction() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveClassName() const { return "HeavisideTimeFunction"; }
    virtual const char *giveInputRecordName() const { return _IFT_HeavisideTimeFunction_Name; }

    virtual double evaluateAtTime(double);
    virtual double evaluateVelocityAtTime(double t) { return 0.; }
    virtual double evaluateAccelerationAtTime(double t) { return 0.; }
};
} // end namespace oofem
#endif // heavisidetimefunction_h
