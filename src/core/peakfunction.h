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

#ifndef peakfunction_h
#define peakfunction_h

#include "function.h"

///@name Input fields for PeakFunction
//@{
#define _IFT_PeakFunction_Name "peakfunction"
#define _IFT_PeakFunction_t "t"
#define _IFT_PeakFunction_ft "f(t)"
//@}

namespace oofem {
/**
 * This class implements a function that is 0 everywhere, except in a single
 * point.
 */
class OOFEM_EXPORT PeakFunction : public Function
{
private:
    /// Specific time when function is nonzero.
    double t;
    /// Value of function at nonzero time.
    double value;

public:
    PeakFunction(int i, Domain * d) : Function(i, d)
    {
        t = 0.0;
        value = 0.0;
    }
    virtual ~PeakFunction() { }

    void initializeFrom(InputRecord &ir) override;
    const char *giveClassName() const override { return "PeakFunction"; }
    const char *giveInputRecordName() const override { return _IFT_PeakFunction_Name; }

    double evaluateAtTime(double) override;
    double evaluateVelocityAtTime(double t) override { return 0.; }
    double evaluateAccelerationAtTime(double t) override { return 0.; }
};
} // end namespace oofem
#endif // peakfunction_h
