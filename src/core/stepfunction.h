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

#ifndef stepfunction_h
#define stepfunction_h

#include "floatarray.h"
#include "function.h"

///@name Input fields for StepFunction
//@{
#define _IFT_StepFunction_Name "stepfunction"
#define _IFT_StepFunction_t "t"
#define _IFT_StepFunction_ft "f(t)"
#define _IFT_StepFunction_dataFile "datafile"
//@}

namespace oofem {
/**
 * This class implements a step function, which is a piecewise constant function having only finitely many pieces.
 * The function is defined by 'numberOfPoints' points. 'dates' and 'values'.
 * The function is defined as:
 * @f[
 * f(t) = \begin{cases}
 * v_1, & t < t_2 \\
 * v_i, & t_i \leq t < t_{i+1}, i= 2, \ldots, n-1 \\
 * v_n, & t \geq t_n
 * \end{cases}
 * @f]
 * where @f$ (t_i, v_i), i=1,\ldots,n @f$ are the points stored in 'dates' and 'values'.
 */

class OOFEM_EXPORT StepFunction : public Function
{
protected:
    std::vector<std::tuple<double,double>> datevalues;

public:
    StepFunction(int i, Domain * d);
    StepFunction(int i, Domain * d, const std::vector<std::tuple<double,double>> &datevalues) :
        Function(i, d), datevalues(datevalues) {}
    virtual ~StepFunction() { }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveClassName() const override { return "StepFunction"; }
    const char *giveInputRecordName() const override { return _IFT_StepFunction_Name; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double evaluateAtTime(double t) override;
    double evaluateVelocityAtTime(double t) override;
    double evaluateAccelerationAtTime(double t) override { return 0.; }
};
} // end namespace oofem
#endif // stepfunction_h
