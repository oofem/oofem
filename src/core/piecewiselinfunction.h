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

#ifndef piecewiselinfunction_h
#define piecewiselinfunction_h

#include "floatarray.h"
#include "function.h"

///@name Input fields for PiecewiseLinFunction
//@{
#define _IFT_PiecewiseLinFunction_Name "piecewiselinfunction"
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
class OOFEM_EXPORT PiecewiseLinFunction : public Function
{
protected:
    FloatArray dates;
    FloatArray values;

public:
    PiecewiseLinFunction(int i, Domain * d);
    virtual ~PiecewiseLinFunction() { }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveClassName() const override { return "PiecewiseLinFunction"; }
    const char *giveInputRecordName() const override { return _IFT_PiecewiseLinFunction_Name; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double evaluateAtTime(double t) override;
    double evaluateVelocityAtTime(double t) override;
    double evaluateAccelerationAtTime(double t) override { return 0.; }
};
} // end namespace oofem
#endif // piecewiselinfunction_h
