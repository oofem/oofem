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

#ifndef constantfunction_h
#define constantfunction_h

#include "function.h"

#define _IFT_ConstantFunction_Name "constantfunction"
#define _IFT_ConstantFunction_f "f(t)" ///@todo Rename this to just "f"

namespace oofem {
/**
 * Class implementing time function that is constant in time; @f$ f(t) = C @f$.
 */
class OOFEM_EXPORT ConstantFunction : public Function
{
private:
    /// Value of receiver.
    double value;

public:
    /**
     * Constructor. Creates constant load time function with given number, belonging to given domain.
     * @param i Load time function number.
     * @param d Domain to which new object will belongs.
     */
    ConstantFunction(int i, Domain * d) : Function(i, d) {
        value = 0;
    }
    /// Destructor.
    virtual ~ConstantFunction() { }

    /// @return Value of receiver.
    double giveValue() { return value; }

    void evaluate(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict, GaussPoint *gp=nullptr, double param=0.) override { answer = FloatArray{this->giveValue()}; }
    double evaluateAtTime(double t) override { return this->giveValue(); }
    double evaluateVelocityAtTime(double t) override { return 0.; }
    double evaluateAccelerationAtTime(double t) override { return 0.; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    const char *giveClassName() const override { return "ConstantFunction"; }
    const char *giveInputRecordName() const override { return _IFT_ConstantFunction_Name; }
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};
} // end namespace oofem
#endif // constantfunction_h
