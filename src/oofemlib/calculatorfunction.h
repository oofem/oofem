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

#ifndef calculatorfunction_h
#define calculatorfunction_h

#include "function.h"

///@name Input fields for CalculatorFunction
//@{
#define _IFT_CalculatorFunction_Name "usrdefltf"
///@todo These aren't limited to just functions of "t" anymore; rename these.
#define _IFT_CalculatorFunction_f "f(t)"
#define _IFT_CalculatorFunction_dfdt "dfdt(t)"
#define _IFT_CalculatorFunction_d2fdt2 "d2fdt2(t)"
//@}

namespace oofem {
/**
 * Class representing user defined load time function. User input is function expression.
 * Uses Parser class to parse given expression. Slow but useful.
 * Load time function typically belongs to domain and is
 * attribute of one or more loads. Generally load time function is real function of time (@f$ y=f(t) @f$).
 */
class OOFEM_EXPORT CalculatorFunction : public Function
{
private:
    /// Expression for the function value.
    std :: string fExpression;
    /// Expression for first time derivative.
    std :: string dfdtExpression;
    /// Expression for second time derivative.
    std :: string d2fdt2Expression;

public:
    /**
     * Constructor. Creates load time function with given number, belonging to given domain.
     * @param n Load time function number.
     * @param d Domain to which new object will belongs..
     */
    CalculatorFunction(int n, Domain * d);
    /// Destructor.
    virtual ~CalculatorFunction() { }

    /**
     * Reads the fields
     * - f(t) (required)
     * - dfdt(t) (optional)
     * - d2fdt2(t) (optional)
     */
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &ir) override;

    void evaluate(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict, GaussPoint *gp=nullptr, double param=0.) override;
    double evaluateAtTime(double t) override;
    double evaluateVelocityAtTime(double t) override;
    double evaluateAccelerationAtTime(double t) override;

    const char *giveClassName() const override { return "CalculatorFunction"; }
    const char *giveInputRecordName() const override { return _IFT_CalculatorFunction_Name; }
};
} // end namespace oofem
#endif // calculatorfunction_h
