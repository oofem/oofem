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

#ifndef pythonexpression_h
#define pythonexpression_h

#include "function.h"

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

///@name Input fields for PythonExpression
//@{
#define _IFT_PythonExpression_Name "pythonexpression"
#define _IFT_PythonExpression_f "f" ///< Expression with return variable named "ret"
#define _IFT_PythonExpression_dfdt "dfdt" ///< Velocity with return variable named "ret"
#define _IFT_PythonExpression_d2fdt2 "d2fdt2" ///< Acceleration with return variable named "ret"
//@}


namespace oofem {
/**
 * Class representing user defined functions as Python expressions
 */
class OOFEM_EXPORT PythonExpression : public Function
{
private:
    /// Expression for the function value.
    std :: string fExpression;
    /// Expression for first time derivative.
    std :: string dfdtExpression;
    /// Expression for second time derivative.
    std :: string d2fdt2Expression;

    PyObject *f;
    PyObject *dfdt;
    PyObject *d2fdt2;

    PyObject *main_dict;

    /// Helper function to convert the std::map to a Python dictionary.
    PyObject *getDict(const std :: map< std :: string, FunctionArgument > &valDict);
    /// Helper function to run given function for given value dictionary.
    void getArray(FloatArray &answer, PyObject *func, const std :: map< std :: string, FunctionArgument > &valDict);
    /// Helper function to run given function for given time
    double getScalar(PyObject *func, double time);

public:
    /**
     * Constructor. Creates load time function with given number, belonging to given domain.
     * @param n Load time function number.
     * @param d Domain to which new object will belongs..
     */
    PythonExpression(int n, Domain * d);
    /// Destructor.
    virtual ~PythonExpression();

    /**
     * Reads the fields
     * - f(t) (required)
     * - dfdt(t) (optional)
     * - d2fdt2(t) (optional)
     */
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &ir) override;

    void evaluate(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict, GaussPoint *gp=nullptr, double param=0.) override;
    void evaluateVelocity(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict);
    void evaluateAcceleration(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict);

    double evaluateAtTime(double t) override;
    double evaluateVelocityAtTime(double t) override;
    double evaluateAccelerationAtTime(double t) override;

    const char *giveClassName() const override { return "PythonExpression"; }
    const char *giveInputRecordName() const override { return _IFT_PythonExpression_Name; }
};
} // end namespace oofem
#endif // pythonexpression_h
