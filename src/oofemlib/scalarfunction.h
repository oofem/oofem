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

#ifndef scalarfunction_h
#define scalarfunction_h

#include "parser.h"
#include "error.h"

#include <string>
#include <map>

namespace oofem {
class FunctionArgument;
class Domain;
class FloatArray;

/**
 * Implementation of Scalar function. The scalar function can be defined as
 * (i)   simple double (constant) value,
 * (ii)  simple expression, that is evaluated by internal parser and that can depend on any number of variables,
 * which are defined by calling context (time finctions depend on time variable 't', variable loads can depend
 * on spatial position described by 'x', 'y', 'z' values, etc.).
 * (iii)  reference to a function, defined in input file and maintained by corresponding Domain.
 *
 * The sclar functions can replace constant variables in many places. The advantage is that they are naturally supported by
 * input readers (simple values are formatted as usual, simple expressions are enclosed in '$', references to external functions
 * are formatted as @i, where i is the number of corresponding external function).
 */
class ScalarFunction
{
    /// Constant, double value.
    double dValue;
    /// Simple expression (evaluated by internal parser)
    std :: string eValue;
    /// Reference to external function
    int fReference;

    /// Enum value determining the dataValue type.
    enum { DV_Undefined, DV_ValueType, DV_SimpleExpressionType, DV_FunctionReferenceType } dvType;

public:
    /**
     * Creates empty constant scalar function.
     */
    ScalarFunction();
    /**
     * Creates constant scalar function defined by given value.
     * @param val Defines the constant value
     */
    ScalarFunction(double val);
    /**
     * Creates scalar funtion defined by given simple expression.
     * @param val String with simple expression
     */
    ScalarFunction(std :: string & val);
    /**
     * Constructor of scalar function defined using external function.
     * @param val External function number
     */
    ScalarFunction(int val);

    ~ScalarFunction();
    /**
     * Sets receiver to be a constant scalar function defined by given value.
     * @param val Defines the constant value
     */
    void setValue(double val);

    /**
     * Sets receiver to be a scalar funtion defined by given simple expression
     * @param val String with simple expression
     */
    void setSimpleExpression(std :: string &val);

    /**
     * Sets receiver to be a scalar function defined using external function
     * @param val External function number
     */
    void setReference(int val);

    /**
     * Evaluates the receiver.
     * @param valDict map defining input parameters in the form  (name, value) pairs
     * @param d domain managing external functions
     */
    double eval(std :: map< std :: string, FunctionArgument >valDict, Domain *d) const;

    /**
     * Evaluates the receiver (convenience function that directly takes just the time as a parameter).
     * @param time Time to evaluate for.
     * @param d domain managing external functions
     */
    double eval(double time, Domain *d) const;

    ///@todo Support FloatArray inputs in the value dictionary, then implement this convenience function;
    /*
     * Evaluates the receiver (convenience function that directly takes just the time and coordinates as parameters).
     * @param time Time to evaluate for.
     * @param coords Coordinates to evaluate for.
     * @param d domain managing external functions
     */
    //double eval(double time, const FloatArray &coords, Domain *d) const;
    
    /**
     * True if receiver is defined.
     */
    bool isDefined() const;

    friend std :: ostream &operator << ( std :: ostream & out, const ScalarFunction & s );
};
} // end namespace OOFEM

#endif // scalarfunction_h
