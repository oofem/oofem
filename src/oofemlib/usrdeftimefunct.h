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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef usrdeftimefunct_h
#define usrdeftimefunct_h

#include "loadtime.h"

namespace oofem {

/**
 * Class representing user defined load time function. User input is function expression.
 * Uses Parser class to parse given expression. Slow but useful.
 * Load time function typically belongs to domain and is
 * attribute of one or more loads. Generally load time function is real function of time (@f$y=f(t)@f$).
 */
class UserDefinedLoadTimeFunction : public LoadTimeFunction
{
private:
    /// Expression for the function value.
    std::string ftExpression;
    /// Expression for first time derivative.
    std::string dfdtExpression;
    /// Expression for second time derivative.
    std::string d2fdt2Expression;

public:
    /**
     * Constructor. Creates load time function with given number, belonging to given domain.
     * @param n Load time function number.
     * @param d Domain to which new object will belongs..
     */
    UserDefinedLoadTimeFunction(int n, Domain *d);
    /// Destructor.
    virtual ~UserDefinedLoadTimeFunction() { }
    
    /**
     * Reads the fields
     * - f(t) (required)
     * - dfdt(t) (optional)
     * - d2fdt2(t) (optional)
     */
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double __at(double t);
    virtual double __derAt(double t);
    virtual double __accelAt(double t);

    virtual classType giveClassID() const { return UserDefinedLoadTimeFunctionClass; }
    virtual const char *giveClassName() const { return "UserDefinedLoadTimeFunction"; }
};
} // end namespace oofem
#endif // usrdeftimefunct_h
