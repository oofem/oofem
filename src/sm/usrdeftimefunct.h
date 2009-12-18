/* $Header: /home/cvs/bp/oofem/sm/src/usrdeftimefunct.h,v 1.5 2003/04/06 14:08:32 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *********************************************
//   *** CLASS USER DEFINED LOAD-TIME FUNCTION ***
//   *********************************************


#ifndef usrdeftimefunct_h
#define usrdeftimefunct_h

#include "loadtime.h"
#include "domain.h"
#include "parser.h"

namespace oofem {

#define UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH 200

/**
 * Class representing usr defined load time function. User input is function expression.
 * Uses Parser class to parse given expresion. Slow but usefull.
 * Load time function typically belongs to domain and is
 * attribute of one or more loads. Generally load time function is real function of time (\f$y=f(t)\f$).
 */
class UserDefinedLoadTimeFunction : public LoadTimeFunction
{
private:
    Parser myParser;
    char ftExpression [ UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH ];
    char dfdtExpression [ UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH ]; // first time derivative
    char d2fdt2Expression [ UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH ]; // second time derivative

public:

    /**
     * Constructor. Creates load time function with given number, belonging to given domain.
     * @param n load time function number
     * @param d domain to which new object will belongs.
     */
    UserDefinedLoadTimeFunction(int i, Domain *d) : LoadTimeFunction(i, d), myParser()
    { dfdtExpression [ 0 ] = d2fdt2Expression [ 0 ] = '\0'; }
    /// Destructor
    virtual ~UserDefinedLoadTimeFunction()  { }

    /// Returns classType id of receiver.
    classType   giveClassID() const { return UserDefinedLoadTimeFunctionClass; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "UserDefinedLoadTimeFunction"; }
    /**
     * Initializes receiver acording to object description stored in input record.
     * Must be implemented in derived classes
     */
    IRResultType initializeFrom(InputRecord *ir);

    /**
     * Returns the value of load time function at given time. Abstract service.
     * Must be implemented by derived classes.
     * @param t time
     * @return load time function value
     */
    virtual double     __at(double);
    /**
     * Returns the first time derivative of load time function at given time. Abstract service.
     * Must be implemented by derived classes.
     * @param t time
     * @return load time function value
     */
    virtual double    __derAt(double);
    /**
     * Returns the second time derivative of load time function at given time. Abstract service.
     * Must be implemented by derived classes.
     * @param t time
     * @return load time function value
     */
    virtual double    __accelAt(double);
};

} // end namespace oofem
#endif // usrdeftimefunct_h
