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

#ifndef usrdeftempfield_h
#define usrdeftempfield_h

#include "structtemperatureload.h"
#include "parser.h"

///@name Input fields for UserDefinedTemperatureField
//@{
#define _IFT_UserDefinedTemperatureField_Name "usrdeftempfield"
#define _IFT_UserDefinedTemperatureField_size "size"
#define _IFT_UserDefinedTemperatureField_t1 "t1(txyz)"
#define _IFT_UserDefinedTemperatureField_t2 "t2(txyz)"
#define _IFT_UserDefinedTemperatureField_t3 "t3(txyz)"
//@}

namespace oofem {
/**
 * Class representing user defined temperature field. User input is function expression,
 * as a function of global x,y and z coordinates and time t.
 *
 * The load time function is not used here, the function provided is
 * supposed to be function of time and coordinates.
 *
 * Uses Parser class to parse given expression. Slow but useful.
 * Temperature load as body load is typically attribute of  domain and is
 * attribute of one or more elements.
 */
class UserDefinedTemperatureField : public StructuralTemperatureLoad
{
private:
    Parser myParser;
    int size;
    std :: string ftExpression [ 3 ];

public:
    /**
     * Constructor. Creates temperature load function with given number, belonging to given domain.
     * @param n Load time function number
     * @param d Domain to which new object will belongs.
     */
    UserDefinedTemperatureField(int n, Domain * d) : StructuralTemperatureLoad(n, d), myParser() { }
    /// Destructor
    virtual ~UserDefinedTemperatureField() { }

    /**
     * Computes components values of temperature field at given point (coordinates given in Global c.s.).
     * taking into account corresponding load time function value respecting load response mode.
     * @param answer Component values at given point and time.
     * @param tStep Time step representing time.
     * @param coords Global coordinates, which are used to evaluate components values.
     * @param mode Determines response mode.
     */
    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);

    virtual const char *giveInputRecordName() const { return _IFT_UserDefinedTemperatureField_Name; }
    virtual const char *giveClassName() const { return "UserDefinedTemperatureField"; }

    virtual IRResultType initializeFrom(InputRecord *ir);
};
} // end namespace oofem
#endif // usrdeftempfield_h
