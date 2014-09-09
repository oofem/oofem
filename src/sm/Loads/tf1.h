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

#ifndef tf1_h
#define tf1_h

#include "structtemperatureload.h"

#define _IFT_TF1_Name "tf1"

namespace oofem {
/**
 * Class representing user defined temperature field.
 * No user input. The expression is hard - coded in the class body
 * as a function of global x,y and z coordinates and time t.
 *
 * The load time function is not used here, the function provided is
 * supposed to be function of time and coordinates.
 */
class TF1 : public StructuralTemperatureLoad
{
public:
    /**
     * Constructor. Creates temperature load function with given number, belonging to given domain.
     * @param n Load time function number.
     * @param d Domain to which new object will belongs.
     */
    TF1(int n, Domain * d) : StructuralTemperatureLoad(n, d) { }
    /// Destructor
    virtual ~TF1()  { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);

    virtual const char *giveInputRecordName() const { return _IFT_TF1_Name; }
    virtual const char *giveClassName() const { return "TF1"; }
};
} // end namespace oofem
#endif // tf1_h
