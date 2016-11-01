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

#ifndef foreigntempfield_h
#define foreigntempfield_h

#include <memory> // std::shared_ptr
#include "structtemperatureload.h"

///@name Input fields for ForeignTemperatureFieldLoad
//@{
#define _IFT_ForeignTemperatureFieldLoad_Name "foreigntempfieldload"
//@}

namespace oofem {

class Field;
/**
 * Class representing foreign temperature field, which asks a field object to return
 * temperature at given point. 
 *
 * The load time function is not used here, the field has to be updated by the user
 * if it is non-constant in time.
 *
 * The field (stored in foreignField) cannot be set in the input file (there
 * is only one mesh defined), it is typically assigned using Python bindings.
 *
 */
class ForeignTemperatureFieldLoad : public StructuralTemperatureLoad
{
private:
public:
    // make public so that it can be simply set from python
    std::shared_ptr<Field> foreignField;
    /**
     * Constructor. Creates temperature load function with given number, belonging to given domain.
     * @param n Load time function number
     * @param d Domain to which new object will belongs.
     */
    ForeignTemperatureFieldLoad(int n, Domain * d) : StructuralTemperatureLoad(n, d) { }
    /// Destructor
    virtual ~ForeignTemperatureFieldLoad() { }

    /**
     * Computes components values of temperature field at given point (coordinates given in Global c.s.).
     * taking into account corresponding load time function value respecting load response mode.
     * @param answer Component values at given point and time.
     * @param tStep Time step representing time.
     * @param coords Global coordinates, which are used to evaluate components values.
     * @param mode Determines response mode.
     */
    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode);

    virtual const char *giveInputRecordName() const { return _IFT_ForeignTemperatureFieldLoad_Name; }
    virtual const char *giveClassName() const { return "ForeignTemperatureFieldLoad"; }

    virtual IRResultType initializeFrom(InputRecord *ir);
};
} // end namespace oofem
#endif // foreigntempfield_h
