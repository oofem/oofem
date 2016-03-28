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

#include "foreigntempfield.h"
#include "timestep.h"
#include "classfactory.h"
#include "field.h"

namespace oofem {
REGISTER_BoundaryCondition(ForeignTemperatureField);

void
ForeignTemperatureField :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
// Returns the value of the receiver at time and given position respecting the mode.
{
    if(!foreignField){
        OOFEM_ERROR("ForeignTemperatureField: foreignField must be assigned from python (is NULL).")  }

    if ( ( mode != VM_Incremental ) && ( mode != VM_Total ) ) {
        OOFEM_ERROR("unknown mode (%s)", __ValueModeTypeToString(mode) );
    }

	 // FIXME: copy coords, since evaluateAt need non-const coords
	 // should be changed in Field (and all derived classes)
	 FloatArray coords2(coords);

    if(foreignField->evaluateAt(answer, coords2, mode, tStep)){
        OOFEM_ERROR("ForeignTemperatureField::foreignField.evaluateAt(...) failed.");
    }
}

IRResultType
ForeignTemperatureField :: initializeFrom(InputRecord *ir)
{
    return IRRT_OK;
}
} // end namespace oofem
