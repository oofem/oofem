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

#include "foreigntempfieldload.h"
#include "timestep.h"
#include "classfactory.h"
#include "field.h"

namespace oofem {
REGISTER_BoundaryCondition(ForeignTemperatureFieldLoad);

void
ForeignTemperatureFieldLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
// Returns the value of the receiver at time and given position respecting the mode.
{
    if(!foreignField){
        OOFEM_ERROR("ForeignTemperatureFieldLoad: foreignField must be assigned from python (is NULL).")  }

    if ( ( mode != VM_Incremental ) && ( mode != VM_Total ) ) {
        OOFEM_ERROR("unknown mode (%s)", __ValueModeTypeToString(mode) );
    }

    if(foreignField->evaluateAt(answer, coords, mode, tStep)){
        OOFEM_ERROR("ForeignTemperatureFieldLoad::foreignField.evaluateAt(...) failed.");
    }
}

void
ForeignTemperatureFieldLoad :: initializeFrom(InputRecord &ir)
{
}
} // end namespace oofem
