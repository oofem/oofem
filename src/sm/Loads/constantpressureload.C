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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "constantpressureload.h"
#include "dynamicinputrecord.h"
#include "function.h"
#include "floatarray.h"
#include "timestep.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(ConstantPressureLoad);

ConstantPressureLoad :: ConstantPressureLoad(int i, Domain *d) : SurfaceLoad(i, d)
{
    this->loadOffset = 0.0;
}

void
ConstantPressureLoad :: initializeFrom(InputRecord &ir)
{
    BoundaryLoad :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, this->loadOffset, _IFT_ConstantPressureLoad_LoadOffset);
#if 1
    if ( ir.hasField(_IFT_GeneralBoundaryCondition_dofs) ) {
        throw ValueInputException(ir, _IFT_GeneralBoundaryCondition_dofs, "Constant pressure load should not specify DOFs");
    }
#endif
}


void
ConstantPressureLoad :: giveInputRecord(DynamicInputRecord &input)
{
    BoundaryLoad :: giveInputRecord(input);
    input.setField(this->loadOffset, _IFT_ConstantPressureLoad_LoadOffset);
}

void
ConstantPressureLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    // we overload general implementation on the boundary load level due
    // to implementation efficiency

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("mode not supported");
    }

    // ask time distribution
    double factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer = componentArray;
    answer.times(factor);
}

void
ConstantPressureLoad :: computeValues(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, const IntArray &dofids, ValueModeType mode)
{
    OOFEM_ERROR("Should not be called for constant pressure loads");
}

} // end namespace oofem
