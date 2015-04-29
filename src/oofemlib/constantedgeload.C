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

#include "constantedgeload.h"
#include "function.h"
#include "floatarray.h"
#include "timestep.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(ConstantEdgeLoad);

IRResultType
ConstantEdgeLoad :: initializeFrom(InputRecord *ir)
{
    return BoundaryLoad :: initializeFrom(ir);
}

void
ConstantEdgeLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("mode not supported");
    }

    if ( !isImposed(tStep) ) {
        answer.zero();
    } else {
        double factor = this->giveTimeFunction()->evaluate(tStep, mode);
        answer.beScaled(factor, componentArray);
    }
}
} // end namespace oofem
