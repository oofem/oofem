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

#include "boundarycondition.h"
#include "timestep.h"
#include "function.h"
#include "verbose.h"
#include "dynamicinputrecord.h"
#include "dof.h"
#include "classfactory.h"
#include "contextioerr.h"
#include "error.h"

namespace oofem {
REGISTER_BoundaryCondition(BoundaryCondition);

double BoundaryCondition :: give(Dof *dof, ValueModeType mode, TimeStep *tStep)
{
    if ( mode == VM_Incremental ) {
        return this->give(dof, VM_Total, tStep->giveTargetTime()) - this->give(dof, VM_Total, tStep->giveTargetTime() - tStep->giveTimeIncrement());
    } else {
        return this->give(dof, mode, tStep->giveIntrinsicTime());
    }
}


double BoundaryCondition :: give(Dof *dof, ValueModeType mode, double time)
{
    double factor = 0;
    if ( (mode == VM_Total) || (mode == VM_TotalIntrinsic)) {
        factor = this->giveTimeFunction()->evaluateAtTime(time);
    } else if ( mode == VM_Velocity ) {
        factor = this->giveTimeFunction()->evaluateVelocityAtTime(time);
    } else if ( mode == VM_Acceleration ) {
        factor = this->giveTimeFunction()->evaluateAccelerationAtTime(time);
    } else {
        OOFEM_ERROR("Should not be called for value mode type then total, velocity, or acceleration.");
    }
    int index = this->dofs.findFirstIndexOf( dof->giveDofID() );
    if ( !index ) {
        index = 1;
    }
    double prescribedValue = this->values.at(index);
    return prescribedValue * factor;
}


void
BoundaryCondition :: initializeFrom(InputRecord &ir)
{
    GeneralBoundaryCondition :: initializeFrom(ir);

    if ( ir.hasField(_IFT_BoundaryCondition_values) ) {
        IR_GIVE_FIELD(ir, values, _IFT_BoundaryCondition_values);
    } else {
        double prescribedValue;
        if ( ir.hasField(_IFT_BoundaryCondition_PrescribedValue) ) {
            IR_GIVE_FIELD(ir, prescribedValue, _IFT_BoundaryCondition_PrescribedValue);
        } else {
            IR_GIVE_FIELD(ir, prescribedValue, _IFT_BoundaryCondition_PrescribedValue_d);
        }
        // Backwards compatibility with old input method:
        if ( this->dofs.giveSize() ) {
            values.resize( this->dofs.giveSize() );
        } else {
            values.resize(1);
        }
        values.zero();
        values.add(prescribedValue);
    }
}


void
BoundaryCondition :: giveInputRecord(DynamicInputRecord &input)
{
    GeneralBoundaryCondition :: giveInputRecord(input);
    input.setField(this->values, _IFT_BoundaryCondition_values);
}


void
BoundaryCondition :: setPrescribedValue(double s)
{
    values.zero();
    values.add(s);
}


void
BoundaryCondition :: scale(double s)
{
    values.times(s);
}


void
BoundaryCondition :: saveContext(DataStream &stream, ContextMode mode)
{
    GeneralBoundaryCondition :: saveContext(stream, mode);

    if ( mode & CM_Definition ) {
        contextIOResultType iores;
        if ( (iores = values.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
BoundaryCondition :: restoreContext(DataStream &stream, ContextMode mode)
{
    GeneralBoundaryCondition :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        contextIOResultType iores;
        if ( (iores = values.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }
}
} // end namespace oofem
