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

#include "boundarycondition.h"
#include "timestep.h"
#include "loadtimefunction.h"
#include "verbose.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(BoundaryCondition);

double BoundaryCondition :: give(Dof *dof, ValueModeType mode, TimeStep *stepN)
// Returns the value at stepN of the prescribed value of the kinematic
// unknown 'u'. Returns 0 if 'u' has no prescribed value.
{
    double factor = this->giveLoadTimeFunction()->evaluate(stepN, mode);
    int index = this->dofs.findFirstIndexOf( dof->giveDofID() );
    if ( !index ) {
        index = 1;
    }
    double prescribedValue = this->values.at(index);
    return prescribedValue * factor;
}


IRResultType
BoundaryCondition :: initializeFrom(InputRecord *ir)
// Sets up the dictionary where the receiver stores the conditions it
// imposes.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    GeneralBoundaryCondition :: initializeFrom(ir);

    if ( ir->hasField(_IFT_BoundaryCondition_values) ) {
        IR_GIVE_FIELD(ir, values, _IFT_BoundaryCondition_values);
    } else {
        double prescribedValue;
        if ( ir->hasField(_IFT_BoundaryCondition_PrescribedValue) ) {
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

    return IRRT_OK;
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
} // end namespace oofem
