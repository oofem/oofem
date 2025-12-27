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

#include "pulsefunction.h"
#include "mathfem.h"
#include "timestep.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Function(PulseFunction);


double 
PulseFunction :: evaluate(TimeStep *tStep, ValueModeType mode)
{
    if (this->mode == 0)  {
        return this->evaluateAtTime( this->time );
    } else if (this->mode == 1) {
        if  (mode == VM_Total ) {
            double ts = tStep->giveTargetTime()-tStep->giveTimeIncrement();
            double te = tStep->giveTargetTime();
            if ((this->time >= ts) && (this->time <= te)) {
                return value;
            } else {
                return 0.;
            }
        } else {
            OOFEM_ERROR("unsupported mode(%d)", mode);
        }
    }
    return 0.0;
}

double
PulseFunction :: evaluateAtTime(double time)
// Returns the value of the receiver at time 'time'.
{
    if (this->mode == 0)  {
        if ( mode == VM_Total ) {
            if ( (time >= tmin) && (time <= tmax) ) {
                return value;
            } else {
                return 0.;
            }
        } else {
            OOFEM_ERROR("unsupported mode(%d)", mode);
        }
    } else if (this->mode == 1) {
            OOFEM_ERROR("unsupported mode(%d)", mode);
    }
    return 0.0;
}

void
PulseFunction :: initializeFrom(InputRecord &ir)
{
    Function :: initializeFrom(ir);

    if (ir.hasField(_IFT_PulseFunction_tsteptime)) {
        mode = 1;
        IR_GIVE_FIELD(ir, time, _IFT_PulseFunction_tsteptime);
    } else {
        mode = 0;
        IR_GIVE_FIELD(ir, tmin, _IFT_PulseFunction_tmin);
        IR_GIVE_FIELD(ir, tmax, _IFT_PulseFunction_tmax);
    }
    this->value = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_PulseFunction_value);
}

} // end namespace oofem
