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

#include "heavisidetimefunction.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Function(HeavisideTimeFunction);

double
HeavisideTimeFunction :: evaluateAtTime(double time)
{
    double relTime = time - this->origin;
    if ( relTime <= 0. ) {
        return 0.;
    }

    return value;
}


IRResultType
HeavisideTimeFunction :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, origin, _IFT_HeavisideTimeFunction_origin);
    IR_GIVE_FIELD(ir, value, _IFT_HeavisideTimeFunction_value);

    return Function :: initializeFrom(ir);
}


void HeavisideTimeFunction :: giveInputRecord(DynamicInputRecord &input)
{
    Function :: giveInputRecord(input);
    input.setField(this->origin, _IFT_HeavisideTimeFunction_origin);
    input.setField(this->value, _IFT_HeavisideTimeFunction_value);
}
} // end namespace oofem
