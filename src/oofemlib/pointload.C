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

#include "pointload.h"
#include "function.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_BoundaryCondition(PointLoad);

void
PointLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    double factor;
    // returns component array for elements which use direct formulae
    Load :: computeComponentArrayAt(answer, tStep, mode);

    // time distribution
    factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer.times(factor);
}

IRResultType
PointLoad :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int dummy;
    IR_GIVE_FIELD(ir, dummy, "ndofs");
    IR_GIVE_FIELD(ir, coords, _IFT_PointLoad_coords);

    int value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_PointLoad_loadtype);
    lType = ( bcType ) value;

    value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_PointLoad_cstype);
    coordSystemType = ( CoordSystType ) value;

    return Load :: initializeFrom(ir);
}


void
PointLoad :: giveInputRecord(DynamicInputRecord &input)
{
    Load :: giveInputRecord(input);
    input.setField(this->lType, _IFT_PointLoad_loadtype);
    input.setField(this->coordSystemType, _IFT_PointLoad_cstype);
    input.setField(this->coords, _IFT_PointLoad_coords);
}
} // end namespace oofem
