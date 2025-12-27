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


#include "field.h"
#include "expressionfield.h"
#include "floatarray.h"
#include "timestep.h"
#include "parser.h"

namespace oofem {
REGISTER_Field(ExpressionField);


/// Constructor.
ExpressionField :: ExpressionField(void) : Field(FieldType::FT_Unknown)
{
}

void ExpressionField :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, this->expression, "f");
    int val;
    IR_GIVE_FIELD(ir, val, "type");
    this->type = static_cast<FieldType>(val);
}


int ExpressionField :: evaluateAt(FloatArray &answer, const FloatArray &coords, ValueModeType mode, TimeStep *tStep)
{
    if (mode == VM_Total) {
        Parser p;
        int err;
        p.setVariableValue("x", 0, coords.at(1));
        p.setVariableValue("y", 0, (coords.giveSize()>1)?coords.at(2):0);
        p.setVariableValue("z", 0, (coords.giveSize()>2)?coords.at(3):0);
        p.setVariableValue("t", 0, tStep->giveIntrinsicTime());
        p.eval(this->expression.c_str(), answer, "f", err); // evaluate the expression; return value of "f" array
        return (err==0);
    } else {
        OOFEM_ERROR("Unsupported mode");
        return 1;
    }
}

int ExpressionField :: evaluateAt(FloatArray &answer, DofManager *dman, ValueModeType mode, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented");
    return 1;
}



} // end namespace oofem
