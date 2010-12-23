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
 *               Copyright (C) 1993 - 2010   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "constant.h"

namespace oofem {
double ConstantFunction :: giveValue()
// Returns the constant value of the receiver. Reads 'value' in the data
// file if it hasn't been done yet.
{
    return value;
}

IRResultType
ConstantFunction :: initializeFrom(InputRecord *ir)
{
    //
    // instanciates receiver according to input record
    //
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LoadTimeFunction :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, value, IFT_LoadTimeFunction_ft, "f(t)"); // Macro

    return IRRT_OK;
}


int
ConstantFunction :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    LoadTimeFunction :: giveInputRecordString(str, keyword);
    sprintf(buff, " f(t) %e", this->value);
    str += buff;

    return 1;
}
} // end namespace oofem
