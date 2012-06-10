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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "heavisideltf.h"

namespace oofem {
double
HeavisideLTF :: __at(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    double relTime = time - this->origin;
    if ( relTime <= 0. ) {
        return 0.;
    }

    return value;
}


IRResultType
HeavisideLTF :: initializeFrom(InputRecord *ir)
//
// initializes according to input record
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LoadTimeFunction :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, origin, IFT_HeavisideLTF_origin, "origin"); // Macro
    IR_GIVE_FIELD(ir, value, IFT_HeavisideLTF_value, "value"); // Macro

    return IRRT_OK;
}

int
HeavisideLTF :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    LoadTimeFunction :: giveInputRecordString(str, keyword);
    sprintf(buff, " origin %e value %e", this->origin, this->value);
    str += buff;

    return 1;
}
} // end namespace oofem
