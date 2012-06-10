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

#include "pointload.h"
#include "loadtime.h"
#include "flotarry.h"

namespace oofem {
void
PointLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode)
{
    double factor;
    // returns component array for elements which use direct formulae
    Load :: computeComponentArrayAt(answer, tStep, mode);

    // time distribution
    factor = this->giveLoadTimeFunction()->evaluate(tStep, mode);
    answer.times(factor);
}

IRResultType
PointLoad :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = Load :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        IR_IOERR(giveClassName(), __proc, IFT_Unknown, "", ir, result);
    }

    IR_GIVE_FIELD(ir, nDofs, IFT_PointLoad_ndofs, "ndofs"); // Macro
    IR_GIVE_FIELD(ir, coords, IFT_PointLoad_coords, "coords"); // Macro

    int value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_PointLoad_loadtype, "loadtype"); // Macro
    lType = ( bcType ) value;

    value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_PointLoad_cstype, "cstype"); // Macro
    coordSystemType = ( PL_CoordSystType ) value;

    return IRRT_OK;
}


int
PointLoad :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    Load :: giveInputRecordString(str, keyword);
    sprintf( buff, " ndofs %d loadtype %d cstype %d coords %d", this->nDofs, ( int ) this->lType,
            ( int ) this->coordSystemType, coords.giveSize() );
    str += buff;

    for ( int i = 1; i <= this->coords.giveSize(); i++ ) {
        sprintf( buff, " %e", this->coords.at(i) );
        str += buff;
    }

    return 1;
}
} // end namespace oofem
