/* $Header: /home/cvs/bp/oofem/sm/src/usrdeftempfield.C,v 1.4 2003/04/06 14:08:32 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#include "usrdeftempfield.h"
#include "timestep.h"
#ifndef __MAKEDEPEND
 #include <math.h>
#endif

namespace oofem {
void
UserDefinedTemperatureField :: computeValueAt(FloatArray &answer, TimeStep *stepN, FloatArray &coords, ValueModeType mode)
// Returns the value of the receiver at time and given position respecting the mode.
{
    int i, err;
    double result;
    FloatArray c(3);

    if ( ( mode != VM_Incremental ) && ( mode != VM_Total ) ) {
        _error2( "computeComponentArrayAt: unknown mode (%s)", __ValueModeTypeToString(mode) );
    }

    for ( i = 1; i <= coords.giveSize(); i++ ) {
        c.at(i) = coords.at(i);
    }

    answer.resize(this->size);
    char buff [ UserDefinedTemperatureField_MAX_EXPR_LENGTH + 80 ];
    for ( i = 1; i <= size; i++ ) {
        sprintf(buff, "x=%e;y=%e;z=%e;t=%e;%s", c.at(1), c.at(2), c.at(3), stepN->giveTargetTime(), ftExpression [ i - 1 ]);
        result = myParser.eval(buff, err);
        if ( err ) {
            _error("computeValueAt: parser syntax error");
        }

        answer.at(i) = result;

        if ( ( mode == VM_Incremental ) && ( !stepN->isTheFirstStep() ) ) {
            sprintf(buff, "x=%e;y=%e;z=%e;t=%e;%s", c.at(1), c.at(2), c.at(3),
                    stepN->giveTargetTime() - stepN->giveTimeIncrement(), ftExpression [ i - 1 ]);

            result = myParser.eval(buff, err);
            if ( err ) {
                _error("computeValueAt: parser syntax error");
            }

            answer.at(i) -= result;
        }
    }

    return;
}

IRResultType
UserDefinedTemperatureField :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    // const char *__keyword;
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, size, IFT_UserDefinedTemperatureField_size, "size"); // Macro
    if ( size > 3 ) {
        size = 3;
    }

    if ( size > 0 ) {
        IR_GIVE_FIELD2(ir, ftExpression [ 0 ], IFT_UserDefinedTemperatureField_t1, "t1(txyz)", UserDefinedTemperatureField_MAX_EXPR_LENGTH);
        /*
         * result = ir->giveField(ftExpression[0], UserDefinedTemperatureField_MAX_EXPR_LENGTH, IFT_UserDefinedTemperatureField_t1, "t1(txyz)");
         * if (result != IRRT_OK) IR_IOERR (giveClassName(), __proc, IFT_UserDefinedTemperatureField_t1, __keyword, ir, result);
         */
    }

    if ( size > 1 ) {
        IR_GIVE_FIELD2(ir, ftExpression [ 1 ], IFT_UserDefinedTemperatureField_t1, "t2(txyz)", UserDefinedTemperatureField_MAX_EXPR_LENGTH);
        /*
         * result = ir->giveField(ftExpression[1], UserDefinedTemperatureField_MAX_EXPR_LENGTH, IFT_UserDefinedTemperatureField_t2, "t2(txyz)");
         * if (result != IRRT_OK) IR_IOERR (giveClassName(), __proc, IFT_UserDefinedTemperatureField_t2, __keyword, ir, result);
         */
    }

    if ( size > 2 ) {
        IR_GIVE_FIELD2(ir, ftExpression [ 2 ], IFT_UserDefinedTemperatureField_t1, "t3(txyz)", UserDefinedTemperatureField_MAX_EXPR_LENGTH);
        /*
         * result = ir->giveField(ftExpression[2], UserDefinedTemperatureField_MAX_EXPR_LENGTH, IFT_UserDefinedTemperatureField_t3, "t3(txyz)");
         * if (result != IRRT_OK) IR_IOERR (giveClassName(), __proc, IFT_UserDefinedTemperatureField_t3, __keyword, ir, result);
         */
    }

    return IRRT_OK;
}
} // end namespace oofem
