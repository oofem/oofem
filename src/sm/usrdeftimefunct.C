/* $Header: /home/cvs/bp/oofem/sm/src/usrdeftimefunct.C,v 1.3 2003/04/06 14:08:32 bp Exp $ */
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

#include "usrdeftimefunct.h"
#include <math.h>

double UserDefinedLoadTimeFunction :: __at(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    int err;
    double result;

    char buff [ UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH + 20 ];
    sprintf(buff, "t=%e;%s", time, ftExpression);
    result = myParser.eval(buff, err);
    if ( err ) {
        _error("at: parser syntax error");
    }

    return result;
}

double UserDefinedLoadTimeFunction :: __derAt(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    int err;
    double result;

    if ( dfdtExpression [ 0 ] == '\0' ) {
        _error("derAt: derivative not provided");
        return 0.;
    }

    char buff [ UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH + 20 ];
    sprintf(buff, "t=%e;%s", time, dfdtExpression);
    result = myParser.eval(buff, err);
    if ( err ) {
        _error("derAt: parser syntax error");
    }

    return result;
}


double UserDefinedLoadTimeFunction :: __accelAt(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    int err;
    double result;

    if ( d2fdt2Expression [ 0 ] == '\0' ) {
        _error("derAt: derivative not provided");
        return 0.;
    }

    char buff [ UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH + 20 ];
    sprintf(buff, "t=%e;%s", time, d2fdt2Expression);
    result = myParser.eval(buff, err);
    if ( err ) {
        _error("accelAt: parser syntax error");
    }

    return result;
}


IRResultType
UserDefinedLoadTimeFunction :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LoadTimeFunction::initializeFrom(ir);

    result = ir->giveField(ftExpression, UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH, IFT_UserDefinedLoadTimeFunction_ft, "f(t)");
    if ( result != IRRT_OK ) {
        IR_IOERR(giveClassName(), __proc, IFT_UserDefinedLoadTimeFunction_ft, "f(t)", ir, result);
    }

    result = ir->giveOptionalField(dfdtExpression, UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH, IFT_UserDefinedLoadTimeFunction_ft, "dfdt(t)");
    result = ir->giveOptionalField(d2fdt2Expression, UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH, IFT_UserDefinedLoadTimeFunction_ft, "d2fdt2(t)");

    //this->readQuotedString (initString, "f(t)", ftExpression, UserDefinedLoadTimeFunction_MAX_EXPR_LENGTH);

    return IRRT_OK;
}
