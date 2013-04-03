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

#include "userdefinedloadtimefunction.h"
#include "parser.h"
#include <sstream>

namespace oofem {
UserDefinedLoadTimeFunction :: UserDefinedLoadTimeFunction(int n, Domain *d) : LoadTimeFunction(n, d) { }

IRResultType
UserDefinedLoadTimeFunction :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    IR_GIVE_FIELD(ir, ftExpression, _IFT_UserDefinedLoadTimeFunction_ft);
    IR_GIVE_OPTIONAL_FIELD(ir, dfdtExpression, _IFT_UserDefinedLoadTimeFunction_dfdt);
    IR_GIVE_OPTIONAL_FIELD(ir, d2fdt2Expression, _IFT_UserDefinedLoadTimeFunction_d2fdt2);

    return LoadTimeFunction :: initializeFrom(ir);
}
    
double UserDefinedLoadTimeFunction :: __at(double time)
{
    Parser myParser;
    int err;
    double result;

    std::ostringstream buff;
    buff << "t=" << time << ";" << ftExpression;
    result = myParser.eval(buff.str().c_str(), err);
    if ( err ) {
        _error("at: parser syntax error");
    }

    return result;
}

double UserDefinedLoadTimeFunction :: __derAt(double time)
{
    Parser myParser;
    int err;
    double result;

    if ( dfdtExpression.size() == 0 ) {
        _error("derAt: derivative not provided");
        return 0.;
    }

    std::ostringstream buff;
    buff << "t=" << time << ";" << dfdtExpression;
    result = myParser.eval(buff.str().c_str(), err);
    if ( err ) {
        _error("derAt: parser syntax error");
    }

    return result;
}


double UserDefinedLoadTimeFunction :: __accelAt(double time)
{
    Parser myParser;
    int err;
    double result;

    if ( d2fdt2Expression.size() == 0 ) {
        _error("derAt: derivative not provided");
        return 0.;
    }

    std::ostringstream buff;
    buff << "t=" << time << ";" << d2fdt2Expression;
    result = myParser.eval(buff.str().c_str(), err);
    if ( err ) {
        _error("accelAt: parser syntax error");
    }

    return result;
}

} // end namespace oofem
