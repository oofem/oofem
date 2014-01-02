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

#include "calculatorfunction.h"
#include "parser.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

#include <sstream>

namespace oofem {
REGISTER_Function(CalculatorFunction);

CalculatorFunction :: CalculatorFunction(int n, Domain *d) : Function(n, d) { }

IRResultType
CalculatorFunction :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    IR_GIVE_FIELD(ir, fExpression, _IFT_CalculatorFunction_f);
    IR_GIVE_OPTIONAL_FIELD(ir, dfdtExpression, _IFT_CalculatorFunction_dfdt);
    IR_GIVE_OPTIONAL_FIELD(ir, d2fdt2Expression, _IFT_CalculatorFunction_d2fdt2);

    return Function :: initializeFrom(ir);
}


void
CalculatorFunction :: giveInputRecord(DynamicInputRecord &input)
{
    Function :: giveInputRecord(input);
    input.setField(this->fExpression, _IFT_CalculatorFunction_f);
    input.setField(this->dfdtExpression, _IFT_CalculatorFunction_dfdt);
    input.setField(this->d2fdt2Expression, _IFT_CalculatorFunction_d2fdt2);
}


FloatArray 
CalculatorFunction :: evaluate(std :: map< std :: string, double > &valDict)
{
    Parser myParser;
    int err;
    FloatArray result(1);

    std :: ostringstream buff;
    //for (auto val : valDict) {
    for (std :: map< std :: string, double > :: iterator val = valDict.begin(); val != valDict.end(); ++val) {
        buff << val->first << "=" << val->second << ";" << fExpression;
    }
    result.at(1) = myParser.eval(buff.str().c_str(), err);
    if ( err ) {
        _error("at: parser syntax error");
    }

    //return {result};
    return result;
}


double CalculatorFunction :: evaluateAtTime(double time)
{
    Parser myParser;
    int err;
    double result;

    std :: ostringstream buff;
    buff << "t=" << time << ";" << fExpression;
    result = myParser.eval(buff.str().c_str(), err);
    if ( err ) {
        _error("at: parser syntax error");
    }

    return result;
}

double CalculatorFunction :: evaluateVelocityAtTime(double time)
{
    Parser myParser;
    int err;
    double result;

    if ( dfdtExpression.size() == 0 ) {
        _error("derAt: derivative not provided");
        return 0.;
    }

    std :: ostringstream buff;
    buff << "t=" << time << ";" << dfdtExpression;
    result = myParser.eval(buff.str().c_str(), err);
    if ( err ) {
        _error("derAt: parser syntax error");
    }

    return result;
}


double CalculatorFunction :: evaluateAccelerationAtTime(double time)
{
    Parser myParser;
    int err;
    double result;

    if ( d2fdt2Expression.size() == 0 ) {
        _error("derAt: derivative not provided");
        return 0.;
    }

    std :: ostringstream buff;
    buff << "t=" << time << ";" << d2fdt2Expression;
    result = myParser.eval(buff.str().c_str(), err);
    if ( err ) {
        _error("accelAt: parser syntax error");
    }

    return result;
}
} // end namespace oofem
