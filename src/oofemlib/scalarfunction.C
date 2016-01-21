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

#include "scalarfunction.h"
#include "floatarray.h"
#include "domain.h"
#include "function.h"
#include "parser.h"
#include "error.h"

#include <map>
#include <string>
#include <sstream>

namespace oofem {

ScalarFunction :: ScalarFunction()
{
    this->dvType = DV_Undefined;
}

ScalarFunction :: ScalarFunction(double val)
{
        this->dvType = DV_ValueType;
        this->setValue(val);
}


ScalarFunction :: ScalarFunction(std :: string &val)
{
    // initialize type to DV_ValueType to prevent clear() method to clear unitialized memory. will be set correctly by set method.
    this->dvType = DV_ValueType;
    this->setSimpleExpression(val);
}

ScalarFunction :: ScalarFunction(int val)
{
    // initialize type to DV_ValueType to prevent clear() method to clear unitialized memory. will be set correctly by set method.
    this->dvType = DV_ValueType;
    this->setReference(val);
}

ScalarFunction :: ~ScalarFunction() { }


void
ScalarFunction :: setValue(double val)
{
    this->dvType = DV_ValueType;
    this->dValue = val;
}


void
ScalarFunction :: setSimpleExpression(std :: string &val)
{
    this->dvType = DV_SimpleExpressionType;
    this->eValue = val;
}


void
ScalarFunction :: setReference(int val)
{
    this->dvType = DV_FunctionReferenceType;
    this->fReference = val;
}


double
ScalarFunction :: eval(std :: map< std :: string, FunctionArgument >valDict, Domain *d) const
{
    if ( this->dvType == DV_ValueType ) {
        return this->dValue;
    } else if ( this->dvType == DV_SimpleExpressionType ) {
        std :: ostringstream buff;
        Parser p;
        int err;
        // process valDict and call internal parser
        for ( const auto &named_arg : valDict ) {
            const FunctionArgument &arg = named_arg.second;
            if ( arg.type == FunctionArgument :: FAT_double ) {
                buff << named_arg.first << "=" << arg.val0 << ";";
            } else if ( arg.type == FunctionArgument :: FAT_FloatArray ) {
                for ( int i = 1; i <= arg.val1.giveSize(); ++i ) {
                    buff << named_arg.first << i << "=" << arg.val1.at(i) << ";";
                }
            } else if ( arg.type == FunctionArgument :: FAT_int ) {
                buff << named_arg.first << "=" << arg.val2 << ";";
            } else if ( arg.type == FunctionArgument :: FAT_IntArray ) {
                for ( int i = 1; i <= arg.val3.giveSize(); ++i ) {
                    buff << named_arg.first << i << "=" << arg.val3.at(i) << ";";
                }
            }
        }
        buff << this->eValue;
        //printf("string is %s\n", buff.str().c_str());
        // evaluate the expression
        double value = p.eval(buff.str().c_str(), err);
        if ( err ) {
            OOFEM_ERROR( "parser syntax error (expr=\"%s\")", buff.str().c_str() );
        }
        return value;
    } else if ( this->dvType == DV_FunctionReferenceType ) {
        FloatArray val;
        d->giveFunction(this->fReference)->evaluate(val, valDict);
        if ( val.giveSize() != 1 ) {
            OOFEM_ERROR( "Function @%d did not return a scalar (size = %d)", this->fReference, val.giveSize() );
        }
        return val.at(1);
    }
    return 0.;
}


double
ScalarFunction :: eval(double time, Domain *d) const
{
    if ( this->dvType == DV_ValueType ) {
        return this->dValue;
    } else {
        std :: map< std :: string, FunctionArgument >valDict;
        valDict.insert( std :: make_pair("t", time) );
        return this->eval(valDict, d);
    }
}

bool
ScalarFunction :: isDefined() const
{
    if ( this->dvType == DV_Undefined ){
        return 0;
    }
    return 1;
}

std :: ostream &operator << ( std :: ostream & out, const ScalarFunction & s )
{
    if ( s.dvType == ScalarFunction :: DV_ValueType ) {
        out << s.dValue;
    } else if ( s.dvType == ScalarFunction :: DV_SimpleExpressionType ) {
        out << '$' << s.eValue << '$';
    } else {
        out << '@' << s.fReference;
    }
    return out;
}
} // end namespace OOFEM
