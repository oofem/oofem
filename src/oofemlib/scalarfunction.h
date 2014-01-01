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

#include "floatarray.h"
#include "domain.h"
#include "parser.h"
#include "error.h"
#include <string>
#include <map>
#include <string.h>
#include <sstream>

namespace oofem {
/**
 *  Implementation of Scalar function. The scalar function can be defined as
 * (i)   simple double (constant) value,
 * (ii)  simple expression, that is evaluated by internal parser and that can depend on any number of variables,
 * which are defined by calling context (time finctions depend on time variable 't', variable loads can depend
 * on spatial position described by 'x', 'y', 'z' values, etc.).
 * (iii)  reference to a function, defined in input file and maintained by corresponding Domain.
 *
 * The sclar functions can replace constant variables in many places. The advantage is that they are naturally supported by
 * input readers (simple values are formatted as usual, simple expressions are enclosed in '$', references to external functions
 * are formatted as @i, where i is the number of corresponding external function).
 */
class ScalarFunction
{
    /// dataValue attribute defining the possible, supported variants
    union {
        /// constant, double value
        double dValue;
        /// simple expression (evaluated by internal parser); note: std::string not allowed in union
        char *eValue;
        /// reference to external function
        int fReference;
    } dataValue;

    /// enum value determining the dataValue type
    enum { DV_ValueType, DV_SimpleExpressionType, DV_FunctionReferenceType } dvType;

public:
    /**
     *  Constructor. Creates constant scalar function defined by given value.
     *  @param val defines the constant value
     */
    ScalarFunction(double val = 0) {
        this->dvType = DV_ValueType;
        this->setValue(val);
    }
    /**
     *  Constructor. Creates scalar funtion defined by given simple expression
     *  @param val string with simple expression
     */
    ScalarFunction(std :: string &val) {
        // initialize type to DV_ValueType to prevent clear() method to clear unitialized memory. will be set correctly by set method.
        this->dvType = DV_ValueType;
        this->setSimpleExpression(val);
    }
    /**
     * Constructor of sclar function defined using external function
     * @param val external function number
     */
    ScalarFunction(int val) {
        // initialize type to DV_ValueType to prevent clear() method to clear unitialized memory. will be set correctly by set method.
        this->dvType = DV_ValueType;
        this->setReference(val);
    }

    ~ScalarFunction() { this->clear(); }
    /**
     *  Sets receiver to be a constant scalar function defined by given value.
     *  @param val defines the constant value
     */
    void setValue(double val) {
        this->clear();
        this->dvType = DV_ValueType;
        this->dataValue.dValue = 0.0;
    }

    /**
     *  Sets receiver to be a scalar funtion defined by given simple expression
     *  @param val string with simple expression
     */
    void setSimpleExpression(std :: string &val) {
        this->clear();
        this->dvType = DV_SimpleExpressionType;
        this->dataValue.eValue  = new char [ val.size() ];
        strcpy( this->dataValue.eValue, val.c_str() );
    };

    /**
     * Sets receiver to be a scalar function defined using external function
     * @param val external function number
     */
    void setReference(int val) {
        this->clear();
        this->dvType = DV_FunctionReferenceType;
        this->dataValue.fReference = val;
    };

    /**
     * Evaluates the receiver.
     * @param valDict map defining input parameters in the form  (name, value) pairs
     * @param d domain managing external functions
     */
    double eval(std :: map< std :: string, double >valDict, Domain *d) const {
        if ( this->dvType == DV_ValueType ) {
            return this->dataValue.dValue;
        } else if ( this->dvType == DV_SimpleExpressionType ) {
            std :: ostringstream buff;
            Parser p;
            int err;
            // process valDict and call internal parser
            std :: map< std :: string, double > :: iterator it;
            for ( it = valDict.begin(); it != valDict.end(); it++ ) {
                buff << it->first << "=" << it->second << ";";
            }
            buff << this->dataValue.eValue;
            // evaluate the expression
            double value = p.eval(buff.str().c_str(), err);
            if ( err ) {
                OOFEM_ERROR2( "ScalarFunction::eval parser syntax error (expr=\"%s\")", buff.str().c_str() );
            }
            return value;
        } else if ( this->dvType == DV_FunctionReferenceType ) {
            // d->giveLoadTimeFunction(this->dataValue.fReference)->evaluate (valDict);
            return 0.0;
        }
        return 0.;
    }

    /**
     * Cleans up before changing receiver
     */
    void clear() {
        if ( ( this->dvType == DV_SimpleExpressionType ) && this->dataValue.eValue ) {
            delete[] this->dataValue.eValue;
            this->dataValue.eValue = NULL;
        }
    }
}; // end class ScalarFunction
} // end namespace OOFEM
