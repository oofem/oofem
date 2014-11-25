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

#include "inputrecord.h"

namespace oofem {
InputRecord :: InputRecord() { }

InputRecord :: InputRecord(const InputRecord &src) {  }

InputRecord &
InputRecord :: operator = ( const InputRecord & src )
{
    return * this;
}

const char *
InputRecord :: strerror(IRResultType rt)
{
    switch ( rt ) {
    case IRRT_NOTFOUND:
        return "Missing Keyword"; // string literal is statically allocated, return safe

    case IRRT_BAD_FORMAT:
        return "Bad format";

    default:
        return "Unknown error";
    }
}

IRResultType
InputRecord :: giveOptionalField(int &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(double &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(bool &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(std :: string &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(FloatArray &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(IntArray &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(FloatMatrix &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(std :: vector< std :: string > &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(Dictionary &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(std :: list< Range > &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}

IRResultType
InputRecord :: giveOptionalField(ScalarFunction &answer, InputFieldType id)
{
    IRResultType r = this->giveField(answer, id);
    if ( r == IRRT_NOTFOUND ) {
        return IRRT_OK;
    } else {
        return r;
    }
}
} // end namespace oofem
