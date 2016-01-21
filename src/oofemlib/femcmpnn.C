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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#include "femcmpnn.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"

#include <cstdarg>


namespace oofem {
contextIOResultType
FEMComponent :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    if ( mode & CM_Definition ) {
        if ( !stream.write(number) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}


contextIOResultType
FEMComponent :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    if ( mode & CM_Definition ) {
        if ( !stream.read(number) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}


void
FEMComponent :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField( this->giveInputRecordName(), this->giveNumber() );
}


std :: string
FEMComponent :: errorInfo(const char *func) const
{
    return std :: string(this->giveClassName()) + "::" + func + ", number: " + std::to_string(this->giveNumber());
}

IRResultType FEMComponent :: initializeFrom(InputRecord* ir)
{
    return IRRT_OK;
}

int FEMComponent :: checkConsistency()
{
    return 1;
}

} // end namespace oofem
