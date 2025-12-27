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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "constantfunction.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
REGISTER_Function(ConstantFunction);

void
ConstantFunction :: initializeFrom(InputRecord &ir)
{
    Function :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, value, _IFT_ConstantFunction_f);
}


void
ConstantFunction :: giveInputRecord(DynamicInputRecord &input)
{
    Function :: giveInputRecord(input);
    input.setField(this->value, _IFT_ConstantFunction_f);
}

void
ConstantFunction :: saveContext(DataStream &stream, ContextMode mode)
{
    Function :: saveContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.write(value) ) {
          THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
ConstantFunction :: restoreContext(DataStream &stream, ContextMode mode)
{
    Function :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.read(value) ) {
          THROW_CIOERR(CIO_IOERR);
        }
    }
}

} // end namespace oofem
