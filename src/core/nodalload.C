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

#include "nodalload.h"
#include "classfactory.h"
#include "datastream.h"
#include "dynamicinputrecord.h"
#include "contextioerr.h"

namespace oofem {
REGISTER_BoundaryCondition(NodalLoad);

void
NodalLoad :: initializeFrom(InputRecord &ir)
{
    Load :: initializeFrom(ir);
    int value = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_NodalLoad_cstype);
    coordSystemType = ( CoordSystType ) value;
}


void NodalLoad :: giveInputRecord(DynamicInputRecord &input)
{
    Load :: giveInputRecord(input);
    input.setField(this->coordSystemType, _IFT_NodalLoad_cstype);
}


void
NodalLoad :: saveContext(DataStream &stream, ContextMode mode)
{
    Load :: saveContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.write(coordSystemType) ) {
          THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
NodalLoad :: restoreContext(DataStream &stream, ContextMode mode)
{
    Load :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        int _val;
        if ( !stream.read(_val) ) {
          THROW_CIOERR(CIO_IOERR);
        }
        coordSystemType = (CoordSystType) _val;
    }
}


} // end namespace oofem
