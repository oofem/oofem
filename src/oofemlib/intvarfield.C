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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "intvarfield.h"
#include "usrdefsub.h"

namespace oofem {
InternalVariableField :: InternalVariableField(InternalStateType ist, FieldType ft, MaterialMappingAlgorithmType mma_type, Domain *d) :
    Field(ft)
{
    this->type = ist;
    this->mma = CreateUsrDefMaterialMappingAlgorithm(mma_type);
    this->domain = d;
}

InternalVariableField :: ~InternalVariableField()
{
    if ( mma ) {
        delete mma;
    }
}

int
InternalVariableField :: evaluateAt(FloatArray &answer, FloatArray &coords, ValueModeType mode, TimeStep *atTime)
{
    IntArray types(1);
    types.at(1) = this->type;
    /// Use MaterialMappingAlgorithm classes to do the job
    this->mma->__init(domain, types, coords, -1, atTime);
    this->mma->__mapVariable(answer, coords, this->type, atTime);

    return 0; // ok
}

int
InternalVariableField::evaluateAt(FloatArray &answer, DofManager* dman, ValueModeType mode, TimeStep *atTime)
{
    if (dman->hasCoordinates()) {
        return this->evaluateAt(answer, *(dman->giveCoordinates()), mode, atTime);
    } else {
        return 1; // failed -> dman without coordinates
    }
}

contextIOResultType
InternalVariableField :: saveContext(DataStream *stream, ContextMode mode)
{
    // int i, type_id = InternalVariableFieldClass;
    // contextIOResultType iores;
    // write class header
    // if (!stream->write(&type_id,1)) return CIO_IOERR;

    return CIO_OK;
}

contextIOResultType
InternalVariableField :: restoreContext(DataStream *stream, ContextMode mode)
{
    // int i, class_id;
    // contextIOResultType iores;
    // read class header
    // if (!stream->read(&class_id,1)) return CIO_IOERR;
    // if (class_id != InternalVariableField) return CIO_BADVERSION;

    return CIO_OK;
}
} // end namespace oofem
