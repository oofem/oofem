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

#include "intvarfield.h"
#include "classfactory.h"
#include "dofmanager.h"
#include "materialmappingalgorithm.h"

namespace oofem {
InternalVariableField :: InternalVariableField(InternalStateType ist, FieldType ft, MaterialMappingAlgorithmType mma_type, Domain *d) :
    Field(ft)
{
    this->type = ist;
    this->mma = classFactory.createMaterialMappingAlgorithm(mma_type);
    this->domain = d;
}

InternalVariableField :: ~InternalVariableField()
{
    if ( mma ) {
        delete mma;
    }
}

int
InternalVariableField :: evaluateAt(FloatArray &answer, FloatArray &coords, ValueModeType mode, TimeStep *tStep)
{
    IntArray types(1);
    types.at(1) = this->type;
    /// Use MaterialMappingAlgorithm classes to do the job
    Set eset(0, domain);
    eset.addAllElements();
    this->mma->__init(domain, types, coords, eset, tStep);
    this->mma->__mapVariable(answer, coords, this->type, tStep);

    return 0; // ok
}

int
InternalVariableField :: evaluateAt(FloatArray &answer, DofManager *dman, ValueModeType mode, TimeStep *tStep)
{
    if ( dman->hasCoordinates() ) {
        return this->evaluateAt(answer, * ( dman->giveCoordinates() ), mode, tStep);
    } else {
        return 1; // failed -> dman without coordinates
    }
}

contextIOResultType
InternalVariableField :: saveContext(DataStream &stream, ContextMode mode)
{
    return CIO_OK;
}

contextIOResultType
InternalVariableField :: restoreContext(DataStream &stream, ContextMode mode)
{
    return CIO_OK;
}
} // end namespace oofem
