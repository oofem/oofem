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

#include "intvarfield.h"
#include "classfactory.h"
#include "dofmanager.h"
#include "materialmappingalgorithm.h"

namespace oofem {
InternalVariableField :: InternalVariableField(InternalStateType ist, FieldType ft, MaterialMappingAlgorithmType mma_type, Domain *d) :
    Field(ft),
    mma(classFactory.createMaterialMappingAlgorithm(mma_type)),
    type(ist),
    domain(d)
{}

int
InternalVariableField :: evaluateAt(FloatArray &answer, const FloatArray &coords, ValueModeType mode, TimeStep *tStep)
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
    return this->evaluateAt(answer, dman->giveCoordinates(), mode, tStep);
}

void
InternalVariableField :: saveContext(DataStream &stream)
{}

void
InternalVariableField :: restoreContext(DataStream &stream)
{}
} // end namespace oofem
