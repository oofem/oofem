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

#ifndef intvarfield_h
#define intvarfield_h

#include "field.h"
#include "materialmappingalgorithmtype.h"
#include "internalstatetype.h"

namespace oofem {
class Domain;
class MaterialMappingAlgorithm;

/**
 * Abstract class representing a field of an internal variable.
 * Field represent the spatial distribution of certain variable and is able to
 * evaluate its value at any point of interest.
 * The field is usually associated to the specific domain.
 *
 * It uses MaterialMappingAlgorithm interface to perform interpolation.
 * Note, that some classes implementing MaterialMappingAlgorithm may require
 * that elements implement corresponding interface.
 */
class OOFEM_EXPORT InternalVariableField : public Field
{
protected:
    /// Material mapping algorithm used.
    MaterialMappingAlgorithm *mma;
    /// InternalStateType.
    InternalStateType type;
    /// Source domain.
    Domain *domain;

public:
    /**
     * Constructor. Creates a internal variable field of given type associated to given domain.
     * @param ist Physical meaning of field.
     * @param b Field type.
     * @param mma_type Algorithm used to map materials.
     * @param d Domain which field belongs to.
     */
    InternalVariableField(InternalStateType ist, FieldType b, MaterialMappingAlgorithmType mma_type, Domain * d);
    virtual ~InternalVariableField();

    virtual int evaluateAt(FloatArray &answer, FloatArray &coords, ValueModeType mode, TimeStep *tStep);
    virtual int evaluateAt(FloatArray &answer, DofManager *dman, ValueModeType mode, TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode);
    virtual const char *giveClassName() const { return "InternalVariableField"; }
};
} // end namespace oofem
#endif // intvarfield_h
