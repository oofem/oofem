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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef intvarfield_h
#define intvarfield_h

#include "domain.h"
#include "field.h"
#include "materialmappingalgorithm.h"
#include "materialmappingalgorithmtype.h"

namespace oofem {
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
class InternalVariableField : public Field
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
    InternalVariableField(InternalStateType ist, FieldType b, MaterialMappingAlgorithmType mma_type, Domain *d);
    ~InternalVariableField();

    virtual int evaluateAt(FloatArray &answer, FloatArray &coords, IntArray &dofId, ValueModeType mode, TimeStep *atTime);
    InternalStateType giveType() { return type; }
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode);
    virtual const char *giveClassName() const { return "InternalVariableField"; }
};
} // end namespace oofem
#endif // intvarfield_h
