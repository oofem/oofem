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

#ifndef maskedprimaryfield_h
#define maskedprimaryfield_h

#include "primaryfield.h"
#include "floatarray.h"
#include "intarray.h"
#include "valuemodetype.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#include <memory>

namespace oofem {
/**
 * Abstract class representing subset of DOFs (identified by DofId mask) of primary field.
 * As the PrimaryField stores the state directly in solution vectors that are usually directly
 * updated by EngngModel, it may contain a mix of different fields (this is especially true for
 * strongly coupled problems). Then masked primary field can be used to select only certain DOFs
 * (based on DofID) from its master PrimaryField.
 */
class OOFEM_EXPORT MaskedPrimaryField : public Field
{
protected:
    PrimaryField *master;
    IntArray mask;
public:
    MaskedPrimaryField(FieldType b, PrimaryField * m, IntArray dofIdMask) : Field(b), mask(std::move(dofIdMask)) {
        master = m;
    }

    virtual int evaluateAt(FloatArray &answer, FloatArray &coords,
                           ValueModeType mode, TimeStep *tStep);
    virtual int evaluateAt(FloatArray &answer, DofManager *dman,
                           ValueModeType mode, TimeStep *tStep);

    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode) { return CIO_OK; }
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode) { return CIO_OK; }

    virtual const char *giveClassName() const { return "MaskedPrimaryField"; }
};
} // end namespace oofem
#endif // maskedprimaryfield_h
