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

#ifndef field_h
#define field_h

#include "oofemcfg.h"
#include "valuemodetype.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "enumitem.h"
#include <string>

namespace oofem {
///@todo FieldType and UnknownType basically determine the same thing. Should be possible to stick to one. Combinations of fields should be possible with logical bitfields.
#define FieldType_DEF \
    ENUM_ITEM_WITH_VALUE(FT_Unknown, 0) \
    ENUM_ITEM_WITH_VALUE(FT_Velocity, 1) \
    ENUM_ITEM_WITH_VALUE(FT_Displacements, 2) \
    ENUM_ITEM_WITH_VALUE(FT_VelocityPressure, 3) \
    ENUM_ITEM_WITH_VALUE(FT_Pressure, 4) \
    ENUM_ITEM_WITH_VALUE(FT_Temperature, 5) \
    ENUM_ITEM_WITH_VALUE(FT_HumidityConcentration, 6) \
    ENUM_ITEM_WITH_VALUE(FT_TransportProblemUnknowns, 7)

/// Physical type of field.
enum FieldType {
    FieldType_DEF
};
#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h

class TimeStep;
class FloatArray;
class DofManager;
class DataStream;

/**
 * Abstract class representing field. Field represent the spatial distribution of certain variable.
 * Field is able to evaluate its value at any point of interest. The field is usually associated to
 * specific domain.
 */
class OOFEM_EXPORT Field
{
protected:
    FieldType type;

public:
    /**
     * Constructor. Creates a field of given type associated to given domain.
     */
    Field(FieldType b) : type(b) { }
    virtual ~Field() { }
    /**
     * Evaluates the field at given point.
     * @param coords Coordinates of the point of interest
     * @param answer Field evaluated at coordinate.
     * @param tStep Time step to evaluate for.
     * @param mode Mode of value (total, velocity,...).
     * @return Zero if ok, otherwise nonzero.
     */
    virtual int evaluateAt(FloatArray &answer, FloatArray &coords,
                           ValueModeType mode, TimeStep *tStep) = 0;

    /**
     * Evaluates the field at given DofManager. This potentially can be resolved quickly, as
     * receiver data may be described using values at dofManagers. Here an additional issue
     * exists: one needs to make sure, that passed dman is from the same domain, so that its
     * number can be used to perform suggested quick evaluation.
     *
     * If this is not the case (the field is described differently),
     * the response can be evaluated using dofman coordinates in a standard way.
     * @param[out] answer Evaluated field for dman.
     * @param dman Reference to dofManager.
     * @param mode Mode of value (total, velocity,...).
     * @param tStep Time step to evaluate for.
     * @return Zero if ok, nonzero Error code (0-ok, 1-failed)
     */
    virtual int evaluateAt(FloatArray &answer, DofManager *dman,
                           ValueModeType mode, TimeStep *tStep) = 0;

    /// Returns the type of receiver
    FieldType giveType() { return type; }

    /**
     * Stores receiver state to output stream.
     * Writes the FEMComponent class-id in order to allow test whether correct data are then restored.
     * @param stream Output stream.
     * @param mode Determines amount of info in stream (state, definition,...).
     * @return contextIOResultType.
     * @exception Throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode) = 0;
    /**
     * Restores the receiver state previously written in stream.
     * Reads the FEMComponent class-id in order to allow test consistency.
     * @param stream Input stream.
     * @param mode Determines amount of info in stream (state, definition,...).
     * @return contextIOResultType.
     * @exception Throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode) = 0;


    /// Returns string for prepending output (used by error reporting macros).
    std :: string errorInfo(const char *func) const;

    /// @return Class name of the receiver.
    virtual const char *giveClassName() const = 0;
};
} // end namespace oofem
#endif // field_h
