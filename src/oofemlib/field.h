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

#ifndef field_h
#define field_h

#include "domain.h"
#include "valuemodetype.h"
#include "contextioresulttype.h"
#include "contextmode.h"

namespace oofem {

/// Physical type of field.
enum FieldType {
    FT_Unknown,

    FT_Velocity,
    FT_Displacements,
    FT_VelocityPressure,

    FT_Pressure,
    FT_Temperature,
    FT_HumidityConcentration,
    FT_TransportProblemUnknowns

};


/**
 * Abstract class representing field. Field represent the spatial distribution of certain variable.
 * Field is able to evaluate its value at any point of interest. The field is usually associated to
 * specific domain.
 */
class Field
{
protected:
    FieldType type;

public:
    /** 
     * Constructor. Creates a field of given type associated to given domain.
     */
    Field(FieldType b) { type = b; }
    virtual ~Field() { }
    /**
     * Evaluates the field at given point.
     * @param coords Coordinates of the point of interest
     * @param answer Field evaluated at coordinate.
     * @param atTime Time step to evaluate for.
     * @param mode Mode of value (total, velocity,...).
     * @return Zero if ok, otherwise nonzero.
     */
    virtual int evaluateAt(FloatArray &answer, FloatArray &coords,
                           ValueModeType mode, TimeStep *atTime) = 0;

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
     * @param atTime Time step to evaluate for.
     * @return Zero if ok, nonzero Error code (0-ok, 1-failed)
     */
    virtual int evaluateAt(FloatArray &answer, DofManager* dman,
                           ValueModeType mode, TimeStep *atTime) = 0;
    
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
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode) = 0;
    /**
     * Restores the receiver state previously written in stream.
     * Reads the FEMComponent class-id in order to allow test consistency.
     * @param stream Input stream.
     * @param mode Determines amount of info in stream (state, definition,...).
     * @return contextIOResultType.
     * @exception Throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode) = 0;

    /**
     * @name Error and warning reporting methods
     * These methods will print error (or warning) message using oofem default loggers.
     * Do not use these methods directly, to avoid specify file and line parameters.
     * More preferably, use these methods via corresponding OOFEM_CLASS_ERROR and OOFEM_CLASS_WARNING macros,
     * that will include file and line parameters automatically.
     *
     * Uses variable number of arguments, so a format string followed by optional arguments is expected
     * (according to printf conventions).
     *
     * @param file Source file name, where error encountered (where error* function called)
     * @param line Source file line number, where error encountered
     *
     */
    //@{
    /// Prints error message and exits
    void error(const char *file, int line, const char *format, ...) const;
    /// Prints warning message
    void warning(const char *file, int line, const char *format, ...) const;
    //@}

    /// @return Class name of the receiver.
    virtual const char *giveClassName() const { return "Field"; }
};
} // end namespace oofem
#endif // field_h
