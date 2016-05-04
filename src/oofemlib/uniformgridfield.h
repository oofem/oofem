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

#ifndef uniformgridfield_h
#define uniformgridfield_h

#include"field.h"
#include"floatarray.h"
#include"intarray.h"

namespace oofem {
/**
 * Field defined by values in uniform grid nodes, with linear interpolation for points inside the grid, or interpolation for the closest point within the grid for points outside.
 */
class OOFEM_EXPORT UniformGridField: public Field
{
protected:
    // points defining axis-aligned box for the grid
    FloatArray lo;
    FloatArray hi;
    // number of divisions along each axis (thus, each axis has div[axis]+1 points!)
    IntArray div;
    // values in all grid points, as linear array; its size must be prod(div)
    // the array is c-ordered, the last index varying the fastest (FIXME: check)
    FloatArray values;
    #if 0
        // precompute useful interally used values (such as cell size)
        void precomputeInternal(); 
    #endif

    void xyz2ijk(const FloatArray& xyz, IntArray& ijk, FloatArray& normXyz) const;
    double nodeValue2d(int i, int j);
    double nodeValue3d(int i, int j, int k);
public:
    /**
     * Constructor. Creates a field, with unspecified field values.
     */
    UniformGridField() : Field(FieldType::FT_Unknown) { }
    virtual ~UniformGridField() { }

    /** Shorthand for defining geometry, with consistency checks. Used primarily from python */
    void setGeometry(const FloatArray& lo_, const FloatArray& hi_, const IntArray& div_);
    /** Accessor for setting nodal values; checks size of the array for correctness. */
    void setValues(const FloatArray& vv);

    /**
        Implementation of Field::evaluateAt for coordinates.
     */
    int evaluateAt(FloatArray &answer, const FloatArray &coords,
                           ValueModeType mode, TimeStep *tStep) override;

    /**
        Implementaton of Field::evaluateAt for DofManager.
     */
    int evaluateAt(FloatArray &answer, DofManager *dman,
                           ValueModeType mode, TimeStep *tStep) override;

    /**
     * Stores receiver state to output stream.
     * Writes the FEMComponent class-id in order to allow test whether correct data are then restored.
     * @param stream Output stream.
     * @param mode Determines amount of info in stream (state, definition,...).
     * @return contextIOResultType.
     * @exception Throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType saveContext(DataStream &stream, ContextMode mode) override { return CIO_OK; }
    /**
     * Restores the receiver state previously written in stream.
     * Reads the FEMComponent class-id in order to allow test consistency.
     * @param stream Input stream.
     * @param mode Determines amount of info in stream (state, definition,...).
     * @return contextIOResultType.
     * @exception Throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType restoreContext(DataStream &stream, ContextMode mode) override { return CIO_OK; }

    /// @return Class name of the receiver.
    const char *giveClassName() const override { return "UniformGridField"; }
};
} // end namespace oofem
#endif // field_h
