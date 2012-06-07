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

#ifndef dofmanvalfield_h
#define dofmanvalfield_h

#include "domain.h"
#include "field.h"

namespace oofem {

/**
 * Class representing field defined by nodal values associated to given domain.
 * Field represent the spatial distribution of certain variable.
 * The implementation allows to set individual dofMan values;
 * However, in the current implementation doe not allow to specify values for different time steps.
 */
class DofManValueField : public Field
{
protected:
    /// Associated domain (need its elements to interpolate)
    Domain* domain;
    /// Array of dofman values
    AList <FloatArray> dmanvallist;

public:
    /**
     * Constructor. Creates an empty field of given type associated to given domain.
     */
    DofManValueField(FieldType b, Domain *d);
    ~DofManValueField() {}
    /**
     * Evaluates the field at given point.
     * @param coords Coordinates of the point of interest
     * @param answer Field evaluated at coordinate.
     * @param atTime Time step to evaluate for.
     * @param mode Mode of value (total, velocity,...).
     * @return Zero if ok, otherwise nonzero.
     */
    virtual int evaluateAt(FloatArray &answer, FloatArray &coords,
                           ValueModeType mode, TimeStep *atTime);

    /**
     * Evaluates the field at given DofManager. This potentially can be resolved quickly, as
     * receiver data may be described using values at dof managers. Here an additional issue
     * exists: one needs to make sure, that passed dman is from the same domain, so that its
     * number can be used to perform suggested quick evaluation.
     *
     * If this is not the case (the field is described differently),
     * the response can be evaluated using dofman coordinates in a standard way.
     * @param[out] answer Evaluated field for dman.
     * @param dman Reference to dof manager.
     * @param mode Mode of value (total, velocity,...).
     * @param atTime Time step to evaluate for.
     * @return Zero if ok, nonzero Error code (0-ok, 1-failed)
     */
    virtual int evaluateAt(FloatArray &answer, DofManager* dman,
                           ValueModeType mode, TimeStep *atTime) ;

    /**
     * Stores receiver state to output stream.
     * Writes the FEMComponent class-id in order to allow test whether correct data are then restored.
     * @param stream Output stream.
     * @param mode Determines amount of info in stream (state, definition,...).
     * @return contextIOResultType.
     * @exception Throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode) ;
    /**
     * Restores the receiver state previously written in stream.
     * Reads the FEMComponent class-id in order to allow test consistency.
     * @param stream Input stream.
     * @param mode Determines amount of info in stream (state, definition,...).
     * @return contextIOResultType.
     * @exception Throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode) ;


    /**
     * Sets the value associated to given dofManager
     */
    void setDofManValue (int dofMan, const FloatArray &value);

    /// @return Class name of the receiver.
    virtual const char *giveClassName() const { return "DofManValueField"; }
};
} // end namespace oofem
#endif // dofmanvalfield_h
