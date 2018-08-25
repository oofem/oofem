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

#ifndef dofmanvalfield_h
#define dofmanvalfield_h

#include "field.h"

#include <vector>

namespace oofem {
class Domain;

/**
 * Class representing field defined by nodal values associated to given domain.
 * Field represent the spatial distribution of certain variable.
 * The implementation allows to set individual dofMan values;
 * However, in the current implementation doe not allow to specify values for different time steps.
 */
class OOFEM_EXPORT DofManValueField : public Field
{
protected:
    /// Associated domain (need its elements to interpolate)
    Domain *domain;
    /// Array of dofman values
    std::vector< FloatArray > dmanvallist;

public:
    /**
     * Constructor. Creates an empty field of given type associated to given domain.
     */
    DofManValueField(FieldType b, Domain * d);
    virtual ~DofManValueField() { }

    int evaluateAt(FloatArray &answer, const FloatArray &coords, ValueModeType mode, TimeStep *tStep) override;

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
     * @param tStep Time step to evaluate for.
     * @return Zero if ok, nonzero Error code (0-ok, 1-failed)
     */
    int evaluateAt(FloatArray &answer, DofManager *dman, ValueModeType mode, TimeStep *tStep) override;

    void saveContext(DataStream &stream) override;
    void restoreContext(DataStream &stream) override;

    /**
     * Sets the value associated to given dofManager
     */
    void setDofManValue(int dofMan, FloatArray value);

    const char *giveClassName() const override { return "DofManValueField"; }
};
} // end namespace oofem
#endif // dofmanvalfield_h
