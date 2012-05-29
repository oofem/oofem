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

#ifndef load_h
#define load_h

#include "generalbc.h"
#include "domain.h"
#include "flotarry.h"
#include "dictionr.h"
#include "valuemodetype.h"

namespace oofem {
/**
 * Load is base abstract class for all loads.
 * Load is an attribute of the domain that is belongs to, and also of several elements, nodes,
 * which are subjected to loading type boundary condition.
 *
 * The value (or the components) of a load will be
 * the product of its value (stored in componentArray) by the value of
 * the associated load time function at given time step.
 */
class Load : public GeneralBoundaryCondition
{
protected:
    /// Components of boundary condition.
    FloatArray componentArray;
    /**
     * The load is specified for all dofs of object to which is associated.
     * For some types of boundary conditions the zero value of load does not mean
     * that the load is not applied (newton's type of bc, for example). Then
     * some mask, which allows to exclude specific dofs is necessary.
     * The dofMask attribute is introduced to allow this.
     * By default it is of the same size as componentArray, filled with zeroes.
     * If some value of dofExcludeMask is set to nonzero, then the corresponding componentArray
     * is set to zero.
     */
    IntArray dofExcludeMask;

public:
    /**
     * Constructor. Creates boundary condition with given number, belonging to given domain.
     * @param n Boundary condition number
     * @param d Domain to which new object will belongs.
     */
    Load(int n, Domain *d);
    /// Destructor.
    virtual ~Load() { }

    /**
     * Computes boundary condition value - its components values at given time.
     * Default implementation returns as the answer its component array multiplied
     * with load time function value (load response mode is taken in to account)
     * @param answer Computed boundary conditions components.
     * @param tStep Time step, for which components are computed.
     * @param mode Determines response mode.
     */
    virtual void computeComponentArrayAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    /**
     * Computes components values of load at given point - global coordinates (coordinates given).
     * Default implementation computes product of approximation matrix and
     * with "vertex" value array attribute and the result is then multiplied by
     * corresponding load time function value respecting load response mode.
     * @param answer Component values at given point and time.
     * @param tStep Time step representing time.
     * @param coords Global (or local) problem coordinates, which are used to
     * evaluate components values.
     * @param mode Determines response mode.
     */
    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, FloatArray &coords, ValueModeType mode) = 0;
    /**
     * Returns the value of dofExcludeMask corresponding to given index.
     * See the description of dofExcludeMask attribute for more details.
     * @param index Index to check at.
     * @return Nonzero if excluded, zero otherwise.
     */
    int isDofExcluded(int index);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual classType giveClassID() const { return LoadClass; }
    virtual const char *giveClassName() const { return "Load"; }

protected:
    /**
     * @return Pointer to receiver component array, where component values of boundary condition are stored.
     */
    FloatArray &giveComponentArray();

public:
    void setComponentArray(FloatArray &arry) { componentArray = arry; }
	 FloatArray giveCopyOfComponentArray() { FloatArray answer = componentArray; return answer; }
};
} // end namespace oofem
#endif // load_h

