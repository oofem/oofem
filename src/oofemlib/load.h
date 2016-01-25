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

#ifndef load_h
#define load_h

#include "generalboundarycondition.h"
#include "intarray.h"
#include "floatarray.h"
#include "valuemodetype.h"
#include "error.h"

///@name Input fields for Load
//@{
#define _IFT_Load_components "components"
#define _IFT_Load_dofexcludemask "dofexcludemask"
//@}

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
class OOFEM_EXPORT Load : public GeneralBoundaryCondition
{
public:
    /**
     * Load coordinate system type. Variable of this type can have following values CST_Global
     * (indicates that load given in global coordinate system) or CST_Local
     * (entity dependent local coordinate system will be  used).
     */
    enum CoordSystType {
        CST_Global, ///< Load is specified in global c.s.
        CST_Local, ///< Load is specified in local element c.s.
        CST_UpdatedGlobal, ///< Load is specified in global c.s. but follows the deformation
    };

    /**
     * Type determining the type of formulation (entity local or global one).
     */
    enum FormulationType {
        FT_Entity,
        FT_Global,
    };

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
    Load(int n, Domain * d);
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
     * @param answer Component values at given point and time.
     * @param tStep Time step representing time.
     * @param coords Global (or local) problem coordinates, which are used to evaluate components values.
     * @param mode Determines response mode.
     */
    virtual void computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode) = 0;
    /**
     * Computes components values for specified dof ids. If a dof id is not defined for the load, zero value is inserted.
     * Typically boundary conditions would specify the full, or a subset of the dofids.
     * @param answer Component values at given point and time.
     * @param tStep Time step representing time.
     * @param coords Global (or local) problem coordinates, which are used to evaluate components values.
     * @param dofids List of DOF IDs to evaluate for.
     * @param mode Determines response mode.
     */    
    virtual void computeValues(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, const IntArray &dofids, ValueModeType mode);
    /**
     * Returns the value of dofExcludeMask corresponding to given index.
     * See the description of dofExcludeMask attribute for more details.
     * @param index Index to check at.
     * @return Nonzero if excluded, zero otherwise.
     */
    int isDofExcluded(int index);

    /**
     * Returns receiver's coordinate system.
     */
    virtual CoordSystType giveCoordSystMode() { return CST_Global; }
    /**
     * Specifies is load should take local or global coordinates.
     */
    virtual FormulationType giveFormulationType() { return FT_Entity; }
    /**
     * @return Approximation order of load geometry.
     */
    virtual int giveApproxOrder() { return 0; }
    /**
     * Returns the value of a property 'aProperty'. Property must be identified
     * by unique integer id.
     * @param aProperty id of property requested
     * @return property value
     */
    virtual double giveProperty(int aProperty, TimeStep *tStep) {
        OOFEM_ERROR("Not supported for this boundary condition.");
        return 0;
    }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    /**
     * @return Pointer to receiver component array, where component values of boundary condition are stored.
     */
    const FloatArray &giveComponentArray() const;
    void setComponentArray(FloatArray &arry) { componentArray = arry; }
};
} // end namespace oofem
#endif // load_h
