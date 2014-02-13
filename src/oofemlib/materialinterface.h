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

#ifndef materialinterface_h
#define materialinterface_h

#include "femcmpnn.h"

namespace oofem {
/**
 * Abstract base class representing (moving) material interfaces.
 * Its typical use to model moving interface (such as free surface)
 * in a fixed-grid methods (as typically used in CFD).
 * The basic tasks are representation of interface and its updating.
 */
class OOFEM_EXPORT MaterialInterface : public FEMComponent
{
public:
    /**
     * Constructor. Takes two two arguments. Creates
     * MaterialInterface instance with given number and belonging to given domain.
     * @param n Component number in particular domain. For instance, can represent
     * node number in particular domain.
     * @param d Domain to which component belongs to.
     */
    MaterialInterface(int n, Domain * d) : FEMComponent(n, d) { }

    virtual const char *giveInputRecordName() const { return NULL; }

    /**
     *  Initializes receiver
     */
    virtual void initialize() { }
    /**
     * Updates the position of interface according to state reached in given solution step.
     */
    virtual void updatePosition(TimeStep *tStep) = 0;
    /**
     * Updates element state after equilibrium in time step has been reached.
     * All temporary history variables,
     * which now describe equilibrium state should be copied into equilibrium ones.
     * The existing internal state is used for update.
     */
    virtual void updateYourself(TimeStep *tStep) = 0;
    /**
     * Computes critical time step induced by receiver integration algorithm
     */
    virtual double computeCriticalTimeStep(TimeStep *tStep) = 0;

    /**
     * Returns relative material contents at given point. Usually only one material is presented in given point,
     * but some smoothing may be applied close to material interface to make transition smooth
     */
    virtual void giveMaterialMixtureAt(FloatArray &answer, FloatArray &position) = 0;
    /**
     * Returns volumetric (or other based measure) of relative material contents in given element.
     */
    virtual void giveElementMaterialMixture(FloatArray &answer, int ielem) = 0;
    /** Returns scalar value representation of material Interface at given point. For visualization */
    virtual double giveNodalScalarRepresentation(int) = 0;
};
} // end namespace oofem
#endif // materialinterface_h
