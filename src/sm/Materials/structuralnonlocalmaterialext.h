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

#ifndef structuralnonlocalmaterialext_h
#define structuralnonlocalmaterialext_h

#include "nonlocalmaterialext.h"
#include "matstatus.h"

namespace oofem {
/**
 * Base class for all nonlocal structural material statuses.
 */
class StructuralNonlocalMaterialStatusExtensionInterface : public NonlocalMaterialStatusExtensionInterface
{
public:
    StructuralNonlocalMaterialStatusExtensionInterface() : NonlocalMaterialStatusExtensionInterface() { }
    ~StructuralNonlocalMaterialStatusExtensionInterface() { }
};


/**
 * Abstract base class for all nonlocal structural materials. Nonlocal in the sense that the response at a particular
 * point depends not only on the state of that point, but also takes into account the state of surrounding
 * points. The response typically depends on some nonlocal quantity obtained as weighted average over
 * some characteristic volume.
 * This class declares the necessary interface for all nonlocal structural constitutive  models.
 * It also stores the parameters needed for nonlocal averaging and provides services for evaluation of the weight function.
 */
class StructuralNonlocalMaterialExtensionInterface : public NonlocalMaterialExtensionInterface
{
public:
    /**
     * Constructor.
     * @param d Domain to which new material will belong.
     */
    StructuralNonlocalMaterialExtensionInterface(Domain * d) : NonlocalMaterialExtensionInterface(d)
    { }
    /// Destructor.
    ~StructuralNonlocalMaterialExtensionInterface() { }

    /**
     * Declares the service updating local variables in given integration points,
     * which take part in nonlocal average process.
     * Because value of single integration point influences nonlocal variables in several near
     * integration points, it is suitable to compute these variables only once. These should be stored
     * in integration point associated statuses.
     * The implementation is left on derived classes.
     * Provide material local strain increment - as is provided to computeRealStresVector.
     * This allows to update internal vars to be averaged to new state
     * @param strainVector Total strain vector in given integration point.
     * @param gp Integration point to update.
     * @param tStep Solution step indicating time of update.
     */
    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) = 0;
};
} // end namespace oofem
#endif // structuralnonlocalmaterialext_h
