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

#ifndef materialmapperinterface_h
#define materialmapperinterface_h

#include "interface.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

/**
 * The class representing the general material model adaptive mapping interface.
 * The basic task is to define the algorithm for mapping of the internal material model
 * variables from one (old) mesh to the given IP of new mesh.
 */
class OOFEM_EXPORT MaterialModelMapperInterface : public Interface
{
protected:

public:
    /// Constructor.
    MaterialModelMapperInterface() : Interface() { }
    /// Destructor.
    virtual ~MaterialModelMapperInterface() { }

    /**
     * @name The interface required by MaterialModelMapperInterface
     */
    //@{
    /**
     * Maps the required internal state variables from
     * old mesh oldd to given ip. The result is stored in gp status.
     * @param gp Integration point belonging to new domain which values will be mapped
     * @param oldd Old mesh reference.
     * @param tStep Time step.
     * @return Nonzero if o.k.
     */
    virtual int MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep) = 0;
    /**
     * Updates the required internal state variables from previously mapped values.
     * The result is stored in gp status. This map and update splitting is necessary,
     * for example for nonlocal models that local quantity to be averaged must be mapped in all integration points
     * and then update can happen, because it may depend on nonlocal variable, which is computed
     * from local values.
     * @param gp Integration point belonging to new domain which values will be mapped.
     * @param tStep Time step.
     * @param elemGPVec Vector passed to MMI_update at material level, probably computed from primary unknowns
     * (for structural elements this represent strain vector).
     * @return Nonzero if o.k.
     */
    virtual int MMI_update(GaussPoint *gp, TimeStep *tStep, FloatArray *elemGPVec = NULL) = 0;
    /**
     * Finishes the mapping for given time step. Used to perform cleanup.
     * Typically some mappers require to compute some global mesh data related to
     * current step, which are valid for example to all IPs - so they are computed only once for
     * all IPs, stored and they need to be deallocated. These mappers are typically class variables,
     * but their finish is invoked by all members.
     */
    virtual int MMI_finish(TimeStep *tStep) = 0;
    //@}
};
} // end namespace oofem
#endif // materialmapperinterface_h
