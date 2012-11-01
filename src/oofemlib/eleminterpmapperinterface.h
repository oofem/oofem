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

#ifndef eleminterpmapperinterface_h
#define eleminterpmapperinterface_h

#include "valuemodetype.h"
#include "interface.h"
#include "flotarry.h"
#include "intarray.h"
#include "error.h"

namespace oofem {

class TimeStep;

/**
 * The element interface class related to Element Interpolation Mappers.
 */
class EIPrimaryUnknownMapperInterface : public Interface
{
public:
    EIPrimaryUnknownMapperInterface() : Interface() { }

    /**
     * Computes the element vector of primary unknowns at given point. Similar to computeVectorOf,
     * but the interpolation from element DOFs to given point using element shape function is done.
     * The method should work also for point outside the volume of element (adaptivity mapping).
     * @param mode    Identifies mode of unknown (eg. total value or velocity of unknown).
     * @param tStep   Time step, when vector of unknowns is requested.
     * @param gcoords Global coordinates of point of interest.
     * @param answer  Vector of unknowns.
     * @return Nonzero if given point is in receiver volume otherwise zero
     */
    virtual int EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                                 TimeStep *tStep, const FloatArray &gcoords,
                                                                 FloatArray &answer) = 0;
    /**
     * Computes the element vector of primary unknowns at given point in the local coordinate system.
     * @see EIPrimaryUnknownMI_computePrimaryUnknownVectorAt
     * @param mode    Identifies mode of unknown (eg. total value or velocity of unknown).
     * @param tStep   Time step, when vector of unknowns is requested.
     * @param lcoords Local coordinates of point of interest.
     * @param answer  Vector of unknowns.
     * @return Nonzero if given point is in receiver volume otherwise zero
     */
    virtual void EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
            TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
    { OOFEM_ERROR("Not implemented\n"); }

    /**
     * Returns the dof meaning of element vector of primary unknowns.
     * @param answer Values of DofIDItem type that identify physical meaning of DOFs.
     */
    virtual void EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer) = 0;
};
} // end namespace oofem
#endif // eleminterpmapperinterface_h
