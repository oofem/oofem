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

#ifndef eleminterpmapperinterface_h
#define eleminterpmapperinterface_h

#include "valuemodetype.h"
#include "interface.h"
#include "floatarray.h"
#include "intarray.h"
#include "error.h"

namespace oofem {
class TimeStep;

/**
 * The element interface class related to Element Interpolation Mappers.
 */
class OOFEM_EXPORT EIPrimaryUnknownMapperInterface : public Interface
{
public:
    EIPrimaryUnknownMapperInterface() : Interface() { }

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
    { OOFEM_ERROR("Not implemented"); }
};
} // end namespace oofem
#endif // eleminterpmapperinterface_h
