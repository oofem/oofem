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

#include "materialmappingalgorithm.h"
#include "gausspoint.h"
#include "element.h"

namespace oofem {
void
MaterialMappingAlgorithm :: init(Domain *dold, IntArray &type, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray coords;
    if ( gp->giveElement()->computeGlobalCoordinates( coords, * ( gp->giveCoordinates() ) ) == 0 ) {
        OOFEM_ERROR("MaterialMappingAlgorithm::init: computeGlobalCoordinates failed");
    }

    this->__init(dold, type, coords, gp->giveElement()->giveRegionNumber(), tStep);
}

int
MaterialMappingAlgorithm :: mapVariable(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatArray coords;
    if ( gp->giveElement()->computeGlobalCoordinates( coords, * ( gp->giveCoordinates() ) ) == 0 ) {
        OOFEM_ERROR("MaterialMappingAlgorithm::init : computeGlobalCoordinates failed");
    }

    return this->__mapVariable(answer, coords, type, tStep);
}
} // end namespace oofem
