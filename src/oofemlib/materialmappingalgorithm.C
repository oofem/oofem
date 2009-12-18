/* $Header: /home/cvs/bp/oofem/oofemlib/src/materialmappingalgorithm.C,v 1.1.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#include "materialmappingalgorithm.h"
#include "gausspnt.h"
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
