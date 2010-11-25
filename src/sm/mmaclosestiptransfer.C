/* $Header: /home/cvs/bp/oofem/sm/src/mmaclosestiptransfer.C,v 1.6.4.1 2004/04/05 15:19:47 bp Exp $ */
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
#include "mmaclosestiptransfer.h"
#include "spatiallocalizer.h"
#include "domain.h"
#include "element.h"
#include "material.h"
#include "gausspnt.h"

namespace oofem {
/*
 * MMAClosestIPTransfer::MMAClosestIPTransfer(GaussPoint* gp, Domain* oldd) : MaterialMappingAlgorithm ()
 * {
 * FloatArray coords;
 * SpatialLocalizer* sl = oldd -> giveSpatialLocalizer();
 * gp->giveElement()->computeGlobalCoordinates (coords, *(gp->giveCoordinates()));
 *
 * this->source = sl->giveClosestIP (coords, gp->giveElement()->giveRegionNumber());
 * }
 *
 * int
 * MMAClosestIPTransfer::mapVariable (FloatArray& answer, InternalStateType type, TimeStep* tStep)
 * {
 *
 * if (source) {
 * source->giveMaterial()->giveIPValue (answer, source, type, tStep);
 * return 1;
 * } else {
 * printf ("MMAClosestIPTransfer::mapVariable: no suitable source found");
 * return 0;
 * }
 * }
 */


MMAClosestIPTransfer :: MMAClosestIPTransfer() : MaterialMappingAlgorithm()
{ }

void
MMAClosestIPTransfer :: __init(Domain *dold, IntArray &type, FloatArray &coords, int region, TimeStep *tStep)
{
    SpatialLocalizer *sl = dold->giveSpatialLocalizer();
    this->source = sl->giveClosestIP(coords, region);

    if ( !source ) {
        OOFEM_ERROR("MMAClosestIPTransfer::__init : no suitable source found");
    }
}

int
MMAClosestIPTransfer :: __mapVariable(FloatArray &answer, FloatArray &coords,
                                      InternalStateType type, TimeStep *tStep)
{
    if ( source ) {
        source->giveMaterial()->giveIPValue(answer, source, type, tStep);
        return 1;
    }

    return 0;
}
} // end namespace oofem
