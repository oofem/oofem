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

#include "mmacontainingelementprojection.h"
#include "spatiallocalizer.h"
#include "domain.h"
#include "element.h"
#include "material.h"
#include "integrationrule.h"
#include "gausspoint.h"

namespace oofem {
MMAContainingElementProjection :: MMAContainingElementProjection() : MaterialMappingAlgorithm()
{ }

void
MMAContainingElementProjection :: __init(Domain *dold, IntArray &type, FloatArray &coords, Set &elemSet, TimeStep *tStep, bool iCohesiveZoneGP)
{
    SpatialLocalizer *sl = dold->giveSpatialLocalizer();
    FloatArray jGpCoords;
    double distance, minDist = 1.e6;
    Element *srcElem;

    if ( ( srcElem = sl->giveElementContainingPoint(coords, elemSet) ) ) {
        this->source = NULL;
        for ( GaussPoint *jGp: *srcElem->giveDefaultIntegrationRulePtr() ) {
            if ( srcElem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                distance = coords.distance(jGpCoords);
                if ( distance < minDist ) {
                    minDist = distance;
                    this->source = jGp;
                }
            }
        }

        if ( !source ) {
            OOFEM_ERROR("no suitable source found");
        }
    } else {
        OOFEM_ERROR("No suitable element found");
    }
}

int
MMAContainingElementProjection :: __mapVariable(FloatArray &answer, FloatArray &coords,
                                                InternalStateType type, TimeStep *tStep)
{
    if ( source ) {
        source->giveMaterial()->giveIPValue(answer, source, type, tStep);
        return 1;
    }

    return 0;
}

int
MMAContainingElementProjection :: mapStatus(MaterialStatus &oStatus) const
{
    OOFEM_ERROR("mapStatus() is not implemented yet.")

    return 0;
}
} // end namespace oofem
