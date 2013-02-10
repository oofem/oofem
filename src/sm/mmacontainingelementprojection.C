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

#include "mmacontainingelementprojection.h"
#include "spatiallocalizer.h"
#include "domain.h"
#include "element.h"
#include "material.h"
#include "integrationrule.h"
#include "gausspnt.h"

namespace oofem {
MMAContainingElementProjection :: MMAContainingElementProjection() : MaterialMappingAlgorithm()
{ }

void
MMAContainingElementProjection :: __init(Domain *dold, IntArray &type, FloatArray &coords, int region, TimeStep *tStep)
{
    SpatialLocalizer *sl = dold->giveSpatialLocalizer();
    IntArray regionList(1);
    regionList.at(1) = region;
    GaussPoint *jGp;
    FloatArray jGpCoords;
    double distance, minDist = 1.e6;
    int j;
    IntegrationRule *iRule;
    Element *srcElem;

    if ( ( srcElem = sl->giveElementContainingPoint(coords, & regionList) ) ) {
        iRule = srcElem->giveDefaultIntegrationRulePtr();

        this->source = NULL;
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            jGp = iRule->getIntegrationPoint(j);
            if ( srcElem->computeGlobalCoordinates( jGpCoords, * ( jGp->giveCoordinates() ) ) ) {
                distance = coords.distance(jGpCoords);
                if ( distance < minDist ) {
                    minDist = distance;
                    this->source = jGp;
                }
            }
        }

        if ( !source ) {
            OOFEM_ERROR("MMAContainingElementProjection::__init : no suitable source found");
        }
    } else {
        OOFEM_ERROR("MMAContainingElementProjection: No suitable element found");
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
} // end namespace oofem
