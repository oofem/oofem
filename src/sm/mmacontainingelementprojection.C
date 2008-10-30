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
#include "mmacontainingelementprojection.h"
#include "spatiallocalizer.h"
#include "domain.h"
#include "element.h"
#include "material.h"
#include "integrationrule.h"
#include "gausspnt.h"

MMAContainingElementProjection :: MMAContainingElementProjection() : MaterialMappingAlgorithm()
{ }

void
MMAContainingElementProjection :: __init(Domain *dold, IntArray &type, FloatArray &coords, int region, TimeStep *tStep)
{
    SpatialLocalizer *sl = dold->giveSpatialLocalizer();	
    IntArray regionList(1); regionList.at(1)=region;
    Element *srcElem = sl->giveElementContainingPoint (coords, &regionList);
    IntegrationRule *iRule = srcElem->giveDefaultIntegrationRulePtr();
    GaussPoint *jGp;
    FloatArray jGpCoords;
    double distance, minDist = 1.e6;
    int j;
    
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













