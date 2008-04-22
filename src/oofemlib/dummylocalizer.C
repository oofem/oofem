/* $Header: /home/cvs/bp/oofem/oofemlib/src/dummylocalizer.C,v 1.6 2003/04/06 14:08:23 bp Exp $ */
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

#include "dummylocalizer.h"
#include "element.h"
#include "domain.h"
#include "integrationrule.h"
#include "gausspnt.h"

Element *
DummySpatialLocalizer :: giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList)
{
    // Dummy implementation here, to be replaced ASAP
    int ielem, nelems = this->giveDomain()->giveNumberOfElements();
    Element *ielemptr;
    SpatialLocalizerInterface *interface;

    for ( ielem = 1; ielem <= nelems; ielem++ ) {
        ielemptr = this->giveDomain()->giveElement(ielem);
        interface = ( SpatialLocalizerInterface * ) ielemptr->giveInterface(SpatialLocalizerInterfaceType);
        if ( interface ) {
            if ( regionList && ( regionList->findFirstIndexOf( ielemptr->giveRegionNumber() ) == 0 ) ) {
                continue;
            }

            if ( interface->SpatialLocalizerI_BBoxContainsPoint(coords) == 0 ) {
                continue;
            }

            if ( interface->SpatialLocalizerI_containsPoint(coords) ) {
                return ielemptr;
            }
        }
    }

    return NULL;
}


Element *
DummySpatialLocalizer :: giveElementCloseToPoint(const FloatArray &coords, const IntArray *regionList)
{
    // Dummy implementation here, to be replaced ASAP
    int ielem, nelems = this->giveDomain()->giveNumberOfElements();
    Element *ielemptr, *answer = NULL;
    SpatialLocalizerInterface *interface;
    double dist = 0.0, currDist;

    for ( ielem = 1; ielem <= nelems; ielem++ ) {
        ielemptr = this->giveDomain()->giveElement(ielem);
        interface = ( SpatialLocalizerInterface * ) ielemptr->giveInterface(SpatialLocalizerInterfaceType);
        if ( interface ) {
            if ( regionList && ( regionList->findFirstIndexOf( ielemptr->giveRegionNumber() ) == 0 ) ) {
                continue;
            }

            currDist = interface->SpatialLocalizerI_giveDistanceFromParametricCenter(coords);
            if ( answer == NULL || currDist < dist ) {
                answer = ielemptr;
                if ( ( dist = currDist ) == 0.0 ) {
                    break;
                }
            }
        }
    }

    return NULL;
}


GaussPoint *
DummySpatialLocalizer :: giveClosestIP(const FloatArray &coords, int region)
{
    int nelem, i, j;
    double minDist = 1.e6, distance;
    Element *ielem;
    GaussPoint *jGp, *answer = NULL;
    IntegrationRule *iRule;
    FloatArray jGpCoords;

    nelem = this->giveDomain()->giveNumberOfElements();


    for ( i = 1; i <= nelem; i++ ) {
        ielem = this->giveDomain()->giveElement(i);
        if ( ( region < 0 ) || ( region == ielem->giveRegionNumber() ) ) {
            iRule = ielem->giveDefaultIntegrationRulePtr();
            for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
                jGp = iRule->getIntegrationPoint(j);
                if ( ielem->computeGlobalCoordinates( jGpCoords, * ( jGp->giveCoordinates() ) ) ) {
                    distance = coords.distance(jGpCoords);
                    if ( distance < minDist ) {
                        minDist = distance;
                        answer = jGp;
                    }
                }
            }
        }
    } // loop over elements

    return answer;
}



void
DummySpatialLocalizer :: giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords, const double radius)
{
    int i, j, nelem;
    double currDist;
    FloatArray jGpCoords;
    Element *ielem;
    IntegrationRule *iRule;
    GaussPoint *jGp;


    nelem = this->giveDomain()->giveNumberOfElements();
    for ( i = 1; i <= nelem; i++ ) {
        ielem = this->giveDomain()->giveElement(i);
        iRule = ielem->giveDefaultIntegrationRulePtr();
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            jGp = iRule->getIntegrationPoint(j);
            if ( ielem->computeGlobalCoordinates( jGpCoords, * ( jGp->giveCoordinates() ) ) ) {
                currDist = coords.distance(jGpCoords);
                if ( currDist <= radius ) {
                    elemSet.insert(i);
                }
            }
        }
    } // end element loop

}
