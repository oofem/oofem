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

#include "dummylocalizer.h"
#include "element.h"
#include "domain.h"
#include "integrationrule.h"
#include "gausspnt.h"
#include "node.h"

namespace oofem {

int
DummySpatialLocalizer :: init(bool force)
{
    if (!force && this->initialized) {
        return true;
    }

    // Count the elements in each cross section.
    int r;
    int nregion = this->domain->giveNumberOfRegions();
    int nelem = this->domain->giveNumberOfElements();

    IntArray region_nelem(nregion);
    region_nelem.zero();
    for (int i = 1; i <= nelem; i++) {
        Element *e = this->domain->giveElement(i);
        r = e->giveRegionNumber();
        region_nelem.at(r)++;
    }

    this->region_elements.resize(nregion);
    // Creates a new int array of correct size for each region
    for (int i = 1; i <= nregion; i++) {
        this->region_elements[i-1].resize(region_nelem.at(i));
    }
    // Add the numbers into the list.
    IntArray c(nregion);
    c.zero();
    for (int i = 1; i <= nelem; i++) {
        Element *e = this->domain->giveElement(i);
        r = e->giveRegionNumber();
        c.at(r)++;
        this->region_elements[r-1].at(c.at(r)) = i;
    }
    return this->initialized = true;
}

Element *
DummySpatialLocalizer :: giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList)
{
    // Dummy implementation here, to be replaced ASAP
    int nelems = this->giveDomain()->giveNumberOfElements();
    Element *ielemptr;
    SpatialLocalizerInterface *interface;

    for ( int ielem = 1; ielem <= nelems; ielem++ ) {
        ielemptr = this->giveDomain()->giveElement(ielem);
        interface = static_cast< SpatialLocalizerInterface * >( ielemptr->giveInterface(SpatialLocalizerInterfaceType) );
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
    int nelems = this->giveDomain()->giveNumberOfElements();
    Element *ielemptr, *answer = NULL;
    SpatialLocalizerInterface *interface;
    double dist = 0.0, currDist;

    for ( int ielem = 1; ielem <= nelems; ielem++ ) {
        ielemptr = this->giveDomain()->giveElement(ielem);
        interface = static_cast< SpatialLocalizerInterface * >( ielemptr->giveInterface(SpatialLocalizerInterfaceType) );
        if ( interface ) {
            if ( regionList && ( regionList->findFirstIndexOf( ielemptr->giveRegionNumber() ) == 0 ) ) {
                continue;
            }

            currDist = interface->SpatialLocalizerI_giveDistanceFromParametricCenter(coords);
            if ( answer == NULL || currDist < dist ) {
                answer = ielemptr;
                dist = currDist;
                if ( dist == 0.0 ) {
                    break;
                }
            }
        }
    }

    return answer;
}


Element *
DummySpatialLocalizer :: giveElementClosestToPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &coords, int region)
{
    int nelems;
    Element *ielemptr, *answer = NULL;
    SpatialLocalizerInterface *interface;
    double dist = 0.0, currDist;
    FloatArray el_coords, el_lcoords;
    lcoords.resize(0);
    closest.resize(0);

    if ( region > 0 ) {
        IntArray &elems = this->region_elements[region-1];
        for ( int ielem = 1; ielem <= elems.giveSize(); ielem++ ) {
            ielemptr = this->domain->giveElement(elems.at(ielem));
            interface = static_cast< SpatialLocalizerInterface * >( ielemptr->giveInterface(SpatialLocalizerInterfaceType) );
            if ( interface ) {
                currDist = interface->SpatialLocalizerI_giveClosestPoint(el_lcoords, el_coords, coords);
                if ( answer == NULL || (currDist < dist && currDist >= 0.0) ) {
                    answer = ielemptr;
                    lcoords = el_lcoords;
                    closest = el_coords;
                    dist = currDist;
                    if ( dist == 0.0 ) {
                        break;
                    }
                }
            }
        }
    } else { // Check them all;
        nelems = this->giveDomain()->giveNumberOfElements();
        for ( int ielem = 1; ielem <= nelems; ielem++ ) {
            ielemptr = this->giveDomain()->giveElement(ielem);
            interface = static_cast< SpatialLocalizerInterface * >( ielemptr->giveInterface(SpatialLocalizerInterfaceType) );
            if ( interface ) {
                if ( region > 0 && region == ielemptr->giveRegionNumber() ) {
                    continue;
                }

                currDist = interface->SpatialLocalizerI_giveClosestPoint(el_lcoords, el_coords, coords);
                if ( answer == NULL || (currDist < dist && currDist >= 0.0) ) {
                    answer = ielemptr;
                    lcoords = el_lcoords;
                    closest = el_coords;
                    dist = currDist;
                    if ( dist == 0.0 ) {
                        break;
                    }
                }
            }
        }
    }

    return answer;
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


void
DummySpatialLocalizer :: giveAllNodesWithinBox(nodeContainerType &nodeSet, const FloatArray &coords, const double radius)
{
    int i, nnode;
    DofManager *idofman;
    Node *inode;

    nnode = this->giveDomain()->giveNumberOfDofManagers();
    for ( i = 1; i <= nnode; i++ ) {
        idofman = this->giveDomain()->giveDofManager(i);
        if ( ( inode = dynamic_cast< Node * >(idofman) ) != NULL ) {
            if ( coords.distance( inode->giveCoordinates() ) <= radius ) {
                nodeSet.push_back(i);
            }
        }
    }
}
} // end namespace oofem
