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

#include "dummylocalizer.h"
#include "element.h"
#include "domain.h"
#include "integrationrule.h"
#include "gausspoint.h"
#include "node.h"

namespace oofem {
int
DummySpatialLocalizer :: init(bool force)
{
    if ( !force && this->initialized ) {
        return true;
    }

    // Count the elements in each cross section.
    int r;
    int nregion = this->domain->giveNumberOfRegions();
    int nelem = this->domain->giveNumberOfElements();

    IntArray region_nelem(nregion);
    region_nelem.zero();
    for ( auto &elem : domain->giveElements() ) {
        r = elem->giveRegionNumber();
        region_nelem.at(r)++;
    }

    this->region_elements.resize(nregion);
    // Creates a new int array of correct size for each region
    for ( int i = 1; i <= nregion; i++ ) {
        this->region_elements [ i - 1 ].resize( region_nelem.at(i) );
    }
    // Add the numbers into the list.
    IntArray c(nregion);
    c.zero();
    for ( int i = 1; i <= nelem; i++ ) {
        Element *e = this->domain->giveElement(i);
        r = e->giveRegionNumber();
        c.at(r)++;
        this->region_elements [ r - 1 ].at( c.at(r) ) = i;
    }
    return this->initialized = true;
}

Element *
DummySpatialLocalizer :: giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList)
{
    for ( auto &elem : domain->giveElements() ) {
        SpatialLocalizerInterface *interface = static_cast< SpatialLocalizerInterface * >( elem->giveInterface(SpatialLocalizerInterfaceType) );
        if ( interface ) {
            if ( regionList && ( regionList->findFirstIndexOf( elem->giveRegionNumber() ) == 0 ) ) {
                continue;
            }

            if ( interface->SpatialLocalizerI_BBoxContainsPoint(coords) == 0 ) {
                continue;
            }

            if ( interface->SpatialLocalizerI_containsPoint(coords) ) {
                return elem.get();
            }
        }
    }

    return NULL;
}


Element *
DummySpatialLocalizer :: giveElementClosestToPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &coords, int region)
{
    int nelems;
    Element *answer = NULL;
    double dist = 0.0;
    FloatArray el_coords, el_lcoords;
    lcoords.clear();
    closest.clear();

    if ( region > 0 ) {
        IntArray &elems = this->region_elements [ region - 1 ];
        for ( int ielem = 1; ielem <= elems.giveSize(); ielem++ ) {
            Element *ielemptr = this->domain->giveElement( elems.at(ielem) );
            SpatialLocalizerInterface *interface = static_cast< SpatialLocalizerInterface * >( ielemptr->giveInterface(SpatialLocalizerInterfaceType) );
            if ( interface ) {
                double currDist = interface->SpatialLocalizerI_giveClosestPoint(el_lcoords, el_coords, coords);
                if ( answer == NULL || ( currDist < dist && currDist >= 0.0 ) ) {
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
            Element *ielemptr = this->giveDomain()->giveElement(ielem);
            SpatialLocalizerInterface *interface = static_cast< SpatialLocalizerInterface * >( ielemptr->giveInterface(SpatialLocalizerInterfaceType) );
            if ( interface ) {
                if ( region > 0 && region == ielemptr->giveRegionNumber() ) {
                    continue;
                }

                double currDist = interface->SpatialLocalizerI_giveClosestPoint(el_lcoords, el_coords, coords);
                if ( answer == NULL || ( currDist < dist && currDist >= 0.0 ) ) {
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
DummySpatialLocalizer :: giveClosestIP(const FloatArray &coords, int region, bool iCohesiveZoneGP)
{
    double minDist = 0.0;
    GaussPoint *answer = NULL;
    FloatArray jGpCoords;

    for ( auto &elem : this->giveDomain()->giveElements() ) {
        if ( ( region < 0 ) || ( region == elem->giveRegionNumber() ) ) {
            for ( GaussPoint *jGp: *elem->giveDefaultIntegrationRulePtr() ) {
                if ( elem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                    double distance = coords.distance(jGpCoords);
                    if ( answer == NULL || distance < minDist ) {
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
    int nelem;
    FloatArray jGpCoords;

    nelem = this->giveDomain()->giveNumberOfElements();
    for ( int i = 1; i <= nelem; i++ ) {
        Element *ielem = this->giveDomain()->giveElement(i);
        for ( GaussPoint *jGp: *ielem->giveDefaultIntegrationRulePtr() ) {
            if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                double currDist = coords.distance(jGpCoords);
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
    int nnode = this->giveDomain()->giveNumberOfDofManagers();
    for ( int i = 1; i <= nnode; i++ ) {
        DofManager *idofman = this->giveDomain()->giveDofManager(i);
        Node *inode = dynamic_cast< Node * >(idofman);
        if ( inode != NULL ) {
            if ( coords.distance( inode->giveCoordinates() ) <= radius ) {
                nodeSet.push_back(i);
            }
        }
    }
}


Node *
DummySpatialLocalizer :: giveNodeClosestToPoint(const FloatArray &coords, double maxDist)
{
    Node *closest = NULL;
    double maxdist = 0.;
    for ( auto &dman : this->giveDomain()->giveDofManagers() ) {
        Node *node = dynamic_cast< Node* >( dman.get() );
        if ( node ) {
            double dist = coords.distance( node->giveCoordinates() );
            if ( closest == NULL || dist < maxdist ) {
                closest = node;
                maxdist = dist;
            }
        }
    }
    return maxdist > maxDist ? closest : NULL;
}

} // end namespace oofem
