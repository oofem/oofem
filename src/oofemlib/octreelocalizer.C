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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "octreelocalizer.h"
#include "element.h"
#include "domain.h"
#include "integrationrule.h"
#include "gausspnt.h"
#include "dofmanager.h"
#include "node.h"
#include "conTable.h"
#include "alist.h"
#include "mathfem.h"
#include "clock.h"

namespace oofem {
OctantRec :: OctantRec(OctreeSpatialLocalizer *loc, OctantRec *parent, FloatArray &origin, double halfWidth)
{
    this->localizer = loc;
    this->parent = parent;
    this->origin = origin;
    this->halfWidth = halfWidth;
    this->nodeList = NULL;
    this->elementIPList = NULL;
    this->depth = parent ? parent->giveCellDepth()+1 : 0;

    for ( int i = 0; i <= 1; i++ ) {
        for ( int j = 0; j <= 1; j++ ) {
            for ( int k = 0; k <= 1; k++ ) {
                this->child [ i ] [ j ] [ k ] = NULL;
            }
        }
    }
}

OctantRec :: ~OctantRec()
{
    // destructor - delete all receivers children recursively
    // delete children first, then itself
    for ( int i = 0; i <= 1; i++ ) {
        for ( int j = 0; j <= 1; j++ ) {
            for ( int k = 0; k <= 1; k++ ) {
                if ( this->child [ i ] [ j ] [ k ] ) {
                    delete this->child [ i ] [ j ] [ k ];
                }
            }
        }
    }

    if ( nodeList ) {
        delete nodeList;
    }

    if ( elementIPList ) {
        delete elementIPList;
    }
    elementList.clear();
}

std :: list< int > *
OctantRec :: giveNodeList()
{
    if ( nodeList == NULL ) {
        nodeList = new std :: list< int >;
    }

    return nodeList;
}

std :: set< int > *
OctantRec :: giveIPElementList()
{
    if ( elementIPList == NULL ) {
        elementIPList = new std :: set< int >;
    }

    return elementIPList;
}

std :: list< int > *
OctantRec :: giveElementList(int region)
{
    if ( elementList.giveSize() < region + 1 ) {
        elementList.growTo(region + 1);
    }
    std :: list < int > *list = elementList.at(region+1);
    if ( list == NULL ) {
        std :: list < int > *newlist = new std :: list < int >();
        elementList.put(region+1, newlist);
        return newlist;
    } else {
        return list;
    }
}


OctantRec *
OctantRec :: giveChild(int xi, int yi, int zi)
{
    if ( ( xi >= 0 ) && ( xi < 2 ) && ( yi >= 0 ) && ( yi < 2 ) && ( zi >= 0 ) && ( zi < 2 ) ) {
        return this->child [ xi ] [ yi ] [ zi ];
    } else {
        OOFEM_ERROR4("OctantRec::giveChild invalid child index (%d,%d,%d)", xi, yi, zi);
    }

    return NULL;
}


OctantRec :: ChildStatus
OctantRec :: giveChildContainingPoint(OctantRec **child, const FloatArray &coords)
{
    IntArray ind(3);

    if ( this->isTerminalOctant() ) {
        * child = NULL;
        return CS_NoChild;
    }

    for ( int i = 1; i <= coords.giveSize(); ++i ) {
        ind.at(i) = localizer->giveOctreeMaskValue(i) && coords.at(i) > this->origin.at(i);
    }

    * child = this->child [ ind.at(1) ] [ ind.at(2) ] [ ind.at(3) ];
    return CS_ChildFound;
}


bool
OctantRec :: isTerminalOctant()
{
    // This implementation is relying on fact, that child [0][0][0]
    // is created even for all degenerated trees
    if ( this->child [ 0 ] [ 0 ] [ 0 ] ) {
        return false;
    }

    return true;
}


void
OctantRec :: divideLocally(int level, const IntArray &octantMask)
{
    if ( this->isTerminalOctant() ) {
        // create corresponding child octants
        FloatArray childOrigin(3);

        for ( int i = 0; i <= octantMask.at(1); i++ ) {
            for ( int j = 0; j <= octantMask.at(2); j++ ) {
                for ( int k = 0; k <= octantMask.at(3); k++ ) {
                    childOrigin.at(1) = this->origin.at(1) + (i - 0.5) * this->halfWidth * octantMask.at(1);
                    childOrigin.at(2) = this->origin.at(2) + (j - 0.5) * this->halfWidth * octantMask.at(2);
                    childOrigin.at(3) = this->origin.at(3) + (k - 0.5) * this->halfWidth * octantMask.at(3);
                    this->child [ i ] [ j ] [ k ] = new OctantRec(localizer, this, childOrigin, this->halfWidth * 0.5);
                }
            }
        }
    }

    int newLevel = level - 1;
    if ( newLevel > 0 ) {
        // propagate message to children recursively with level param decreased
        for ( int i = 0; i <= octantMask.at(1); i++ ) {
            for ( int j = 0; j <= octantMask.at(2); j++ ) {
                for ( int k = 0; k <= octantMask.at(3); k++ ) {
                    if ( this->child [ i ] [ j ] [ k ] ) {
                        this->child [ i ] [ j ] [ k ]->divideLocally(newLevel, octantMask);
                    }
                }
            }
        }
    }
}


OctantRec :: BoundingBoxStatus
OctantRec :: testBoundingBox(const FloatArray &coords, double radius)
{
    int size = coords.giveSize();
    double bb0, bb1;
    double oct0, oct1;
    bool bbInside = true;

    for ( int i = 1; i <= size; i++ ) {
        if ( localizer->giveOctreeMaskValue(i) ) {
            bb0 = coords.at(i) - radius;
            bb1 = coords.at(i) + radius;
            oct0 = this->origin.at(i) - this->halfWidth;
            oct1 = this->origin.at(i) + this->halfWidth;
            
            if ( oct1 < bb0 || oct0 > bb1 ) { // Then its definitely outside, no need to go on
                return BBS_OutsideCell;
            }
            // Other cases; We first check if bb is completely inside the cell, in every spatial dimension
            bbInside = bbInside && oct0 < bb0 && oct1 > bb1;
        }
    }
    return bbInside ? BBS_InsideCell : BBS_ContainsCell;
}


void OctantRec :: printYourself()
{
    if (this->isTerminalOctant()) {
        printf(" center = {%e, %e, %e} size = %f, nodes = %d, elem_ips = %d\n", 
                this->origin.at(1), this->origin.at(2), this->origin.at(3), 
                this->halfWidth*2., this->nodeList->size(), this->elementIPList->size());
    } else {
        if (this->depth == 0) {
            printf("*ROOTCELL*");
        }
        printf("\n");
        for ( int i = 0; i <= 1; i++ ) {
            for ( int j = 0; j <= 1; j++ ) {
                for ( int k = 0; k <= 1; k++ ) {
                    if ( this->child [ i ] [ j ] [ k ] ) {
                        for (int q = 0; q < this->depth - 1; q++) printf("  ");
                        printf("+");
                        this->child [ i ] [ j ] [ k ]->printYourself();
                    }
                }
            }
        }
    }
}


OctantRec *
OctreeSpatialLocalizer :: findTerminalContaining(OctantRec *startCell, const FloatArray &coords)
{
    OctantRec :: ChildStatus result;
    OctantRec *currCell = startCell;
    // found terminal octant containing node
    while ( !currCell->isTerminalOctant() ) {
        result = currCell->giveChildContainingPoint(& currCell, coords);
        if ( result == OctantRec :: CS_NoChild ) {
            _error("findTerminalContaining: internal error - octree inconsistency");
        }
    }
    return currCell;
}


bool
OctreeSpatialLocalizer :: buildOctreeDataStructure()
{
    OOFEM_LOG_INFO("Initializing Octree structure\n");
    int init = 1, nnode = this->domain->giveNumberOfDofManagers();
    double rootSize, resolutionLimit;
    FloatArray minc(3), maxc(3), * coords;
    DofManager *dman;

    // test if tree already built
    if ( rootCell ) {
        return true;
    }

    this->elementListsInitialized.resize(this->domain->giveNumberOfRegions()+1);
    this->elementListsInitialized.zero();
    this->elementIPListsInitialized = false;

    // measure time consumed by octree build phase
    oofem_timeval ut, tstart;
    getUtime(tstart);

    // first determine domain extends (bounding box), and check for degenerated domain type
    for ( int i = 1; i <= nnode; i++ ) {
        dman = domain->giveDofManager(i);
        if ( ( dman->giveClassID() == NodeClass ) || ( dman->giveClassID() == RigidArmNodeClass ) ) {
            coords = ( ( ( Node * ) dman )->giveCoordinates() );
            if ( init ) {
                init = 0;
                for ( int j = 1; j <= coords->giveSize(); j++ ) {
                    minc.at(j) = maxc.at(j) = coords->at(j);
                }
            } else {
                for ( int j = 1; j <= coords->giveSize(); j++ ) {
                    if ( coords->at(j) < minc.at(j) ) {
                        minc.at(j) = coords->at(j);
                    }

                    if ( coords->at(j) > maxc.at(j) ) {
                        maxc.at(j) = coords->at(j);
                    }
                }
            }
        }
    } // end loop over nodes

    // determine root size
    rootSize = 0.0;
    for ( int i = 1; i <= 3; i++ ) {
        rootSize = 1.000001 * max( rootSize, maxc.at(i) - minc.at(i) );
    }

    // check for degenerated domain
    resolutionLimit = min(1.e-3, rootSize / 1.e6);
    for ( int i = 1; i <= 3; i++ ) {
        if ( ( maxc.at(i) - minc.at(i) ) > resolutionLimit ) {
            this->octreeMask.at(i) = 1;
        } else {
            this->octreeMask.at(i) = 0;
        }
    }

    // Create root Octant
    FloatArray center = minc;
    center.add(maxc);
    center.times(0.5);
    this->rootCell = new OctantRec(this, NULL, center, rootSize*0.5);

    // Build octree tree
    if ( nnode > OCTREE_MAX_NODES_LIMIT ) {
        this->rootCell->divideLocally(1, this->octreeMask);
    }

    // insert domain nodes into tree
    for ( int i = 1; i <= nnode; i++ ) {
        dman = domain->giveDofManager(i);
        if ( ( dman->giveClassID() == NodeClass ) || ( dman->giveClassID() == RigidArmNodeClass ) ) {
            coords = ( ( ( Node * ) dman )->giveCoordinates() );
            this->insertNodeIntoOctree(this->rootCell, i, * coords);
        }
    }

    getRelativeUtime(ut, tstart);

    // compute max. tree depth
    int treeDepth = 0;
    this->giveMaxTreeDepthFrom(this->rootCell, treeDepth);
    // compute processor time used by the program
    double nsec = ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM );
    OOFEM_LOG_DEBUG("Octree init [depth %d in %.2fs]\n", treeDepth, nsec);

    return true;
}


void
OctreeSpatialLocalizer :: initElementIPDataStructure()
{
    // Original implementation
    //
    int nip;
    int nelems = this->domain->giveNumberOfElements();
    Element *ielem;
    IntegrationRule *iRule;
    GaussPoint *jGp;
    FloatArray jGpCoords;

    if ( this->elementIPListsInitialized ) {
        return;
    }

    // insert IP records into tree (the tree topology is determined by nodes)
    for ( int i = 1; i <= nelems; i++ ) {
        // only default IP are taken into account
        ielem = this->giveDomain()->giveElement(i);
        if (ielem->giveNumberOfIntegrationRules() > 0) {
            iRule = ielem->giveDefaultIntegrationRulePtr();
            nip = iRule->getNumberOfIntegrationPoints();
            if ( nip ) {
                for ( int j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
                    jGp = iRule->getIntegrationPoint(j);
                    if ( ielem->computeGlobalCoordinates( jGpCoords, * ( jGp->giveCoordinates() ) ) ) {
                        this->insertIPElementIntoOctree(this->rootCell, i, jGpCoords);
                    } else {
                        _error("initElementIPDataStructure: computeGlobalCoordinates failed");
                    }
                }
                continue;
            }
        }
        // there are no IP (belonging to default integration rule of an element)
        // but the element should be present in octree data structure
        // this is needed by some services (giveElementContainingPoint, for example)
        FloatArray *nc;
        for ( int j = 1; j <= ielem->giveNumberOfNodes(); j++ ) {
            nc = ielem->giveNode(j)->giveCoordinates();
            this->insertIPElementIntoOctree(this->rootCell, i, * nc);
        }
    }

    // Initializes the element lists  in octree data structure.
    // This implementation requires that the list of nodes in terminate cells exists
    // simply all shared elements to nodes in terminal cell are added.
    // If this is added to existing implementation based on adding elements only if integration point is in the cell
    // This leads to more complete element list in terminal cell.
    // Can improve the searching for background element, but will deteriorate
    // the performance of closest IP search and search of elements in given volume.

    // Note: since in general, the integration point of an element may fall into
    // an octant, where are not the element nodes, the original algorithm is
    // necessary.

    //this->insertElementsUsingNodalConnectivitiesIntoOctree (this->rootCell);
    this->elementIPListsInitialized = true;
}


void
OctreeSpatialLocalizer :: insertIPElementIntoOctree(OctantRec *rootCell, int elemNum, const FloatArray &coords)
{
    OctantRec *currCell;
    // found terminal octant containing IP
    currCell = this->findTerminalContaining(rootCell, coords);
    currCell->addElementIP(elemNum);
}


void
OctreeSpatialLocalizer :: initElementDataStructure(int region)
{
    Element *ielem;
    SpatialLocalizerInterface *interface;
    FloatArray b0, b1;

    this->init();
    if ( this->elementListsInitialized.giveSize() >= region + 1 && this->elementListsInitialized(region) ) {
        return;
    }

    for ( int i = 1; i <= this->domain->giveNumberOfElements(); i++ ) {
        ielem = this->giveDomain()->giveElement(i);
        if ( ielem->giveRegionNumber() == region || region == 0 ) {
            interface = ( SpatialLocalizerInterface * ) ielem->giveInterface(SpatialLocalizerInterfaceType);
            if (interface) {
                interface->SpatialLocalizerI_giveBBox(b0, b1);
                this->insertElementIntoOctree(this->rootCell, region, i, b0, b1);
            }
        }
    }
    this->elementListsInitialized(region) = true;
}


void
OctreeSpatialLocalizer :: insertElementIntoOctree(OctantRec *rootCell, int region, int elemNum, const FloatArray &b0, const FloatArray &b1)
{
    // Compare the bounding box corners to the center to determine which region is overlaps
    // Checks: b0 <= center, b1 >= center for each entry.
    // This is bundled into an array for more convenient code in the loop. 
    IntArray bbc[2] = { IntArray(3), IntArray(3) };
    FloatArray origin;
    rootCell->giveOrigin(origin);
    for ( int i = 1; i <= b0.giveSize(); i++ ) {
        if ( this->octreeMask.at(i) ) {
            bbc[0].at(i) = b0.at(i) <= origin.at(i);
            bbc[1].at(i) = b1.at(i) >= origin.at(i);
        } else {
            bbc[0].at(i) = bbc[1].at(i) = true;
        }
    }
    for ( int i = b0.giveSize()+1; i <= 3; i++ ) {
        bbc[0].at(i) = bbc[1].at(i) = true;
    }
    // Check terminal or recurse
    if ( rootCell->isTerminalOctant() ) {
        rootCell->addElement( region, elemNum );
    } else {
        // The for loops below check if the bounding box overlaps the region before or after the cell center 
        // e.g. when i = 0, then we check for the children that have x-coordinate <= than root center
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            if ( bbc[i].at(1) ) {
                for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                    if ( bbc[j].at(2) ) {
                        for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                            if ( bbc[k].at(3) && rootCell->giveChild(i, j, k) ) {
                                this->insertElementIntoOctree( rootCell->giveChild(i, j, k), region, elemNum, b0, b1 );
                            }
                        }
                    }
                }
            }
        }
    }
}


void
OctreeSpatialLocalizer :: insertElementsUsingNodalConnectivitiesIntoOctree(OctantRec *currCell)
{
    if ( currCell->isTerminalOctant() ) {
        nodeContainerType *cellNodes = currCell->giveNodeList();
        nodeContainerType :: iterator pos;
        const IntArray *dofmanConnectivity;
        int size;

        if ( !cellNodes->empty() ) {
            for ( pos = cellNodes->begin(); pos != cellNodes->end(); ++pos ) {
                // loop over cell nodes and ask connectivity table for shared elements
                dofmanConnectivity = domain->giveConnectivityTable()->giveDofManConnectivityArray(* pos);
                if ( dofmanConnectivity ) {
                    size = dofmanConnectivity->giveSize();
                    for ( int i = 1; i <= size; i++ ) {
                        currCell->addElementIP( dofmanConnectivity->at(i) );
                    }
                }
            }
        }
    } else {
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    if ( currCell->giveChild(i, j, k) ) {
                        this->insertElementsUsingNodalConnectivitiesIntoOctree( currCell->giveChild(i, j, k) );
                    }
                }
            }
        }
    }
}


void
OctreeSpatialLocalizer :: insertNodeIntoOctree(OctantRec *rootCell, int nodeNum, const FloatArray &coords)
{
    int nCellItems, cellDepth;
    OctantRec :: ChildStatus result;
    OctantRec *currCell;
    nodeContainerType *cellNodeList;
    nodeContainerType :: const_iterator pos;
    DofManager *dman;
    FloatArray *nodeCoords;

    // found terminal octant containing node
    currCell = this->findTerminalContaining(rootCell, coords);
    // request cell node list
    cellNodeList = currCell->giveNodeList();
    nCellItems = cellNodeList->size();
    cellDepth  = currCell->giveCellDepth();
    // check for refinement criteria
    // should also include max refinement level criteria
    //
    if ( ( nCellItems > OCTREE_MAX_NODES_LIMIT ) && ( cellDepth <= OCTREE_MAX_DEPTH ) ) {
        // refine tree one level
        currCell->divideLocally(1, this->octreeMask);
        // propagate all nodes already assigned to currCell to children
        for ( pos = cellNodeList->begin(); pos != cellNodeList->end(); ++pos ) {
            dman = domain->giveDofManager(* pos);
            nodeCoords = ( ( ( Node * ) dman )->giveCoordinates() );
            this->insertNodeIntoOctree(currCell, * pos, * nodeCoords);
        }

        // remove node list at relative root cell
        currCell->deleteNodeList();
        // now the new node record is saved simply to one of created children.
        // Generally, one has to check if child is candidate for further refinement
        // (case, when all parent nodes are belonging to particular child), refine it
        // and check children again. This can be implemented by recursive procedure, but
        // probability of this case is relatively small.
        // Current implementation simply insert node to child created in first refinement,
        // which can lead ne node violation of refinement criteria.
        // But implementation is simple and effective.
        // Please note, that the next insertion in particular cell cause the refinement.

        // find child containing new node
        result = currCell->giveChildContainingPoint(& currCell, coords);
        if ( result != OctantRec :: CS_ChildFound ) {
            _error("insertNodeIntoOctree: internal error - octree inconsistency");
        }
    }

    currCell->addNode(nodeNum);
}


Element *
OctreeSpatialLocalizer :: giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList)
{
    Element *answer;
    //SpatialLocalizerInterface* interface;
    OctantRec *currCell, *childCell = NULL;

    this->init();
    this->initElementIPDataStructure();

    // found terminal octant containing point
    currCell = this->findTerminalContaining(rootCell, coords);

    while ( currCell != NULL ) {
        // loop over all elements in currCell, skip search on child cell already scanned
        answer = this->giveElementContainingPoint(currCell, coords, childCell, regionList);
        if ( answer ) {
            return answer;
        }

        childCell = currCell;
        // terminal cell does not contain node and its connected elements containing the given point
        // search at parent level
        currCell = currCell->giveParent();
    }

    // there isn't any element containing point.
    return NULL;
}


Element *
OctreeSpatialLocalizer :: giveElementContainingPoint(OctantRec *cell, const FloatArray &coords,
                                                     OctantRec *scannedChild, const IntArray *regionList)
{
    elementContainerType *elementList = cell->giveIPElementList();

    // recursive implementation
    if ( cell->isTerminalOctant() && ( !elementList->empty() ) ) {
        int ielem;
        Element *ielemptr;
        SpatialLocalizerInterface *interface;
        elementContainerType :: iterator pos;

        for ( pos = elementList->begin(); pos !=  elementList->end(); ++pos ) {
            ielem = * pos;
            ielemptr = this->giveDomain()->giveElement(ielem);

            /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
            if ( ielemptr->giveParallelMode() == Element_remote ) {
                continue;
            }
#endif

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
    } else {
        Element *answer;
        // receiver is not terminal octant -> call same service on childs
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    if ( cell->giveChild(i, j, k) ) {
                        if ( scannedChild == cell->giveChild(i, j, k) ) {
                            continue;
                        }

                        answer = this->giveElementContainingPoint(cell->giveChild(i, j, k), coords, NULL, regionList);
                        if ( answer ) {
                            return answer;
                        }
                    }
                }
            }
        }
    }

    return NULL;
}


//
//  WARNING !! LIMITED IMPLEMENTATION HERE !!
//
//  a close element in the same region (or in other if check using region list is not applied)
//  but topologically of large distance from the point may be returned;
//  maybe nonlocal barrier should be used
//

Element *
OctreeSpatialLocalizer :: giveElementCloseToPoint(const FloatArray &coords, const IntArray *regionList)
{
    Element *ielemptr, *answer = NULL;
    double currDist, minDist = 1.1 * rootCell->giveWidth();
    elementContainerType :: iterator pos;
    elementContainerType *elementList;
    OctantRec *currCell;

    this->init();
    this->initElementIPDataStructure();

    // found terminal octant containing point
    currCell = this->findTerminalContaining(rootCell, coords);

    elementList = currCell->giveIPElementList();
    if ( !elementList->empty() ) {
        SpatialLocalizerInterface *interface;

        for ( pos = elementList->begin(); pos !=  elementList->end(); ++pos ) {
            ielemptr = this->giveDomain()->giveElement(* pos);

            /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
            if ( ielemptr->giveParallelMode() == Element_remote ) {
                continue;
            }
#endif

            interface = ( SpatialLocalizerInterface * ) ielemptr->giveInterface(SpatialLocalizerInterfaceType);
            if ( interface ) {
                if ( regionList && ( regionList->findFirstIndexOf( ielemptr->giveRegionNumber() ) == 0 ) ) {
                    continue;
                }

                currDist = interface->SpatialLocalizerI_giveDistanceFromParametricCenter(coords);
                if ( currDist < minDist ) {
                    answer = ielemptr;
                    minDist = currDist;
                }
            } else {
                _error("giveElementCloseToPoint: no interface");
            }
        }

        if ( minDist == 0.0 || currCell == rootCell ) {
            return answer;
        }
    }

    OctantRec :: BoundingBoxStatus BBStatus;
    FloatArray bborigin = coords;

    // construct bounding box and test its position within currCell
    BBStatus = currCell->testBoundingBox(bborigin, minDist);
    if ( BBStatus == OctantRec :: BBS_InsideCell ) {
        return answer;
    } else if ( BBStatus == OctantRec :: BBS_ContainsCell ) {
        std :: list< OctantRec * >cellList;
        std :: list< OctantRec * > :: iterator cellListIt;

        // found terminal octant containing point
        OctantRec *startCell = this->findTerminalContaining(rootCell, coords);
        if ( startCell == NULL ) {
            startCell = rootCell;
        }

        // go up, until cell containing bbox is found
        if ( startCell != rootCell ) {
            while ( startCell->testBoundingBox(coords, minDist) != OctantRec :: BBS_InsideCell ) {
                startCell = startCell->giveParent();
                if ( startCell == rootCell ) {
                    break;
                }
            }
        }

        this->giveListOfTerminalCellsInBoundingBox(cellList, bborigin, minDist, 0, startCell);

        for ( cellListIt = cellList.begin(); cellListIt != cellList.end(); ++cellListIt ) {
            if ( currCell == ( * cellListIt ) ) {
                continue;
            }

            this->giveElementCloseToPointWithinOctant( ( * cellListIt ), coords, minDist, & answer, regionList );
        }
    }

    return answer;
}


Element *
OctreeSpatialLocalizer :: giveElementClosestToPoint(FloatArray &lcoords, FloatArray &closest,
                                                    const FloatArray &gcoords, int region)
{
    Element *answer = NULL;
    std :: list< OctantRec * > cellList;
    OctantRec *currCell;
    FloatArray currLcoords;
    FloatArray currClosest;
    double radius, prevRadius;
    FloatArray c;

    this->rootCell->giveOrigin(c);
    // Maximum distance given coordinate and furthest terminal cell ( center_distance + width/2*sqrt(3) )
    double minDist = c.distance(gcoords) + this->rootCell->giveWidth() * 0.87;

    this->initElementDataStructure(region);

    // found terminal octant containing point
    currCell = this->findTerminalContaining(rootCell, gcoords);

    // Look in center, then expand.
    this->giveElementClosestToPointWithinOctant(currCell, gcoords, minDist, lcoords, closest, answer, region);
    prevRadius = 0.;
    radius = currCell->giveWidth();
    while (radius < minDist) {
        this->giveListOfTerminalCellsInBoundingBox(cellList, gcoords, radius, prevRadius, this->rootCell);
        for (std :: list< OctantRec * > :: iterator it = cellList.begin(); it != cellList.end(); ++it) {
            this->giveElementClosestToPointWithinOctant(*it, gcoords, minDist, lcoords, closest, answer, region);
        }
        prevRadius = radius;
        radius *= 2.; // Keep expanding the scope until the radius is larger than the entire root cell, then give up (because we have checked every possible cell)
    }
    return answer;
}


void
OctreeSpatialLocalizer :: giveElementClosestToPointWithinOctant(OctantRec *currCell, const FloatArray &gcoords,
        double &minDist, FloatArray &lcoords, FloatArray &closest, Element *&answer, int region)
{
    SpatialLocalizerInterface *interface;
    Element *ielemptr;

    double currDist;
    std :: list< int > :: iterator pos;
    std :: list< int > *elementList;
    FloatArray currLcoords;
    FloatArray currClosest;

    elementList = currCell->giveElementList(region);
    if ( !elementList->empty() ) {
        for ( pos = elementList->begin(); pos !=  elementList->end(); ++pos ) {
            ielemptr = this->giveDomain()->giveElement(* pos);

#ifdef __PARALLEL_MODE
            if ( ielemptr->giveParallelMode() == Element_remote ) {
                continue;
            }
#endif
            interface = ( SpatialLocalizerInterface * ) ielemptr->giveInterface(SpatialLocalizerInterfaceType);
            if ( region > 0 && ielemptr->giveRegionNumber() != region ) {
                continue;
            }
            currDist = interface->SpatialLocalizerI_giveClosestPoint(currLcoords, currClosest, gcoords);
            if ( currDist < minDist ) {
                lcoords = currLcoords;
                closest = currClosest;
                answer = ielemptr;
                minDist = currDist;
            }
        }
        if ( minDist == 0.0 ) {
            return;
        }
    }
}


void
OctreeSpatialLocalizer :: giveElementCloseToPointWithinOctant(OctantRec *cell, const FloatArray &coords,
                                                              double &minDist, Element **answer,
                                                              const IntArray *regionList)
{
    if ( cell->isTerminalOctant() ) {
        Element *ielemptr;
        SpatialLocalizerInterface *interface;
        elementContainerType *elementList = cell->giveIPElementList();
        elementContainerType :: iterator pos;
        double currDist;

        if ( !elementList->empty() ) {
            for ( pos = elementList->begin(); pos !=  elementList->end(); ++pos ) {
                ielemptr = this->giveDomain()->giveElement(* pos);

                /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                if ( ielemptr->giveParallelMode() == Element_remote ) {
                    continue;
                }

#endif

                interface = ( SpatialLocalizerInterface * ) ielemptr->giveInterface(SpatialLocalizerInterfaceType);
                if ( interface ) {
                    if ( regionList && ( regionList->findFirstIndexOf( ielemptr->giveRegionNumber() ) == 0 ) ) {
                        continue;
                    }

                    currDist = interface->SpatialLocalizerI_giveDistanceFromParametricCenter(coords);
                    if ( currDist < minDist ) {
                        * answer = ielemptr;
                        minDist = currDist;
                    }
                } else {
                    _error("giveElementCloseToPointWithinOctant: no interface");
                }
            }
        }

        return;
    } else {
        OctantRec :: BoundingBoxStatus BBStatus;

        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    if ( cell->giveChild(i, j, k) ) {
                        // test if box hits the cell
                        BBStatus = cell->giveChild(i, j, k)->testBoundingBox(coords, minDist);
                        if ( ( BBStatus == OctantRec :: BBS_InsideCell ) || ( BBStatus == OctantRec :: BBS_ContainsCell ) ) {
                            // if yes call this method for such cell
                            this->giveElementCloseToPointWithinOctant(cell->giveChild(i, j, k), coords, minDist, answer, regionList);
                        }
                    }
                }
            }
        }
    }
}


GaussPoint *
OctreeSpatialLocalizer :: giveClosestIP(const FloatArray &coords, int region)
{
    int j;
    double dist, minDist;
    OctantRec *currCell;
    GaussPoint *nearestGp, *jGp;
    OctantRec :: BoundingBoxStatus BBStatus;
    elementContainerType :: iterator pos;
    elementContainerType *elementList;
    Element *ielem;
    FloatArray jGpCoords;
    IntegrationRule *iRule;

    this->init();
    this->initElementIPDataStructure();

    // list of already tested elements
    // elementContainerType visitedElems;
    minDist = 1.1 * rootCell->giveWidth();
    // found terminal octant containing point
    currCell = this->findTerminalContaining(rootCell, coords);
    elementList = currCell->giveIPElementList();
    // find nearest ip in this terminal cell
    if ( !elementList->empty() ) {
        for ( pos = elementList->begin(); pos !=  elementList->end(); ++pos ) {
            ielem = domain->giveElement(* pos);

            /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
            if ( ielem->giveParallelMode() == Element_remote ) {
                continue;
            }

#endif

            //igp   = ipContainer[i].giveGp();
            if ( ( region > 0 ) && ( region != ielem->giveRegionNumber() ) ) {
                continue;
            }

            // test if element already visited
            // if (!visitedElems.insert(*pos).second) continue;
            iRule = ielem->giveDefaultIntegrationRulePtr();
            for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
                jGp = iRule->getIntegrationPoint(j);
                if ( ielem->computeGlobalCoordinates( jGpCoords, * ( jGp->giveCoordinates() ) ) ) {
                    // compute distance
                    dist = coords.distance(jGpCoords);
                    if ( dist < minDist ) {
                        minDist   = dist;
                        nearestGp = jGp;
                    }
                } else {
                    _error("giveClosestIP: computeGlobalCoordinates failed");
                }
            }
        }
    }

    FloatArray bborigin = coords;
    // all cell element ip's scanned
    // construct bounding box and test its position within currCell
    BBStatus = currCell->testBoundingBox(bborigin, minDist);
    if ( BBStatus == OctantRec :: BBS_InsideCell ) {
        return nearestGp;
    } else if ( BBStatus == OctantRec :: BBS_ContainsCell ) {
        std :: list< OctantRec * >cellList;
        std :: list< OctantRec * > :: iterator cellListIt;

        // found terminal octant containing point
        OctantRec *startCell = this->findTerminalContaining(rootCell, coords);
        // go up, until cell containing bbox is found
        if ( startCell != rootCell ) {
            while ( startCell->testBoundingBox(coords, minDist) != OctantRec :: BBS_InsideCell ) {
                startCell = startCell->giveParent();
                if ( startCell == rootCell ) {
                    break;
                }
            }
        }

        this->giveListOfTerminalCellsInBoundingBox(cellList, bborigin, minDist, 0, startCell);

        for ( cellListIt = cellList.begin(); cellListIt != cellList.end(); ++cellListIt ) {
            if ( currCell == ( * cellListIt ) ) {
                continue;
            }

            this->giveClosestIPWithinOctant( ( * cellListIt ), coords, region, minDist, & nearestGp );
        }

        return nearestGp;

#if 0
        // found terminal octant containing point
        OctantRec* startCell = this->findTerminalContaining (rootCell, coords);
        // go up, until cell containing bbox is found
        if (startCell != rootCell) {
            while (startCell->testBoundingBox (coords, minDist) != OctantRec::BBS_InsideCell) {
                startCell = startCell->giveParent();
                if (startCell == rootCell) break;
            }
        }
        // loop over all child (if any) and found all nodes meeting the criteria
        // this->giveClosestIPWithinOctant (startCell, visitedElems, coords, region, minDist, &nearestGp);
        this->giveClosestIPWithinOctant (startCell, coords, region, minDist, &nearestGp);
        return nearestGp;
#endif


#if 0
        set <int> nearElementList;
        // find all elements within bbox and find nearest ip
        this->giveAllElementsWithIpWithinBox (nearElementList, bborigin, minDist);
        // loop over those elements and find nearest ip and return it
        // check for region also
        if (!nearElementList.empty()) {
            for (pos=nearElementList.begin(); pos !=  nearElementList.end(); ++pos) {
                ielem = domain->giveElement(*pos);
                //igp   = ipContainer[i].giveGp();
                if (region == ielem->giveRegionNumber()) {
                    iRule = ielem->giveDefaultIntegrationRulePtr ();
                    for (j=0 ; j < iRule->getNumberOfIntegrationPoints() ; j++) {
                        jGp = iRule->getIntegrationPoint(j) ;
                        if (ielem->computeGlobalCoordinates (jGpCoords, *(jGp->giveCoordinates()))) {
                            // compute distance
                            dist = coords.distance(jGpCoords);
                            if (dist < minDist) {
                                minDist   = dist;
                                nearestGp = jGp;
                            }
                        } else _error ("giveClosestIP: computeGlobalCoordinates failed");
                    }
                }
            }
        }
#endif
    } else {
        _error("giveClosestIP: octree inconsistency found");
    }

    return NULL;
}


void
OctreeSpatialLocalizer :: giveClosestIPWithinOctant(OctantRec *currentCell, //elementContainerType& visitedElems,
                                                    const FloatArray &coords,
                                                    int region, double &dist, GaussPoint **answer)
{
    if ( currentCell->isTerminalOctant() ) {
        double currDist;
        elementContainerType :: iterator pos;
        elementContainerType *elementList;
        //FloatArray *ipCoords;
        Element *ielem;
        IntegrationRule *iRule;
        FloatArray jGpCoords;

        // loop over cell elements and check if they meet the criteria
        elementList = currentCell->giveIPElementList();
        if ( !elementList->empty() ) {
            for ( pos = elementList->begin(); pos !=  elementList->end(); ++pos ) {
                // ask for element
                ielem = domain->giveElement(* pos);

                /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                if ( ielem->giveParallelMode() == Element_remote ) {
                    continue;
                }

#endif

                if ( ( region > 0 ) && ( region != ielem->giveRegionNumber() ) ) {
                    continue;
                }

                // test if element already visited
                // if (!visitedElems.insert(*pos).second) continue;
                // is one of his ip's  within given bbox -> inset it into elemSet
                iRule = ielem->giveDefaultIntegrationRulePtr();
                for ( int j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
                    if ( ielem->computeGlobalCoordinates( jGpCoords, * ( iRule->getIntegrationPoint(j)->giveCoordinates() ) ) ) {
                        currDist = coords.distance(jGpCoords);
                        // multiple insertion are handled by STL set implementation
                        if ( currDist <= dist ) {
                            dist = currDist;
                            * answer = iRule->getIntegrationPoint(j);
                        }
                    } else {
                        _error("giveClosestIPWithinOctant: computeGlobalCoordinates failed");
                    }
                }
            }
        }

        return;
    } else {
        OctantRec :: BoundingBoxStatus BBStatus;

        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    if ( currentCell->giveChild(i, j, k) ) {
                        // test if box hits the cell
                        BBStatus = currentCell->giveChild(i, j, k)->testBoundingBox(coords, dist);
                        if ( ( BBStatus == OctantRec :: BBS_InsideCell ) || ( BBStatus == OctantRec :: BBS_ContainsCell ) ) {
                            // if yes call this method for such cell
                            //this->giveClosestIPWithinOctant (currentCell->giveChild(i,j,k), visitedElems, coords, region, dist, answer);
                            this->giveClosestIPWithinOctant(currentCell->giveChild(i, j, k), coords, region, dist, answer);
                        }
                    }
                }
            }
        }
    }
}


void
OctreeSpatialLocalizer :: giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords,
                                                         const double radius)
{
    this->init();
    this->initElementIPDataStructure();
    // found terminal octant containing point
    OctantRec *currCell = this->findTerminalContaining(rootCell, coords);
    // go up, until cell containing bbox is found
    if ( currCell != rootCell ) {
        while ( currCell->testBoundingBox(coords, radius) != OctantRec :: BBS_InsideCell ) {
            currCell = currCell->giveParent();
            if ( currCell == rootCell ) {
                break;
            }
        }
    }

    // loop over all child (if any) and found all nodes meeting the criteria
    this->giveElementsWithIPWithinBox(elemSet, currCell, coords, radius);
    if ( elemSet.empty() ) {
        _error("OctreeSpatialLocalizer::giveAllElementsWithIpWithinBox empty set found");
    }
}


void
OctreeSpatialLocalizer :: giveElementsWithIPWithinBox(elementContainerType &elemSet, OctantRec *currentCell,
                                                      const FloatArray &coords, const double radius)
{
    if ( currentCell->isTerminalOctant() ) {
        double currDist;
        elementContainerType :: iterator pos;
        elementContainerType *elementList;
        Element *ielem;
        IntegrationRule *iRule;
        FloatArray jGpCoords;

        // loop over cell elements and check if they meet the criteria
        elementList = currentCell->giveIPElementList();
        if ( !elementList->empty() ) {
            for ( pos = elementList->begin(); pos !=  elementList->end(); ++pos ) {
                // test if element is already present
                if ( elemSet.find(* pos) != elemSet.end() ) {
                    continue;
                }

                // ask for element
                ielem = domain->giveElement(* pos);

                /* #ifdef __PARALLEL_MODE
                 *       if(ielem -> giveParallelMode() == Element_remote)continue;
                 *#endif
                 */
                // is one of his ip's  within given bbox -> inset it into elemSet
                iRule = ielem->giveDefaultIntegrationRulePtr();
                for ( int j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
                    if ( ielem->computeGlobalCoordinates( jGpCoords, * ( iRule->getIntegrationPoint(j)->giveCoordinates() ) ) ) {
                        currDist = coords.distance(jGpCoords);
                        // multiple insertion are handled by STL set implementation
                        if ( currDist <= radius ) {
                            elemSet.insert(* pos);
                        }
                    } else {
                        _error("giveElementsWithIPWithinBox: computeGlobalCoordinates failed");
                    }
                }
            }
        }
    } else {
        OctantRec :: BoundingBoxStatus BBStatus;


        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    if ( currentCell->giveChild(i, j, k) ) {
                        // test if box hits the cell
                        BBStatus = currentCell->giveChild(i, j, k)->testBoundingBox(coords, radius);
                        if ( ( BBStatus == OctantRec :: BBS_InsideCell ) || ( BBStatus == OctantRec :: BBS_ContainsCell ) ) {
                            // if yes call this method for such cell
                            this->giveElementsWithIPWithinBox(elemSet, currentCell->giveChild(i, j, k), coords, radius);
                        }
                    }
                }
            }
        }
    }
}


void
OctreeSpatialLocalizer :: giveAllNodesWithinBox(nodeContainerType &nodeSet, const FloatArray &coords, const double radius)
{
    this->init();
    // found terminal octant containing point
    OctantRec *currCell = this->findTerminalContaining(rootCell, coords);
    // go up, until cell containing bbox is found
    if ( currCell != rootCell ) {
        while ( currCell->testBoundingBox(coords, radius) != OctantRec :: BBS_InsideCell ) {
            currCell = currCell->giveParent();
            if ( currCell == rootCell ) {
                break;
            }
        }
    }

    // loop over all child (if any) and found all nodes meeting the criteria
    this->giveNodesWithinBox(nodeSet, currCell, coords, radius);
}


void
OctreeSpatialLocalizer :: giveNodesWithinBox(nodeContainerType &nodeList, OctantRec *currentCell,
                                             const FloatArray &coords, const double radius)
{
    nodeContainerType *cellNodes = currentCell->giveNodeList();
    nodeContainerType :: iterator pos;

    if ( currentCell->isTerminalOctant() ) {
        FloatArray *nodeCoords;
        if ( !cellNodes->empty() ) {
            for ( pos = cellNodes->begin(); pos != cellNodes->end(); ++pos ) {
                // loop over cell nodes and check if they meet the criteria
                nodeCoords = domain->giveNode(* pos)->giveCoordinates();
                // is node within bbox
                if ( nodeCoords->distance(coords) <= radius ) {
                    // if yes, append them into set
                    nodeList.push_back(* pos);
                }
            }
        }
    } else {
        OctantRec :: BoundingBoxStatus BBStatus;

        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    if ( currentCell->giveChild(i, j, k) ) {
                        // test if box hits the cell
                        BBStatus = currentCell->giveChild(i, j, k)->testBoundingBox(coords, radius);
                        if ( ( BBStatus == OctantRec :: BBS_InsideCell ) || ( BBStatus == OctantRec :: BBS_ContainsCell ) ) {
                            // if yes call this method for such cell
                            this->giveNodesWithinBox(nodeList, currentCell->giveChild(i, j, k), coords, radius);
                        }
                    }
                }
            }
        }
    }
}


void
OctreeSpatialLocalizer :: giveMaxTreeDepthFrom(OctantRec *root, int &maxDepth)
{
    int depth = root->giveCellDepth();
    maxDepth = max(maxDepth, depth);

    for ( int i = 0; i <= octreeMask.at(1); i++ ) {
        for ( int j = 0; j <= octreeMask.at(2); j++ ) {
            for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                if ( root->giveChild(i, j, k) ) {
                    this->giveMaxTreeDepthFrom(root->giveChild(i, j, k), maxDepth);
                }
            }
        }
    }
}


void
OctreeSpatialLocalizer :: giveListOfTerminalCellsInBoundingBox(std :: list< OctantRec * > &cellList, const FloatArray &coords,
                                                               double radius, double innerRadius, OctantRec *currentCell)
{
    OctantRec :: BoundingBoxStatus BBStatus;

    BBStatus = currentCell->testBoundingBox(coords, radius);
    if ( BBStatus != OctantRec :: BBS_OutsideCell ) {
        if ( currentCell->isTerminalOctant() ) {
            if ( currentCell->testBoundingBox(coords, innerRadius) == OctantRec :: BBS_OutsideCell ) {
                cellList.push_back(currentCell);
            }
        } else {
            for ( int i = 0; i <= octreeMask.at(1); i++ ) {
                for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                    for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                        if ( currentCell->giveChild(i, j, k) ) {
                            this->giveListOfTerminalCellsInBoundingBox( cellList, coords, radius, innerRadius, currentCell->giveChild(i, j, k) );
                        }
                    }
                }
            }
        }
    }
}


int
OctreeSpatialLocalizer :: init(bool force)
{
    if ( force ) {
        if ( rootCell ) {
            delete rootCell;
        }
        rootCell = NULL;
        elementIPListsInitialized = false;
        elementListsInitialized.zero();
    }

    if ( !rootCell ) {
        return this->buildOctreeDataStructure();
    } else {
        return 0;
    }
}
} // end namespace oofem
