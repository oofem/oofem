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

#include "octreelocalizer.h"
#include "element.h"
#include "domain.h"
#include "integrationrule.h"
#include "gausspoint.h"
#include "dofmanager.h"
#include "node.h"
#include "connectivitytable.h"
#include "mathfem.h"
#include "timer.h"
#include "error.h"
#include "xfem/xfemelementinterface.h"

#include <iostream>

namespace oofem {
OctantRec :: OctantRec(OctantRec *parent, FloatArray origin, double halfWidth) :
    parent(parent),
    origin(std::move(origin)),
    halfWidth(halfWidth)
{
    this->depth = parent ? parent->giveCellDepth() + 1 : 0;
}

std :: list< int > &
OctantRec :: giveNodeList()
{
    return nodeList;
}

IntArray &
OctantRec :: giveIPElementList()
{
    return elementIPList;
}

std :: list< int > &
OctantRec :: giveElementList(int region)
{
    if ( (int)elementList.size() < region + 1 ) {
        elementList.resize(region + 1);
    }
    return elementList[region];
}


OctantRec *
OctantRec :: giveChild(int xi, int yi, int zi)
{
    if ( ( xi >= 0 ) && ( xi < 2 ) && ( yi >= 0 ) && ( yi < 2 ) && ( zi >= 0 ) && ( zi < 2 ) ) {
        return this->child [ xi ] [ yi ] [ zi ].get();
        ///@todo We could inline the entire block if we do: std::unique_ptr<OctantRec[2][2][2]> child;
        //return &(*child) [ xi ] [ yi ] [ zi ];
    } else {
        OOFEM_ERROR("invalid child index (%d,%d,%d)", xi, yi, zi);
    }

    return nullptr;
}


OctantRec :: ChildStatus
OctantRec :: giveChildContainingPoint(OctantRec *&child, const FloatArray &coords, const IntArray &mask)
{
    if ( this->isTerminalOctant() ) {
        child = nullptr;
        return CS_NoChild;
    }

    IntArray ind(3);
    for ( int i = 0; i < coords.giveSize(); ++i ) {
        ind[i] = mask[i] && coords[i] > this->origin[i];
    }

    child = this->child [ ind[0] ] [ ind[1] ] [ ind[2] ].get();
    return CS_ChildFound;
}


bool
OctantRec :: isTerminalOctant()
{
    // This implementation is relying on fact, that child [0][0][0]
    // is created even for all degenerated trees
    return ! this->child [ 0 ] [ 0 ] [ 0 ];
}


void
OctantRec :: divideLocally(int level, const IntArray &mask)
{
    if ( this->isTerminalOctant() ) {
        // create corresponding child octants
        for ( int i = 0; i <= mask.at(1); i++ ) {
            for ( int j = 0; j <= mask.at(2); j++ ) {
                for ( int k = 0; k <= mask.at(3); k++ ) {
                    FloatArray childOrigin = {
                        this->origin.at(1) + ( i - 0.5 ) * this->halfWidth * mask.at(1),
                        this->origin.at(2) + ( j - 0.5 ) * this->halfWidth * mask.at(2),
                        this->origin.at(3) + ( k - 0.5 ) * this->halfWidth * mask.at(3)
                    };
                    this->child [ i ] [ j ] [ k ] = std::make_unique<OctantRec>(this, std::move(childOrigin), this->halfWidth * 0.5);
                }
            }
        }
    }

    int newLevel = level - 1;
    if ( newLevel > 0 ) {
        // propagate message to children recursively with level param decreased
        for ( int i = 0; i <= mask.at(1); i++ ) {
            for ( int j = 0; j <= mask.at(2); j++ ) {
                for ( int k = 0; k <= mask.at(3); k++ ) {
                    if ( this->child [ i ] [ j ] [ k ] ) {
                        this->child [ i ] [ j ] [ k ]->divideLocally(newLevel, mask);
                    }
                }
            }
        }
    }
}


OctantRec :: BoundingBoxStatus
OctantRec :: testBoundingBox(const FloatArray &coords, double radius, const IntArray &mask)
{
    bool bbInside = true;

    for ( int i = 1; i <= coords.giveSize(); i++ ) {
        if ( mask.at(i) ) {
            double bb0 = coords.at(i) - radius;
            double bb1 = coords.at(i) + radius;
            double oct0 = this->origin.at(i) - this->halfWidth;
            double oct1 = this->origin.at(i) + this->halfWidth;

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
    if ( this->isTerminalOctant() ) {
        std :: cout << " center = {" << this->origin.at(1) << "," << this->origin.at(2) << "," << this->origin.at(3)
                    << "} size = " << ( this->halfWidth * 2. ) << " nodes = " << this->nodeList.size() << " elem_ips = " << this->elementIPList.giveSize() << "\n";
    } else {
        if ( this->depth == 0 ) {
            printf("*ROOTCELL*");
        }
        printf("\n");
        for ( int i = 0; i <= 1; i++ ) {
            for ( int j = 0; j <= 1; j++ ) {
                for ( int k = 0; k <= 1; k++ ) {
                    if ( this->child [ i ] [ j ] [ k ] ) {
                        for ( int q = 0; q < this->depth - 1; q++ ) {
                            printf("  ");
                        }
                        printf("+");
                        this->child [ i ] [ j ] [ k ]->printYourself();
                    }
                }
            }
        }
    }
}


OctreeSpatialLocalizer :: OctreeSpatialLocalizer(Domain* d) : SpatialLocalizer(d),
    octreeMask(3),
    elementIPListsInitialized(false)
{
#ifdef _OPENMP
    omp_init_lock(&ElementIPDataStructureLock);
#endif
}


OctantRec *
OctreeSpatialLocalizer :: findTerminalContaining(OctantRec &startCell, const FloatArray &coords)
{
    OctantRec *currCell = &startCell;
    // found terminal octant containing node
    while ( !currCell->isTerminalOctant() ) {
        auto result = currCell->giveChildContainingPoint(currCell, coords, octreeMask);
        if ( result == OctantRec :: CS_NoChild ) {
            OOFEM_ERROR("internal error - octree inconsistency");
        }
    }
    return currCell;
}


bool
OctreeSpatialLocalizer :: buildOctreeDataStructure()
{
//    OOFEM_LOG_INFO("Initializing Octree structure\n");
    int init = 1, nnode = this->domain->giveNumberOfDofManagers();
    double rootSize, resolutionLimit;
    FloatArray minc(3), maxc(3);

    // test if tree already built
    if ( rootCell ) {
        return true;
    }

    this->elementListsInitialized.resize(this->domain->giveNumberOfRegions() + 1);
    this->elementListsInitialized.zero();
    this->elementIPListsInitialized = false;

    // measure time consumed by octree build phase
    Timer timer;
    timer.startTimer();

    // first determine domain extends (bounding box), and check for degenerated domain type
    for ( int i = 1; i <= nnode; i++ ) {
        Node *node = domain->giveNode(i);
        if ( node ) {
            const auto &coords = node->giveCoordinates();
            if ( init ) {
                init = 0;
                for ( int j = 1; j <= coords.giveSize(); j++ ) {
                    minc.at(j) = maxc.at(j) = coords.at(j);
                }
            } else {
                for ( int j = 1; j <= coords.giveSize(); j++ ) {
                    if ( coords.at(j) < minc.at(j) ) {
                        minc.at(j) = coords.at(j);
                    }

                    if ( coords.at(j) > maxc.at(j) ) {
                        maxc.at(j) = coords.at(j);
                    }
                }
            }
        }
    } // end loop over nodes

    // determine root size
    rootSize = 0.0;
    for ( int i = 1; i <= 3; i++ ) {
        rootSize = max( rootSize, 1.000001 * (maxc.at(i) - minc.at(i)) );
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
    this->rootCell = std::make_unique<OctantRec>(nullptr, center, rootSize * 0.5);

    // Build octree tree
    if ( nnode > OCTREE_MAX_NODES_LIMIT ) {
        this->rootCell->divideLocally(1, this->octreeMask);
    }

    // insert domain nodes into tree
    for ( int i = 1; i <= nnode; i++ ) {
        Node *node = domain->giveNode(i);
        if ( node ) {
            const auto &coords = node->giveCoordinates();
            this->insertNodeIntoOctree(*this->rootCell, i, coords);
        }
    }

    timer.stopTimer();

    // compute max. tree depth
    int treeDepth = this->giveMaxTreeDepthFrom(*this->rootCell);
    OOFEM_LOG_DEBUG( "Octree init [depth %d in %.2fs]\n", treeDepth, timer.getUtime() );

    return true;
}


void
OctreeSpatialLocalizer :: initElementIPDataStructure()
{
    // Original implementation
    //
    int nelems = this->domain->giveNumberOfElements();
    FloatArray jGpCoords;
    if ( this->elementIPListsInitialized ) {
        return;
    }
#ifdef _OPENMP
    omp_set_lock(&ElementIPDataStructureLock); // if not initialized yet; one thread can proceed with init; others have to wait until init completed
    if ( this->elementIPListsInitialized ) {
        omp_unset_lock(&ElementIPDataStructureLock); 
        return;
    }
#endif   
    // insert IP records into tree (the tree topology is determined by nodes)
    for ( int i = 1; i <= nelems; i++ ) {
        // only default IP are taken into account
        Element *ielem = this->giveDomain()->giveElement(i);
        if ( ielem->giveNumberOfIntegrationRules() > 0 ) {
            for ( GaussPoint *jGp: *ielem->giveDefaultIntegrationRulePtr() ) {
                if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                    this->insertIPElementIntoOctree(*this->rootCell, i, jGpCoords);
                } else {
                    OOFEM_ERROR("computeGlobalCoordinates failed");
                }
            }
        }
        // there are no IP (belonging to default integration rule of an element)
        // but the element should be present in octree data structure
        // this is needed by some services (giveElementContainingPoint, for example)
        for ( int j = 1; j <= ielem->giveNumberOfNodes(); j++ ) {
            const auto &nc = ielem->giveNode(j)->giveCoordinates();
            this->insertIPElementIntoOctree(*this->rootCell, i, nc);
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
#ifdef _OPENMP
    omp_unset_lock(&ElementIPDataStructureLock);
#endif

}


void
OctreeSpatialLocalizer :: insertIPElementIntoOctree(OctantRec &rootCell, int elemNum, const FloatArray &coords)
{
    // found terminal octant containing IP
    OctantRec *currCell = this->findTerminalContaining(rootCell, coords);
    currCell->addElementIP(elemNum);
}


void
OctreeSpatialLocalizer :: initElementDataStructure(int region)
{
    FloatArray b0, b1;

    this->init();
    if ( this->elementListsInitialized.giveSize() >= region + 1 && this->elementListsInitialized[region] ) {
        return;
    }

    for ( int i = 1; i <= this->domain->giveNumberOfElements(); i++ ) {
        Element *ielem = this->giveDomain()->giveElement(i);
        if ( ielem->giveRegionNumber() == region || region == 0 ) {
            SpatialLocalizerInterface *interface = static_cast< SpatialLocalizerInterface * >( ielem->giveInterface(SpatialLocalizerInterfaceType) );
            if ( interface ) {
                interface->SpatialLocalizerI_giveBBox(b0, b1);
                this->insertElementIntoOctree(*this->rootCell, region, i, b0, b1);
            }
        }
    }
    this->elementListsInitialized[region] = true;
}


void
OctreeSpatialLocalizer :: insertElementIntoOctree(OctantRec &rootCell, int region, int elemNum, const FloatArray &b0, const FloatArray &b1)
{
    // Compare the bounding box corners to the center to determine which region is overlaps
    // Checks: b0 <= center, b1 >= center for each entry.
    // This is bundled into an array for more convenient code in the loop.
    IntArray bbc [ 2 ] = {
        IntArray(3), IntArray(3)
    };
    FloatArray origin = rootCell.giveOrigin();
    for ( int i = 1; i <= b0.giveSize(); i++ ) {
        if ( this->octreeMask.at(i) ) {
            bbc [ 0 ].at(i) = b0.at(i) <= origin.at(i);
            bbc [ 1 ].at(i) = b1.at(i) >= origin.at(i);
        } else {
            bbc [ 0 ].at(i) = bbc [ 1 ].at(i) = true;
        }
    }
    for ( int i = b0.giveSize() + 1; i <= 3; i++ ) {
        bbc [ 0 ].at(i) = bbc [ 1 ].at(i) = true;
    }
    // Check terminal or recurse
    if ( rootCell.isTerminalOctant() ) {
        rootCell.addElement(region, elemNum);
    } else {
        // The for loops below check if the bounding box overlaps the region before or after the cell center
        // e.g. when i = 0, then we check for the children that have x-coordinate <= than root center
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            if ( bbc [ i ].at(1) ) {
                for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                    if ( bbc [ j ].at(2) ) {
                        for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                            if ( bbc [ k ].at(3) && rootCell.giveChild(i, j, k) ) {
                                this->insertElementIntoOctree(*rootCell.giveChild(i, j, k), region, elemNum, b0, b1);
                            }
                        }
                    }
                }
            }
        }
    }
}


void
OctreeSpatialLocalizer :: insertElementsUsingNodalConnectivitiesIntoOctree(OctantRec &currCell)
{
    if ( currCell.isTerminalOctant() ) {
        for ( int inod: currCell.giveNodeList() ) {
            // loop over cell nodes and ask connectivity table for shared elements
            const IntArray *dofmanConnectivity = domain->giveConnectivityTable()->giveDofManConnectivityArray(inod);
            if ( dofmanConnectivity ) {
                for ( int d: *dofmanConnectivity ) {
                    currCell.addElementIP( d );
                }
            }
        }
    } else {
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    auto child = currCell.giveChild(i, j, k);
                    if ( child ) {
                        this->insertElementsUsingNodalConnectivitiesIntoOctree(*child);
                    }
                }
            }
        }
    }
}


void
OctreeSpatialLocalizer :: insertNodeIntoOctree(OctantRec &rootCell, int nodeNum, const FloatArray &coords)
{
    // found terminal octant containing node
    OctantRec *currCell = this->findTerminalContaining(rootCell, coords);
    // request cell node list
    nodeContainerType &cellNodeList = currCell->giveNodeList();
    int nCellItems = cellNodeList.size();
    int cellDepth = currCell->giveCellDepth();
    // check for refinement criteria
    // should also include max refinement level criteria
    //
    if ( nCellItems > OCTREE_MAX_NODES_LIMIT && cellDepth <= OCTREE_MAX_DEPTH ) {
        // refine tree one level
        currCell->divideLocally(1, this->octreeMask);
        // propagate all nodes already assigned to currCell to children
        for ( int inod: cellNodeList ) {
            Node *node = domain->giveNode(inod);
            const auto &nodeCoords = node->giveCoordinates();
            this->insertNodeIntoOctree(*currCell, inod, nodeCoords);
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
        auto result = currCell->giveChildContainingPoint(currCell, coords, this->octreeMask);
        if ( result != OctantRec :: CS_ChildFound ) {
            OOFEM_ERROR("internal error - octree inconsistency");
        }
        this->insertNodeIntoOctree(*currCell, nodeNum, coords);
    } else {
        currCell->addNode(nodeNum);
    }
}


Element *
OctreeSpatialLocalizer :: giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList)
{
    OctantRec *currCell, *childCell = nullptr;

    this->init();
    this->initElementIPDataStructure();

    // found terminal octant containing point
    currCell = this->findTerminalContaining(*rootCell, coords);

    while ( currCell ) {
        Element *answer = this->giveElementContainingPoint(*currCell, coords, childCell, regionList);
        if ( answer ) {
            return answer;
        }

        childCell = currCell;
        // terminal cell does not contain node and its connected elements containing the given point
        // search at parent level
        currCell = currCell->giveParent();
    }

    // there isn't any element containing point.
    return nullptr;
}

Element *
OctreeSpatialLocalizer :: giveElementContainingPoint(const FloatArray &coords, const Set &eset)
{
    OctantRec *currCell, *childCell = nullptr;

    this->init();
    this->initElementIPDataStructure();

    // found terminal octant containing point
    currCell = this->findTerminalContaining(*rootCell, coords);

    while ( currCell ) {
        // loop over all elements in currCell, skip search on child cell already scanned
        Element *answer = this->giveElementContainingPoint(*currCell, coords, childCell, & eset);
        if ( answer ) {
            return answer;
        }

        childCell = currCell;
        // terminal cell does not contain node and its connected elements containing the given point
        // search at parent level
        currCell = currCell->giveParent();
    }

    // there isn't any element containing point.
    return nullptr;
}


Element *
OctreeSpatialLocalizer :: giveElementContainingPoint(OctantRec &cell, const FloatArray &coords,
                                                     OctantRec *scannedChild, const IntArray *regionList)
{
    // recursive implementation
    if ( cell.isTerminalOctant() ) {

        for ( int iel: cell.giveIPElementList() ) {
            Element *ielemptr = this->giveDomain()->giveElement(iel);

            if ( ielemptr->giveParallelMode() == Element_remote ) {
                continue;
            }

            SpatialLocalizerInterface *interface = static_cast< SpatialLocalizerInterface * >( ielemptr->giveInterface(SpatialLocalizerInterfaceType) );
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

        return nullptr;
    } else {
        // receiver is not terminal octant -> call same service on childs
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    auto child = cell.giveChild(i, j, k);
                    if ( child ) {
                        if ( scannedChild == child ) {
                            continue;
                        }

                        Element *answer = this->giveElementContainingPoint(*child, coords, nullptr, regionList);
                        if ( answer ) {
                            return answer;
                        }
                    }
                }
            }
        }
    }

    return nullptr;
}

Element *
OctreeSpatialLocalizer :: giveElementContainingPoint(OctantRec &cell, const FloatArray &coords,
                                                     OctantRec *scannedChild, const Set *eset)
{
    // recursive implementation
    if ( cell.isTerminalOctant() ) {

        for ( int iel: cell.giveIPElementList() ) {
            Element *ielemptr = this->giveDomain()->giveElement(iel);

            if ( ielemptr->giveParallelMode() == Element_remote ) {
                continue;
            }

            SpatialLocalizerInterface *interface = static_cast< SpatialLocalizerInterface * >( ielemptr->giveInterface(SpatialLocalizerInterfaceType) );
            if ( interface ) {
                if ( eset && ( !eset->hasElement( ielemptr->giveNumber() ) ) ) {
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

        return nullptr;
    } else {
        // receiver is not terminal octant -> call same service on childs
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    auto child = cell.giveChild(i, j, k);
                    if ( child ) {
                        if ( scannedChild == child ) {
                            continue;
                        }

                        Element *answer = this->giveElementContainingPoint(*child, coords, nullptr, eset);
                        if ( answer ) {
                            return answer;
                        }
                    }
                }
            }
        }
    }

    return nullptr;
}


Element *
OctreeSpatialLocalizer :: giveElementClosestToPoint(FloatArray &lcoords, FloatArray &closest,
                                                    const FloatArray &gcoords, int region)
{
    Element *answer = nullptr;
    std :: list< OctantRec * >cellList;
    OctantRec *currCell;
    double radius, prevRadius;
    const FloatArray &c = this->rootCell->giveOrigin();

    this->initElementDataStructure(region);

    // Maximum distance given coordinate and furthest terminal cell ( center_distance + width/2*sqrt(3) )
    double minDist = distance(c, gcoords) + this->rootCell->giveWidth() * 0.87;

    // found terminal octant containing point
    currCell = this->findTerminalContaining(*rootCell, gcoords);

    // Look in center, then expand.
    this->giveElementClosestToPointWithinOctant(*currCell, gcoords, minDist, lcoords, closest, answer, region);
    prevRadius = 0.;
    radius = currCell->giveWidth();
    while ( radius < minDist ) {
        this->giveListOfTerminalCellsInBoundingBox(cellList, gcoords, radius, prevRadius, *this->rootCell);
        for ( OctantRec *icell: cellList ) {
            this->giveElementClosestToPointWithinOctant(*icell, gcoords, minDist, lcoords, closest, answer, region);
        }
        prevRadius = radius;
        radius *= 2.; // Keep expanding the scope until the radius is larger than the entire root cell, then give up (because we have checked every possible cell)
    }
    return answer;
}


void
OctreeSpatialLocalizer :: giveElementClosestToPointWithinOctant(OctantRec &currCell, const FloatArray &gcoords,
                                                                double &minDist, FloatArray &lcoords, FloatArray &closest, Element * &answer, int region)
{
    FloatArray currLcoords;
    FloatArray currClosest;

    auto &elementList = currCell.giveElementList(region);
    for ( int iel: elementList ) {
        Element *ielemptr = this->giveDomain()->giveElement(iel);

        if ( ielemptr->giveParallelMode() == Element_remote ) {
            continue;
        }

        SpatialLocalizerInterface *interface = static_cast< SpatialLocalizerInterface * >( ielemptr->giveInterface(SpatialLocalizerInterfaceType) );
        if ( region > 0 && ielemptr->giveRegionNumber() != region ) {
            continue;
        }
        double currDist = interface->SpatialLocalizerI_giveClosestPoint(currLcoords, currClosest, gcoords);
        if ( currDist < minDist ) {
            lcoords = currLcoords;
            closest = currClosest;
            answer = ielemptr;
            minDist = currDist;
            if ( minDist == 0.0 ) {
                return;
            }
        }
    }
}


GaussPoint *
OctreeSpatialLocalizer :: giveClosestIP(const FloatArray &coords, int region, bool iCohesiveZoneGP)
{
    GaussPoint *nearestGp = nullptr;
    FloatArray jGpCoords;

    this->init();
    this->initElementIPDataStructure();

    // list of already tested elements
    // elementContainerType visitedElems;
    double minDist = 1.1 * rootCell->giveWidth();
    // found terminal octant containing point
    OctantRec *currCell = this->findTerminalContaining(*rootCell, coords);
    auto &elementList = currCell->giveIPElementList();
    // find nearest ip in this terminal cell

    for ( int iel: elementList ) {
        Element *ielem = domain->giveElement(iel);

        if ( ielem->giveParallelMode() == Element_remote ) {
            continue;
        }

        if ( region > 0 && region != ielem->giveRegionNumber() ) {
            continue;
        }

        if ( !iCohesiveZoneGP ) {
            // test if element already visited
            // if (!visitedElems.insert(*pos).second) continue;
            for ( auto &jGp: *ielem->giveDefaultIntegrationRulePtr() ) {
                if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                    // compute distance
                    double dist = distance(coords, jGpCoords);
                    if ( dist < minDist ) {
                        minDist = dist;
                        nearestGp = jGp;
                    }
                } else {
                    OOFEM_ERROR("computeGlobalCoordinates failed");
                }
            }
        } else {
            ////////////////////////////////
            // Check for cohesive zone Gauss points
            XfemElementInterface *xFemEl = dynamic_cast< XfemElementInterface * >(ielem);

            if ( xFemEl ) {
                size_t numCZRules = xFemEl->mpCZIntegrationRules.size();
                for ( size_t czRuleIndex = 0; czRuleIndex < numCZRules; czRuleIndex++ ) {
                    std :: unique_ptr< IntegrationRule > &iRule = xFemEl->mpCZIntegrationRules [ czRuleIndex ];
                    if ( iRule ) {
                        for ( auto &jGp: *iRule ) {
                            if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                                // compute distance
                                double dist = distance(coords, jGpCoords);
                                //printf("czRuleIndex: %d j: %d dist: %e\n", czRuleIndex, j, dist);
                                if ( dist < minDist ) {
                                    minDist = dist;
                                    nearestGp = jGp;
                                }
                            } else {
                                OOFEM_ERROR("computeGlobalCoordinates failed");
                            }
                        }
                    } else {
                        OOFEM_ERROR("iRule is null");
                    }
                }
            }
        }
    }

    FloatArray bborigin = coords;
    // all cell element ip's scanned
    // construct bounding box and test its position within currCell
    auto BBStatus = currCell->testBoundingBox(bborigin, minDist, this->octreeMask);
    if ( BBStatus == OctantRec :: BBS_InsideCell ) {
        return nearestGp;
    } else if ( BBStatus == OctantRec :: BBS_ContainsCell ) {
        std :: list< OctantRec * >cellList;

        // found terminal octant containing point
        OctantRec *startCell = this->findTerminalContaining(*rootCell, coords);
        // go up, until cell containing bbox is found
        if ( startCell != rootCell.get() ) {
            while ( startCell->testBoundingBox(coords, minDist, this->octreeMask) != OctantRec :: BBS_InsideCell ) {
                startCell = startCell->giveParent();
                if ( startCell == rootCell.get() ) {
                    break;
                }
            }
        }

        this->giveListOfTerminalCellsInBoundingBox(cellList, bborigin, minDist, 0, *startCell);

        for ( OctantRec *icell: cellList ) {
            if ( currCell == icell ) {
                continue;
            }

            this->giveClosestIPWithinOctant(*icell, coords, region, minDist, nearestGp, iCohesiveZoneGP);
        }

        return nearestGp;

#if 0
        // found terminal octant containing point
        OctantRec *startCell = this->findTerminalContaining(rootCell, coords);
        // go up, until cell containing bbox is found
        if ( startCell != rootCell ) {
            while ( startCell->testBoundingBox(coords, minDist, this->octreeMask) != OctantRec :: BBS_InsideCell ) {
                startCell = startCell->giveParent();
                if ( startCell == rootCell ) {
                    break;
                }
            }
        }
        // loop over all child (if any) and found all nodes meeting the criteria
        // this->giveClosestIPWithinOctant (startCell, visitedElems, coords, region, minDist, &nearestGp);
        this->giveClosestIPWithinOctant(startCell, coords, region, minDist, & nearestGp);
        return nearestGp;

#endif


#if 0
        set< int >nearElementList;
        // find all elements within bbox and find nearest ip
        this->giveAllElementsWithIpWithinBox(nearElementList, bborigin, minDist);
        // loop over those elements and find nearest ip and return it
        // check for region also
        if ( !nearElementList.empty() ) {
            for ( int elnum: nearElementList ) {
                ielem = domain->giveElement(elnum);
                //igp   = ipContainer[i].giveGp();
                if ( region == ielem->giveRegionNumber() ) {
                    iRule = ielem->giveDefaultIntegrationRulePtr();
                    for ( auto &jGp: *iRule ) {
                        if ( ielem->computeGlobalCoordinates( jGpCoords, * ( jGp->giveCoordinates() ) ) ) {
                            // compute distance
                            double dist = distance(coords, jGpCoords);
                            if ( dist < minDist ) {
                                minDist   = dist;
                                nearestGp = jGp;
                            }
                        } else {
                            OOFEM_ERROR("computeGlobalCoordinates failed");
                        }
                    }
                }
            }
        }
#endif
    } else {
        coords.printYourself("coords");
        OOFEM_ERROR("octree inconsistency found");
    }

    return nullptr;
}


void
OctreeSpatialLocalizer :: giveClosestIPWithinOctant(OctantRec &currentCell, //elementContainerType& visitedElems,
                                                    const FloatArray &coords,
                                                    int region, double &dist, GaussPoint *&answer, bool iCohesiveZoneGP)
{
    if ( currentCell.isTerminalOctant() ) {
        FloatArray jGpCoords;

        // loop over cell elements and check if they meet the criteria
        for ( int iel: currentCell.giveIPElementList() ) {
            // ask for element
            Element *ielem = domain->giveElement(iel);

            if ( ielem->giveParallelMode() == Element_remote ) {
                continue;
            }

            if ( ( region > 0 ) && ( region != ielem->giveRegionNumber() ) ) {
                continue;
            }

            if ( !iCohesiveZoneGP ) {
                // test if element already visited
                // if (!visitedElems.insert(*pos).second) continue;
                // is one of his ip's  within given bbox -> inset it into elemSet
                for ( auto &gp: *ielem->giveDefaultIntegrationRulePtr() ) {
                    if ( ielem->computeGlobalCoordinates( jGpCoords, gp->giveNaturalCoordinates() ) ) {
                        double currDist = distance(coords, jGpCoords);
                        // multiple insertion are handled by STL set implementation
                        if ( currDist <= dist ) {
                            dist = currDist;
                            answer = gp;
                        }
                    } else {
                        OOFEM_ERROR("computeGlobalCoordinates failed");
                    }
                }
            } else {
                //////////////////////////////////////////////////////////
                // Check for cohesive zone Gauss points
                XfemElementInterface *xFemEl = dynamic_cast< XfemElementInterface * >(ielem);

                if ( xFemEl ) {
                    size_t numCZRules = xFemEl->mpCZIntegrationRules.size();
                    for ( size_t czRuleIndex = 0; czRuleIndex < numCZRules; czRuleIndex++ ) {
                        std :: unique_ptr< IntegrationRule > &iRule = xFemEl->mpCZIntegrationRules [ czRuleIndex ];
                        if ( iRule ) {
                            for ( auto &gp: *iRule ) {
                                if ( ielem->computeGlobalCoordinates( jGpCoords, gp->giveNaturalCoordinates() ) ) {
                                    double currDist = distance(coords, jGpCoords);
                                    // multiple insertion are handled by STL set implementation
                                    if ( currDist <= dist ) {
                                        dist = currDist;
                                        answer = gp;
                                    }
                                } else {
                                    OOFEM_ERROR("computeGlobalCoordinates failed");
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    if ( currentCell.giveChild(i, j, k) ) {
                        // test if box hits the cell
                        auto BBStatus = currentCell.giveChild(i, j, k)->testBoundingBox(coords, dist, this->octreeMask);
                        if ( BBStatus == OctantRec :: BBS_InsideCell || BBStatus == OctantRec :: BBS_ContainsCell ) {
                            // if yes call this method for such cell
                            //this->giveClosestIPWithinOctant(*currentCell.giveChild(i,j,k), visitedElems, coords, region, dist, answer);
                            this->giveClosestIPWithinOctant(*currentCell.giveChild(i, j, k), coords, region, dist, answer, iCohesiveZoneGP);
                        }
                    }
                }
            }
        }
    }
}


GaussPoint *
OctreeSpatialLocalizer :: giveClosestIP(const FloatArray &coords, Set &elementSet, bool iCohesiveZoneGP)
{
    GaussPoint *nearestGp = nullptr;
    FloatArray jGpCoords;

    this->init();
    this->initElementIPDataStructure();

    // list of already tested elements
    // elementContainerType visitedElems;
    double minDist = 1.1 * rootCell->giveWidth();
    // found terminal octant containing point
    OctantRec *currCell = this->findTerminalContaining(*rootCell, coords);
    // find nearest ip in this terminal cell
    for ( int iel: currCell->giveIPElementList()) {
        Element *ielem = domain->giveElement(iel);

        if ( ielem->giveParallelMode() == Element_remote ) {
            continue;
        }

        //igp   = ipContainer[i].giveGp();
        if ( !elementSet.hasElement( ielem->giveNumber() ) ) {
            continue;
        }

        if ( !iCohesiveZoneGP ) {
            // test if element already visited
            // if (!visitedElems.insert(ielem).second) continue;
            for ( auto &jGp: *ielem->giveDefaultIntegrationRulePtr() ) {
                if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                    // compute distance
                    double dist = distance(coords, jGpCoords);
                    if ( dist < minDist ) {
                        minDist   = dist;
                        nearestGp = jGp;
                    }
                } else {
                    OOFEM_ERROR("computeGlobalCoordinates failed");
                }
            }
        } else {
            ////////////////////////////////
            // Check for cohesive zone Gauss points
            XfemElementInterface *xFemEl = dynamic_cast< XfemElementInterface * >(ielem);

            if ( xFemEl ) {
                size_t numCZRules = xFemEl->mpCZIntegrationRules.size();
                for ( size_t czRuleIndex = 0; czRuleIndex < numCZRules; czRuleIndex++ ) {
                    std :: unique_ptr< IntegrationRule > &iRule = xFemEl->mpCZIntegrationRules [ czRuleIndex ];
                    if ( iRule ) {
                        for ( auto &jGp: *iRule ) {
                            if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                                // compute distance
                                double dist = distance(coords, jGpCoords);
                                //printf("czRuleIndex: %d j: %d dist: %e\n", czRuleIndex, j, dist);
                                if ( dist < minDist ) {
                                    minDist   = dist;
                                    nearestGp = jGp;
                                }
                            } else {
                                OOFEM_ERROR("computeGlobalCoordinates failed");
                            }
                        }
                    } else {
                        OOFEM_ERROR("iRule is null");
                    }
                }
            }
        }
    }

    FloatArray bborigin = coords;
    // all cell element ip's scanned
    // construct bounding box and test its position within currCell
    auto BBStatus = currCell->testBoundingBox(bborigin, minDist, this->octreeMask);
    if ( BBStatus == OctantRec :: BBS_InsideCell ) {
        return nearestGp;
    } else if ( BBStatus == OctantRec :: BBS_ContainsCell ) {
        std :: list< OctantRec * >cellList;

        // found terminal octant containing point
        OctantRec *startCell = this->findTerminalContaining(*rootCell, coords);
        // go up, until cell containing bbox is found
        if ( startCell != rootCell.get() ) {
            while ( startCell->testBoundingBox(coords, minDist, this->octreeMask) != OctantRec :: BBS_InsideCell ) {
                startCell = startCell->giveParent();
                if ( startCell == rootCell.get() ) {
                    break;
                }
            }
        }

        this->giveListOfTerminalCellsInBoundingBox(cellList, bborigin, minDist, 0, *startCell);

        for ( OctantRec *icell: cellList ) {
            if ( currCell == icell ) {
                continue;
            }

            this->giveClosestIPWithinOctant(*icell, coords, elementSet, minDist, nearestGp, iCohesiveZoneGP );
        }

        return nearestGp;

#if 0
        // found terminal octant containing point
        OctantRec *startCell = this->findTerminalContaining(rootCell, coords);
        // go up, until cell containing bbox is found
        if ( startCell != rootCell ) {
            while ( startCell->testBoundingBox(coords, minDist, this->octreeMask) != OctantRec :: BBS_InsideCell ) {
                startCell = startCell->giveParent();
                if ( startCell == rootCell ) {
                    break;
                }
            }
        }
        // loop over all child (if any) and found all nodes meeting the criteria
        // this->giveClosestIPWithinOctant(startCell, visitedElems, coords, region, minDist, &nearestGp);
        this->giveClosestIPWithinOctant(startCell, coords, region, minDist, & nearestGp);
        return nearestGp;

#endif


#if 0
        set< int >nearElementList;
        // find all elements within bbox and find nearest ip
        this->giveAllElementsWithIpWithinBox(nearElementList, bborigin, minDist);
        // loop over those elements and find nearest ip and return it
        // check for region also
        if ( !nearElementList.empty() ) {
            for ( int iel: nearElementList ) {
                Element *ielem = domain->giveElement(iel);
                //igp   = ipContainer[i].giveGp();
                if ( region == ielem->giveRegionNumber() ) {
                    for ( auto &jGp: *ielem->giveDefaultIntegrationRulePtr() ) {
                        if ( ielem->computeGlobalCoordinates( jGpCoords, * ( jGp->giveCoordinates() ) ) ) {
                            // compute distance
                            double dist = distance(coords, jGpCoords);
                            if ( dist < minDist ) {
                                minDist   = dist;
                                nearestGp = jGp;
                            }
                        } else {
                            OOFEM_ERROR("computeGlobalCoordinates failed");
                        }
                    }
                }
            }
        }
#endif
    } else {
        printf("coords: ");
        coords.printYourself();
        OOFEM_ERROR("octree inconsistency found");
    }

    return nullptr;
}


void
OctreeSpatialLocalizer :: giveClosestIPWithinOctant(OctantRec &currentCell, //elementContainerType& visitedElems,
                                                    const FloatArray &coords,
                                                    Set &elementSet, double &dist, GaussPoint *&answer, bool iCohesiveZoneGP)
{
    if ( currentCell.isTerminalOctant() ) {
        double currDist;
        FloatArray jGpCoords;

        // loop over cell elements and check if they meet the criteria
        for ( int iel: currentCell.giveIPElementList() ) {
            // ask for element
            Element *ielem = domain->giveElement(iel);

            if ( ielem->giveParallelMode() == Element_remote ) {
                continue;
            }

            if ( !elementSet.hasElement( ielem->giveNumber() ) ) {
                continue;
            }

            if ( !iCohesiveZoneGP ) {
                // test if element already visited
                // if (!visitedElems.insert(iel).second) continue;
                // is one of his ip's  within given bbox -> inset it into elemSet
                for ( auto &gp: *ielem->giveDefaultIntegrationRulePtr() ) {
                    if ( ielem->computeGlobalCoordinates( jGpCoords, gp->giveNaturalCoordinates() ) ) {
                        currDist = distance(coords, jGpCoords);
                        // multiple insertion are handled by STL set implementation
                        if ( currDist <= dist ) {
                            dist = currDist;
                            answer = gp;
                        }
                    } else {
                        OOFEM_ERROR("computeGlobalCoordinates failed");
                    }
                }
            } else {
                //////////////////////////////////////////////////////////
                // Check for cohesive zone Gauss points
                XfemElementInterface *xFemEl = dynamic_cast< XfemElementInterface * >(ielem);

                if ( xFemEl ) {
                    size_t numCZRules = xFemEl->mpCZIntegrationRules.size();
                    for ( size_t czRuleIndex = 0; czRuleIndex < numCZRules; czRuleIndex++ ) {
                        std :: unique_ptr< IntegrationRule > &iRule = xFemEl->mpCZIntegrationRules [ czRuleIndex ];
                        if ( iRule ) {
                            for ( auto &gp: *iRule ) {
                                if ( ielem->computeGlobalCoordinates( jGpCoords, gp->giveNaturalCoordinates() ) ) {
                                    currDist = distance(coords, jGpCoords);
                                    // multiple insertion are handled by STL set implementation
                                    if ( currDist <= dist ) {
                                        dist = currDist;
                                        answer = gp;
                                    }
                                } else {
                                    OOFEM_ERROR("computeGlobalCoordinates failed");
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    if ( currentCell.giveChild(i, j, k) ) {
                        // test if box hits the cell
                        auto BBStatus = currentCell.giveChild(i, j, k)->testBoundingBox(coords, dist, this->octreeMask);
                        if ( BBStatus == OctantRec :: BBS_InsideCell || BBStatus == OctantRec :: BBS_ContainsCell ) {
                            // if yes call this method for such cell
                            //this->giveClosestIPWithinOctant (currentCell->giveChild(i,j,k), visitedElems, coords, region, dist, answer);
                            this->giveClosestIPWithinOctant(*currentCell.giveChild(i, j, k), coords, elementSet, dist, answer, iCohesiveZoneGP);
                        }
                    }
                }
            }
        }
    }
}


void
OctreeSpatialLocalizer :: giveAllElementsWithIpWithinBox_EvenIfEmpty(elementContainerType &elemSet, const FloatArray &coords,
                                                         const double radius, bool iCohesiveZoneGP)
{
    this->init();
    this->initElementIPDataStructure();
    // found terminal octant containing point
    OctantRec *currCell = this->findTerminalContaining(*rootCell, coords);
    // go up, until cell containing bbox is found
    if ( currCell != rootCell.get() ) {
        while ( currCell->testBoundingBox(coords, radius, this->octreeMask) != OctantRec :: BBS_InsideCell ) {
            currCell = currCell->giveParent();
            if ( currCell == rootCell.get() ) {
                break;
            }
        }
    }

    // loop over all child (if any) and found all nodes meeting the criteria
    this->giveElementsWithIPWithinBox(elemSet, *currCell, coords, radius, iCohesiveZoneGP);
}


void
OctreeSpatialLocalizer :: giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords,
                                                         const double radius, bool iCohesiveZoneGP)
{
    this->giveAllElementsWithIpWithinBox_EvenIfEmpty(elemSet, coords, radius, iCohesiveZoneGP);
    if ( elemSet.isEmpty() ) {
        OOFEM_ERROR("empty set found");
    }
}


void
OctreeSpatialLocalizer :: giveElementsWithIPWithinBox(elementContainerType &elemSet, OctantRec &currentCell,
                                                      const FloatArray &coords, const double radius, bool iCohesiveZoneGP)
{
    if ( currentCell.isTerminalOctant() ) {
        double currDist;
        FloatArray jGpCoords;

        // loop over cell elements and check if they meet the criteria
        for ( int iel: currentCell.giveIPElementList() ) {
            // test if element is already present
            if ( elemSet.findSorted(iel) ) {
                continue;
            }

            // ask for element
            Element *ielem = domain->giveElement(iel);

            //if(ielem -> giveParallelMode() == Element_remote)continue;
            if ( !iCohesiveZoneGP ) {
                // is one of his ip's  within given bbox -> inset it into elemSet
                for ( auto &gp: *ielem->giveDefaultIntegrationRulePtr() ) {
                    if ( ielem->computeGlobalCoordinates( jGpCoords, gp->giveNaturalCoordinates() ) ) {
                        currDist = distance(coords, jGpCoords);
                        // multiple insertion are handled by STL set implementation
                        if ( currDist <= radius ) {
                            elemSet.insertSortedOnce(iel);
                            break;
                        }
                    } else {
                        OOFEM_ERROR("computeGlobalCoordinates failed");
                    }
                }
            } else {
                ///////////////////////////////////////////////////

                XfemElementInterface *xFemEl = dynamic_cast< XfemElementInterface * >(ielem);

                if ( xFemEl ) {
                    size_t numCZRules = xFemEl->mpCZIntegrationRules.size();
                    for ( size_t czRuleIndex = 0; czRuleIndex < numCZRules; czRuleIndex++ ) {
                        std :: unique_ptr< IntegrationRule > &iRule = xFemEl->mpCZIntegrationRules [ czRuleIndex ];
                        if ( iRule ) {
                            for ( auto &gp: *iRule ) {
                                if ( ielem->computeGlobalCoordinates( jGpCoords, gp->giveNaturalCoordinates() ) ) {
                                    currDist = distance(coords, jGpCoords);
                                    // multiple insertion are handled by STL set implementation
                                    if ( currDist <= radius ) {
                                        elemSet.insertSortedOnce(iel);
                                        czRuleIndex = numCZRules;
                                        break;
                                    }
                                } else {
                                    OOFEM_ERROR("computeGlobalCoordinates failed");
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    auto child = currentCell.giveChild(i, j, k);
                    if ( child ) {
                        // test if box hits the cell
                        auto BBStatus = child->testBoundingBox(coords, radius, this->octreeMask);
                        if ( BBStatus == OctantRec :: BBS_InsideCell || BBStatus == OctantRec :: BBS_ContainsCell ) {
                            // if yes call this method for such cell
                            this->giveElementsWithIPWithinBox(elemSet, *child, coords, radius, iCohesiveZoneGP);
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
    OctantRec *currCell = this->findTerminalContaining(*rootCell, coords);
    // go up, until cell containing bbox is found
    if ( currCell != rootCell.get() ) {
        while ( currCell->testBoundingBox(coords, radius, this->octreeMask) != OctantRec :: BBS_InsideCell ) {
            currCell = currCell->giveParent();
            if ( currCell == rootCell.get() ) {
                break;
            }
        }
    }

    // loop over all child (if any) and found all nodes meeting the criteria
    this->giveNodesWithinBox(nodeSet, *currCell, coords, radius);
}


Node *
OctreeSpatialLocalizer :: giveNodeClosestToPoint(const FloatArray &gcoords, double maxDist)
{
    Node *answer = nullptr;
    std :: list< OctantRec * > cellList;

    // Maximum distance given coordinate and furthest terminal cell ( center_distance + width/2*sqrt(3) )
    double minDist = maxDist;

    // found terminal octant containing point
    OctantRec *currCell = this->findTerminalContaining(*rootCell, gcoords);

    // Look in center, then expand.
    this->giveNodeClosestToPointWithinOctant(*currCell, gcoords, minDist, answer);
    double prevRadius = 0.;
    double radius = min(currCell->giveWidth(), minDist);
    do {
        this->giveListOfTerminalCellsInBoundingBox(cellList, gcoords, radius, prevRadius, *this->rootCell);
        for ( OctantRec *cell: cellList ) {
            this->giveNodeClosestToPointWithinOctant(*cell, gcoords, minDist, answer);
        }
        prevRadius = radius;
        radius *= 2.; // Keep expanding the scope until the radius is larger than the entire root cell, then give up (because we have checked every possible cell)
    } while ( prevRadius < minDist );
    return answer;
}


void
OctreeSpatialLocalizer :: giveNodeClosestToPointWithinOctant(OctantRec &currCell, const FloatArray &gcoords,
                                                                double &minDist, Node * &answer)
{
    double minDist2 = minDist*minDist;
    for ( auto &inod : currCell.giveNodeList() ) {
        Node *node = domain->giveNode(inod);

        double currDist2 = distance_square(gcoords, node->giveCoordinates());

        if ( currDist2 < minDist2 ) {
            answer = node;
            minDist2 = currDist2;
        }
    }
    minDist = sqrt(minDist2);
}


void
OctreeSpatialLocalizer :: giveNodesWithinBox(nodeContainerType &nodeList, OctantRec &currentCell,
                                             const FloatArray &coords, const double radius)
{
    if ( currentCell.isTerminalOctant() ) {
        nodeContainerType &cellNodes = currentCell.giveNodeList();
        if ( !cellNodes.empty() ) {
            for ( int inod: cellNodes ) {
                // loop over cell nodes and check if they meet the criteria
                const auto &nodeCoords = domain->giveNode(inod)->giveCoordinates();
                // is node within bbox
                if ( distance(nodeCoords, coords) <= radius ) {
                    // if yes, append them into set
                    nodeList.push_back(inod);
                }
            }
        }
    } else {
        for ( int i = 0; i <= octreeMask.at(1); i++ ) {
            for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                    auto child = currentCell.giveChild(i, j, k);
                    if ( child ) {
                        // test if box hits the cell
                        auto BBStatus = child->testBoundingBox(coords, radius, this->octreeMask);
                        if ( BBStatus == OctantRec :: BBS_InsideCell || BBStatus == OctantRec :: BBS_ContainsCell ) {
                            // if yes call this method for such cell
                            this->giveNodesWithinBox(nodeList, *child, coords, radius);
                        }
                    }
                }
            }
        }
    }
}


int
OctreeSpatialLocalizer :: giveMaxTreeDepthFrom(OctantRec &root)
{
    int maxDepth = 0;
    for ( int i = 0; i <= octreeMask.at(1); i++ ) {
        for ( int j = 0; j <= octreeMask.at(2); j++ ) {
            for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                auto child = root.giveChild(i, j, k);
                if ( child ) {
                    maxDepth = std::max(maxDepth, this->giveMaxTreeDepthFrom(*child));
                }
            }
        }
    }
    return maxDepth + 1;
}


void
OctreeSpatialLocalizer :: giveListOfTerminalCellsInBoundingBox(std :: list< OctantRec * > &cellList, const FloatArray &coords,
                                                               double radius, double innerRadius, OctantRec &currentCell)
{
    auto BBStatus = currentCell.testBoundingBox(coords, radius, this->octreeMask);
    if ( BBStatus != OctantRec :: BBS_OutsideCell ) {
        if ( currentCell.isTerminalOctant() ) {
            //if ( currentCell.testBoundingBox(coords, innerRadius, this->octreeMask) == OctantRec :: BBS_OutsideCell ) {
            cellList.push_back(&currentCell);
            //}
        } else {
            for ( int i = 0; i <= octreeMask.at(1); i++ ) {
                for ( int j = 0; j <= octreeMask.at(2); j++ ) {
                    for ( int k = 0; k <= octreeMask.at(3); k++ ) {
                        auto child = currentCell.giveChild(i, j, k);
                        if ( child ) {
                            this->giveListOfTerminalCellsInBoundingBox( cellList, coords, radius, innerRadius, *child );
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
        int ans;
#ifdef _OPENMP
#pragma omp single
#endif
        {
            rootCell = nullptr;
            elementIPListsInitialized = false;
            elementListsInitialized.zero();
            ans = this->buildOctreeDataStructure();
            //this->initElementIPDataStructure();
        }
        return ans;
    }

    if ( !rootCell ) {
        int ans;
#ifdef _OPENMP
#pragma omp single
#endif
        {
            OOFEM_LOG_INFO("OctreeLocalizer: init\n");
            ans = this->buildOctreeDataStructure();
            //this->initElementIPDataStructure();
        }
        return ans;
    } else {
        return 0;
    }
}
} // end namespace oofem
