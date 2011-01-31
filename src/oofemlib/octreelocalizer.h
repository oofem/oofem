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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef octreelocalizer_h
#define octreelocalizer_h

#include "spatiallocalizer.h"
#include "compiler.h"

#ifndef __MAKEDEPEND
 #include <set>
 #include <list>
#endif

namespace oofem {
class Domain;
class Element;
class TimeStep;
class OctreeSpatialLocalizer;
/// Max desired number of nodes per octant
#define OCTREE_MAX_NODES_LIMIT 50
/// Max octree depth
#define OCTREE_MAX_DEPTH 15


/**
 * Class representing the octant of octree.
 * It maintains the link to parent cell or if it it the root cell, this link pointer is set to NULL.
 * Maintains links to possible child octree cells as well as its position and size.
 * Also list of node numbers contained in given octree cell can be maintained if cell is terminal cell.
 */
class oofemOctantRec
{
protected:

    /// Link to octree class.
    OctreeSpatialLocalizer *localizer;
    /// Link to parent cell record.
    oofemOctantRec *parent;
    /// Link to octant children.
    oofemOctantRec *child [ 2 ] [ 2 ] [ 2 ];
    /// Octant origin coordinates.
    FloatArray origin;
    /// Octant size.
    double size;

    /// Octant node list.
    std :: list< int > *nodeList;
    /// Element list, containing all elements having IP in cell.
    std :: set< int > *elementList;

public:
    enum boundingBoxStatus { BBS_OUTSIDECELL, BBS_INSIDECELL, BBS_CONTAINSCELL };
    enum ChildStatus { CS_ChildFound, CS_NoChild, CS_PointOutside };

    /// Constructor.
    oofemOctantRec(OctreeSpatialLocalizer *loc, oofemOctantRec *parent, FloatArray &origin, double size);
    /// Destructor.
    ~oofemOctantRec();

    /// @return Reference to parent; NULL if root.
    oofemOctantRec *giveParent() { return this->parent; }
    /**
     * Gives the cell origin.
     * @param answer Cell origin.
     */
    void giveOrigin(FloatArray &answer) { answer = this->origin; }
    /// @return Cell size.
    double giveSize() { return this->size; }
    /**
     * Returns nonzero if octant contains given point.
     * If not 3-coordinates are given, then missing coordinates are
     * not included in test.
     * @param coords Coordinate.
     */
    int containsPoint(const FloatArray &coords);
    /**
     * Gives the Child at given local indices.
     * @param xi First index.
     * @param yi Second index.
     * @param zi Third index.
     * @return Child cell with given local cell coordinates.
     */
    oofemOctantRec *giveChild(int xi, int yi, int zi);
    /**
     * Returns the child containing given point.
     * If not full 3d coordinates are provided, then only provided coordinates are taken into account,
     * assuming remaining to be same as origin.
     * If point is contained by receiver, corresponding child is set, otherwise.
     * child is set to NULL.
     * @param child
     * @param coords Coordinate which child should contain.
     * @return Child status.
     */
    ChildStatus giveChildContainingPoint(oofemOctantRec **child, const FloatArray &coords);
    /// @return True if octant is terminal (no children).
    bool isTerminalOctant();
    /// @return Reference to node List.
    std :: list< int > *giveNodeList();
    /// @return Reference to IPelement set.
    std :: set< int > *giveIPElementList();

    /**
     * Divide receiver further, creating corresponding children.
     */
    int divideLocally(int level, const IntArray &octantMask);
    /**
     * Test if receiver within bounding box (sphere).
     * @return boundingBoxStatus status.
     */
    boundingBoxStatus testBoundingBox(const FloatArray &coords, double radius);
    /**
     * Adds given element to cell list of elements having IP within this cell.
     * @param elementNum Element number to add.
     */
    void addElementIP(int elementNum) { this->giveIPElementList()->insert(elementNum); }
    /**
     * Adds given Node to node list of nodes contained by receiver.
     * @param nodeNum Node number to add.
     */
    void addNode(int nodeNum) { this->giveNodeList()->push_back(nodeNum); }
    /**
     * Clears and deletes the nodeList.
     */
    void deleteNodeList() { // TODO CHECK THIS. This actually deletes elements, not nodes?!
        if ( elementList ) { delete elementList; }

        elementList = NULL;
    }
};





/**
 * The implementation of spatial localizer based on octree technique.
 * The basic task is to provide spatial information and localization for domain, to which receiver is associated.
 * The octree data structure is build over domain nodes, the element type searches are completed using
 * nodal connectivity informations provided by ConTable.
 * Typical services include searching the closes node to give position, searching of an element containing given point, etc.
 * If special element algorithms required, these should be included using interface concept.
 */
class OctreeSpatialLocalizer : public SpatialLocalizer
{
protected:
    /// Root cell of octree.
    oofemOctantRec *rootCell;
    /// Octree degenerate mask.
    IntArray octreeMask;
    /// Flag indicating elementIP tables are initialized.
    int elementIPListsInitialized;

public:
    /// Constructor
    OctreeSpatialLocalizer(int n, Domain *d) : SpatialLocalizer(n, d), octreeMask(3) {
        rootCell = NULL;
        elementIPListsInitialized = 0;
    }
    /// Destructor - deletes the octree tree
    ~OctreeSpatialLocalizer() { if ( rootCell ) { delete rootCell; } }

    /**
     * Returns the octreeMask value given by the index
     * @param indx Index of mask.
     * @return Mask for given index.
     */
    int giveOctreeMaskValue(int indx) { return octreeMask.at(indx); }

    /**
     * Initialize receiver data structure if not done previously.
     * Current implementation calls and returns the buildOctreeDataStructure service response.
     */
    int init(bool force = false);

    Element *giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList = NULL);

    /**
     * Returns the element containing given point.
     * The search is done only for given cell and its children, skipping the given child from search
     * @param cell Top level cell to search.
     * @param coords Point coordinates.
     * @param scannedChild Child pointer to exclude from search.
     * @param regionList Only elements within given regions are considered, if NULL all regions are considered.
     */
    Element *giveElementContainingPoint(oofemOctantRec *cell, const FloatArray &coords,
                                        oofemOctantRec *scannedChild = NULL, const IntArray *regionList = NULL);

    Element *giveElementCloseToPoint(const FloatArray &coords, const IntArray *regionList = NULL);

    GaussPoint *giveClosestIP(const FloatArray &coords, int region);
    void giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords, const double radius);
    void giveAllNodesWithinBox(nodeContainerType &nodeList, const FloatArray &coords, const double radius);

    const char *giveClassName() const { return "OctreeSpatialLocalizer"; }
    classType giveClassID() const { return OctreeSpatialLocalizerClass; }

protected:
    /**
     * Builds the underlying octree data structure.
     * The desired tree level is determined by following rules:
     * - if number of nodes exceed threshold, cell is subdivided.
     * - there is maximal octree level.
     * - in current implementation, the neighbor cell size difference is allowed to be > 2.
     * @return Nonzero if successful, otherwise zero.
     */
    int buildOctreeDataStructure();
    /**
     * Insert IP records into tree (the tree topology is determined by nodes).
     * @return Nonzero if successful, otherwise zero.
     */
    int initElementIPDataStructure();
    /**
     * Finds the terminal octant containing the given point.
     * @param startCell Cell used to start search.
     * @param coords Coordinates of point of interest.
     * @return Pointer to terminal octant, NULL if point outside startingCell.
     */
    oofemOctantRec *findTerminalContaining(oofemOctantRec *startCell, const FloatArray &coords);
    /**
     * Inserts the given node (identified by its number and position) to the octree structure.
     * The tree is traversed until terminal octant containing given position is found and node is then inserted
     * into octant nodal list. If there is too much nodes per cell, this is subdivided further and
     * its assigned nodes are propagated to children.
     * @param rootCell Starting cell for insertion.
     * @param nodeNum Node number.
     * @param coords Corresponding node coordinates.
     * @return Nonzero if node insertion was successful.
     */
    int insertNodeIntoOctree(oofemOctantRec *rootCell, int nodeNum, const FloatArray &coords);
    /**
     * Inserts the given integration point (or more precisely the element owning it) to the octree data structure.
     * The tree is traversed until terminal octant containing given position (ip coordinates) is found
     * and corresponding entry is then inserted into corresponding octant list.
     * @param rootCell Starting cell for insertion.
     * @param elemNum Element number.
     * @param coords Global IP coordinates.
     * @return Nonzero if insertion was successful.
     */
    int insertIPElementIntoOctree(oofemOctantRec *rootCell, int elemNum, const FloatArray &coords);
    /**
     * Initializes the element lists  in octree data structure.
     * This implementation requires that the list of nodes in terminate cells exists
     * simply all shared elements to nodes in terminal cell are added.
     * If this is added to existing implementation based on adding elements only if integration point is in the cell.
     * this leads to more complete element list in terminal cell.
     * @param rootCell Starting cell for octree transversal.
     */
    void insertElementsUsingNodalConnectivitiesIntoOctree(oofemOctantRec *rootCell);
    /**
     * Returns container (set) of elements having integration point within given box and given root cell.
     * @param elemSet answer containing the list of elements meeting the criteria.
     * @param currentCell The starting cell to be transversed.
     * @param coords Center of box of interest.
     * @param radius Radius of bounding sphere.
     */
    void giveElementsWithIPWithinBox(elementContainerType &elemSet, oofemOctantRec *currentCell,
                                     const FloatArray &coords, const double radius);
    /**
     * Returns container (list) of nodes within given box and given root cell.
     * @param nodeList Answer containing the list of nodes meeting the criteria.
     * @param currentCell The starting cell to be transversed.
     * @param coords Center of box of interest.
     * @param radius Radius of bounding sphere.
     */
    void giveNodesWithinBox(nodeContainerType &nodeList, oofemOctantRec *currentCell,
                            const FloatArray &coords, const double radius);

    /**
     * Returns closest IP to given point contained within given octree cell.
     * @param currentCell Starting cell to search, all children will be searched too
     * @param coords Point coordinates.
     * @param region Region id of elements.
     * @param dist Threshold distance, only update answer param, if distance is smaller, distance is updated too.
     * @param answer Pointer to IP, which has the smallest distance "distance" from given point.
     */
    void giveClosestIPWithinOctant(oofemOctantRec *currentCell, //elementContainerType& visitedElems,
                                   const FloatArray &coords,
                                   int region, double &dist, GaussPoint **answer);
    /**
     * Returns the element close to given point as element which center is closest to the given point.
     * Only elements with IP in given cell are considered.
     * The search is done only for given cell and its children, skipping the given child from search.
     * @param cell Top level cell to search.
     * @param coords Point coordinates.
     * @param minDist Distance from the center of returned element.
     * @param answer Requested element.
     * @param regionList Only elements within given regions are considered, if NULL all regions are considered.
     */
    void giveElementCloseToPointWithinOctant(oofemOctantRec *cell, const FloatArray &coords,
                                             double &minDist, Element **answer,
                                             const IntArray *regionList);
    /**
     * Returns the octree depth for given cell. The depth is not parameter of octree cells, but is
     * computed from cell size and root cell size.
     * @param cell Cell to check depth for.
     * @return Cell depth.
     */
    int giveCellDepth(oofemOctantRec *cell);
    /**
     * Determines the max tree depth computed for given tree cell and its children.
     * To obtain total tree depth, root cell should be supplied.
     * The tree depth is always measured from the root cell.
     * Before call, maxDepth should be set to zero.
     * @param root Root of tree.
     * @param maxDepth The maximum depth in the tree.
     */
    void giveMaxTreeDepthFrom(oofemOctantRec *root, int &maxDepth);
    /**
     * Builds the list of terminal cells contained within given box (coords, radius), starting from given currentCell.
     * @param cellList List of terminal cell pointers contained by bounding box.
     * @param coords Center of box of interest.
     * @param radius radius of bounding sphere.
     * @param currentCell starting cell.
     */
    void giveListOfTerminalCellsInBoundingBox(std :: list< oofemOctantRec * > &cellList, const FloatArray &coords,
                                              const double radius, oofemOctantRec *currentCell);
};
} // end namespace oofem
#endif // octreelocalizer_h
