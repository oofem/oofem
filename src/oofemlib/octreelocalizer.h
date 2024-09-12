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

#ifndef octreelocalizer_h
#define octreelocalizer_h

#include "oofemenv.h"
#include "spatiallocalizer.h"
#include "floatarray.h"
#include "intarray.h"

#include <set>
#include <list>
#include <vector>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace oofem {
class Domain;
class Element;
class TimeStep;
class OctreeSpatialLocalizer;
/// Max desired number of nodes per octant
#define OCTREE_MAX_NODES_LIMIT 10
/// Max octree depth
#define OCTREE_MAX_DEPTH 15


/**
 * Class representing the octant of octree.
 * It maintains the link to parent cell or if it it the root cell, this link pointer is set to NULL.
 * Maintains links to possible child octree cells as well as its position and size.
 * Also list of node numbers contained in given octree cell can be maintained if cell is terminal cell.
 */
class OOFEM_NO_EXPORT OctantRec
{
protected:
    /// Link to parent cell record.
    OctantRec *parent;
    /// Link to octant children.
    std::unique_ptr<OctantRec> child [ 2 ] [ 2 ] [ 2 ];
    //std::unique_ptr<OctantRec[2][2][2]> child2;
    /// Octant origin coordinates (lower corner)
    FloatArray origin;
    /// Octant size.
    double halfWidth;
    /// Tree depth
    int depth;

    /// Octant node list.
    std :: list< int > nodeList;
    /// Element list, containing all elements having IP in cell.
    IntArray elementIPList;
    /// Element list of all elements close to the cell.
    std :: vector< std :: list< int > >elementList;


public:
    enum BoundingBoxStatus { BBS_OutsideCell, BBS_InsideCell, BBS_ContainsCell };
    enum ChildStatus { CS_ChildFound, CS_NoChild };

    /// Constructor.
    OctantRec(OctantRec * parent, FloatArray origin, double halfWidth);
    /// Destructor.
    ~OctantRec() {}

    /// @return Reference to parent; NULL if root.
    OctantRec *giveParent() { return this->parent; }
    /**
     * Gives the cell origin.
     * @param answer Cell origin.
     */
    const FloatArray & giveOrigin() { return this->origin; }
    /// @return Half the cell width.
    double giveWidth() { return 2. * this->halfWidth; }
    /// @return Depth in the tree for this octant.
    int giveCellDepth() { return this->depth; }
    /**
     * Gives the Child at given local indices.
     * @param xi First index.
     * @param yi Second index.
     * @param zi Third index.
     * @return Child cell with given local cell coordinates.
     */
    OctantRec *giveChild(int xi, int yi, int zi);
    /**
     * Returns the child containing given point.
     * If not full 3d coordinates are provided, then only provided coordinates are taken into account,
     * assuming remaining to be same as origin.
     * If point is contained by receiver, corresponding child is set, otherwise.
     * child is set to NULL.
     * @param child
     * @param coords Coordinate which child should contain.
     * @param mask Mask for which dimensions are in used (size 3, 0 or 1 values)
     * @return Child status.
     */
    ChildStatus giveChildContainingPoint(OctantRec *&child, const FloatArray &coords, const IntArray &mask);
    /// @return True if octant is terminal (no children).
    bool isTerminalOctant();
    /// @return Reference to node List.
    std :: list< int > &giveNodeList();
    /// @return Reference to IPelement set.
    IntArray &giveIPElementList();
    /// @return Reference to closeElement list.
    std :: list< int > &giveElementList(int region);

    /**
     * Divide receiver further, creating corresponding children.
     * @param level Depth of tree.
     * @param octantMask Masking of dimensions.
     */
    void divideLocally(int level, const IntArray &octantMask);
    /**
     * Test if receiver within bounding box (sphere).
     * @param coords Center of sphere.
     * @param radius Radius of sphere.
     * @param mask Mask for which dimensions are in used (size 3, 0 or 1 values)
     * @return BoundingBoxStatus status.
     */
    BoundingBoxStatus testBoundingBox(const FloatArray &coords, double radius, const IntArray &mask);
    /**
     * Adds given element to cell list of elements having IP within this cell.
     * @param elementNum Element number to add.
     */
    void addElementIP(int elementNum) { this->giveIPElementList().insertSortedOnce(elementNum); }
    /**
     * Adds given element to cell list of elements having IP within this cell.
     * @param region Element region number (0 for global).
     * @param elementNum Element number to add.
     */
    void addElement(int region, int elementNum) { this->giveElementList(region).push_back(elementNum); }
    /**
     * Adds given Node to node list of nodes contained by receiver.
     * @param nodeNum Node number to add.
     */
    void addNode(int nodeNum) { this->giveNodeList().push_back(nodeNum); }
    /**
     * Clears and deletes the nodeList.
     */
    void deleteNodeList() {
        nodeList.clear();
    }
    /// Recursively prints structure.
    void printYourself();
    /// Error printing helper.
    std :: string errorInfo(const char *func) const { return std :: string("OctantRec") + func; }
};


/**
 * The implementation of spatial localizer based on octree technique.
 * The basic task is to provide spatial information and localization for domain, to which receiver is associated.
 * The octree data structure is build over domain nodes, the element type searches are completed using
 * nodal connectivity informations provided by ConTable.
 * Typical services include searching the closes node to give position, searching of an element containing given point, etc.
 * If special element algorithms required, these should be included using interface concept.
 */
class OOFEM_EXPORT OctreeSpatialLocalizer : public SpatialLocalizer
{
protected:
    /// Root cell of octree.
    std::unique_ptr<OctantRec> rootCell;
    /// Octree degenerate mask.
    IntArray octreeMask;
    /// Flag indicating elementIP tables are initialized.
    bool elementIPListsInitialized;
    IntArray elementListsInitialized;
    bool initialized;
#ifdef _OPENMP
    omp_lock_t initLock;
    omp_lock_t buildOctreeDataStructureLock;
    omp_lock_t ElementIPDataStructureLock;
#endif
public:
    /// Constructor
    OctreeSpatialLocalizer(Domain * d);
    /// Destructor - deletes the octree tree
    virtual ~OctreeSpatialLocalizer() {}

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
    int init(bool force = false) override;

    Element *giveElementContainingPoint(const FloatArray &coords, const IntArray *regionList = nullptr) override;
    Element *giveElementContainingPoint(const FloatArray &coords, const Set &eset) override;
    Element *giveElementClosestToPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &gcoords, int region) override;

    GaussPoint *giveClosestIP(const FloatArray &coords, int region, bool iCohesiveZoneGP = false) override;
    GaussPoint *giveClosestIP(const FloatArray &coords, Set &elemSet, bool iCohesiveZoneGP = false) override;

    void giveAllElementsWithIpWithinBox_EvenIfEmpty(elementContainerType &elemSet, const FloatArray &coords, const double radius) override { giveAllElementsWithIpWithinBox_EvenIfEmpty(elemSet, coords, radius, false); }
    void giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords, const double radius) override { giveAllElementsWithIpWithinBox(elemSet, coords, radius, false); }
    void giveAllElementsWithIpWithinBox_EvenIfEmpty(elementContainerType &elemSet, const FloatArray &coords, const double radius, bool iCohesiveZoneGP);
    void giveAllElementsWithIpWithinBox(elementContainerType &elemSet, const FloatArray &coords, const double radius, bool iCohesiveZoneGP);

    void giveAllNodesWithinBox(nodeContainerType &nodeList, const FloatArray &coords, const double radius) override;
    Node * giveNodeClosestToPoint(const FloatArray &coords, double maxDist) override;

    const char *giveClassName() const override { return "OctreeSpatialLocalizer"; }

protected:
    /**
     * Builds the underlying octree data structure.
     * The desired tree level is determined by following rules:
     * - if number of nodes exceed threshold, cell is subdivided.
     * - there is maximal octree level.
     * - in current implementation, the neighbor cell size difference is allowed to be > 2.
     */
    bool buildOctreeDataStructure();
    /**
     * Insert IP records into tree (the tree topology is determined by nodes).
     * @return Nonzero if successful, otherwise zero.
     */
    void initElementIPDataStructure();
    /**
     * Insert element into tree (the tree topology is determined by nodes).
     */
    void initElementDataStructure(int region = 0);
    /**
     * Finds the terminal octant containing the given point.
     * @param startCell Cell used to start search.
     * @param coords Coordinates of point of interest.
     * @return Pointer to terminal octant, NULL if point outside startingCell.
     */
    OctantRec *findTerminalContaining(OctantRec &startCell, const FloatArray &coords);
    /**
     * Inserts the given node (identified by its number and position) to the octree structure.
     * The tree is traversed until terminal octant containing given position is found and node is then inserted
     * into octant nodal list. If there is too much nodes per cell, this is subdivided further and
     * its assigned nodes are propagated to children.
     * @param rootCell Starting cell for insertion.
     * @param nodeNum Node number.
     * @param coords Corresponding node coordinates.
     */
    void insertNodeIntoOctree(OctantRec &rootCell, int nodeNum, const FloatArray &coords);
    /**
     * Inserts the given integration point (or more precisely the element owning it) to the octree data structure.
     * The tree is traversed until terminal octant containing given position (ip coordinates) is found
     * and corresponding entry is then inserted into corresponding octant list.
     * @param rootCell Starting cell for insertion.
     * @param elemNum Element number.
     * @param coords Global IP coordinates.
     */
    void insertIPElementIntoOctree(OctantRec &rootCell, int elemNum, const FloatArray &coords);
    /**
     * Inserts an element with the given bounding box.
     * @param rootCell Starting cell for insertion.
     * @param region Element region number.
     * @param elemNum Element number.
     * @param b0 Lower bounding box.
     * @param b1 Upper bounding box.
     */
    void insertElementIntoOctree(OctantRec &rootCell, int region, int elemNum, const FloatArray &b0, const FloatArray &b1);
    /**
     * Initializes the element lists  in octree data structure.
     * This implementation requires that the list of nodes in terminate cells exists
     * simply all shared elements to nodes in terminal cell are added.
     * If this is added to existing implementation based on adding elements only if integration point is in the cell.
     * this leads to more complete element list in terminal cell.
     * @param rootCell Starting cell for octree transversal.
     */
    void insertElementsUsingNodalConnectivitiesIntoOctree(OctantRec &rootCell);
    /**
     * Returns container (set) of elements having integration point within given box and given root cell.
     * @param elemSet answer containing the list of elements meeting the criteria.
     * @param currentCell The starting cell to be transversed.
     * @param coords Center of box of interest.
     * @param radius Radius of bounding sphere.
     */
    void giveElementsWithIPWithinBox(elementContainerType &elemSet, OctantRec &currentCell,
                                     const FloatArray &coords, const double radius, bool iCohesiveZoneGP = false);
    /**
     * Returns container (list) of nodes within given box and given root cell.
     * @param nodeList Answer containing the list of nodes meeting the criteria.
     * @param currentCell The starting cell to be transversed.
     * @param coords Center of box of interest.
     * @param radius Radius of bounding sphere.
     */
    void giveNodesWithinBox(nodeContainerType &nodeList, OctantRec &currentCell,
                            const FloatArray &coords, const double radius);

    /**
     * Returns closest IP to given point contained within given octree cell.
     * @param currentCell Starting cell to search, all children will be searched too
     * @param coords Point coordinates.
     * @param region Region id of elements.
     * @param dist Threshold distance, only update answer param, if distance is smaller, distance is updated too.
     * @param answer Pointer to IP, which has the smallest distance "distance" from given point.
     */
    void giveClosestIPWithinOctant(OctantRec &currentCell, //elementContainerType& visitedElems,
                                   const FloatArray &coords,
                                   int region, double &dist, GaussPoint *&answer, bool iCohesiveZoneGP);
    /**
     * Returns closest IP to given point contained within given octree cell.
     * @param currentCell Starting cell to search, all children will be searched too
     * @param coords Point coordinates.
     * @param elemSet set of considered elements.
     * @param dist Threshold distance, only update answer param, if distance is smaller, distance is updated too.
     * @param answer Pointer to IP, which has the smallest distance "distance" from given point.
     */
    void giveClosestIPWithinOctant(OctantRec &currentCell, //elementContainerType& visitedElems,
                                   const FloatArray &coords,
                                   Set &elemSet, double &dist, GaussPoint *&answer, bool iCohesiveZoneGP);
    /**
     * Returns the element containing given point.
     * The search is done only for given cell and its children, skipping the given child from search
     * @param cell Top level cell to search.
     * @param coords Point coordinates.
     * @param scannedChild Child pointer to exclude from search.
     * @param regionList Only elements within given regions are considered, if NULL all regions are considered.
     * @note regions depreceted, use sets insteed
     */
    Element *giveElementContainingPoint(OctantRec &cell, const FloatArray &coords,
                                        OctantRec *scannedChild = nullptr, const IntArray *regionList = nullptr);
    /**
     * Returns the element containing given point.
     * The search is done only for given cell and its children, skipping the given child from search
     * @param cell Top level cell to search.
     * @param coords Point coordinates.
     * @param scannedChild Child pointer to exclude from search.
     * @param elset Only elements in gibven set are considered, if NULL all regions are considered.
     */
    Element *giveElementContainingPoint(OctantRec &cell, const FloatArray &coords,
                                        OctantRec *scannedChild = nullptr, const Set *elset = nullptr);
    /**
     * Returns the element closest to the given point within the cell.
     * @param currCell Terminal cell to look in.
     * @param gcoords Point coordinates.
     * @param lcoords Local coordinates in closest element.
     * @param[in,out] minDist Distance from the center of returned element.
     * @param closest Coordinates in closest element.
     * @param answer Requested element.
     * @param region Region to consider.
     */
    void giveElementClosestToPointWithinOctant(OctantRec &currCell, const FloatArray &gcoords,
                                               double &minDist, FloatArray &lcoords, FloatArray &closest, Element * &answer, int region);
    /**
     * Returns the node closest to the given point within the cell.
     * @param currCell Terminal cell to look in.
     * @param gcoords Point coordinates.
     * @param[in,out] minDist Distance from the center of returned element.
     * @param answer Requested node.
     */
    void giveNodeClosestToPointWithinOctant(OctantRec &cell, const FloatArray &gcoords, double &minDist, Node * &answer);
    /**
     * Determines the max tree depth computed for given tree cell and its children.
     * To obtain total tree depth, root cell should be supplied.
     * The tree depth is always measured from the root cell.
     * Before call, maxDepth should be set to zero.
     * @param root Root of tree.
     * @return The maximum depth in the tree from the given root
     */
    int giveMaxTreeDepthFrom(OctantRec &root);
    /**
     * Builds the list of terminal cells contained within given box (coords, radius), starting from given currentCell.
     * @param cellList List of terminal cell pointers contained by bounding box.
     * @param coords Center of box of interest.
     * @param radius Radius of bounding sphere.
     * @param innerRadius Inner radius of bounding sphere.
     * @param currentCell Starting cell.
     */
    void giveListOfTerminalCellsInBoundingBox(std :: list< OctantRec * > &cellList, const FloatArray &coords,
                                              const double radius, double innerRadius, OctantRec &currentCell);
};
} // end namespace oofem
#endif // octreelocalizer_h
