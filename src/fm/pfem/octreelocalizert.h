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

//   **************************************
//   *** CLASS OCTREE SPATIAL LOCALIZER ***
//   **************************************
// general Octree localizer - using template

#ifndef octreelocalizert_h
#define octreelocalizert_h

#include "spatiallocalizer.h"
#include "octreelocalizer.h"
#include "delaunaytriangle.h"
#include "dofmanager.h"
#include "node.h"
#include "element.h"
#include "mathfem.h"
#include "timer.h"

#include <set>
#include <list>

namespace oofem {
class Domain;
class Element;
class TimeStep;
class DofManager;
class IntArray;

template< class T >class OctreeSpatialLocalizerT;


#define TEMPLATED_OCTREE_MAX_NODES_LIMIT 300

#define TEMPLATED_OCTREE_MAX_DEPTH 15

///
enum boundingSphereStatus { SphereOutsideCell, SphereInsideCell, SphereContainsCell };

/**
 * Squared bounding box for templated octree localizer
 *
 * @author David Krybus
 */
class BoundingBox
{
protected:
    /// Starting point
    FloatArray origin;
    /// Bounding box size length
    double size;
    /// Spatial dimension mask
    IntArray spatialMask;

public:
    /// Constructor
    BoundingBox() : spatialMask(3) { }

    /// Destructor
    ~BoundingBox() { }

    /// Sets the origin of the bounding box
    void setOrigin(FloatArray &coords)
    {
        origin.resize( coords.giveSize() );
        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            origin.at(i) = coords.at(i);
        }
    }

    /// Sets the size of the bounding box (all sides are equal)
    void setSize(double s) { size = s; }
    /// Sets the spatial mask
    void setMask(int i, int mask) { spatialMask.at(i) = mask; }
    /**
     * Returns the bounding box origin
     * @param answer the bounding box origin
     */
    void giveOrigin(FloatArray &answer) { answer = this->origin; }
    /// Gives the size of the bounding box
    double giveSize() { return size; }
    /**
     * Gives the spatial mask of the bounding box
     * @param answer the mask of the bounding box
     */
    void giveMask(IntArray &answer) { answer = this->spatialMask; }
};

/**
 * Templated octree cell containing data of T type
 *
 * @author David Krybus
 */
template< class T >
class OctantRecT
{
public:
    typedef OctreeSpatialLocalizerT< T > *LocalizerPtrType;
    typedef OctantRecT< T > *CellPtrType;
    typedef typename std :: list< T > :: iterator listIteratorType;

protected:
    /// Octree localizer whose part of is the cell
    LocalizerPtrType localizer;
    /// Parent octree cell
    CellPtrType parent;
    /// Child cells
    CellPtrType child [ 2 ] [ 2 ] [ 2 ];
    /// Origin of the cell
    FloatArray origin;
    /// Size of the cell
    double size;

    /// Octant cell member list
    std :: list< T > *dataList;

public:

    /// Constructor
    OctantRecT(LocalizerPtrType loc, CellPtrType parent, FloatArray &origin, double size)
    {
        this->localizer = loc;
        this->parent = parent;
        this->origin = origin;
        this->size   = size;
        this->dataList = NULL;


        for ( int i = 0; i <= 1; i++ ) {
            for ( int j = 0; j <= 1; j++ ) {
                for ( int k = 0; k <= 1; k++ ) {
                    this->child [ i ] [ j ] [ k ] = NULL;
                }
            }
        }
    }
    /// Destructor
    ~OctantRecT()
    {
        for ( int i = 0; i <= 1; i++ ) {
            for ( int j = 0; j <= 1; j++ ) {
                for ( int k = 0; k <= 1; k++ ) {
                    if ( this->child [ i ] [ j ] [ k ] ) {
                        delete this->child [ i ] [ j ] [ k ];
                    }
                }
            }
        }

        if ( dataList ) {
            delete dataList;
        }
    }

    /// @return Reference to parent; NULL if root.
    CellPtrType giveParent() { return this->parent; }
    /**
     * Gives the cell origin.
     * @param answer Cell origin.
     */
    void giveOrigin(FloatArray &answer) { answer = this->origin; }
    /// Gives the size of the cell
    double giveSize() { return this->size; }
    /**
     * Gives 1 if a given point is contained in the cell, 0 otherwise.
     * @param coords coordinates of tested point.
     */
    int containsPoint(const FloatArray &coords) {
        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            if ( localizer->giveOctreeMaskValue(i) ) {
                if ( coords.at(i) < this->origin.at(i) ) {
                    return 0;
                }

                if ( coords.at(i) > ( this->origin.at(i) + this->size ) ) {
                    return 0;
                }
            }
        }

        return 1;
    }
    /**
     * Gives the Child at given local indices.
     * @param xi First index.
     * @param yi Second index.
     * @param zi Third index.
     * @return Child cell with given local cell coordinates.
     */
    CellPtrType giveChild(int xi, int yi, int zi) {
        if ( ( xi >= 0 ) && ( xi < 2 ) && ( yi >= 0 ) && ( yi < 2 ) && ( zi >= 0 ) && ( zi < 2 ) ) {
            return this->child [ xi ] [ yi ] [ zi ];
        } else {
            OOFEM_ERROR("OctantRecT::giveChild invalid child index (%d,%d,%d)", xi, yi, zi);
        }
        return NULL;
    }

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
    int giveChildContainingPoint(CellPtrType *child, const FloatArray &coords) {
        IntArray ind(3);
        if ( this->containsPoint(coords) ) {
            if ( this->isTerminalOctant() ) {
                * child = NULL;
                return -1;
            }

            for ( int i = 1; i <= coords.giveSize(); i++ ) {
                if ( localizer->giveOctreeMaskValue(i) && ( coords.at(i) > ( this->origin.at(i) + this->size / 2. ) ) ) {
                    ind.at(i) = 1;
                } else {
                    ind.at(i) = 0;
                }
            }

            * child = this->child [ ind.at(1) ] [ ind.at(2) ] [ ind.at(3) ];
            return 1;
        } else {
            * child = NULL;
            return -2;
        }
    }

    /// @return True if octant is terminal (no children).
    int isTerminalOctant()
    {
        if ( this->child [ 0 ] [ 0 ] [ 0 ] ) {
            return 0;
        }

        return 1;
    }

    /**
     * Divide receiver further, creating corresponding children.
     * @param level Depth of tree.
     * @param octantMask Masking of dimensions.
     */
    int divideLocally(int level, const IntArray &octantMask)
    {
        int i, j, k, result = 1;

        if ( this->isTerminalOctant() ) {
            // create corresponding child octants
            int i, j, k;
            FloatArray childOrigin(3);

            for ( i = 0; i <= octantMask.at(1); i++ ) {
                for ( j = 0; j <= octantMask.at(2); j++ ) {
                    for ( k = 0; k <= octantMask.at(3); k++ ) {
                        childOrigin.at(1) = this->origin.at(1) + i * ( this->size / 2. );
                        childOrigin.at(2) = this->origin.at(2) + j * ( this->size / 2. );
                        childOrigin.at(3) = this->origin.at(3) + k * ( this->size / 2. );

                        this->child [ i ] [ j ] [ k ] = new OctantRecT(localizer, this, childOrigin, this->size / 2.0);
                    }
                }
            }
        }

        int newLevel = level - 1;
        if ( newLevel > 0 ) {
            // propagate message to childs recursivelly with level param decreased
            for ( i = 0; i <= octantMask.at(1); i++ ) {
                for ( j = 0; j <= octantMask.at(2); j++ ) {
                    for ( k = 0; k <= octantMask.at(3); k++ ) {
                        if ( this->child [ i ] [ j ] [ k ] ) {
                            result &= this->child [ i ] [ j ] [ k ]->divideLocally(newLevel, octantMask);
                        }
                    }
                }
            }
        }

        return result;
    }

    /**
     * Tests the position of a bounding box in relation to the octant cell
     * @param testedBBX Tested bounding box
     * @returns BoundingBoxStatus status
     */
    OctantRec :: BoundingBoxStatus testBoundingBox(BoundingBox &testedBBX)
    {
        int i, test = 0;
        double BBXSize = testedBBX.giveSize();
        FloatArray BBXOrigin;
        testedBBX.giveOrigin(BBXOrigin);
        int BBXOriginInside = this->containsPoint(BBXOrigin);
        if ( BBXOriginInside ) {
            for ( i = 1; i <= 3; i++ ) {
                if ( localizer->giveOctreeMaskValue(i) ) {
                    if ( this->size > ( BBXSize + ( BBXOrigin.at(i) - this->origin.at(i) ) ) ) {
                        test = 1;
                    }
                }
            }

            if ( test ) {
                return OctantRec :: BBS_InsideCell;
            } else {
                return OctantRec :: BBS_ContainsCell;
            }
        } else {
            for ( i = 1; i <= 3; i++ ) {
                if ( localizer->giveOctreeMaskValue(i) ) {
                    if ( BBXOrigin.at(i) > ( this->size + this->origin.at(i) ) ) {
                        return OctantRec :: BBS_OutsideCell;
                    } else {
                        if ( BBXOrigin.at(i) + BBXSize > this->origin.at(i) ) {
                            return OctantRec :: BBS_ContainsCell;
                        } else {
                            return OctantRec :: BBS_OutsideCell;
                        }
                    }
                }
            }
        }
    }

    /**
     *  Tests the position of a bounding spere in relation to the octant cell
     * @param coords Coordinates of the spere center
     * @param radius Radius of the sphere
     * @returns boundingSpereStatus status
     */
    boundingSphereStatus testBoundingSphere(const FloatArray &coords, double radius)
    {
        int i, test = 1, nsd, size = coords.giveSize();
        double dist;

        nsd = localizer->giveOctreeMaskValue(1) + localizer->giveOctreeMaskValue(2) + localizer->giveOctreeMaskValue(3);

        FloatArray cellCenter = this->origin;
        double cellRadius = sqrt( nsd * ( this->size / 2.0 ) * ( this->size / 2. ) );
        for ( i = 1; i < 4; i++ ) {
            cellCenter.at(i) += this->size / 2.0;
        }

        for ( i = 1, dist = 0.0; i <= size; i++ ) {
            if ( localizer->giveOctreeMaskValue(i) ) {
                dist += ( cellCenter.at(i) - coords.at(i) ) * ( cellCenter.at(i) - coords.at(i) );
            }
        }

        dist = sqrt(dist);
        if ( dist > ( cellRadius + radius ) ) {
            return SphereOutsideCell;
        }

        int centerInside = this->containsPoint(coords);
        if ( centerInside ) {                 // test if whole bsphere inside
            for ( i = 1; i <= size; i++ ) {
                if ( localizer->giveOctreeMaskValue(i) ) {
                    if ( ( this->origin.at(i) > ( coords.at(i) - radius ) )  ||  ( ( this->origin.at(i) + this->size ) < ( coords.at(i) + radius ) ) ) {
                        test = 0;
                    }
                }
            }

            if ( test ) {
                return SphereInsideCell;
            } else {
                return SphereContainsCell;
            }
        } else {                 // BSphere center not inside cell, but may hit cell
                                 // test if bounding sphere hits boundary surface
            int inBounds = 1;
            for ( i = 1; i <= size; i++ ) {
                if ( localizer->giveOctreeMaskValue(i) && ( fabs( cellCenter.at(i) - coords.at(i) ) > ( this->size / 2. + radius ) ) ) {
                    inBounds = 0;
                }
            }

            if ( inBounds ) {
                return SphereContainsCell;
            }
        }

        return SphereOutsideCell;
    }


    /// Return reference to member List
    std :: list< T > *giveDataList() {
        if ( dataList == NULL ) {
            dataList = new std :: list< T >;
        }

        return dataList;
    }

    /// Removes the data list
    void deleteDataList() {
        if ( dataList ) {
            delete dataList;
            dataList = NULL;
        }
    }

    /// Adds the object in to the data list returning iterator position to the object in the list
    listIteratorType addMember(T &member)
    {
        this->giveDataList()->push_front(member);
        return this->giveDataList()->begin();
    }

    /// Removes member from data list
    void removeMember(T &member) { this->giveDataList()->remove(member); }

    friend class OctreeSpatialLocalizerT< T >;
};

/**
 * Help class for storing pointer to octant cell and position of the member in the data list
 *
 * @author David Krybus
 */
template< class T >
class LocalInsertionData
{
public:
    /// Constructor
    LocalInsertionData() { }
    typedef OctantRecT< T > *CellPtrType;
    typedef typename std :: list< T > :: iterator listIteratorType;

    /// Octant cell containing object
    OctantRecT< T > *containedInCell;
    /// Iterator position in the list of cell objects
    listIteratorType posInCellDataList;
};



/**
 * Functor base class responsible for insertion of members into the octree cell
 *
 * @author David Krybus
 */
template< class T >class SL_Insertion_Functor
{
protected:
    typedef typename std :: list< T > :: iterator listIteratorType;

public:
    /// Evaluates wether the member should be stored in the octant cell
    virtual bool evaluate(T &member, OctantRecT< T > *cell) = 0;
    /// Stores LocalInsertionData on the member
    virtual void registerInsertion(T &member, LocalInsertionData< T >LIdata) = 0;
    /// Returns list of LocalInsertionData stored on the member
    virtual std :: list< LocalInsertionData< T > > *giveInsertionList(T &member) = 0;
};


/**
 * Functor for storing nodes in the octree
 *
 * @author David Krybus
 */
class InsertNode : public SL_Insertion_Functor< int >
{
protected:
    Domain *domain;

public:
    /// Constuctor
    InsertNode(Domain *d) {
        domain = d;
    }
    /// Destructor
    ~InsertNode() { }

    /**
     * Evaluates the position of a node
     * @param nodeNr Number of the node in the domain
     * @param cell Octant cell to be tested
     * @returns true if the node is contained in the cell, false otherwise
     */
    bool evaluate(int &nodeNr, OctantRecT< int > *cell)
    {
        FloatArray coords(3);
        for ( int i = 1; i <= 3; i++ ) {
            coords.at(i) = this->domain->giveNode(nodeNr)->giveCoordinate(i);
        }

        if ( cell->containsPoint(coords) ) {
            return true;
        } else {
            return false;
        }
    };

    void registerInsertion(int &nodeNr, LocalInsertionData< int >LIdata) { }

    std :: list< LocalInsertionData< int > > *giveInsertionList(int &nodeNr)
    {
        return NULL;
    }
};

/**
 * Functor for storing triangles in the octree according to theirs circumscribed circles
 *
 * @author David Krybus
 */
class InsertTriangleBasedOnCircumcircle : public SL_Insertion_Functor< DelaunayTriangle * >
{
protected:
    Domain *domain;

public:
    /// Constructor
    InsertTriangleBasedOnCircumcircle(Domain *d) {
        domain = d;
    }
    /// Destructor
    ~InsertTriangleBasedOnCircumcircle() { }

    /**
     * Evaluates the position of a triangle
     * @param DTptr Delaunay triangle
     * @param cell Octant cell to be tested
     * @returns true if the cirumscribed circle of the circle lies inside or contains a part of the cell, false otherwise
     */
    bool evaluate(DelaunayTriangle * &DTptr, OctantRecT< DelaunayTriangle * > *cell)
    {
        double radius = DTptr->giveCircumRadius();
        FloatArray centerCoords(3);
        centerCoords.at(1) = DTptr->giveXCenterCoordinate();
        centerCoords.at(2) = DTptr->giveYCenterCoordinate();
        centerCoords.at(3) = 0.0;

        if ( cell->testBoundingSphere(centerCoords, radius) == SphereOutsideCell ) {
            return false;
        } else {
            return true;
        }
    }
    void registerInsertion(DelaunayTriangle * &TEptr, LocalInsertionData< DelaunayTriangle * >LIdata)
    {
        ( TEptr->giveListOfCellsAndPosition() )->push_back(LIdata);
    }
    std :: list< LocalInsertionData< DelaunayTriangle * > > *giveInsertionList(DelaunayTriangle * &DTptr)
    {
        return DTptr->giveListOfCellsAndPosition();
    }
};

/**
 * Functor base class for evaluating search tasks on the octree according given condition
 *
 * @author David Krybus
 */
template< class T >class SL_Evaluation_Functor
{
public:
    /// Evaluates wether the search condition is accomplished or not
    virtual bool evaluate(T &obj) = 0;
    /// Gives the starting position of the search
    virtual void giveStartingPosition(FloatArray &answer) = 0;
    /// Gives a container with found objects
    virtual void giveResult(std :: list< T > &answer) = 0;

    /// Stage1 means, we are looking for objects in a distance given by some boundingBox
    /// (e.g. IPs around some given coordinates)
    virtual bool isBBXStage1Defined(BoundingBox &BBXStage1) = 0;

    /// Stage2BBX is given by results of a prior search.
    /// e.g. We found a closest point to another point in an octant, BBX is defined by the distance from each other and the found closest point.
    /// Now we have to check surrounding octants whithin this bounding, which may contain points closer to starting point.
    virtual bool isBBXStage2Defined(BoundingBox &BBXStage2) = 0;
};


/**
 * Functor for closest node search
 *
 * @author David Krybus
 */
class ClosestNode : public SL_Evaluation_Functor< int >
{
protected:
    FloatArray *startingPosition;
    Domain *domain;
    int CNindex;
    double distanceToClosestNode;
    bool initFlag;

    std :: list< int >closestNodeIndices;

public:
    /**
     * Constructor
     * @param pos Starting position of the search
     * @param d Domain containing nodes
     */
    ClosestNode(FloatArray *pos, Domain *d) {
        initFlag = false;
        startingPosition = pos;
        domain = d;
    }
    /// Destructor
    ~ClosestNode() { }

    /**
     * Gives the starting position of the search
     * @param position startingPosition
     */
    void giveStartingPosition(FloatArray &position)
    {
        position = * startingPosition;
    }

    /**
     * Evaluates a node. The closest nodes are stored in a container, if their distance to starting position is the same.
     * Is the distance smaller than previous one, the container is emptied and new node is added.
     * @param nodeNr Number of the node in the domain list
     * @returns true after evaluation is processed.
     */
    bool evaluate(int &nodeNr)
    {
        if ( initFlag ) {
            double distance = startingPosition->distance( this->domain->giveNode(nodeNr)->giveCoordinates() );

            if ( ( distance - distanceToClosestNode ) <= distance * 0.001 ) {
                if ( ( distance - distanceToClosestNode ) >= -0.001 * distance ) {
                    closestNodeIndices.push_back(nodeNr);
                } else {
                    closestNodeIndices.clear();
                    closestNodeIndices.push_back(nodeNr);
                    distanceToClosestNode = distance;
                }
            }
        } else {
            closestNodeIndices.push_back(nodeNr);
            distanceToClosestNode = startingPosition->distance( this->domain->giveNode( * ( closestNodeIndices.begin() ) )->giveCoordinates() );
            initFlag = true;
        }

        return true;
    }

    /**
     * Gives the closest nodes
     * @param answer List containing numbers of the closest nodes
     */
    void giveResult(std :: list< int > &answer) {
        answer = closestNodeIndices;
    }

    bool isBBXStage1Defined(BoundingBox &BBXStage1)
    {
        return false;
    }

    bool isBBXStage2Defined(BoundingBox &BBXStage2)
    {
        if ( initFlag ) {
            BBXStage2.setOrigin(* startingPosition);
            BBXStage2.setSize(distanceToClosestNode);
            return true;
        }                 else {
            return false;
        }
    }
};


/**
 * Functor for finding triangles whose circumscribed circles contains given node
 *
 * @author David Krybus
 */
class ElementCircumCirclesContainingNode : public SL_Evaluation_Functor< DelaunayTriangle * >
{
protected:
    FloatArray *startingPosition;
    Domain *domain;
    std :: list< DelaunayTriangle * >result;

public:
    /**
     * Constructor
     * @param pos Starting position of the search
     * @param d Domain containing nodes the triangles are defined by
     */
    ElementCircumCirclesContainingNode(FloatArray *pos, Domain *d)
    {
        startingPosition = pos;
        domain = d;
    }
    ~ElementCircumCirclesContainingNode() { }

    /**
     * Evaluates a triangle upon its circumscribed cricle.
     * @param DTptr Delaunay triangle. nodeNr Number of the node in the domain list
     * @returns true if the circumscribed circle of the Delaunay triangle contains the point, false otherwise
     */
    bool evaluate(DelaunayTriangle * &DTptr)
    {
        double radius = DTptr->giveCircumRadius();
        FloatArray centerCoords(3);
        centerCoords.at(1) = DTptr->giveXCenterCoordinate();
        centerCoords.at(2) = DTptr->giveYCenterCoordinate();
        centerCoords.at(3) = 0.0;

        if ( ( startingPosition->distance(centerCoords) ) < radius ) {
            result.push_back(DTptr);
            return true;
        } else {
            return false;
        }
    }

    /**
     * Gives the starting position of the search
     * @param position startingPosition
     */
    void giveStartingPosition(FloatArray &answer)
    {
        answer = * startingPosition;
    }

    /**
     * Gives the triangles containing the node
     * @param answer List containing Delaunay triangles
     */
    void giveResult(std :: list< DelaunayTriangle * > &answer)
    {
        answer = result;
    }

    // neither stage1, nor stage 2 are defined
    bool isBBXStage1Defined(BoundingBox &BBXStage1) { return false; }
    bool isBBXStage2Defined(BoundingBox &BBXStage2) { return false; }
};

/**
 * Templated octree spatial localizer
 *
 * @author David Krybus
 */
template< class T >class OctreeSpatialLocalizerT
{
protected:

    typedef OctantRecT< T > *CellPtrType;
    typedef std :: list< T >dataContainerType;
    typedef typename std :: list< T > :: iterator listIteratorType;
    typedef typename std :: list< T > :: const_iterator listConstIteratorType;
    IntArray octreeMask;
    CellPtrType rootCell;
    Domain *domain;
    int maxDepthReached;

public:
    /// Constructor
    OctreeSpatialLocalizerT(int n, Domain *d) {
        rootCell = NULL;
        domain = d;
        maxDepthReached = 0;
    }
    /// Destructor
    ~OctreeSpatialLocalizerT() {
        if ( rootCell ) {
            delete rootCell;
        }
    }

    /// Initilizes the octree structure
    int init(BoundingBox &BBX, int initialDivision = 0) {
        if ( !rootCell ) {
            FloatArray origin;
            BBX.giveOrigin(origin);
            this->rootCell = new OctantRecT< T >( this, NULL, origin, BBX.giveSize() );
            BBX.giveMask(octreeMask);
            if ( initialDivision ) {
                this->rootCell->divideLocally(initialDivision, this->octreeMask);
            }

            return 1;
        } else {
            return 0;
        }
    }

    /// Calculates the bounding box base on the domain's nodes
    void computeBBXBasedOnNodeData(BoundingBox &BBX)
    {
        int i, j, init = 1, nnode = this->domain->giveNumberOfDofManagers();
        double rootSize, resolutionLimit;
        FloatArray minc(3), maxc(3), * coords;
        DofManager *dman;

        // first determine domain extends (bounding box), and check for degenerated domain type
        for ( i = 1; i <= nnode; i++ ) {
            dman = domain->giveDofManager(i);
            coords = static_cast< Node * >( dman )->giveCoordinates();
            if ( init ) {
                init = 0;
                for ( j = 1; j <= coords->giveSize(); j++ ) {
                    minc.at(j) = maxc.at(j) = coords->at(j);
                }
            } else {
                for ( j = 1; j <= coords->giveSize(); j++ ) {
                    if ( coords->at(j) < minc.at(j) ) {
                        minc.at(j) = coords->at(j);
                    }

                    if ( coords->at(j) > maxc.at(j) ) {
                        maxc.at(j) = coords->at(j);
                    }
                }
            }
        }                 // end loop over nodes

        BBX.setOrigin(minc);

        // determine root size
        rootSize = 0.0;
        for ( i = 1; i <= 3; i++ ) {
            rootSize = 1.000001 * max( rootSize, maxc.at(i) - minc.at(i) );
        }

        BBX.setSize(rootSize);

        // check for degenerated domain
        resolutionLimit = min(1.e-3, rootSize / 1.e6);
        for ( i = 1; i <= 3; i++ ) {
            if ( ( maxc.at(i) - minc.at(i) ) > resolutionLimit ) {
                BBX.setMask(i, 1);
            } else {
                BBX.setMask(i, 0);
            }
        }
    }

    /// Builds the octree from domain's nodes - NOT IN USE
    int buildOctreeDataStructureFromNodes() {
        int i, j, init = 1, nnode = this->domain->giveNumberOfDofManagers();
        DofManager *dman;

        // measure time consumed by octree build phase
        clock_t sc = :: clock();
        clock_t ec;

        ec = :: clock();

        // compute max. tree depth
        int treeDepth = 0;
        this->giveMaxTreeDepthFrom(this->rootCell, treeDepth);
        // compute processor time used by the program
        long nsec = ( ec - sc ) / CLOCKS_PER_SEC;
        OOFEM_LOG_INFO("Octree init [depth %d in %lds]\n", treeDepth, nsec);

        return 1;
    }

    /// Inserts member into the octree using functor for the evaluation
    int insertMemberIntoOctree(T &memberID, SL_Insertion_Functor< T > &functor)
    {
        int result = 0;

        if ( functor.evaluate(memberID, rootCell) ) {
            this->insertMemberIntoCell(memberID, functor, rootCell);
            result = 1;
        }

        return result;
    }

    /// Removes member from octree using insertion functor - NOT IN USE
    int removeMemberFromOctree(T &memberID, SL_Insertion_Functor< T > &functor)
    {
        this->removeMemberFromCell(memberID, functor, rootCell);
        return 1;
    }

    /// Evalutes the search accoring used functor a fills the list with results - NOT IN USE
    void giveDataOnFilter(std :: list< T > &answer, SL_Evaluation_Functor< T > &filter)
    {
        typename std :: list< T > *cellDataList;

        listIteratorType pos;
        BoundingBox BBS1, BBS2;
        typename std :: list< CellPtrType > *cellListPostSearch;
        typename std :: list< CellPtrType > :: iterator cellListPos;


        if ( filter.isBBXStage1Defined(BBS1) ) { }

        CellPtrType terminal;
        FloatArray startPos;
        filter.giveStartingPosition(startPos);
        terminal = this->findTerminalContaining(rootCell, startPos);

        cellDataList = terminal->giveDataList();
        for ( pos = cellDataList->begin(); pos != cellDataList->end(); ++pos ) {
            filter.evaluate(* pos);
        }

        if ( filter.isBBXStage2Defined(BBS2) ) {
            OctantRec :: BoundingBoxStatus BBStatus = terminal->testBoundingBox(BBS2);
            if ( BBStatus == OctantRec :: BBS_ContainsCell ) {
                giveListOfTerminalCellsInBoundingBox(* cellListPostSearch, BBS2, rootCell);
                for ( cellListPos = cellListPostSearch->begin(); cellListPos != cellListPostSearch->end(); ++cellListPos ) {
                    if ( * cellListPos != terminal ) {
                        cellDataList = ( * cellListPos )->giveDataList();
                        for ( pos = cellDataList->begin(); pos != cellDataList->end(); ++pos ) {
                            filter.evaluate(* pos);
                        }
                    }
                }
            }
        }

        filter.giveResult(answer);
    }

    /// Applies the evaluation functor, fills the answer list with results and removes found object from octree
    /// USED BY DELAUNAY TRIANGULATOR
    void proceedDataOnFilterAndRemoveFromOctree(std :: list< T > &answer, SL_Evaluation_Functor< T > &filter, SL_Insertion_Functor< T > &insertor, Timer &searchingTimer)
    {
        typename std :: list< T > *cellDataList;
        listIteratorType pos, positionInDataList;
        BoundingBox BBS1, BBS2;
        std :: list< LocalInsertionData< T > > *insertionDataList;
        typedef typename std :: list< LocalInsertionData< T > > :: iterator LIDiterator;

        if ( filter.isBBXStage1Defined(BBS1) ) { }

        CellPtrType terminal, cell;
        FloatArray startingPosition;
        filter.giveStartingPosition(startingPosition);
        terminal = this->findTerminalContaining(rootCell, startingPosition);

        cellDataList = terminal->giveDataList();
        for ( pos = cellDataList->begin(); pos != cellDataList->end(); ) {
            T _obj = * pos;
            bool ok = filter.evaluate(_obj);
            if ( ok ) {
                insertionDataList = insertor.giveInsertionList(_obj);
                if ( insertionDataList ) {
                    for ( LIDiterator insDataIT = insertionDataList->begin(); insDataIT != insertionDataList->end(); insDataIT++ ) {
                        cell = ( * insDataIT ).containedInCell;
                        positionInDataList = ( * insDataIT ).posInCellDataList;
                        cell->giveDataList()->erase(positionInDataList);
                    }
                } else {
                    // should not happen
                    removeMemberFromCell(_obj, insertor, terminal);
                }

                // actual iterator position might not be valid, so start from beginning
                pos = cellDataList->begin();
            } else {
                pos++;
            }
        }

        if ( filter.isBBXStage2Defined(BBS2) ) { }

        filter.giveResult(answer);
    }

    const char *giveClassName() const { return "OctreeSpatialLocalizerT"; }
    int giveOctreeMaskValue(int indx) { return octreeMask.at(indx); }
    void giveOctreeMask(IntArray &answer) { answer = this->octreeMask; }
    CellPtrType giveRootCell() { return rootCell; }
    void giveReport()
    {
        // compute max. tree depth
        int treeDepth = 0;
        this->giveMaxTreeDepthFrom(rootCell, treeDepth);

        OOFEM_LOG_INFO("Octree structure [depth %d]\n", treeDepth);
    }

protected:
    /// Returns terminal octant cell containing node with coords
    CellPtrType findTerminalContaining(CellPtrType startCell, const FloatArray &coords) {
        int result;
        CellPtrType currCell = startCell;
        if ( startCell->containsPoint(coords) ) {
            // found terminal octant containing node
            while ( !currCell->isTerminalOctant() ) {
                result = currCell->giveChildContainingPoint(& currCell, coords);
                if ( result == -2 ) {
                    // report some error
                }
            }

            return currCell;
        } else {
            return NULL;
        }
    }

    /// Returns the depth of the cell
    int giveCellDepth(CellPtrType cell) {
        return ( int ) ( log( this->rootCell->giveSize() / cell->giveSize() ) / M_LN2 );
    }

    /// Inserts member into the cell by evaluating the insertion functor
    /// Method is called recursivelly until terminal cell is reached
    int insertMemberIntoCell(T &memberID, SL_Insertion_Functor< T > &functor, CellPtrType cell)
    {
        int i, j, k;
        int nCellItems, cellDepth;
        dataContainerType *cellDataList;
        listIteratorType pos;
        listIteratorType insertedPosition;
        LocalInsertionData< T >LIdata;
        std :: list< LocalInsertionData< T > > *insData;
        typedef typename std :: list< LocalInsertionData< T > > :: iterator LIDiterator;

        if ( cell == NULL ) {
            return 0;
        }

        if ( cell->isTerminalOctant() ) {
            cellDataList = cell->giveDataList();
            nCellItems = cellDataList->size();
            cellDepth  = this->giveCellDepth(cell);
            if ( cellDepth > maxDepthReached ) {
                maxDepthReached = cellDepth;
                //printf("Reached cell depth: %i \n", maxDepthReached);
            }

            if ( ( nCellItems > TEMPLATED_OCTREE_MAX_NODES_LIMIT ) && ( cellDepth <= TEMPLATED_OCTREE_MAX_DEPTH ) ) {
                cell->divideLocally(1, this->octreeMask);
                for ( pos = cellDataList->begin(); pos != cellDataList->end(); ++pos ) {
                    this->insertMemberIntoCell(* pos, functor, cell);
                    insData = functor.giveInsertionList(* pos);
                    if ( insData ) {
                        for ( LIDiterator insDataIT = insData->begin(); insDataIT != insData->end(); ) {
                            if ( ( * insDataIT ).containedInCell == cell ) {
                                insDataIT = insData->erase(insDataIT);
                            } else {
                                insDataIT++;
                            }
                        }
                    }
                }

                cell->deleteDataList();
                this->insertMemberIntoCell(memberID, functor, cell);
            } else {                     // insertion without refinement
                insertedPosition = cell->addMember(memberID);
                LIdata.containedInCell = cell;
                LIdata.posInCellDataList = insertedPosition;

                functor.registerInsertion(memberID, LIdata);
                return 1;
            }
        } else {                 // not a terminal octant, visit children
            CellPtrType cptr;
            for ( i = 0; i <= 1; i++ ) {
                for ( j = 0; j <= 1; j++ ) {
                    for ( k = 0; k <= 1; k++ ) {
                        cptr = cell->giveChild(i, j, k);
                        if ( cptr ) {
                            if ( functor.evaluate( memberID, cell->giveChild(i, j, k) ) ) {
                                this->insertMemberIntoCell( memberID, functor, cell->giveChild(i, j, k) );
                            }
                        }
                    }
                }
            }
        }

        return 1;
    }

    /// Removes member from cell using insertion functor to ensure member is contained in
    int removeMemberFromCell(T &memberID, SL_Insertion_Functor< T > &functor, CellPtrType cell)
    {
        int i, j, k;
        listIteratorType pos;

        if ( cell ) {
            if ( cell->isTerminalOctant() ) {
                if ( functor.evaluate(memberID, cell) ) {
                    cell->removeMember(memberID);
                    return 1;
                } else {
                    return 0;
                }
            } else {
                for ( i = 0; i <= 1; i++ ) {
                    for ( j = 0; j <= 1; j++ ) {
                        for ( k = 0; k <= 1; k++ ) {
                            this->removeMemberFromCell( memberID, functor, cell->giveChild(i, j, k) );
                        }
                    }
                }
            }
        }

        return 1;
    }

    /// Gives the maximal tree depth from given cell
    void giveMaxTreeDepthFrom(CellPtrType root, int &maxDepth)
    {
        int i, j, k, depth = this->giveCellDepth(root);
        maxDepth = max(maxDepth, depth);

        for ( i = 0; i <= octreeMask.at(1); i++ ) {
            for ( j = 0; j <= octreeMask.at(2); j++ ) {
                for ( k = 0; k <= octreeMask.at(3); k++ ) {
                    if ( root->giveChild(i, j, k) ) {
                        this->giveMaxTreeDepthFrom(root->giveChild(i, j, k), maxDepth);
                    }
                }
            }
        }
    }

    /// Gives a list of terminal cells in a bounding box
    void giveListOfTerminalCellsInBoundingBox(std :: list< CellPtrType > &cellList, BoundingBox &BBX, CellPtrType currentCell)
    {
        int i, j, k;
        OctantRec :: BoundingBoxStatus BBStatus = currentCell->testBoundingBox(BBX);
        if ( ( BBStatus == OctantRec :: BBS_InsideCell ) || ( BBStatus == OctantRec :: BBS_ContainsCell ) ) {
            if ( currentCell->isTerminalOctant() ) {
                cellList.push_back(currentCell);
            } else {
                for ( i = 0; i <= octreeMask.at(1); i++ ) {
                    for ( j = 0; j <= octreeMask.at(2); j++ ) {
                        for ( k = 0; k <= octreeMask.at(3); k++ ) {
                            if ( currentCell->giveChild(i, j, k) ) {
                                this->giveListOfTerminalCellsInBoundingBox( cellList, BBX, currentCell->giveChild(i, j, k) );
                            }
                        }
                    }
                }
            }
        }
    }

    friend class OctantRecT< T >;
};
} // end namespace oofem
#endif // octreelocalizer_h
