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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
//
//   Generic templated Octree localizer

#ifndef octreelocalizert_h
#define octreelocalizert_h

#include "spatiallocalizer.h"
#include "octreelocalizer.h"
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


#define TEMPLATED_OCTREE_MAX_NODES_LIMIT 500

#define TEMPLATED_OCTREE_MAX_DEPTH 8

//#define TEMPLATED_OCTREE_DEBUG 

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
    BoundingBox() : size(0), spatialMask(3)  { }

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
    /**
     * Sets all BBOx parameters in ince
     */
    void init (FloatArray& _origin, double _size, IntArray &mask)
    {
      this->origin = _origin;
      this->size   = _size;
      this->spatialMask = mask;
    }
    bool contains(const FloatArray& coords) const
    {
        for ( int i = 1; i <= coords.giveSize(); i++ ) {
            if ( spatialMask.at(i) ) {
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
    OctantRec :: BoundingBoxStatus testBoundingBox(BoundingBox &A)
    {
        double aSize = A.giveSize();
        FloatArray sA,eA;
        A.giveOrigin(sA);
        eA = sA;
        eA.add(aSize);

        bool sAInside = this->containsPoint(sA);
        bool eAInside = this->containsPoint(eA);
        if (sAInside && eAInside) {
            return OctantRec :: BBS_InsideCell;
        } else if (sAInside || eAInside) {
            return OctantRec :: BBS_ContainsCell;
        } else {
            return OctantRec :: BBS_OutsideCell;
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
    //Domain *domain;
    int maxDepthReached;

public:
    /// Constructor
    OctreeSpatialLocalizerT() {
        rootCell = NULL;
        maxDepthReached = 0;
    }
    /// Destructor
    ~OctreeSpatialLocalizerT() {
        if ( rootCell ) {
            delete rootCell;
        }
    }

    void clear () {
      if (rootCell) delete rootCell;
      rootCell = NULL;
      octreeMask.zero();
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
        typename std :: list< CellPtrType > *cellListPostSearch = NULL;
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
#ifdef TEMPLATED_OCTREE_DEBUG
                printf("Reached cell depth: %i \n", maxDepthReached);
#endif
            }

            if ( ( nCellItems > TEMPLATED_OCTREE_MAX_NODES_LIMIT ) && ( cellDepth < TEMPLATED_OCTREE_MAX_DEPTH ) ) {
                cell->divideLocally(1, this->octreeMask);


#if 1 // more memory efficient implementation
                while (!cellDataList->empty())
                  {
                    this->insertMemberIntoCell(cellDataList->back(), functor, cell);
                    insData = functor.giveInsertionList(cellDataList->back());
                    if ( insData ) {
                        for ( LIDiterator insDataIT = insData->begin(); insDataIT != insData->end(); ) {
                            if ( ( * insDataIT ).containedInCell == cell ) {
                                insDataIT = insData->erase(insDataIT);
                            } else {
                                insDataIT++;
                            }
                        }
                    }
                    cellDataList->pop_back();
                  }
#else
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
#endif
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
