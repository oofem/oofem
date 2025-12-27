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

//   ****************************************************
//   *** UTILITY CLASSES FOR OCTREE SPATIAL LOCALIZER ***
//   ****************************************************

#ifndef octreelocalizertutil_h
#define octreelocalizertutil_h

#include "spatiallocalizer.h"
#include "octreelocalizer.h"
#include "octreelocalizert.h"
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
    InsertNode(Domain *d) :
        domain(d)
    { }
    /// Destructor
    ~InsertNode() { }

    /**
     * Evaluates the position of a node
     * @param nodeNr Number of the node in the domain
     * @param cell Octant cell to be tested
     * @returns true if the node is contained in the cell, false otherwise
     */
    bool evaluate(int &nodeNr, OctantRecT< int > *cell) override
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

    void registerInsertion(int &nodeNr, LocalInsertionData< int >LIdata) override { }

    std :: list< LocalInsertionData< int > > *giveInsertionList(int &nodeNr) override
    {
        return nullptr;
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
    InsertTriangleBasedOnCircumcircle(Domain *d) :
        domain(d)
    { }
    /// Destructor
    ~InsertTriangleBasedOnCircumcircle() { }

    /**
     * Evaluates the position of a triangle
     * @param DTptr Delaunay triangle
     * @param cell Octant cell to be tested
     * @returns true if the cirumscribed circle of the circle lies inside or contains a part of the cell, false otherwise
     */
    bool evaluate(DelaunayTriangle * &DTptr, OctantRecT< DelaunayTriangle * > *cell) override
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
    void registerInsertion(DelaunayTriangle * &TEptr, LocalInsertionData< DelaunayTriangle * >LIdata) override
    {
        ( TEptr->giveListOfCellsAndPosition() )->push_back(LIdata);
    }
    std :: list< LocalInsertionData< DelaunayTriangle * > > *giveInsertionList(DelaunayTriangle * &DTptr) override
    {
        return DTptr->giveListOfCellsAndPosition();
    }
};

/**
 * Functor for closest node search
 *
 * @author David Krybus
 */
class ClosestNode : public SL_Evaluation_Functor< int >
{
protected:
    FloatArray startingPosition;
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
    ClosestNode(const FloatArray &pos, Domain *d) :
        startingPosition(pos),
        domain(d),
        initFlag(false)
    { }
    /// Destructor
    ~ClosestNode() { }

    /**
     * Gives the starting position of the search
     * @param position startingPosition
     */
    void giveStartingPosition(FloatArray &position) override
    {
        position = startingPosition;
    }

    /**
     * Evaluates a node. The closest nodes are stored in a container, if their distance to starting position is the same.
     * Is the distance smaller than previous one, the container is emptied and new node is added.
     * @param nodeNr Number of the node in the domain list
     * @returns true after evaluation is processed.
     */
    bool evaluate(int &nodeNr) override
    {
        if ( initFlag ) {
            double dist = distance(startingPosition, this->domain->giveNode(nodeNr)->giveCoordinates());

            if ( ( dist - distanceToClosestNode ) <= dist * 0.001 ) {
                if ( ( dist - distanceToClosestNode ) >= -0.001 * dist ) {
                    closestNodeIndices.push_back(nodeNr);
                } else {
                    closestNodeIndices.clear();
                    closestNodeIndices.push_back(nodeNr);
                    distanceToClosestNode = dist;
                }
            }
        } else {
            closestNodeIndices.push_back(nodeNr);
            distanceToClosestNode = distance(startingPosition, this->domain->giveNode( * ( closestNodeIndices.begin() ) )->giveCoordinates());
            initFlag = true;
        }

        return true;
    }

    /**
     * Gives the closest nodes
     * @param answer List containing numbers of the closest nodes
     */
    void giveResult(std :: list< int > &answer) override {
        answer = closestNodeIndices;
    }

    bool isBBXStage1Defined(BoundingBox &BBXStage1) override
    {
        return false;
    }

    bool isBBXStage2Defined(BoundingBox &BBXStage2) override
    {
        if ( initFlag ) {
            BBXStage2.setOrigin(startingPosition);
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
    FloatArray startingPosition;
    Domain *domain;
    std :: list< DelaunayTriangle * >result;

public:
    /**
     * Constructor
     * @param pos Starting position of the search
     * @param d Domain containing nodes the triangles are defined by
     */
    ElementCircumCirclesContainingNode(FloatArray pos, Domain *d) :
        startingPosition(std::move(pos)),
        domain(d)
    { }
    ~ElementCircumCirclesContainingNode() { }

    /**
     * Evaluates a triangle upon its circumscribed cricle.
     * @param DTptr Delaunay triangle. nodeNr Number of the node in the domain list
     * @returns true if the circumscribed circle of the Delaunay triangle contains the point, false otherwise
     */
    bool evaluate(DelaunayTriangle * &DTptr) override
    {
        double radius = DTptr->giveCircumRadius();
        FloatArray centerCoords(3);
        centerCoords.at(1) = DTptr->giveXCenterCoordinate();
        centerCoords.at(2) = DTptr->giveYCenterCoordinate();
        centerCoords.at(3) = 0.0;

        if ( distance(startingPosition, centerCoords) < radius ) {
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
    void giveStartingPosition(FloatArray &answer) override
    {
        answer = startingPosition;
    }

    /**
     * Gives the triangles containing the node
     * @param answer List containing Delaunay triangles
     */
    void giveResult(std :: list< DelaunayTriangle * > &answer) override
    {
        answer = result;
    }

    // neither stage1, nor stage 2 are defined
    bool isBBXStage1Defined(BoundingBox &BBXStage1) override { return false; }
    bool isBBXStage2Defined(BoundingBox &BBXStage2) override { return false; }
};

} // end namespace oofem
#endif // octreelocalizertutil_h
