/* $Header: /home/cvs/bp/oofem/oofemlib/src/sloangraphnode.h,v 1.5 2003/04/06 14:08:25 bp Exp $ */
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

//   ************************************************
//   *** CLASS SLOAN PROFILE OPTIMIZER GRAPH NODE ***
//   ************************************************


/* Modified and optimized by: Borek Patzak */
/* Author: Milan Jirasek */

#ifndef sloangraphnode_h
#define sloangraphnode_h

#include "dynalist.h"

class SloanGraph;

/**
 * Class representing node in undirected graph, used by Sloan profile optimizer.
 * Node kepps its old and new number, status, degree, distance and corresponding priority.
 * The list of neigbouring nodes is also kept.
 */
class SloanGraphNode
{
public:
    /// Status type definition
    enum SloanGraphNode_StatusType { Inactive, Preactive, Active, Postactive };

private:

    /// Associated graph structure, to which node belongs
    SloanGraph *graph;
    /// Old (original) number
    int NumberOld;
    /// New (Optimized) number
    int NumberNew;
    /// Status of node
    SloanGraphNode_StatusType nodeStatus;
    /// Node degree (number of adjacent edges)
    int Degree;
    /// Node distance
    int Distance;
    /// Node priority
    int Priority;
    /// List of neighbouring nodes (represent graph edges)
    dynaList< int >neighborList;

public:

    /// Constructor. Creates node belonging to given graph with given old number.
    SloanGraphNode(SloanGraph *graph, int numOld);
    /// Destructor
    ~SloanGraphNode();

    /**
     * Add neighbouring node to corresponding list.
     * The entry is added only if not exist before,
     * the Degree member is updated accordingly.
     */
    void  addNeighbor(int neighbor);

    /// Returns new number of receiver
    int   giveNewNumber() { return NumberNew; }
    /// Returns old number of receiver
    int   giveOldNumber() { return NumberOld; }
    /// Returns receiver status
    SloanGraphNode_StatusType giveStatus() { return nodeStatus; }
    /// Return the receiver's degree
    int  giveDegree() { return Degree; }
    /// Returns distance of receiver
    int  giveDistance() { return Distance; }
    /// Returns priority of receiver
    int  givePriority() { return Priority; }
    /// Returns the neighbor list of receiver
    dynaList< int > *giveNeighborList()  { return & neighborList; }
    /// sets new number equal to old one
    void assignOldNumber() { NumberNew = NumberOld; }

    /// Sets the receiver distance to given number
    void  setDistance(int d) { Distance = d; }
    /// Sets the receiver priority to given value
    void  setPriority(int p) { Priority = p; }
    /// Sets the status of receiver to given value
    void  setStatus(SloanGraphNode_StatusType s) { nodeStatus = s; }
    /// Sets the new number of receiver
    void  setNewNumber(int n) { NumberNew = n; }
    /// Increases the priority of receiver by given value
    void  increasePriorityBy(int p) { Priority += p; }
    /**
     * Computes the profile height corresponding to receiver from
     * current new and old numbers.
     */
    int   computeProfileHeight();
};

#endif // sloangraphnode_h
