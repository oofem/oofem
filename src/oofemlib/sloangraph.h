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

#ifndef sloangraph_h
#define sloangraph_h

#include "oofemcfg.h"
#include "sloangraphnode.h"
#include "sloanlevelstruct.h"
#include "dofmanager.h"
#include "intarray.h"

#include <list>

namespace oofem {
class TimeStep;

#define SLOAN_TIME_CHUNK 60
/**
 * If defined, the support for Multiple subdomains is turned on.
 * Then the algorithm properly handles the isolated nodes as well as
 * the case of multiple subdomains that are not connected
 */
#define MDC

class Domain;

/**
 * Graph representing the undirected graph used for Sloan algorithm for symmetric matrix profile
 * reduction. The Sloan algorithm is described by following papers:
 * Sloan, S., W: An algorithm for profile and wavefront reduction of sparse matrices, IJNME, vol. 23, 239-251, 1986
 * Sloan, S., W: A Fortran program for profile and wavefront reduction, IJNME, vol. 28, 2651-2679, 1989.
 * The current implementation undergoes the description of algorithm given in second paper.
 *
 * The rigid arm nodes and slave dofs are supported. A fictitious element (graph edge) connecting corresponding slave and
 * master nodes is added to reflect the connection.
 *
 * The current implementation allows to select quality.
 * The recommended value for quality is Good, causing Sloan pseudo-peripheral node selection algorithm is used.
 * The Best causes the level structure to be generated for all nodes and selecting the best one.
 * However, this is extremely time consuming.
 *
 * @author Milan Jirasek
 * @author Borek Patzak
 */
class OOFEM_EXPORT SloanGraph
{
public:
    /// Quality type definition.
    enum SpineQualityType { Good, Best };
private:
    /// Domain asoociated to graph
    Domain *domain;

    /// List of graph nodes.
    std::vector< SloanGraphNode >nodes;
    /// List of dof managers corresponding to nodes.
    std::vector< DofManager* >dmans;
    /// Start peripheral node.
    int startNode;
    /// End peripheral node.
    int endNode;
    /// Priority queue of active or preactive nodes.
    std :: list< int >queue;
    /// Integer distance weight.
    int WeightDistance;
    /// Integer degree weight.
    int WeightDegree;
    SpineQualityType SpineQuality;
    /// Minimal profile size obtained.
    int MinimalProfileSize;
    /// Optimal degree weight.
    int OptimalWeightDegree;
    /// Optimal distance weight.
    int OptimalWeightDistance;
    /// Flag indicating that node distances from endNode were already computed.
    int nodeDistancesFlag;

    /**
     * Inverse renumbering table. Contains at i-th position the old node number corresponding to
     * newly generated i-th node. Used directly in equation numbering, the nodes (still referenced by old numbers)
     * are numbered in the order given by OptimalRenumberingTable values.
     */
    IntArray OptimalRenumberingTable;

public:
    /// Constructor. Creates the graph associated to given domain.
    SloanGraph(Domain * d);
    /// Destructor
    ~SloanGraph();

    /// Returns associated domain
    Domain *giveDomain() { return this->domain; }
    /// Initialize graph from domain description
    void  initialize();
    /// Resets the receiver state. Clears the startNode, endNode and nodeDistancesFlag values.
    void resetAll() { startNode = endNode = nodeDistancesFlag = 0; }

    /// Return graph node
    SloanGraphNode &giveNode(int num);

    /// Finds the peripheral nodes (rooted in optimal start node) according to receiver quality and current weights.
    void findPeripheralNodes();

    int computeTrueDiameter();
    int findBestRoot();
    int giveFullProfileSize();
    /// Returns the optimal profile found.
    int giveOptimalProfileSize() { return MinimalProfileSize; }
    /// Returns the optimal density of mesh.
    double giveOptimalProfileDensity();

    /**
     * Assigns the New numbers by node labeling algorithm
     * (old numbers are used when both weights are zero)
     * and returns the profile size.
     */
    int computeProfileSize();
    /// Numbers all the DOFs according to the optimal renumbering found.
    void askNewOptimalNumbering(TimeStep *tStep);
    /// Returns the optimal reverse renumbering table
    IntArray &giveOptimalRenumberingTable() { return OptimalRenumberingTable; }
    void writeRenumberingTable(FILE *file);
    int writeOptimalRenumberingTable(FILE *file);

    /// Sets weight distance to given value.
    void setWeightDistance(int w) {
        if ( w >= 0 ) {
            WeightDistance = w;
        }
    }
    /// Sets weight degree to given value.
    void setWeightDegree(int w) {
        if ( w >= 0 ) {
            WeightDegree = w;
        }
    }
    /// Select spine quality generation.
    void  setSpineQuality(SpineQualityType q) {
        SpineQuality = q;
        resetAll();
    }

    /// Prints actual parameters
    void printParameters();
    /// Sets weight degee and weight dist to given values
    void setParameters(int wdeg, int wdis);
    /**
     * Generates the new nodal numbering based on given parameters.
     * Can be used multiple times, the best result is stored in OptimalRenumberingTable,
     * MinimalProfileSize, OptimalWeightDegree and OptimalWeightDistance attributes.
     */
    void tryParameters(int wdeg, int wdis);

private:
    /// Returns graph node number with minimal degree.
    int giveNodeWithMinDegree();
    /**
     * Extract candidates from given level structure.
     * The list of candidates contains only one node of each degree
     * of last level of active spine.
     */
    void extractCandidates(std :: list< int > &candidates, SloanLevelStructure &Spine);
    /// Initializes statuses and priority of nodes of receiver
    void initStatusAndPriority();
    /// Evaluates the nodal distances from backSpine. The backSpine is generated if not available.
    void evaluateNodeDistances();
    /// Assigns old node numbers as new ones. Used to compute the profile of existing old numbering.
    void assignOldNumbers();
    /// Implementation of node labeling algorithm.
    void assignNewNumbers();
    /// Inserts inactive neighbours of given node as preactive ones.
    void insertNeigborsOf(int);
    /// Modifies the priority around node with max priority.
    void modifyPriorityAround(int);
    /// Finds node with highest priority in queue, removes its entry and returns its number.
    int findTopPriorityInQueue();
#ifdef MDC
    /// Numbers isolated nodes (those with degree equal to 0).
    void numberIsolatedNodes(int &NextNumber, int &labeledNodes);
#endif
};

class OOFEM_EXPORT SloanNodalDegreeOrderingCrit
{
    SloanGraph *graph;
public:
    SloanNodalDegreeOrderingCrit(SloanGraph * g) {
        graph = g;
    }
    int operator() (const int n1, const int n2) {
        return graph->giveNode(n1).giveDegree() - graph->giveNode(n2).giveDegree();
    }
};
} // end namespace oofem
#endif // sloangraph_h
