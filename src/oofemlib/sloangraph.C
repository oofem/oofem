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

/* Modified and optimized by: Borek Patzak */
/* Author: Milan Jirasek */

#include "sloangraph.h"
#include "domain.h"
#include "element.h"
#include "dof.h"
#include "dofmanager.h"
#include "simpleslavedof.h"
#include "intarray.h"
#include "generalboundarycondition.h"

#include <memory>
#include <ctime>
#include <set>

namespace oofem {
SloanGraph :: SloanGraph(Domain *d)  : nodes(), queue(), OptimalRenumberingTable()
{
    domain  = d;
    WeightDistance = 1;
    WeightDegree = 2;
    SpineQuality  = Good;
    MinimalProfileSize    = 0;
    OptimalWeightDegree   = 0;
    OptimalWeightDistance = 0;
    startNode = endNode   = 0;
    nodeDistancesFlag     = 0;
}

SloanGraph :: ~SloanGraph()
{
}

void SloanGraph :: initialize()
{
    int nnodes = domain->giveNumberOfDofManagers();

    this->nodes.reserve(nnodes);
    this->dmans.reserve(nnodes);
    std :: map< DofManager*, int > dman2map;

    // Add dof managers.
    int count = 0;
    for ( auto &dman : domain->giveDofManagers() ) {
        nodes.emplace_back(this, ++count);
        dmans.push_back(dman.get());
        dman2map.insert({dman.get(), count});
    }
    // Add element internal dof managers
    for ( auto &elem : domain->giveElements() ) {
        this->nodes.reserve( nodes.size() + elem->giveNumberOfInternalDofManagers() );
        this->dmans.reserve( dmans.size() + elem->giveNumberOfInternalDofManagers() );
        for ( int j = 1; j <= elem->giveNumberOfInternalDofManagers(); ++j ) {
            nodes.emplace_back(this, ++count);
            dmans.push_back( elem->giveInternalDofManager(j) );
            dman2map.insert({elem->giveInternalDofManager(j), count});
        }
    }
    // Add boundary condition internal dof managers
    for ( auto &bc : domain->giveBcs() ) {
        this->nodes.reserve( nodes.size() + bc->giveNumberOfInternalDofManagers() );
        this->dmans.reserve( dmans.size() + bc->giveNumberOfInternalDofManagers() );
        for ( int j = 1; j <= bc->giveNumberOfInternalDofManagers(); ++j ) {
            nodes.emplace_back(this, ++count);
            dmans.push_back( bc->giveInternalDofManager(j) );
            dman2map.insert({bc->giveInternalDofManager(j), count});
        }
    }

    IntArray connections;
    for ( auto &elem : domain->giveElements() ) {
        int ielemnodes = elem->giveNumberOfDofManagers();
        int ielemintdmans = elem->giveNumberOfInternalDofManagers();
        int ndofmans = ielemnodes + ielemintdmans;
        connections.resize(ndofmans);
        for ( int j = 1; j <= ielemnodes; j++ ) {
            connections.at(j) = dman2map[ elem->giveDofManager(j) ];
        }
        for ( int j = 1; j <= ielemintdmans; j++ ) {
            connections.at(ielemnodes + j) = dman2map[ elem->giveInternalDofManager(j) ];
        }
        for ( int j = 1; j <= ndofmans; j++ ) {
            for ( int k = j + 1; k <= ndofmans; k++ ) {
                // Connect both ways
                this->giveNode( connections.at(j) ).addNeighbor( connections.at(k) );
                this->giveNode( connections.at(k) ).addNeighbor( connections.at(j) );
            }
        }
    }
    ///@todo Add connections from dof managers to boundary condition internal dof managers.

    IntArray dofMasters;
    for ( auto &dman : dmans ) {
        // according to forum discussion Peter & Mikael
        count = 0;
        ++count;
        if ( dman->hasAnySlaveDofs() ) {
            std :: set< int >masters;
            // slave dofs are present in dofManager
            // first - ask for masters, these may be different for each dof
            for ( Dof *dof: *dman ) {
                if ( !dof->isPrimaryDof() ) {
                    ///@todo FIXME This is inherently limited; It assumes that dofmanagers cannot be slaves to internal dofmanagers in BCs or Elements.
                    /// If we were able to get the masters dofmanagers are pointers we could fix this.
                    dof->giveMasterDofManArray(dofMasters);
                    for ( int m : dofMasters ) {
                        masters.insert( m );
                    }
                }
            }

            for ( int connection: masters ) {
                this->giveNode(count).addNeighbor(connection);
                this->giveNode(connection).addNeighbor(count);
            }
        }
    } // end dof man loop
}


SloanGraphNode &
SloanGraph :: giveNode(int num)
{
    return nodes[num-1];
}

int
SloanGraph :: giveNodeWithMinDegree()
{
    int nnodes = (int)nodes.size();
    int node_min = 0;
    int min = nnodes + 1;

    for ( int i = 1; i <= nnodes; i++ ) {
        int deg = this->giveNode(i).giveDegree();
        if ( deg < min && deg > 0 && ( this->giveNode(i).giveNewNumber() == 0 ) ) {
            min      = deg;
            node_min = i;
        }
    }

    return node_min;
}

void
SloanGraph :: findPeripheralNodes()
{
    //if (startNode && endNode) return;

    int InitialRoot;
    std :: unique_ptr< SloanLevelStructure > Spine;

    // clock_t time_0 = ::clock();
    if ( SpineQuality == Best ) {
        InitialRoot = findBestRoot();
    } else if ( SpineQuality == Good ) {
        InitialRoot = giveNodeWithMinDegree();
    } else {
        OOFEM_WARNING("Unsupported value of SpineQuality, using (Good)");
        InitialRoot = giveNodeWithMinDegree();
    }

    // initial spine
    Spine.reset( new SloanLevelStructure(this, InitialRoot) );
    this->startNode = InitialRoot;

    int CurrentDiameter = Spine->giveDepth();
    int TrialDepth, TrialWidth;
    std :: list< int >candidates;
    int newStartNode = 1;

    while ( newStartNode ) {
        newStartNode = 0;

        this->extractCandidates(candidates, *Spine);
        int MinimumWidth = domain->giveNumberOfDofManagers();

        for ( int Root: candidates ) {
            std :: unique_ptr< SloanLevelStructure > TrialSpine( new SloanLevelStructure(this, Root) );
            // abort level struc assembly if TrialWidth>MinimumWidth.
            if ( TrialSpine->formYourself(MinimumWidth) == 0 ) {
                continue;
            }

            TrialDepth = TrialSpine->giveDepth();
            TrialWidth = TrialSpine->giveWidth();
            if ( TrialWidth < MinimumWidth && TrialDepth > CurrentDiameter ) {
                Spine = std :: move( TrialSpine );
                this->startNode = Root;
                CurrentDiameter = Spine->giveDepth();
                newStartNode = 1;
                break; // exit for loop
            }

            if ( TrialWidth < MinimumWidth ) {
                MinimumWidth = TrialWidth;
                this->endNode = Root;
            }
        } // end loop over candidates
    }

    //clock_t time_1 = ::clock();

    //printf("Time needed to select the starting node %10lds\n",(time_1-time_0)/ CLOCKS_PER_SEC);
}


void
SloanGraph :: extractCandidates(std :: list< int > &candidates, SloanLevelStructure &Spine)
{
    int NumberOfLevels = Spine.giveDepth();
    IntArray LastLevel = Spine.giveLevel(NumberOfLevels);
    sort( LastLevel, SloanNodalDegreeOrderingCrit(this) );

    candidates.clear();
    // shrink candidates to contain only one node of each degree
    int lastDegree = 0;
    for ( int node: LastLevel ) {
        if ( lastDegree != this->giveNode( node ).giveDegree() ) {
            lastDegree = this->giveNode( node ).giveDegree();
            candidates.push_back( node );
        }
    }
}


int
SloanGraph :: computeTrueDiameter()
{
    SloanLevelStructure LSC( this, findBestRoot() );
    return LSC.giveDepth();
}

int
SloanGraph :: findBestRoot()
/* this is a very expensive method */
{
    int BestRoot = 0;
    int Diameter = 0;
    int nnodes = (int)nodes.size();
    clock_t time_1, time_0 = :: clock();
    for ( int i = 1; i <= nnodes; i++ ) {
        SloanLevelStructure LSC(this, i);
        int Depth = LSC.giveDepth();
        if ( Depth > Diameter ) {
            Diameter = Depth;
            BestRoot = i;
        }

        ///@todo Use the timer functionalities from OOFEM here, instead of clock() directly
        time_1 = :: clock();
        if ( ( time_1 - time_0 ) / CLOCKS_PER_SEC > SLOAN_TIME_CHUNK ) {
            OOFEM_LOG_INFO("%d roots (%5.1f per cent) checked: largest pseudo-diameter = %d\n", i, float ( 100 * i ) / nnodes, Diameter);
            //fflush(stdout);
            //if (get_yes_or_no() == NO) break;
            time_0 = time_1;
        }
    }

    return BestRoot;
}


double
SloanGraph :: giveOptimalProfileDensity()
{
    int nnodes = (int)nodes.size();
    double d = 0;
    if ( nnodes > 0 ) {
        d = 100.0 * MinimalProfileSize;
        d /= ( nnodes * ( nnodes + 1 ) );
    }

    return d;
}

void
SloanGraph :: initStatusAndPriority()
{
    this->findPeripheralNodes();
#ifndef MDC
    this->evaluateNodeDistances();
    int nnodes = domain->giveNumberOfDofManagers();

    for ( auto &node: nodes ) {
        int Distance = node->giveDistance();
        int Degree   = node->giveDegree();
        int Priority = WeightDistance * Distance - WeightDegree * ( Degree + 1 );
        node->setPriority(Priority);
        node->setStatus(SloanGraphNode :: Inactive);
    }

#else
    int NumLevels;
    int End = endNode;

    int Distance, Degree, Priority;

    SloanLevelStructure BackSpine(this, End);
    NumLevels = BackSpine.giveDepth();

    for ( int i = 1; i <= NumLevels; i++ ) {
        for ( int nodeNum: BackSpine.giveLevel(i) ) {
            Distance = i - 1;
            this->giveNode(nodeNum).setDistance(Distance);
            Degree   = this->giveNode(nodeNum).giveDegree();
            Priority = WeightDistance * Distance - WeightDegree * ( Degree + 1 );
            this->giveNode(nodeNum).setPriority(Priority);
            this->giveNode(nodeNum).setStatus(SloanGraphNode :: Inactive);
        }
    }

#endif
}


void
SloanGraph :: evaluateNodeDistances()
{
    int NumLevels;

    if ( this->nodeDistancesFlag ) {
        return;
    }

    int End = endNode;
    SloanLevelStructure BackSpine(this, End);
    NumLevels = BackSpine.giveDepth();

    for ( int i = 1; i <= NumLevels; i++ ) {
        for ( int nodeNum: BackSpine.giveLevel(i) ) {
            this->giveNode(nodeNum).setDistance(i - 1);
        }
    }

    this->nodeDistancesFlag = 1;
}

void
SloanGraph :: assignOldNumbers()
{
    for ( auto &node: nodes ) {
        node.assignOldNumber();
    }
}

#ifdef MDC
void
SloanGraph :: numberIsolatedNodes(int &NextNumber, int &labeledNodes)
{
    for ( auto &node: nodes ) {
        if ( node.giveDegree() == 0 ) {
            node.setNewNumber(++NextNumber);
            labeledNodes++;
        }
    }
}
#endif

void
SloanGraph :: assignNewNumbers()
{
    int Start, inext, NextNumber = 0;
    int labeledNodes = 0;

    for ( auto &node: nodes ) {
        node.setNewNumber(0);
    }

#ifdef MDC
    this->numberIsolatedNodes(NextNumber, labeledNodes);
#endif

    this->initStatusAndPriority();
    //std::list<int>::iterator next;

    Start = this->startNode;
    this->queue.clear();
    this->queue.push_back(Start);
    this->giveNode(Start).setStatus(SloanGraphNode :: Preactive);

#ifdef MDC
    for ( ; ; ) {
#endif
    while ( !this->queue.empty() ) {
        // finds top priority, returns the corresponding node and deletes the entry
        inext = findTopPriorityInQueue();
        // this->queue.erase(next); - done by findTopPriority
        SloanGraphNode &nextNode = this->giveNode(inext);
        if ( nextNode.giveStatus() == SloanGraphNode :: Preactive ) {
            this->insertNeigborsOf(inext);
        }

        nextNode.setNewNumber(++NextNumber);
        nextNode.setStatus(SloanGraphNode :: Postactive);
        modifyPriorityAround(inext);
        labeledNodes++;
    }

#ifdef MDC
    if ( labeledNodes == domain->giveNumberOfDofManagers() ) {
        break;
    }

    // proceed next subdomain
    this->initStatusAndPriority();
    Start = this->startNode;
    this->queue.clear();
    this->queue.push_back(Start);
    this->giveNode(Start).setStatus(SloanGraphNode :: Preactive);
}
#else
    if ( labeledNodes != domain->giveNumberOfDofManagers() ) {
        OOFEM_ERROR("Internal error:\n%s", "Isolated nodes or separated sub-domains exist");
    }

#endif
}


void
SloanGraph :: insertNeigborsOf(int next)
{
    SloanGraphNode &node = this->giveNode(next);

    for ( int nodeNum: node.giveNeighborList() ) {
        if ( this->giveNode(nodeNum).giveStatus() == SloanGraphNode :: Inactive ) {
            this->giveNode(nodeNum).setStatus(SloanGraphNode :: Preactive);
            this->queue.push_front(nodeNum);
        }
    }
}

void
SloanGraph :: modifyPriorityAround(int next)
{
    SloanGraphNode :: SloanGraphNode_StatusType status;
    SloanGraphNode &nextNode = this->giveNode(next);

    for ( int nodeNum: nextNode.giveNeighborList() ) {
        SloanGraphNode &neighborNode = this->giveNode(nodeNum);
        if ( neighborNode.giveStatus() == SloanGraphNode :: Preactive ) {
            neighborNode.increasePriorityBy(WeightDegree);
            neighborNode.setStatus(SloanGraphNode :: Active);
            for ( int nnodeNum: neighborNode.giveNeighborList() ) {
                SloanGraphNode &neighborOfNeighbor = this->giveNode(nnodeNum);
                status = neighborOfNeighbor.giveStatus();
                if ( status == SloanGraphNode :: Active || status == SloanGraphNode :: Preactive ) {
                    neighborOfNeighbor.increasePriorityBy(WeightDegree);
                } else if ( status == SloanGraphNode :: Inactive ) {
                    neighborOfNeighbor.increasePriorityBy(WeightDegree);
                    this->queue.push_front(nnodeNum);
                    neighborOfNeighbor.setStatus(SloanGraphNode :: Preactive);
                }
            }
        }
    }
}

int
SloanGraph :: findTopPriorityInQueue()
{
    int candidate = 0, priority, pmax = -WeightDegree * ( domain->giveNumberOfDofManagers() + 1 );
    std :: list< int > :: iterator toDel;

    for ( auto pos = queue.begin(); pos != queue.end(); ++pos ) {
        priority = this->giveNode(* pos).givePriority();
        if ( pmax < priority ) {
            pmax = priority;
            toDel = pos;
            candidate = * pos;
        }
    }

    if ( candidate ) {
        this->queue.erase(toDel);
    }

    return candidate;
}


int
SloanGraph :: computeProfileSize()
{
    int ProfSize = 0;
    //clock_t time_1, time_0 = ::clock();

    if ( WeightDistance || WeightDegree ) {
        assignNewNumbers();
    } else {
        assignOldNumbers();
    }

    for ( auto &node: nodes ) {
        ProfSize += node.computeProfileHeight();
    }

    //time_1 = ::clock();
    //printf("Time needed to renumber the nodes       %10lds\n",(time_1-time_0)/CLOCKS_PER_SEC);

    return ProfSize;
}


void
SloanGraph :: askNewOptimalNumbering(TimeStep *tStep)
{
    for ( int dmanNum: OptimalRenumberingTable ) {
        dmans[ dmanNum - 1 ]->askNewEquationNumbers(tStep);
    }
}


void
SloanGraph :: writeRenumberingTable(FILE *file)
{
    int nnodes = (int)nodes.size();

    for ( int i = 1; i <= nnodes; i++ ) {
        int inew = this->giveNode(i).giveNewNumber();
        fprintf(file, "%8i %8i\n", i, inew);
    }
}

int
SloanGraph :: writeOptimalRenumberingTable(FILE *OutputFile)
{
    if ( OptimalRenumberingTable.isEmpty() ) {
        return 0;
    }

    int nnodes = OptimalRenumberingTable.giveSize();
    for ( int i = 1; i <= nnodes; i++ ) {
        fprintf( OutputFile, "%8i %8i\n", i, OptimalRenumberingTable.at(i) );
    }

    //fclose(OutputFile);
    return 1;
}

void
SloanGraph :: printParameters()
{
    printf("\nCurrent parameter values:\n");
    printf("  1) weight of degree    = %d\n", WeightDegree);
    printf("  2) weight of distance  = %d\n", WeightDistance);
    printf("  3) diameter quality    = %d\n", SpineQuality);
}

void
SloanGraph :: setParameters(int wdeg, int wdis)
{
    if ( wdeg >= 0 ) {
        setWeightDegree(wdeg);
    }

    if ( wdis >= 0 ) {
        setWeightDistance(wdis);
    }
}

void
SloanGraph :: tryParameters(int wdeg, int wdis)
{

    setParameters(wdeg, wdis);
    int psize = computeProfileSize();
#ifndef SILENT
    // if (wdeg || wdis)
    //  printf("\nRenumbering with weights %2d / %2d",wdeg,wdis);
    // else
    //  printf("\nExisting node numbering         ");
    // printf(" profile size %d",psize);
#endif
    if ( psize < MinimalProfileSize || MinimalProfileSize == 0 ) {
        int nnodes = (int)nodes.size();
        MinimalProfileSize    = psize;
        OptimalWeightDegree   = wdeg;
        OptimalWeightDistance = wdis;
        OptimalRenumberingTable.resize(nnodes);
        for ( int i = 1; i <= nnodes; i++ ) {
            OptimalRenumberingTable.at( this->giveNode(i).giveNewNumber() ) = i;
        }
    }
}


int
SloanGraph :: giveFullProfileSize()
{
    int n = (int)nodes.size();
    return n * ( n + 1 );
}
} // end namespace oofem
