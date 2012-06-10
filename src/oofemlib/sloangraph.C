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

/* Modified and optimized by: Borek Patzak */
/* Author: Milan Jirasek */

#include "sloangraph.h"
#include "domain.h"
#include "element.h"
#include "dof.h"
#include "dofmanager.h"
#include "simpleslavedof.h"
#include "intarray.h"
#include "generalbc.h"

#ifndef __MAKEDEPEND
 #include <time.h>
 #include <set>
#endif

namespace oofem {
SloanGraph :: SloanGraph(Domain *d)  : nodes(0), queue(), OptimalRenumberingTable()
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
    dmans.clear(false); // Otherwise AList will delete the dof managers!
}

void SloanGraph :: initialize()
{
    int i, j, k, ielemnodes, ielemintdmans, ndofmans;
    int nnodes = domain->giveNumberOfDofManagers();
    int nelems = domain->giveNumberOfElements();
    int nbcs = domain->giveNumberOfBoundaryConditions();
    Element *ielem;
    GeneralBoundaryCondition *ibc;

    ///@todo Use dynaList for this first part instead (suboptimization?)
    this->nodes.growTo(nnodes);

    // Add dof managers.
    for ( i = 1; i <= nnodes; i++ ) {
        SloanGraphNode *node = new SloanGraphNode(this, i);
        nodes.put(i, node);
        dmans.put(i, domain->giveDofManager(i) );
    }
    k = nnodes;
    // Add element internal dof managers
    for ( i = 1; i <= nelems; i++ ) {
        ielem = domain->giveElement(i);
        this->nodes.growTo(k+ielem->giveNumberOfInternalDofManagers());
        for ( j = 1; j <= ielem->giveNumberOfInternalDofManagers(); ++j ) {
            SloanGraphNode *node = new SloanGraphNode(this, i);
            nodes.put(++k, node);
            dmans.put(++k, ielem->giveInternalDofManager(i) );
        }
    }
    // Add boundary condition internal dof managers
    for ( i = 1; i <= nbcs; i++ ) {
        ibc = domain->giveBc(i);
        if (ibc) {
            this->nodes.growTo(k+ibc->giveNumberOfInternalDofManagers());
            for ( j = 1; j <= ibc->giveNumberOfInternalDofManagers(); ++j ) {
                SloanGraphNode *node = new SloanGraphNode(this, i);
                nodes.put(++k, node);
                dmans.put(++k, ibc->giveInternalDofManager(i) );
            }
        }
    }

    IntArray connections;
    for ( i = 1; i <= nelems; i++ ) {
        ielem = domain->giveElement(i);
        ielemnodes = ielem->giveNumberOfDofManagers();
        ielemintdmans = ielem->giveNumberOfInternalDofManagers();
        ndofmans = ielemnodes + ielemintdmans;
        connections.resize(ndofmans);
        for ( j = 1; j <= ielemnodes; j++ ) {
            connections.at(j) = ielem->giveDofManager(j)->giveNumber();
        }
        for ( j = 1; j <= ielemintdmans; j++ ) {
            connections.at(ielemnodes+j) = ielem->giveInternalDofManager(j)->giveNumber();
        }
        for ( j = 1; j <= ndofmans; j++ ) {
            for ( k = j + 1; k <= ndofmans; k++ ) {
                // Connect both ways
                this->giveNode( connections.at(j) )->addNeighbor( connections.at(k) );
                this->giveNode( connections.at(k) )->addNeighbor( connections.at(j) );
            }
        }
    }
    ///@todo Add connections from dof managers to boundary condition internal dof managers.

    // loop over dof managers and test if there are some "slave" or rigidArm connection
    // if yes, such dependency is reflected in the graph by introducing additional
    // graph edges between slaves and corresponding masters
    /*
     * DofManager* iDofMan;
     * for (i=1; i <= nnodes; i++){
     * if (domain->giveDofManager (i)->hasAnySlaveDofs()) {
     * iDofMan = domain->giveDofManager (i);
     * if (iDofMan->giveClassID() == RigidArmNodeClass) {
     *   // rigid arm node -> has only one master
     *   int master = ((RigidArmNode*)iDofMan)->giveMasterDofMngr()->giveNumber();
     *   // add edge
     *   this->giveNode(i)->addNeighbor (master);
     *   this->giveNode(master)->addNeighbor(i);
     *
     * } else {
     * // slave dofs are present in dofManager
     * // first - ask for masters, these may be different for each dof
     *   int j;
     *   for (j=1; j<=iDofMan->giveNumberOfDofs(); j++)
     *     if (iDofMan->giveDof (j)->giveClassID() == SimpleSlaveDofClass) {
     *       int master = ((SimpleSlaveDof*) iDofMan->giveDof (j))->giveMasterDofManagerNum();
     *       // add edge
     *       this->giveNode(i)->addNeighbor (master);
     *       this->giveNode(master)->addNeighbor(i);
     *
     *     }
     * }
     * }
     * } // end dof man loop */

    std :: set< int, std :: less< int > >masters;
    std :: set< int, std :: less< int > > :: iterator it;

    IntArray dofMasters;
    DofManager *iDofMan;
    for ( i = 1; i <= nnodes; i++ ) {
        if ( domain->giveDofManager(i)->hasAnySlaveDofs() ) {
            // slave dofs are present in dofManager
            // first - ask for masters, these may be different for each dof
            masters.clear();
            iDofMan = domain->giveDofManager(i);
            int j, k;
            for ( j = 1; j <= iDofMan->giveNumberOfDofs(); j++ ) {
                if ( !iDofMan->giveDof(j)->isPrimaryDof() ) {
                    iDofMan->giveDof(j)->giveMasterDofManArray(dofMasters);
                    for ( k = 1; k <= dofMasters.giveSize(); k++ ) {
                        masters.insert( dofMasters.at(k) );
                    }
                }
            }

            for ( it = masters.begin(); it != masters.end(); ++it ) {
                this->giveNode(i)->addNeighbor( * ( it ) );
                this->giveNode( * ( it ) )->addNeighbor(i);
            }
        }
    } // end dof man loop

}


SloanGraphNode *
SloanGraph :: giveNode(int num)
{
    SloanGraphNode *node = nodes.at(num);
    return node;
}

int
SloanGraph :: giveNodeWithMinDegree()
{
    int nnodes = domain->giveNumberOfDofManagers();
    int node_min = 0;
    int min = nnodes + 1;
    int i, deg;

    for ( i = 1; i <= nnodes; i++ ) {
        deg = this->giveNode(i)->giveDegree();
        if ( deg < min && deg > 0 && ( this->giveNode(i)->giveNewNumber() == 0 ) ) {
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
    SloanLevelStructure *Spine;

    // clock_t time_0 = ::clock();
    if ( SpineQuality == Best ) {
        InitialRoot = findBestRoot();
    } else if ( SpineQuality == Good ) {
        InitialRoot = giveNodeWithMinDegree();
    } else {
        OOFEM_WARNING("SloanGraph::findPeripheralNodes : Unsupported value of SpineQuality, using (Good)");
        InitialRoot = giveNodeWithMinDegree();
    }

    // initial spine
    Spine = new SloanLevelStructure(this, InitialRoot);
    this->startNode = InitialRoot;

    int MinimumWidth = domain->giveNumberOfDofManagers();
    int CurrentDiameter = Spine->giveDepth();
    int TrialDepth, TrialWidth, Root;
    dynaList< int >candidates;
    dynaList< int > :: iterator pos;
    int newStartNode = 1;

    while ( newStartNode ) {
        newStartNode = 0;

        this->extractCandidates(candidates, Spine);
        MinimumWidth = domain->giveNumberOfDofManagers();

        for ( pos = candidates.begin(); pos != candidates.end(); ++pos ) {
            Root = * pos;
            SloanLevelStructure *TrialSpine = new SloanLevelStructure(this, Root);
            // abort level struc assembly if TrialWidth>MinimumWidth.
            if ( TrialSpine->formYourself(MinimumWidth) == 0 ) {
                delete TrialSpine;
                continue;
            }

            TrialDepth = TrialSpine->giveDepth();
            TrialWidth = TrialSpine->giveWidth();
            if ( TrialWidth < MinimumWidth && TrialDepth > CurrentDiameter ) {
                delete Spine;
                Spine = TrialSpine;
                this->startNode = Root;
                CurrentDiameter = Spine->giveDepth();
                newStartNode = 1;
                break; // exit for loop
            }

            if ( TrialWidth < MinimumWidth ) {
                MinimumWidth = TrialWidth;
                this->endNode = Root;
            }

            delete TrialSpine;
        } // end loop over candidates

    }

    delete Spine;

    //clock_t time_1 = ::clock();

    //printf("Time needed to select the starting node %10lds\n",(time_1-time_0)/ CLOCKS_PER_SEC);
}


void
SloanGraph :: extractCandidates(dynaList< int > &candidates, SloanLevelStructure *Spine)
{
    if ( !Spine ) {
        OOFEM_ERROR("SloanGraph::extractCandidates : Invalid spine");
    }

    int i, NumberOfLevels = Spine->giveDepth();
    IntArray *LastLevel = Spine->giveLevel(NumberOfLevels);
    sort( * LastLevel, SloanNodalDegreeOrderingCrit(this) );

    candidates.clear();
    // shrink candidates to contain only one node of each degree
    int lastDegree = 0;
    for ( i = 1; i <= LastLevel->giveSize(); i++ ) {
        if ( lastDegree != this->giveNode( LastLevel->at(i) )->giveDegree() ) {
            lastDegree = this->giveNode( LastLevel->at(i) )->giveDegree();
            candidates.pushBack( LastLevel->at(i) );
        }
    }
}


int
SloanGraph :: computeTrueDiameter()
{
    SloanLevelStructure *LSC = new SloanLevelStructure( this, findBestRoot() );
    int Diameter = LSC->giveDepth();
    delete LSC;
    return Diameter;
}

int
SloanGraph :: findBestRoot()
/* this is a very expensive method */
{
    int BestRoot = 0;
    int Diameter = 0;
    int i, nnodes = domain->giveNumberOfDofManagers();
    SloanLevelStructure *LSC;
    clock_t time_1, time_0 = :: clock();
    for ( i = 1; i <= nnodes; i++ ) {
        LSC = new SloanLevelStructure(this, i);
        int Depth = LSC->giveDepth();
        if ( Depth > Diameter ) {
            Diameter = Depth;
            BestRoot = i;
        }

        delete LSC;
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
    int nnodes = domain->giveNumberOfDofManagers();
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
    int i;

    this->findPeripheralNodes();
#ifndef MDC
    this->evaluateNodeDistances();
    int nnodes = domain->giveNumberOfDofManagers();

    for ( i = 1; i <= nnodes; i++ ) {
        int Distance = giveNode(i)->giveDistance();
        int Degree   = giveNode(i)->giveDegree();
        int Priority = WeightDistance * Distance - WeightDegree * ( Degree + 1 );
        giveNode(i)->setPriority(Priority);
        giveNode(i)->setStatus(SloanGraphNode :: Inactive);
    }

#else
    int NumLevels;
    int j;
    IntArray *Level;
    int End = endNode;

    int Distance, Degree, Priority;

    SloanLevelStructure BackSpine(this, End);
    NumLevels = BackSpine.giveDepth();

    for ( i = 1; i <= NumLevels; i++ ) {
        Level = BackSpine.giveLevel(i);
        for ( j = 1; j <= Level->giveSize(); j++ ) {
            Distance = i - 1;
            this->giveNode( Level->at(j) )->setDistance(Distance);
            Degree   = this->giveNode( Level->at(j) )->giveDegree();
            Priority = WeightDistance * Distance - WeightDegree * ( Degree + 1 );
            this->giveNode( Level->at(j) )->setPriority(Priority);
            this->giveNode( Level->at(j) )->setStatus(SloanGraphNode :: Inactive);
        }
    }

#endif
}


void
SloanGraph :: evaluateNodeDistances()
{
    int NumLevels;
    int i, j;
    IntArray *Level;

    if ( this->nodeDistancesFlag ) {
        return;
    }

    int End = endNode;
    SloanLevelStructure BackSpine(this, End);
    NumLevels = BackSpine.giveDepth();

    for ( i = 1; i <= NumLevels; i++ ) {
        Level = BackSpine.giveLevel(i);
        for ( j = 1; j <= Level->giveSize(); j++ ) {
            this->giveNode( Level->at(j) )->setDistance(i - 1);
        }
    }

    this->nodeDistancesFlag = 1;
}

void
SloanGraph :: assignOldNumbers()
{
    int i, nnodes = domain->giveNumberOfDofManagers();
    for ( i = 1; i <= nnodes; i++ ) {
        giveNode(i)->assignOldNumber();
    }
}

#ifdef MDC
void
SloanGraph :: numberIsolatedNodes(int &NextNumber, int &labeledNodes)
{
    int i, nnodes = domain->giveNumberOfDofManagers();
    for ( i = 1; i <= nnodes; i++ ) {
        if ( giveNode(i)->giveDegree() == 0 ) {
            giveNode(i)->setNewNumber(++NextNumber);
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

    int i, nnodes = domain->giveNumberOfDofManagers();
    for ( i = 1; i <= nnodes; i++ ) {
        giveNode(i)->setNewNumber(0);
    }

#ifdef MDC
    this->numberIsolatedNodes(NextNumber, labeledNodes);
#endif

    this->initStatusAndPriority();
    SloanGraphNode *nextNode;
    //dynaList<int>::iterator next;

    Start = this->startNode;
    this->queue.clear();
    this->queue.pushBack(Start);
    this->giveNode(Start)->setStatus(SloanGraphNode :: Preactive);

#ifdef MDC
    for ( ; ; ) {
#endif
    while ( !this->queue.isEmpty() ) {
        // finds top priority, returns the corresponding node and deletes the entry
        inext = findTopPriorityInQueue();
        // this->queue.erase(next); - done by findTopPriority
        nextNode = this->giveNode(inext);
        if ( nextNode->giveStatus() == SloanGraphNode :: Preactive ) {
            this->insertNeigborsOf(inext);
        }

        nextNode->setNewNumber(++NextNumber);
        nextNode->setStatus(SloanGraphNode :: Postactive);
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
    this->queue.pushBack(Start);
    this->giveNode(Start)->setStatus(SloanGraphNode :: Preactive);
}
#else
    if ( labeledNodes != domain->giveNumberOfDofManagers() ) {
        OOFEM_ERROR2("SloanGraph :: assignNewNumbers: Internal error:\n%s", "Isolated nodes or separated sub-domains exist");
    }

#endif
}


void
SloanGraph :: insertNeigborsOf(int next)
{
    int neighborNode;
    SloanGraphNode *node = this->giveNode(next);
    dynaList< int > *neighbors = node->giveNeighborList();
    dynaList< int > :: iterator pos;

    for ( pos = neighbors->begin(); pos != neighbors->end(); ++pos ) {
        neighborNode = * pos;
        if ( this->giveNode(neighborNode)->giveStatus() == SloanGraphNode :: Inactive ) {
            this->giveNode(neighborNode)->setStatus(SloanGraphNode :: Preactive);
            this->queue.pushFront(neighborNode);
        }
    }
}

void
SloanGraph :: modifyPriorityAround(int next)
{
    SloanGraphNode :: SloanGraphNode_StatusType status;
    SloanGraphNode *neighborNode, *neighborOfNeighbor, *nextNode = this->giveNode(next);
    dynaList< int > *neighborsOfneighbors, *neighbors = nextNode->giveNeighborList();
    dynaList< int > :: iterator pos, npos;

    for ( pos = neighbors->begin(); pos != neighbors->end(); ++pos ) {
        neighborNode = this->giveNode(* pos);
        if ( neighborNode->giveStatus() == SloanGraphNode :: Preactive ) {
            neighborNode->increasePriorityBy(WeightDegree);
            neighborNode->setStatus(SloanGraphNode :: Active);
            neighborsOfneighbors =  neighborNode->giveNeighborList();
            for ( npos = neighborsOfneighbors->begin(); npos != neighborsOfneighbors->end(); ++npos ) {
                neighborOfNeighbor = this->giveNode(* npos);
                status = neighborOfNeighbor->giveStatus();
                if ( status == SloanGraphNode :: Active || status == SloanGraphNode :: Preactive ) {
                    neighborOfNeighbor->increasePriorityBy(WeightDegree);
                } else if ( status == SloanGraphNode :: Inactive ) {
                    neighborOfNeighbor->increasePriorityBy(WeightDegree);
                    this->queue.pushFront(* npos);
                    neighborOfNeighbor->setStatus(SloanGraphNode :: Preactive);
                }
            }
        }
    }
}

int
SloanGraph :: findTopPriorityInQueue()
{
    int candidate = 0, priority, pmax = -WeightDegree * ( domain->giveNumberOfDofManagers() + 1 );
    dynaList< int > :: iterator pos, toDel;

    for ( pos = queue.begin(); pos != queue.end(); ++pos ) {
        priority = this->giveNode(* pos)->givePriority();
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
    int nnodes = this->domain->giveNumberOfDofManagers();
    //clock_t time_1, time_0 = ::clock();

    if ( WeightDistance || WeightDegree ) {
        assignNewNumbers();
    } else {
        assignOldNumbers();
    }

    int i;
    for ( i = 1; i <= nnodes; i++ ) {
        ProfSize += this->giveNode(i)->computeProfileHeight();
    }

    //time_1 = ::clock();
    //printf("Time needed to renumber the nodes       %10lds\n",(time_1-time_0)/CLOCKS_PER_SEC);

    return ProfSize;
}


void
SloanGraph :: askNewOptimalNumbering(TimeStep *tStep)
{
    int ndmans = dmans.giveSize();
    for ( int i = 1; i <= ndmans; ++i ) {
        dmans.at( OptimalRenumberingTable.at(i) )->askNewEquationNumbers(tStep);
    }
}


void
SloanGraph :: writeRenumberingTable(FILE *file)
{
    int i, inew, nnodes = this->domain->giveNumberOfDofManagers();

    for ( i = 1; i <= nnodes; i++ ) {
        inew = this->giveNode(i)->giveNewNumber();
        fprintf(file, "%8i %8i\n", i, inew);
    }

    //fclose(file);
}

int
SloanGraph :: writeOptimalRenumberingTable(FILE *OutputFile)
{
    if ( OptimalRenumberingTable.isEmpty() ) {
        return 0;
    }

    int i, nnodes = domain->giveNumberOfDofManagers();
    for ( i = 1; i <= nnodes; i++ ) {
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
    int i, nnodes = domain->giveNumberOfDofManagers();

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
        MinimalProfileSize    = psize;
        OptimalWeightDegree   = wdeg;
        OptimalWeightDistance = wdis;
        OptimalRenumberingTable.resize(nnodes);
        for ( i = 1; i <= nnodes; i++ ) {
            OptimalRenumberingTable.at( this->giveNode(i)->giveNewNumber() ) = i;
        }
    }
}


int
SloanGraph ::  giveFullProfileSize()
{
    int n = this->domain->giveNumberOfDofManagers();
    return n * ( n + 1 );
}
} // end namespace oofem
