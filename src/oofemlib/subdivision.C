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

#include "subdivision.h"
#include "material.h"
#include "errorestimator.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "mathfem.h"
#include "masterdof.h"
#include "simpleslavedof.h"
#include "nonlocalbarrier.h"
#include "initialcondition.h"
#include "classfactory.h"
#include "outputmanager.h"
#include "crosssection.h"
#include "function.h"
#include "timer.h"
#include "remeshingcrit.h"
#include "dynamicinputrecord.h"
#include "engngm.h"

#include <queue>
#include <set>

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

#ifdef __PARALLEL_MODE
 #include "parallel.h"
 #include "problemcomm.h"
 #include "communicator.h"
 #include "datastream.h"
 #include "domaintransactionmanager.h"
#endif

namespace oofem {
//#define __VERBOSE_PARALLEL
#define DEBUG_CHECK
//#define DEBUG_INFO
//#define DEBUG_SMOOTHING

// QUICK_HACK enables comparison of sequential and parallel version when smoothing is applied;
// sequential nodes on parallel interface are marked as boundary (in input file);
// important: there is a limitation:
// single element is allowed to have only 2 (in 2D) or 3 (in 3D) nodes marked as boundary
// if quick hack is used, standard identification of boundary nodes during bisection is
// replaced by simplified version;
// also smoothing of nodes is driven by coordinates based ordering which should ensure
// that sequential nodes (on individual parallel partitions) are smoothed in the same
// order as parallel nodes

//#define QUICK_HACK

#ifdef __OOFEG
 #define DRAW_IRREGULAR_NODES
 #define DRAW_REMOTE_ELEMENTS
  //#define DRAW_MESH_BEFORE_BISECTION
 #define DRAW_MESH_AFTER_BISECTION
  //#define DRAW_MESH_AFTER_EACH_BISECTION_LEVEL
#endif



//#define NM   // Nooru Mohamed
//#define BRAZIL_2D   // splitting brazil test
//#define THREEPBT_3D   // 3pbt
#define HEADEDSTUD    // headed stud



/* builds connectivity of top level node only however at any level of the subdivision */

int Subdivision :: RS_Node :: buildTopLevelNodeConnectivity(ConnectivityTable *ct)
{
    IntArray me(1), connElems;
    int i, el;
    me.at(1) = this->giveNumber();

    ct->giveNodeNeighbourList(connElems, me);
    this->connectedElements.preallocate( connElems.giveSize() );    // estimate size of the list

    for ( i = 1; i <= connElems.giveSize(); i++ ) {
        el = connElems.at(i);

#ifdef __PARALLEL_MODE
        if ( this->mesh->giveElement(el)->giveParallelMode() != Element_local ) {
            continue;
        }
#endif

        // use nonzero chunk, the estimated size may be not large enough
        this->mesh->giveElement(el)->buildTopLevelNodeConnectivity(this);
    }

    return ( connectedElements.giveSize() );
}


void Subdivision :: RS_Element :: buildTopLevelNodeConnectivity(Subdivision :: RS_Node *node)
{
    int el, i;
    Subdivision :: RS_Element *elem;

    if ( this->isTerminal() ) {
        node->insertConnectedElement( this->giveNumber() );
        return;
    }

    for ( el = 1; el <= this->giveChildren()->giveSize(); el++ ) {
        elem = mesh->giveElement( this->giveChildren()->at(el) );
        for ( i = 1; i <= elem->giveNodes()->giveSize(); i++ ) {
            if ( elem->giveNode(i) == node->giveNumber() ) {
                elem->buildTopLevelNodeConnectivity(node);
                break;
            }
        }
    }
}


#ifdef __PARALLEL_MODE

int
Subdivision :: RS_Node :: importConnectivity(ConnectivityTable *ct)
{
    IntArray me(1), connElems;
    int i, el, cnt = 0;
    me.at(1) = this->giveNumber();

    ct->giveNodeNeighbourList(connElems, me);
    this->connectedElements.preallocate( connElems.giveSize() );    // estimate size of the list

    for ( i = 1; i <= connElems.giveSize(); i++ ) {
        el = connElems.at(i);
        if ( this->mesh->giveElement(el)->giveParallelMode() != Element_local ) {
            continue;
        }

        // use zero chunk, the array is large enough
        this->connectedElements.insertSorted(el, 0);
        cnt++;
    }

    return ( cnt );
}



void
Subdivision :: RS_Node :: numberSharedEdges()
{
    int ie, num = this->giveNumber();
    IntArray connNodes;
    const IntArray *connElems = this->giveConnectedElements();

    for ( ie = 1; ie <= connElems->giveSize(); ie++ ) {
        this->mesh->giveElement( connElems->at(ie) )->numberSharedEdges(num, connNodes);
    }
}


void
Subdivision :: RS_Triangle :: numberSharedEdges(int iNode, IntArray &connNodes)
{
    int eIndex, jnode, jNode, eNum = this->mesh->giveNumberOfEdges();
    Subdivision :: RS_SharedEdge *_edge;

    for ( jnode = 1; jnode <= 3; jnode++ ) {
        jNode = nodes.at(jnode);
        if ( jNode == iNode ) {
            continue;
        }

        if ( this->mesh->giveNode(jNode)->giveParallelMode() != DofManager_shared ) {
            continue;
        }

        // potential shared edge

        if ( connNodes.contains(jNode) ) {
            continue;                                   // edge already processed (from iNode)
        }

        connNodes.followedBy(jNode, 10);

        eIndex = giveEdgeIndex(iNode, jNode);

        if ( shared_edges.giveSize() ) {
            if ( shared_edges.at(eIndex) ) {
                continue;                                     // edge already processed (from jNode)
            }
        }

        // create edge (even if it might not be shared)
        _edge = new Subdivision :: RS_SharedEdge(this->mesh);
        _edge->setEdgeNodes(iNode, jNode);

        this->mesh->addEdge(_edge);
        eNum++;

        // get the tentative partitions
        // they must be confirmed via mutual communication
        IntArray partitions;
        if ( _edge->giveSharedPartitions(partitions) ) {
            _edge->setPartitions(partitions);

            // I do rely on the fact that in 2D shared edge can be incident only to one element !!!
            if ( !shared_edges.giveSize() ) {
                makeSharedEdges();
            }

            shared_edges.at(eIndex) = eNum;
        }
    }
}


void
Subdivision :: RS_Tetra :: numberSharedEdges(int iNode, IntArray &connNodes)
{
    int ie, elems, eIndex, jnode, jNode, eNum = this->mesh->giveNumberOfEdges();
    Subdivision :: RS_SharedEdge *_edge;
    Subdivision :: RS_Element *elem;

    for ( jnode = 1; jnode <= 4; jnode++ ) {
        jNode = nodes.at(jnode);
        if ( jNode == iNode ) {
            continue;
        }

        if ( this->mesh->giveNode(jNode)->giveParallelMode() != DofManager_shared ) {
            continue;
        }

        // potential shared edge

        if ( connNodes.contains(jNode) ) {
            continue;                                   // edge already processed (from iNode)
        }

        connNodes.followedBy(jNode, 20);

        eIndex = giveEdgeIndex(iNode, jNode);

        if ( shared_edges.giveSize() ) {
            if ( shared_edges.at(eIndex) ) {
                continue;                                     // edge already processed (from jNode)
            }
        }

        // create edge (even if it might not be shared)
        _edge = new Subdivision :: RS_SharedEdge(this->mesh);
        _edge->setEdgeNodes(iNode, jNode);

        this->mesh->addEdge(_edge);
        ++eNum;

        // get the tentative partitions
        // they must be confirmed via mutual communication
        IntArray partitions;
        if ( _edge->giveSharedPartitions(partitions) ) {
            _edge->setPartitions(partitions);
            if ( !this->giveSharedEdges()->giveSize() ) {
                this->makeSharedEdges();
            }

            shared_edges.at(eIndex) = eNum;

            // put edge on all elements sharing it
            const IntArray *iElems, *jElems;

            iElems = mesh->giveNode(iNode)->giveConnectedElements();
            jElems = mesh->giveNode(jNode)->giveConnectedElements();

            IntArray common;
            if ( iElems->giveSize() <= jElems->giveSize() ) {
                common.preallocate( iElems->giveSize() );
            } else {
                common.preallocate( jElems->giveSize() );
            }

            // I do rely on the fact that the arrays are ordered !!!
            // I am using zero chunk because common is large enough
            elems = iElems->findCommonValuesSorted(* jElems, common, 0);
 #ifdef DEBUG_CHECK
            if ( !elems ) {
                OOFEM_ERROR("potentionally shared edge %d is not shared by common elements", eNum);
            }

 #endif
            for ( ie = 1; ie <= elems; ie++ ) {
                elem = mesh->giveElement( common.at(ie) );
                if ( elem == this ) {
                    continue;
                }

                eIndex = elem->giveEdgeIndex(iNode, jNode);
                if ( !elem->giveSharedEdges()->giveSize() ) {
                    elem->makeSharedEdges();
                }

                elem->setSharedEdge(eIndex, eNum);
            }
        }
    }
}


// CAUTION: gives the tentative partitions;
// they must be confirmed via mutual communication

int
Subdivision :: RS_SharedEdge :: giveSharedPartitions(IntArray &partitions)
{
    int count = 0;
    int myrank = this->mesh->giveSubdivision()->giveRank();
    const IntArray *iPartitions = mesh->giveNode(iNode)->givePartitions();
    const IntArray *jPartitions = mesh->giveNode(jNode)->givePartitions();

    // find common partitions;
    // this could be done more efficiently (using method findCommonValuesSorted)
    // if array iPartitions and jPartitions were sorted;
    // however they are not sorted in current version
    int ipsize = iPartitions->giveSize();
    for ( int i = 1; i <= ipsize; i++ ) {
        if ( jPartitions->contains( iPartitions->at(i) ) && ( iPartitions->at(i) != myrank ) ) {
            partitions.followedBy(iPartitions->at(i), 2);
            count++;
        }
    }

    return count;
}

#endif


double
Subdivision :: RS_Element :: giveRequiredDensity()
{
    double answer = 0.0;
    for ( int n: nodes ) {
        answer += mesh->giveNode( n )->giveRequiredDensity();
    }

    answer = answer / nodes.giveSize();
    return answer;
}


int
Subdivision :: RS_Element :: giveTopParent()
{
    RS_Element *elem = this;
    while ( elem->giveNumber() != elem->giveParent() ) {
        elem = mesh->giveElement( elem->giveParent() );
    }

    return ( elem->giveNumber() );
}


Subdivision :: RS_Triangle :: RS_Triangle(int number, RS_Mesh *mesh, int parent, IntArray &nodes) :
    Subdivision :: RS_Element(number, mesh, parent, nodes)
{
    irregular_nodes.resize(3);
    irregular_nodes.zero();
    neghbours_base_elements.resize(3);
    neghbours_base_elements.zero();

#ifdef DEBUG_CHECK
    if ( nodes.findFirstIndexOf(0) ) {
        OOFEM_ERROR("0 node of element %d", this->number);
    }

#endif
}


Subdivision :: RS_Tetra :: RS_Tetra(int number, RS_Mesh *mesh, int parent, IntArray &nodes) :
    Subdivision :: RS_Element(number, mesh, parent, nodes)
{
    irregular_nodes.resize(6);
    irregular_nodes.zero();
    side_leIndex.resize(4);
    side_leIndex.zero();
    neghbours_base_elements.resize(4);
    neghbours_base_elements.zero();

#ifdef DEBUG_CHECK
    if ( nodes.findFirstIndexOf(0) ) {
        OOFEM_ERROR("0 node of element %d", this->number);
    }

#endif
}


int
Subdivision :: RS_Triangle :: evaluateLongestEdge()
{
    if ( !this->leIndex ) {     // prevent multiple evaluation
        double elength1, elength2, elength3;

        elength1 = mesh->giveNode( this->nodes.at(1) )->giveCoordinates()->distance_square( * ( mesh->giveNode( this->nodes.at(2) )->giveCoordinates() ) );
        elength2 = mesh->giveNode( this->nodes.at(2) )->giveCoordinates()->distance_square( * ( mesh->giveNode( this->nodes.at(3) )->giveCoordinates() ) );
        elength3 = mesh->giveNode( this->nodes.at(3) )->giveCoordinates()->distance_square( * ( mesh->giveNode( this->nodes.at(1) )->giveCoordinates() ) );

        leIndex = 1;
        if ( elength2 > elength1 ) {
            leIndex = 2;
            if ( elength3 > elength2 ) {
                leIndex = 3;
            }
        } else if ( elength3 > elength1 ) {
            leIndex = 3;
        }
    }

    return ( leIndex );
}


int
Subdivision :: RS_Tetra :: evaluateLongestEdge()
{
    if ( !this->leIndex ) {    // prevent multiple evaluation
        // IMPORTANT: ambiguity must be handled !!!

        // in order to achieve that shared tetra faces are subdivided in compatible way,
        // always the same longest edge must be identified;
        // this may be problem if there are two or more edges of the "same" length
        // then the longest is dependent on the order in which edges are processed;
        // also in order to get always the same distance of two nodes, the nodes should be processes in the same order;

        // however it is not clear whether one may rely that the same distance is always obtained for the same pair of nodes
        // (in proper order) especially if processed on different processors
        // in parallel environment some additional communication may be needed !!! HUHU

        // if the second largest edge is not equal to the largest edge and if the largest and second largest edges
        // are opposite edges, no action is theoretically required;
        // but because it is not clear whether different result would be obtained if edge lengths are evaluated
        // using node ordering, this possibility (of no action) is not considered

        double elength [ 6 ], maxlength, maxlen [ 4 ];
        int i, j, side, l_nd [ 4 ], ind [ 6 ];
        // array ed_side contains face numbers shared by the edge (indexing from 0)
        int nd [ 4 ] = {
            0, 1, 2, 3
        }, ed_side [ 6 ] [ 2 ] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 3, 1 }, { 1, 2 }, { 2, 3 } };
        // array ed contains edge numbers (1 to 6) connecting the node with each node (1 to 4)
        // (indexing from 1, zero means, that node is not connected by an edge with itself)
        int ed [ 4 ] [ 4 ] = { { 0, 1, 3, 4 }, { 1, 0, 2, 5 }, { 3, 2, 0, 6 }, { 4, 5, 6, 0 } };
#ifdef __PARALLEL_MODE
        int g_nd [ 4 ];
#endif
        bool swap;

        // sort node ids
        for ( i = 0; i < 4; i++ ) {
            l_nd [ i ] = nodes.at(i + 1);
#ifdef __PARALLEL_MODE
            g_nd [ i ] = mesh->giveNode(l_nd [ i ])->giveGlobalNumber();
#endif
        }

        swap = true;
        while ( swap ) {
            swap = false;
            for ( i = 0; i < 3; i++ ) {
#ifdef __PARALLEL_MODE
                if ( g_nd [ i ] > g_nd [ i + 1 ] ) {
#else
                if ( l_nd [ i ] > l_nd [ i + 1 ] ) {
#endif
                    int tmp;
#ifdef __PARALLEL_MODE
                    tmp = g_nd [ i ];
                    g_nd [ i ] = g_nd [ i + 1 ];
                    g_nd [ i + 1 ] = tmp;
#endif
                    tmp = l_nd [ i ];
                    l_nd [ i ] = l_nd [ i + 1 ];
                    l_nd [ i + 1 ] = tmp;
                    tmp = nd [ i ];
                    nd [ i ] = nd [ i + 1 ];
                    nd [ i + 1 ] = tmp;
                    swap = true;
                }
            }
        }

        // assign edge indices from smallest node id !!!
        ind [ 0 ] = ed [ nd [ 0 ] ] [ nd [ 1 ] ];
        ind [ 1 ] = ed [ nd [ 0 ] ] [ nd [ 2 ] ];
        ind [ 2 ] = ed [ nd [ 0 ] ] [ nd [ 3 ] ];
        ind [ 3 ] = ed [ nd [ 1 ] ] [ nd [ 2 ] ];
        ind [ 4 ] = ed [ nd [ 1 ] ] [ nd [ 3 ] ];
        ind [ 5 ] = ed [ nd [ 2 ] ] [ nd [ 3 ] ];

        // evaluate edge lengths from smallest node id !!!
        elength [ 0 ] = mesh->giveNode(l_nd [ 0 ])->giveCoordinates()->distance_square( * ( mesh->giveNode(l_nd [ 1 ])->giveCoordinates() ) );
        elength [ 1 ] = mesh->giveNode(l_nd [ 0 ])->giveCoordinates()->distance_square( * ( mesh->giveNode(l_nd [ 2 ])->giveCoordinates() ) );
        elength [ 2 ] = mesh->giveNode(l_nd [ 0 ])->giveCoordinates()->distance_square( * ( mesh->giveNode(l_nd [ 3 ])->giveCoordinates() ) );
        elength [ 3 ] = mesh->giveNode(l_nd [ 1 ])->giveCoordinates()->distance_square( * ( mesh->giveNode(l_nd [ 2 ])->giveCoordinates() ) );
        elength [ 4 ] = mesh->giveNode(l_nd [ 1 ])->giveCoordinates()->distance_square( * ( mesh->giveNode(l_nd [ 3 ])->giveCoordinates() ) );
        elength [ 5 ] = mesh->giveNode(l_nd [ 2 ])->giveCoordinates()->distance_square( * ( mesh->giveNode(l_nd [ 3 ])->giveCoordinates() ) );

        // get absolutely largest edge and the largest edge on individual sides
        // process edges in the same order in which their length was evaluated !!!
        maxlength = maxlen [ 0 ] = maxlen [ 1 ] = maxlen [ 2 ] = maxlen [ 3 ] = 0.0;
        for ( i = 0; i < 6; i++ ) {
            if ( elength [ i ] > maxlength ) {
                maxlength = elength [ i ];
                leIndex = ind [ i ];
                for ( j = 0; j < 2; j++ ) {
                    side = ed_side [ ind [ i ] - 1 ] [ j ];
                    maxlen [ side ] = elength [ i ];
                    side_leIndex.at(side + 1) = ind [ i ];
                }
            } else {
                for ( j = 0; j < 2; j++ ) {
                    side = ed_side [ ind [ i ] - 1 ] [ j ];
                    if ( elength [ i ] > maxlen [ side ] ) {
                        maxlen [ side ] = elength [ i ];
                        side_leIndex.at(side + 1) = ind [ i ];
                    }
                }
            }
        }
    }

    return ( leIndex );
}


int
Subdivision :: RS_Triangle :: giveEdgeIndex(int iNode, int jNode)
{
    /* returns zero, if triangle does not have iNode or jNode */
    int in = 0, jn = 0;

    int k;
    for ( k = 1; k <= 3; k++ ) {
        if ( nodes.at(k) == iNode ) {
            in = k;
        }

        if ( nodes.at(k) == jNode ) {
            jn = k;
        }
    }

    if ( in == 0 || jn == 0 ) {
        OOFEM_ERROR( "there is no edge connecting %d and %d on element %d",
                     iNode, jNode, this->giveNumber() );
        return 0;
    }

    if ( in < jn ) {
        if ( jn == in + 1 ) {
            return ( in );
        }

        return ( 3 );
    } else {
        if ( in == jn + 1 ) {
            return ( jn );
        }

        return ( 3 );
    }
}


int
Subdivision :: RS_Tetra :: giveEdgeIndex(int iNode, int jNode)
{
    /* returns zero, if tetra does not have iNode or jNode */
    int in = 0, jn = 0;

    int k;
    for ( k = 1; k <= 4; k++ ) {
        if ( nodes.at(k) == iNode ) {
            in = k;
        }

        if ( nodes.at(k) == jNode ) {
            jn = k;
        }
    }

    if ( in == 0 || jn == 0 ) {
        OOFEM_ERROR( "there is no edge connecting %d and %d on element %d",
                     iNode, jNode, this->giveNumber() );
        return 0;
    }

    if ( in < jn ) {
        if ( jn == 4 ) {
            return ( in + 3 );
        }

        if ( jn == in + 1 ) {
            return ( in );
        }

        return ( 3 );
    } else {
        if ( in == 4 ) {
            return ( jn + 3 );
        }

        if ( in == jn + 1 ) {
            return ( jn );
        }

        return ( 3 );
    }
}


void
Subdivision :: RS_Triangle :: bisect(std :: queue< int > &subdivqueue, std :: list< int > &sharedIrregularsQueue)
{
    /* this is symbolic bisection - no new elements are added, only irregular nodes are introduced */
    int inode, jnode, iNode, jNode, iNum, eInd;
    double density;
    bool boundary = false;
    FloatArray coords;
    Subdivision :: RS_Element *elem;
#ifdef __PARALLEL_MODE
    Subdivision :: RS_SharedEdge *edge;
#endif

    if ( !irregular_nodes.at(leIndex) ) {
        RS_IrregularNode *irregular;
        // irregular on the longest edge does not exist
        // introduce new irregular node on the longest edge
        inode = leIndex;
        jnode = ( inode < 3 ) ? inode + 1 : 1;

        iNode = nodes.at(inode);
        jNode = nodes.at(jnode);
        // compute coordinates of new irregular
        coords = * ( mesh->giveNode(iNode)->giveCoordinates() );
        coords.add( * mesh->giveNode(jNode)->giveCoordinates() );
        coords.times(0.5);
        // compute required density of a new node
        density = 0.5 * ( mesh->giveNode(iNode)->giveRequiredDensity() +
                         mesh->giveNode(jNode)->giveRequiredDensity() );
        // create new irregular
        iNum = mesh->giveNumberOfNodes() + 1;
        irregular = new Subdivision :: RS_IrregularNode(iNum, mesh, 0, coords, density, boundary);
        mesh->addNode(irregular);
        // add irregular to receiver
        this->irregular_nodes.at(leIndex) = iNum;

#ifdef __OOFEG
 #ifdef DRAW_IRREGULAR_NODES
        irregular->drawGeometry();
 #endif
#endif

#ifdef __PARALLEL_MODE
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_INFO("[%d] RS_Triangle::bisecting %d nodes %d %d %d, leIndex %d, new irregular %d\n", mesh->giveSubdivision()->giveRank(), this->number, nodes.at(1), nodes.at(2), nodes.at(3), leIndex, iNum);
#endif
#endif

#ifdef QUICK_HACK
        if ( mesh->giveNode( nodes.at(1) )->isBoundary() && mesh->giveNode( nodes.at(2) )->isBoundary() && mesh->giveNode( nodes.at(3) )->isBoundary() ) {
            OOFEM_ERROR( "quick hack not applicable due to element %d", this->giveNumber() );
        }

        if ( mesh->giveNode(iNode)->isBoundary() && mesh->giveNode(jNode)->isBoundary() ) {
            boundary = true;
        }

#else
        // check whether new node is boundary
        if ( this->neghbours_base_elements.at(leIndex) ) {
            Domain *dorig = mesh->giveSubdivision()->giveDomain();
            // I rely on tha fact that nodes on intermaterial interface are marked as boundary
            // however this might not be true
            // therefore in smoothing the boundary flag is setuped again if not set boundary from here
            if ( mesh->giveNode(iNode)->isBoundary() || mesh->giveNode(jNode)->isBoundary() ) {
                if ( dorig->giveElement( this->giveTopParent() )->giveRegionNumber() != dorig->giveElement( mesh->giveElement( this->neghbours_base_elements.at(leIndex) )->giveTopParent() )->giveRegionNumber() ) {
                    boundary = true;
                }
            }
        } else {
            boundary = true;
        }

#endif

        if ( boundary ) {
            irregular->setBoundary(true);
        }

        if ( this->neghbours_base_elements.at(leIndex) ) {
            // add irregular to neighbour
            elem = mesh->giveElement( this->neghbours_base_elements.at(leIndex) );
            eInd = elem->giveEdgeIndex(iNode, jNode);
            elem->setIrregular(eInd, iNum);

            if ( !elem->giveQueueFlag() ) {
                // add neighbour to list of elements for subdivision
                subdivqueue.push( this->neghbours_base_elements.at(leIndex) );
                elem->setQueueFlag(true);
            }
        }

#ifdef __PARALLEL_MODE
        else {
            // check if there are (potentionally) shared edges
            if ( shared_edges.giveSize() ) {
                // check if the edge is (really) shared
                if ( shared_edges.at(leIndex) ) {
                    edge = mesh->giveEdge( shared_edges.at(leIndex) );

 #ifdef DEBUG_CHECK
                    if ( !edge->givePartitions()->giveSize() ) {
                        OOFEM_ERROR( "unshared edge %d of element %d is marked as shared",
                                     shared_edges.at(leIndex), this->giveNumber() );
                    }

 #endif

                    // new node is on shared interpartition boundary
                    irregular->setParallelMode(DofManager_shared);
                    // partitions are inherited from shared edge
                    irregular->setPartitions( * ( edge->givePartitions() ) );
                    irregular->setEdgeNodes(iNode, jNode);
                    // put its number into queue of shared irregulars that is later used to inform remote partitions about this fact
                    sharedIrregularsQueue.push_back(iNum);
 #ifdef __VERBOSE_PARALLEL
                    OOFEM_LOG_INFO("RS_Triangle::bisect: Shared irregular detected, number %d nodes %d %d [%d %d], elem %d\n", iNum, iNode, jNode, mesh->giveNode(iNode)->giveGlobalNumber(), mesh->giveNode(jNode)->giveGlobalNumber(), this->number);
 #endif
                }
            }
        }
#endif
    }

    this->setQueueFlag(false);
}


void
Subdivision :: RS_Tetra :: bisect(std :: queue< int > &subdivqueue, std :: list< int > &sharedIrregularsQueue)
{
    /* this is symbolic bisection - no new elements are added, only irregular nodes are introduced */
    int i, j, inode, jnode, iNode, jNode, ngb, eIndex, eInd, reg, iNum, side, cnt = 0, elems;
    // array ed_side contains face numbers NOT shared by the edge (indexing from 1)
    int ed_side [ 6 ] [ 2 ] = { { 3, 4 }, { 4, 2 }, { 2, 3 }, { 1, 3 }, { 1, 4 }, { 1, 2 } }, ed [ 4 ] = {
        0, 0, 0, 0
    }, opp_ed [ 6 ] = {
        6, 4, 5, 2, 3, 1
    };
    // array side_ed contains edge numbers bounding faces (indexing from 1)
    int side_ed [ 4 ] [ 3 ] = { { 1, 2, 3 }, { 1, 5, 4 }, { 2, 6, 5 }, { 3, 4, 6 } };
    double density;
    bool shared, boundary, iboundary, jboundary, opposite = false;
    Subdivision :: RS_Element *elem1, *elem2, *elem;
    FloatArray coords;
    Domain *dorig;
    RS_IrregularNode *irregular;
    const IntArray *iElems, *jElems;
#ifdef __PARALLEL_MODE
    Subdivision :: RS_SharedEdge *edge;
#endif

    dorig = mesh->giveSubdivision()->giveDomain();
    reg = dorig->giveElement( this->giveTopParent() )->giveRegionNumber();

    // first resolve whether there will be inserted irregular on the edge opposite to longest edge
    // this will happen if the opposite edge is longest for a side not shared by the longest edge
    // and simultaneously there is already an irregular on that side
    if ( !irregular_nodes.at(opp_ed [ leIndex - 1 ]) ) {
        for ( i = 0; i < 2; i++ ) {
            // get side not incident to leIndex edge
            side = ed_side [ leIndex - 1 ] [ i ];
            // check whether the longest edge of that side is opposite edge
            if ( side_leIndex.at(side) != opp_ed [ leIndex - 1 ] ) {
                continue;
            }

            for ( j = 0; j < 3; j++ ) {
                // check for irregulars on that side
                if ( irregular_nodes.at(side_ed [ side - 1 ] [ j ]) ) {
                    opposite = true;
                    irregular_nodes.at(opp_ed [ leIndex - 1 ]) = -1;                    // fictitious irregular
                    break;
                }
            }

            if ( opposite ) {
                break;
            }
        }
    }

    // insert irregular on absolutely longest edge first;
    // this controls the primary bisection;
    // if there is an irregular (including a fictitious one) on side not shared by longest edge
    // insert irregular on longest edge of that side;
    ed [ cnt++ ] = leIndex;
    for ( i = 0; i < 2; i++ ) {
        // get side not incident to leIndex edge
        side = ed_side [ leIndex - 1 ] [ i ];
        for ( j = 0; j < 3; j++ ) {
            // check for irregulars on that side
            if ( irregular_nodes.at(side_ed [ side - 1 ] [ j ]) ) {
                ed [ cnt++ ] = side_leIndex.at(side);
                break;
            }
        }
    }

    // remove fictitious irregular
    if ( opposite ) {
        irregular_nodes.at(opp_ed [ leIndex - 1 ]) = 0;
    }

    // check for the case when longest edge of the sides not shared by the absolutely longest edge is the same
    // (it is opposite to absolutely longest edge)
    if ( cnt == 3 ) {
        if ( ed [ 1 ] == ed [ 2 ] ) {
            cnt = 2;
#ifdef DEBUG_CHECK
            if ( ed [ 1 ] != opp_ed [ leIndex - 1 ] ) {
                OOFEM_ERROR( "unexpected situation on element %d", this->giveNumber() );
            }

#endif
        }
    }

#ifdef QUICK_HACK
    if ( mesh->giveNode( nodes.at(1) )->isBoundary() && mesh->giveNode( nodes.at(2) )->isBoundary() &&
        mesh->giveNode( nodes.at(3) )->isBoundary() && mesh->giveNode( nodes.at(4) )->isBoundary() ) {
        OOFEM_ERROR( "quick hack not applicable due to element %d", this->giveNumber() );
    }

#endif

    // insert all relevant irregulars
    for ( i = 0; i < cnt; i++ ) {
        eIndex = ed [ i ];
        if ( !irregular_nodes.at(eIndex) ) {
            // irregular on this edge does not exist;
            // introduce new irregular node on this edge on this element and on all local elements sharing that edge;
            // if the edge is local, neigbours are processed;
            // if the edge is shared, elements sharing simultaneously both end nodes are processed;
            if ( eIndex <= 3 ) {
                inode = eIndex;
                jnode = ( eIndex < 3 ) ? inode + 1 : 1;
                ngb = 1;
            } else {
                inode = eIndex - 3;
                jnode = 4;
                ngb = inode + 1;
            }

            iNode = nodes.at(inode);
            jNode = nodes.at(jnode);
            // compute coordinates of new irregular
            coords = * ( mesh->giveNode(iNode)->giveCoordinates() );
            coords.add( * mesh->giveNode(jNode)->giveCoordinates() );
            coords.times(0.5);
#ifdef HEADEDSTUD
            double dist, rad, rate;
            FloatArray *c;

            c = mesh->giveNode(iNode)->giveCoordinates();
            dist = c->at(1) * c->at(1) + c->at(3) * c->at(3);
            if ( c->at(2) > 69.9999999 ) {
                rad = 7.0;
            } else if ( c->at(2) < 64.5000001 ) {
                rad = 18.0;
            } else {
                rad = 18.0 - 11.0 / 5.5 * ( c->at(2) - 64.5 );
            }

            if ( fabs(dist - rad * rad) < 0.01 ) {            // be very tolerant (geometry is not precise)
                c = mesh->giveNode(jNode)->giveCoordinates();
                dist = c->at(1) * c->at(1) + c->at(3) * c->at(3);
                if ( c->at(2) > 69.9999999 ) {
                    rad = 7.0;
                } else if ( c->at(2) < 64.5000001 ) {
                    rad = 18.0;
                } else {
                    rad = 18.0 - 11.0 / 5.5 * ( c->at(2) - 64.5 );
                }

                if ( fabs(dist - rad * rad) < 0.01 ) {                // be very tolerant (geometry is not precise)
                    dist = coords.at(1) * coords.at(1) + coords.at(3) * coords.at(3);
                    if ( coords.at(2) > 69.9999999 ) {
                        rad = 7.0;
                    } else if ( coords.at(2) < 64.5000001 ) {
                        rad = 18.0;
                    } else {
                        rad = 18.0 - 11.0 / 5.5 * ( coords.at(2) - 64.5 );
                    }

                    rate = rad / sqrt(dist);
                    coords.at(1) *= rate;
                    coords.at(3) *= rate;
                }
            }

#endif
            // compute required density of a new node
            density = 0.5 * ( mesh->giveNode(iNode)->giveRequiredDensity() +
                             mesh->giveNode(jNode)->giveRequiredDensity() );
            // create new irregular
            iNum = mesh->giveNumberOfNodes() + 1;
            irregular = new Subdivision :: RS_IrregularNode(iNum, mesh, 0, coords, density, false);
            mesh->addNode(irregular);
            // add irregular to receiver
            this->irregular_nodes.at(eIndex) = iNum;

#ifdef __OOFEG
 #ifdef DRAW_IRREGULAR_NODES
            irregular->drawGeometry();
 #endif
#endif

#ifdef DEBUG_INFO
 #ifdef __PARALLEL_MODE
            // do not print global numbers of elements because they are not available (they are assigned at once after bisection);
            // do not print global numbers of irregulars as these may not be available yet
            OOFEM_LOG_INFO( "[%d] Irregular %d added on %d [%d] (edge %d, nodes %d %d [%d %d], nds %d %d %d %d [%d %d %d %d], ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                            mesh->giveSubdivision()->giveRank(), iNum, this->number, this->giveGlobalNumber(), eIndex, iNode, jNode,
                           mesh->giveNode(iNode)->giveGlobalNumber(), mesh->giveNode(jNode)->giveGlobalNumber(),
                           nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4),
                           mesh->giveNode( nodes.at(1) )->giveGlobalNumber(), mesh->giveNode( nodes.at(2) )->giveGlobalNumber(),
                           mesh->giveNode( nodes.at(3) )->giveGlobalNumber(), mesh->giveNode( nodes.at(4) )->giveGlobalNumber(),
                           neghbours_base_elements.at(1), neghbours_base_elements.at(2),
                           neghbours_base_elements.at(3), neghbours_base_elements.at(4),
                           irregular_nodes.at(1), irregular_nodes.at(2), irregular_nodes.at(3),
                           irregular_nodes.at(4), irregular_nodes.at(5), irregular_nodes.at(6) );
 #else
            OOFEM_LOG_INFO( "Irregular %d added on %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                           iNum, this->number, eIndex, iNode, jNode,
                           nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4),
                           neghbours_base_elements.at(1), neghbours_base_elements.at(2),
                           neghbours_base_elements.at(3), neghbours_base_elements.at(4),
                           irregular_nodes.at(1), irregular_nodes.at(2), irregular_nodes.at(3),
                           irregular_nodes.at(4), irregular_nodes.at(5), irregular_nodes.at(6) );
 #endif
#endif

            shared = boundary = false;

#ifdef __PARALLEL_MODE
            // check if there are (potentionally) shared edges
            if ( shared_edges.giveSize() ) {
                // check if the edge is (really) shared
                if ( shared_edges.at(eIndex) ) {
                    edge = mesh->giveEdge( shared_edges.at(eIndex) );

 #ifdef DEBUG_CHECK
                    if ( !edge->givePartitions()->giveSize() ) {
                        OOFEM_ERROR( "unshared edge %d of element %d is marked as shared",
                                     shared_edges.at(eIndex), this->giveNumber() );
                    }

 #endif

                    shared = boundary = true;
                    // new node is on shared interpartition boundary
                    irregular->setParallelMode(DofManager_shared);
                    irregular->setPartitions( * ( edge->givePartitions() ) );
                    irregular->setEdgeNodes(iNode, jNode);
                    // put its number into queue of shared irregulars that is later used to inform remote partitions about this fact
                    sharedIrregularsQueue.push_back(iNum);
 #ifdef __VERBOSE_PARALLEL
                    OOFEM_LOG_INFO("RS_Tetra::bisect: Shared irregular detected, number %d nodes %d %d [%d %d], elem %d\n", iNum, iNode, jNode, mesh->giveNode(iNode)->giveGlobalNumber(), mesh->giveNode(jNode)->giveGlobalNumber(), this->number);
 #endif

                    iElems = mesh->giveNode(iNode)->giveConnectedElements();
                    jElems = mesh->giveNode(jNode)->giveConnectedElements();

                    IntArray common;
                    if ( iElems->giveSize() <= jElems->giveSize() ) {
                        common.preallocate( iElems->giveSize() );
                    } else {
                        common.preallocate( jElems->giveSize() );
                    }

                    // I do rely on the fact that the arrays are ordered !!!
                    // I am using zero chunk because common is large enough
                    elems = iElems->findCommonValuesSorted(* jElems, common, 0);
 #ifdef DEBUG_CHECK
                    if ( !elems ) {
                        OOFEM_ERROR( "shared edge %d is not shared by common elements",
                                     shared_edges.at(eIndex) );
                    }

 #endif
                    // after subdivision there will be twice as much of connected local elements
                    irregular->preallocateConnectedElements(elems * 2);
                    // put the new node on appropriate edge of all local elements (except "this") sharing both nodes
                    for ( j = 1; j <= elems; j++ ) {
                        elem = mesh->giveElement( common.at(j) );
                        if ( elem == this ) {
                            continue;
                        }

 #ifdef DEBUG_CHECK
                        if ( !elem->giveSharedEdges()->giveSize() ) {
                            OOFEM_ERROR( "element %d incident to shared edge %d does not have shared edges",
                                         common.at(j), shared_edges.at(eIndex) );
                        }

 #endif

                        eInd = elem->giveEdgeIndex(iNode, jNode);
                        elem->setIrregular(eInd, iNum);

                        if ( !elem->giveQueueFlag() ) {
                            // add elem to list of elements for subdivision
                            subdivqueue.push( elem->giveNumber() );
                            elem->setQueueFlag(true);
                        }

 #ifdef DEBUG_INFO
  #ifdef __PARALLEL_MODE
                        // do not print global numbers of elements because they are not available (they are assigned at once after bisection);
                        // do not print global numbers of irregulars as these may not be available yet
                        OOFEM_LOG_INFO( "[%d] Irregular %d added on %d [%d] (edge %d, nodes %d %d [%d %d], nds %d %d %d %d [%d %d %d %d], ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                                       mesh->giveSubdivision()->giveRank(), iNum,
                                        elem->giveNumber(), this->giveGlobalNumber(), eInd, iNode, jNode,
                                       mesh->giveNode(iNode)->giveGlobalNumber(), mesh->giveNode(jNode)->giveGlobalNumber(),
                                       elem->giveNode(1), elem->giveNode(2), elem->giveNode(3), elem->giveNode(4),
                                       mesh->giveNode( elem->giveNode(1) )->giveGlobalNumber(),
                                       mesh->giveNode( elem->giveNode(2) )->giveGlobalNumber(),
                                       mesh->giveNode( elem->giveNode(3) )->giveGlobalNumber(),
                                       mesh->giveNode( elem->giveNode(4) )->giveGlobalNumber(),
                                       elem->giveNeighbor(1), elem->giveNeighbor(2), elem->giveNeighbor(3), elem->giveNeighbor(4),
                                       elem->giveIrregular(1), elem->giveIrregular(2), elem->giveIrregular(3),
                                       elem->giveIrregular(4), elem->giveIrregular(5), elem->giveIrregular(6) );
  #else
                        OOFEM_LOG_INFO( "Irregular %d added on %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                                       iNum, elem->giveNumber(), eInd, iNode, jNode,
                                       elem->giveNode(1), elem->giveNode(2), elem->giveNode(3), elem->giveNode(4),
                                       elem->giveNeighbor(1), elem->giveNeighbor(2), elem->giveNeighbor(3), elem->giveNeighbor(4),
                                       elem->giveIrregular(1), elem->giveIrregular(2), elem->giveIrregular(3),
                                       elem->giveIrregular(4), elem->giveIrregular(5), elem->giveIrregular(6) );
  #endif
 #endif
                    }
                }
            }

#endif

            if ( !shared ) {
                iboundary = mesh->giveNode(iNode)->isBoundary();
                jboundary = mesh->giveNode(jNode)->isBoundary();
#ifdef QUICK_HACK
                if ( iboundary == true && jboundary == true ) {
                    boundary = true;
                }

#endif

                // traverse neighbours
                elem1 = this;
                elem2 = NULL;
                while ( elem1->giveNeighbor(ngb) ) {
                    elem2 = mesh->giveElement( elem1->giveNeighbor(ngb) );
                    if ( elem2 == this ) {
                        break;
                    }

                    eInd = elem2->giveEdgeIndex(iNode, jNode);
                    elem2->setIrregular(eInd, iNum);

                    if ( !elem2->giveQueueFlag() ) {
                        // add neighbour to list of elements for subdivision
                        subdivqueue.push( elem2->giveNumber() );
                        elem2->setQueueFlag(true);
                    }

#ifndef QUICK_HACK
                    if ( !boundary ) {
                        // I rely on the fact that nodes on intermaterial interface are marked as boundary
                        // however this might not be true
                        // therefore in smoothing the boundary flag is setuped again if not set boundary from here
                        if ( iboundary == true || jboundary == true ) {
                            if ( dorig->giveElement( elem2->giveTopParent() )->giveRegionNumber() != reg ) {
                                boundary = true;
                            }
                        }
                    }

#endif

#ifdef DEBUG_INFO
 #ifdef __PARALLEL_MODE
                    // do not print global numbers of elements because they are not available (they are assigned at once after bisection);
                    // do not print global numbers of irregulars as these may not be available yet
                    OOFEM_LOG_INFO( "[%d] Irregular %d added on %d [%d] (edge %d, nodes %d %d [%d %d], nds %d %d %d %d [%d %d %d %d], ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                                   mesh->giveSubdivision()->giveRank(), iNum,
                                    elem2->giveNumber(), this->giveGlobalNumber(), eInd, iNode, jNode,
                                   mesh->giveNode(iNode)->giveGlobalNumber(), mesh->giveNode(jNode)->giveGlobalNumber(),
                                   elem2->giveNode(1), elem2->giveNode(2), elem2->giveNode(3), elem2->giveNode(4),
                                   mesh->giveNode( elem2->giveNode(1) )->giveGlobalNumber(),
                                   mesh->giveNode( elem2->giveNode(2) )->giveGlobalNumber(),
                                   mesh->giveNode( elem2->giveNode(3) )->giveGlobalNumber(),
                                   mesh->giveNode( elem2->giveNode(4) )->giveGlobalNumber(),
                                   elem2->giveNeighbor(1), elem2->giveNeighbor(2), elem2->giveNeighbor(3), elem2->giveNeighbor(4),
                                   elem2->giveIrregular(1), elem2->giveIrregular(2), elem2->giveIrregular(3),
                                   elem2->giveIrregular(4), elem2->giveIrregular(5), elem2->giveIrregular(6) );
 #else
                    OOFEM_LOG_INFO( "Irregular %d added on %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                                   iNum, elem2->giveNumber(), eInd, iNode, jNode,
                                   elem2->giveNode(1), elem2->giveNode(2), elem2->giveNode(3), elem2->giveNode(4),
                                   elem2->giveNeighbor(1), elem2->giveNeighbor(2), elem2->giveNeighbor(3), elem2->giveNeighbor(4),
                                   elem2->giveIrregular(1), elem2->giveIrregular(2), elem2->giveIrregular(3),
                                   elem2->giveIrregular(4), elem2->giveIrregular(5), elem2->giveIrregular(6) );
 #endif
#endif

                    if ( eInd <= 3 ) {
                        if ( elem2->giveNeighbor(1) == elem1->giveNumber() ) {
                            ngb = eInd + 1;
                        } else {
                            ngb = 1;
                        }
                    } else {
                        if ( elem2->giveNeighbor(eInd - 2) == elem1->giveNumber() ) {
                            ngb = ( eInd > 4 ) ? eInd - 3 : 4;
                        } else {
                            ngb = eInd - 2;
                        }
                    }

                    elem1 = elem2;
                }

                if ( elem2 != this ) {
#ifdef DEBUG_CHECK
 #ifdef THREEPBT_3D
                    if ( coords.at(1) > 0.000001 && coords.at(1) < 1999.99999 &&
                         coords.at(2) > 0.000001 && coords.at(2) < 249.99999 &&
                         coords.at(3) > 0.000001 && coords.at(3) < 499.99999 ) {
                        if ( 987.5 - coords.at(1) > 0.000001 || coords.at(1) - 1012.5 > 0.000001 || 300.0 - coords.at(3) > 0.000001 ) {
                            OOFEM_ERROR("Irregular %d [%d %d] not on boundary", iNum, iNode, jNode);
                        }
                    }

 #endif
#endif
#ifndef QUICK_HACK
                    boundary = true;
#endif
                    // edge is on outer boundary

                    // I do rely on the fact that if the list of connected elements is not availale
                    // then the node is on top level
                    iElems = mesh->giveNode(iNode)->giveConnectedElements();
                    if ( !iElems->giveSize() ) {
                        mesh->giveNode(iNode)->buildTopLevelNodeConnectivity( mesh->giveSubdivision()->giveDomain()->giveConnectivityTable() );
                    }

                    jElems = mesh->giveNode(jNode)->giveConnectedElements();
                    if ( !jElems->giveSize() ) {
                        mesh->giveNode(jNode)->buildTopLevelNodeConnectivity( mesh->giveSubdivision()->giveDomain()->giveConnectivityTable() );
                    }

                    IntArray common;
                    if ( iElems->giveSize() <= jElems->giveSize() ) {
                        common.preallocate( iElems->giveSize() );
                    } else {
                        common.preallocate( jElems->giveSize() );
                    }

                    // I do rely on the fact that the arrays are ordered !!!
                    // I am using zero chunk because common is large enough
                    elems = iElems->findCommonValuesSorted(* jElems, common, 0);
#ifdef DEBUG_CHECK
                    if ( !elems ) {
                        OOFEM_ERROR("local outer edge %d %d is not shared by common elements",
                                     iNode, jNode);
                    }

#endif
                    // after subdivision there will be twice as much of connected local elements
                    irregular->preallocateConnectedElements(elems * 2);
                    irregular->setNumber(-iNum);                                                 // mark local unshared irregular for connectivity setup
                    // put the new node on appropriate edge of all local elements sharing both nodes
                    // (if not yet done during neighbour traversal)
                    for ( j = 1; j <= elems; j++ ) {
                        elem = mesh->giveElement( common.at(j) );
                        if ( elem == this ) {
                            continue;
                        }

                        eInd = elem->giveEdgeIndex(iNode, jNode);
                        if ( !elem->giveIrregular(eInd) ) {
                            elem->setIrregular(eInd, iNum);

#ifdef DEBUG_INFO
 #ifdef __PARALLEL_MODE
                            // do not print global numbers of elements because they are not available (they are assigned at once after bisection);
                            // do not print global numbers of irregulars as these may not be available yet
                            OOFEM_LOG_INFO( "[%d] Irregular %d added on %d [%d] (edge %d, nodes %d %d [%d %d], nds %d %d %d %d [%d %d %d %d], ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                                           mesh->giveSubdivision()->giveRank(), iNum,
                                            elem->giveNumber(), this->giveGlobalNumber(), eInd, iNode, jNode,
                                           mesh->giveNode(iNode)->giveGlobalNumber(), mesh->giveNode(jNode)->giveGlobalNumber(),
                                           elem->giveNode(1), elem->giveNode(2), elem->giveNode(3), elem->giveNode(4),
                                           mesh->giveNode( elem->giveNode(1) )->giveGlobalNumber(),
                                           mesh->giveNode( elem->giveNode(2) )->giveGlobalNumber(),
                                           mesh->giveNode( elem->giveNode(3) )->giveGlobalNumber(),
                                           mesh->giveNode( elem->giveNode(4) )->giveGlobalNumber(),
                                           elem->giveNeighbor(1), elem->giveNeighbor(2), elem->giveNeighbor(3), elem->giveNeighbor(4),
                                           elem->giveIrregular(1), elem->giveIrregular(2), elem->giveIrregular(3),
                                           elem->giveIrregular(4), elem->giveIrregular(5), elem->giveIrregular(6) );
 #else
                            OOFEM_LOG_INFO( "Irregular %d added on %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                                           iNum, elem->giveNumber(), eInd, iNode, jNode,
                                           elem->giveNode(1), elem->giveNode(2), elem->giveNode(3), elem->giveNode(4),
                                           elem->giveNeighbor(1), elem->giveNeighbor(2), elem->giveNeighbor(3), elem->giveNeighbor(4),
                                           elem->giveIrregular(1), elem->giveIrregular(2), elem->giveIrregular(3),
                                           elem->giveIrregular(4), elem->giveIrregular(5), elem->giveIrregular(6) );
 #endif
#endif

                            if ( !elem->giveQueueFlag() ) {
                                // add elem to list of elements for subdivision
                                subdivqueue.push( elem->giveNumber() );
                                elem->setQueueFlag(true);
                            }
                        }
                    }
                }
            }

            if ( boundary ) {
                // OOFEM_LOG_INFO("Irregular %d set boundary\n", abs(irregular->giveNumber()));
                irregular->setBoundary(true);
            }
        }
    }

    this->setQueueFlag(false);
}


void
Subdivision :: RS_Triangle :: generate(std :: list< int > &sharedEdgesQueue)
{
#ifdef DEBUG_CHECK
    if ( this->queue_flag ) {
        OOFEM_ERROR("unexpected queue flag of %d ", this->number);
    }

#endif

    /* generates the children elements of already bisected element and adds them into mesh;
     * also children array of receiver is updated to contain children numbers */
    if ( irregular_nodes.containsOnlyZeroes() ) {
        // no subdivision of receiver required
        children.clear();
    } else {
        int i, ind, nIrregulars = 0;
        int childNum, iedge, jedge, kedge, inode, jnode, knode;
        IntArray _nodes(3);
        Subdivision :: RS_Triangle *child;
        Subdivision :: RS_Element *ngb;
#ifdef __PARALLEL_MODE
        bool ishared = false, jshared = false, kshared = false;
#endif

        for ( i = 1; i <= irregular_nodes.giveSize(); i++ ) {
            if ( irregular_nodes.at(i) ) {
                nIrregulars++;
            }
        }

        this->children.resize(nIrregulars + 1);

        // leIndex determines primary edge to be subdivided
        /*
         *                         inode
         *                           o
         *
         *
         *             kedge x                x jedge
         *
         *
         *          o                x                o
         *        jnode            iedge               knode
         *                        =leIndex
         */
        iedge = leIndex;
        jedge = ( iedge < 3 ) ? iedge + 1 : 1;
        kedge = ( jedge < 3 ) ? jedge + 1 : 1;

        inode = ( iedge > 1 ) ? iedge - 1 : 3;
        jnode = ( inode < 3 ) ? inode + 1 : 1;
        knode = ( jnode < 3 ) ? jnode + 1 : 1;

#ifdef __PARALLEL_MODE
        if ( shared_edges.giveSize() ) {
            if ( shared_edges.at(iedge) ) {
                ishared = true;
            }

            if ( shared_edges.at(jedge) ) {
                jshared = true;
            }

            if ( shared_edges.at(kedge) ) {
                kshared = true;
            }
        }

#endif

        if ( irregular_nodes.at(iedge) && irregular_nodes.at(jedge) && irregular_nodes.at(kedge) ) {
            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = irregular_nodes.at(kedge);
            _nodes.at(3) = irregular_nodes.at(iedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(1) = childNum;
            // set neigbour info
            child->setNeighbor( 1, -this->giveNeighbor(kedge) );
            child->setNeighbor(2, childNum + 2);
            child->setNeighbor(3, childNum + 1);
#ifdef __PARALLEL_MODE
            // set shared info
            if ( kshared ) {
                child->makeSharedEdges();
                child->setSharedEdge( 1, shared_edges.at(kedge) );
            }

#endif

            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = irregular_nodes.at(iedge);
            _nodes.at(3) = irregular_nodes.at(jedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(2) = childNum;
            // set neigbour info
            child->setNeighbor(1, childNum - 1);
            child->setNeighbor(2, childNum + 2);
            child->setNeighbor( 3, -this->giveNeighbor(jedge) );
#ifdef __PARALLEL_MODE
            // set shared info
            if ( jshared ) {
                child->makeSharedEdges();
                child->setSharedEdge( 3, shared_edges.at(jedge) );
            }

#endif

            _nodes.at(1) = nodes.at(jnode);
            _nodes.at(2) = irregular_nodes.at(iedge);
            _nodes.at(3) = irregular_nodes.at(kedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(3) = childNum;
            // set neigbour info
            child->setNeighbor( 1, -this->giveNeighbor(iedge) );
            child->setNeighbor(2, childNum - 2);
            child->setNeighbor( 3, -this->giveNeighbor(kedge) );
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || kshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 1, shared_edges.at(iedge) );
                }

                if ( kshared ) {
                    child->setSharedEdge( 3, shared_edges.at(kedge) );
                }
            }

#endif

            _nodes.at(1) = nodes.at(knode);
            _nodes.at(2) = irregular_nodes.at(jedge);
            _nodes.at(3) = irregular_nodes.at(iedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(4) = childNum;
            // set neigbour info
            child->setNeighbor( 1, -this->giveNeighbor(jedge) );
            child->setNeighbor(2, childNum - 2);
            child->setNeighbor( 3, -this->giveNeighbor(iedge) );

#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || jshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 3, shared_edges.at(iedge) );
                }

                if ( jshared ) {
                    child->setSharedEdge( 1, shared_edges.at(jedge) );
                }
            }

#endif

#ifdef __PARALLEL_MODE
            // set connected elements
            if ( mesh->giveNode( nodes.at(inode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(inode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 2);
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 3);
            }

            if ( mesh->giveNode( nodes.at(jnode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(jnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(knode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(knode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum);
            }

            if ( ishared ) {
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 1);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 2);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 3);
            }

            if ( jshared ) {
                mesh->giveNode( irregular_nodes.at(jedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(jedge) )->insertConnectedElement(childNum - 2);
            }

            if ( kshared ) {
                mesh->giveNode( irregular_nodes.at(kedge) )->insertConnectedElement(childNum - 1);
                mesh->giveNode( irregular_nodes.at(kedge) )->insertConnectedElement(childNum - 3);
            }

#endif
        } else if ( irregular_nodes.at(iedge) && irregular_nodes.at(jedge) ) {
            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = nodes.at(jnode);
            _nodes.at(3) = irregular_nodes.at(iedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(1) = childNum;
            // set neigbour info
            child->setNeighbor( 1, -this->giveNeighbor(kedge) );
            child->setNeighbor( 2, -this->giveNeighbor(iedge) );
            child->setNeighbor(3, childNum + 1);
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || kshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 2, shared_edges.at(iedge) );
                }

                if ( kshared ) {
                    child->setSharedEdge( 1, shared_edges.at(kedge) );
                }
            }

#endif

            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = irregular_nodes.at(iedge);
            _nodes.at(3) = irregular_nodes.at(jedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(2) = childNum;
            // set neigbour info
            child->setNeighbor(1, childNum - 1);
            child->setNeighbor(2, childNum + 1);
            child->setNeighbor( 3, -this->giveNeighbor(jedge) );
#ifdef __PARALLEL_MODE
            // set shared info
            if ( jshared ) {
                child->makeSharedEdges();
                child->setSharedEdge( 3, shared_edges.at(jedge) );
            }

#endif

            _nodes.at(1) = nodes.at(knode);
            _nodes.at(2) = irregular_nodes.at(jedge);
            _nodes.at(3) = irregular_nodes.at(iedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(3) = childNum;
            // set neigbour info
            child->setNeighbor( 1, -this->giveNeighbor(jedge) );
            child->setNeighbor(2, childNum - 1);
            child->setNeighbor( 3, -this->giveNeighbor(iedge) );
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || jshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 3, shared_edges.at(iedge) );
                }

                if ( jshared ) {
                    child->setSharedEdge( 1, shared_edges.at(jedge) );
                }
            }

#endif

#ifdef __PARALLEL_MODE
            // set connected elements
            if ( mesh->giveNode( nodes.at(inode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(inode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 1);
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 2);
            }

            if ( mesh->giveNode( nodes.at(jnode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(jnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum - 2);
            }

            if ( mesh->giveNode( nodes.at(knode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(knode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum);
            }

            if ( ishared ) {
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 1);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 2);
            }

            if ( jshared ) {
                mesh->giveNode( irregular_nodes.at(jedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(jedge) )->insertConnectedElement(childNum - 1);
            }

#endif
        } else if ( irregular_nodes.at(iedge) && irregular_nodes.at(kedge) ) {
            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = irregular_nodes.at(kedge);
            _nodes.at(3) = irregular_nodes.at(iedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(1) = childNum;
            // set neigbour info
            child->setNeighbor( 1, -this->giveNeighbor(kedge) );
            child->setNeighbor(2, childNum + 1);
            child->setNeighbor(3, childNum + 2);
#ifdef __PARALLEL_MODE
            // set shared info
            if ( kshared ) {
                child->makeSharedEdges();
                child->setSharedEdge( 1, shared_edges.at(kedge) );
            }

#endif

            _nodes.at(1) = nodes.at(jnode);
            _nodes.at(2) = irregular_nodes.at(iedge);
            _nodes.at(3) = irregular_nodes.at(kedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(2) = childNum;
            // set neigbour info
            child->setNeighbor( 1, -this->giveNeighbor(iedge) );
            child->setNeighbor(2, childNum - 1);
            child->setNeighbor( 3, -this->giveNeighbor(kedge) );
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || kshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 1, shared_edges.at(iedge) );
                }

                if ( kshared ) {
                    child->setSharedEdge( 3, shared_edges.at(kedge) );
                }
            }

#endif

            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = irregular_nodes.at(iedge);
            _nodes.at(3) = nodes.at(knode);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(3) = childNum;
            // set neigbour info
            child->setNeighbor(1, childNum - 2);
            child->setNeighbor( 2, -this->giveNeighbor(iedge) );
            child->setNeighbor( 3, -this->giveNeighbor(jedge) );
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || jshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 2, shared_edges.at(iedge) );
                }

                if ( jshared ) {
                    child->setSharedEdge( 3, shared_edges.at(jedge) );
                }
            }

#endif

#ifdef __PARALLEL_MODE
            // set connected elements
            if ( mesh->giveNode( nodes.at(inode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(inode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 2);
            }

            if ( mesh->giveNode( nodes.at(jnode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(jnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(knode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(knode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum);
            }

            if ( ishared ) {
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 1);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 2);
            }

            if ( kshared ) {
                mesh->giveNode( irregular_nodes.at(kedge) )->insertConnectedElement(childNum - 1);
                mesh->giveNode( irregular_nodes.at(kedge) )->insertConnectedElement(childNum - 2);
            }

#endif
        } else if ( irregular_nodes.at(iedge) ) {
            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = nodes.at(jnode);
            _nodes.at(3) = irregular_nodes.at(iedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(1) = childNum;
            // set neigbour info
            child->setNeighbor( 1, -this->giveNeighbor(kedge) );
            child->setNeighbor( 2, -this->giveNeighbor(iedge) );
            child->setNeighbor(3, childNum + 1);
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || kshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 2, shared_edges.at(iedge) );
                }

                if ( kshared ) {
                    child->setSharedEdge( 1, shared_edges.at(kedge) );
                }
            }

#endif

            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = irregular_nodes.at(iedge);
            _nodes.at(3) = nodes.at(knode);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Triangle(childNum, mesh, this->number, _nodes); // number, parent, coords
            mesh->addElement(child);
            children.at(2) = childNum;
            // set neigbour info
            child->setNeighbor(1, childNum - 1);
            child->setNeighbor( 2, -this->giveNeighbor(iedge) );
            child->setNeighbor( 3, -this->giveNeighbor(jedge) );
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || jshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 2, shared_edges.at(iedge) );
                }

                if ( jshared ) {
                    child->setSharedEdge( 3, shared_edges.at(jedge) );
                }
            }

#endif

#ifdef __PARALLEL_MODE
            // set connected elements
            if ( mesh->giveNode( nodes.at(inode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(inode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(jnode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(jnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(knode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(knode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum);
            }

            if ( ishared ) {
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 1);
            }

#endif
        } else {
            OOFEM_ERROR("element %d internal data inconsistency", this->number);
        }

        // if there is neighbor of "this" not designated for bisection
        // its neihgbor corresponding to "this" must be made negative to enforce update_neighbors
        for ( i = 1; i <= 3; i++ ) {
            if ( this->neghbours_base_elements.at(i) ) {
#ifdef DEBUG_CHECK
                if ( this->neghbours_base_elements.at(i) < 0 ) {
                    OOFEM_ERROR("negative neighbor of %d not expected", this->number);
                }

#endif

                ngb = mesh->giveElement( this->neghbours_base_elements.at(i) );
                if ( !ngb->hasIrregulars() ) {
                    ind = ngb->giveNeighbors()->findFirstIndexOf(this->number);
                    if ( ind ) {
                        ngb->setNeighbor(ind, -this->number);
                    }
                }
            }
        }
    }
}


void
Subdivision :: RS_Tetra :: generate(std :: list< int > &sharedEdgesQueue)
{
#ifdef DEBUG_CHECK
    if ( this->queue_flag ) {
        OOFEM_ERROR("unexpected queue flag of %d ", this->number);
    }

#endif

    /* generates the children elements of already bisected element and adds them into mesh;
     *             this is done recursively;
     * also children array of receiver is updated to contain children numbers */

    if ( irregular_nodes.containsOnlyZeroes() ) {
        // no subdivision of receiver required
        children.clear();
    } else {
        int irregulars1 = 0, irregulars2 = 0;
        int childNum, iedge, jedge, kedge, iiedge, jjedge, kkedge, inode, jnode, knode, nnode, iside, jside, kside, nside;
        int i, ind, leIndex1 = 0, leIndex2 = 0;
        IntArray _nodes(4);
        Subdivision :: RS_Tetra *child;
        Subdivision :: RS_Element *ngb;
#ifdef __PARALLEL_MODE
        int eNum = this->mesh->giveNumberOfEdges();
        Subdivision :: RS_SharedEdge *_edge;
        bool ishared = false, jshared = false, kshared = false;
        bool iishared = false, jjshared = false, kkshared = false;
#endif

        // leIndex determines primary edge to be subdivided
        // it is identified either with iedge or iiedge
        /*
         *
         *                         inode
         *                           o
         *
         *                        iiedge=leIndex(if>3)
         *                           x
         *                                                     inode, jnode, knode - base nodes
         *             kedge x     nnode      x jedge          nnode - top node
         *                           o                         iiedge, jjedge, kkedge - edges connecting base nodes with top node
         *
         *            jjedge x                x kkedge
         *
         *          o                x                o
         *        jnode            iedge               knode
         *                        =leIndex(if<=3)
         *
         *           children keep all nodes in the same order except for that one which is replaced by irregular on leIndex edge !!!
         */
        if ( leIndex <= 3 ) {
            iedge = leIndex;
        } else {
            iedge = ( leIndex == 6 ) ? 1 : leIndex - 2;
        }

        jedge = ( iedge < 3 ) ? iedge + 1 : 1;
        kedge = ( jedge < 3 ) ? jedge + 1 : 1;

        inode = ( iedge > 1 ) ? iedge - 1 : 3;
        jnode = ( inode < 3 ) ? inode + 1 : 1;
        knode = ( jnode < 3 ) ? jnode + 1 : 1;
        nnode = 4;

        iiedge = inode + 3;
        jjedge = jnode + 3;
        kkedge = knode + 3;

        iside = iedge + 1;
        jside = jedge + 1;
        kside = kedge + 1;
        nside = 1;

        this->children.resize(2);

#ifdef __PARALLEL_MODE
        if ( shared_edges.giveSize() ) {
            if ( shared_edges.at(iedge) ) {
                ishared = true;
            }

            if ( shared_edges.at(jedge) ) {
                jshared = true;
            }

            if ( shared_edges.at(kedge) ) {
                kshared = true;
            }

            if ( shared_edges.at(iiedge) ) {
                iishared = true;
            }

            if ( shared_edges.at(jjedge) ) {
                jjshared = true;
            }

            if ( shared_edges.at(kkedge) ) {
                kkshared = true;
            }
        }

#endif

        if ( leIndex <= 3 ) {
            // count number of irregulars on each child;
            // if there is one irregular only, the set leIndex is valid
            // otherwise it is corrected later
            if ( irregular_nodes.at(iiedge) ) {
                irregulars1++;
                irregulars2++;
                leIndex1 = leIndex2 = 4;
            }

            if ( irregular_nodes.at(kedge) ) {
                irregulars1++;
                leIndex1 = 1;
            }

            if ( irregular_nodes.at(jjedge) ) {
                irregulars1++;
                leIndex1 = 5;
            }

            if ( irregular_nodes.at(jedge) ) {
                irregulars2++;
                leIndex2 = 3;
            }

            if ( irregular_nodes.at(kkedge) ) {
                irregulars2++;
                leIndex2 = 6;
            }

            if ( irregulars1 > 1 ) {
                // use parent information about the longest edge on kside
                leIndex1 = this->side_leIndex.at(kside);

                if ( leIndex1 == jjedge ) {
                    leIndex1 = 5;
                } else if ( leIndex1 == kedge ) {
                    leIndex1 = 1;
                } else {
#ifdef DEBUG_CHECK
                    if ( leIndex1 != iiedge ) {
                        OOFEM_ERROR("side longest edge inconsistency on %d", this->number);
                    }

#endif
                    leIndex1 = 4;
                }
            }

            if ( irregulars2 > 1 ) {
                // use parent information about the longest edge on jside
                leIndex2 = this->side_leIndex.at(jside);

                if ( leIndex2 == kkedge ) {
                    leIndex2 = 6;
                } else if ( leIndex2 == jedge ) {
                    leIndex2 = 3;
                } else {
#ifdef DEBUG_CHECK
                    if ( leIndex2 != iiedge ) {
                        OOFEM_ERROR("side longest edge inconsistency on %d", this->number);
                    }

#endif
                    leIndex2 = 4;
                }
            }

#ifdef __PARALLEL_MODE
            int i_shared_id = 0, n_shared_id = 0;

            // check whether new edges are potentially shared
            if ( ishared ) {
                if ( mesh->giveNode( nodes.at(inode) )->giveParallelMode() == DofManager_shared ) {
                    if ( !this->giveNeighbor(nside) ) {
                        _edge = new Subdivision :: RS_SharedEdge(this->mesh);
                        _edge->setEdgeNodes( irregular_nodes.at(iedge), nodes.at(inode) );

                        this->mesh->addEdge(_edge);
                        eNum++;

                        // get the tentative partitions
                        // they must be confirmed via mutual communication
                        IntArray partitions;
                        if ( _edge->giveSharedPartitions(partitions) ) {
                            i_shared_id = eNum;
                            _edge->setPartitions(partitions);
                            // put edge number into queue of shared edges to resolve remote partitions
                            sharedEdgesQueue.push_back(eNum);
                        }
                    }
                }

                if ( mesh->giveNode( nodes.at(nnode) )->giveParallelMode() == DofManager_shared ) {
                    if ( !this->giveNeighbor(iside) ) {
                        _edge = new Subdivision :: RS_SharedEdge(this->mesh);
                        _edge->setEdgeNodes( irregular_nodes.at(iedge), nodes.at(nnode) );

                        this->mesh->addEdge(_edge);
                        eNum++;

                        // get the tentative partitions
                        // they must be confirmed via mutual communication
                        IntArray partitions;
                        if ( _edge->giveSharedPartitions(partitions) ) {
                            n_shared_id = eNum;
                            _edge->setPartitions(partitions);
                            // put edge number into queue of shared edges to resolve remote partitions
                            sharedEdgesQueue.push_back(eNum);
                        }
                    }
                }
            }

#endif

            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = nodes.at(jnode);
            _nodes.at(3) = irregular_nodes.at(iedge);
            _nodes.at(4) = nodes.at(nnode);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Tetra(childNum, mesh, this->number, _nodes);            // number, parent, coords
            mesh->addElement(child);
            children.at(1) = childNum;

            // set neigbour info
            if ( irregulars1 ) {
                child->setNeighbor( 1, this->giveNeighbor(nside) );
                child->setNeighbor( 2, this->giveNeighbor(kside) );
                child->setNeighbor( 3, this->giveNeighbor(iside) );

                child->setIrregular( 1, irregular_nodes.at(kedge) );
                child->setIrregular( 5, irregular_nodes.at(jjedge) );
                child->setIrregular( 4, irregular_nodes.at(iiedge) );
            } else {
                child->setNeighbor( 1, -this->giveNeighbor(nside) );
                child->setNeighbor( 2, -this->giveNeighbor(kside) );
                child->setNeighbor( 3, -this->giveNeighbor(iside) );
            }

            // neihgbor4 of child1 is changed to negative during subdivision (if any) of child2
            child->setNeighbor(4, childNum + 1);

#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || kshared || iishared || jjshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 2, shared_edges.at(iedge) );
                }

                if ( kshared ) {
                    child->setSharedEdge( 1, shared_edges.at(kedge) );
                }

                if ( iishared ) {
                    child->setSharedEdge( 4, shared_edges.at(iiedge) );
                }

                if ( jjshared ) {
                    child->setSharedEdge( 5, shared_edges.at(jjedge) );
                }

                child->setSharedEdge(3, i_shared_id);
                child->setSharedEdge(6, n_shared_id);
            }

#endif

#ifdef DEBUG_INFO
 #ifdef __PARALLEL_MODE
            OOFEM_LOG_INFO("[%d] Child %d generated on parent %d [%d] (leIndex %d, nds %d %d %d %d [%d %d %d %d], ngbs %d %d %d %d, irr %d %d %d %d %d %d [%d %d %d %d %d %d])\n",
                           mesh->giveSubdivision()->giveRank(), childNum,
                           this->number, this->giveGlobalNumber(), this->leIndex,
                           _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
                           mesh->giveNode( _nodes.at(1) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(2) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(3) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(4) )->giveGlobalNumber(),
                           child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
                           child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3),
                           child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6),
                           ( child->giveIrregular(1) ) ? mesh->giveNode( child->giveIrregular(1) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(2) ) ? mesh->giveNode( child->giveIrregular(2) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(3) ) ? mesh->giveNode( child->giveIrregular(3) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(4) ) ? mesh->giveNode( child->giveIrregular(4) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(5) ) ? mesh->giveNode( child->giveIrregular(5) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(6) ) ? mesh->giveNode( child->giveIrregular(6) )->giveGlobalNumber() : 0);
 #else
            OOFEM_LOG_INFO( "Child %d generated on parent %d (leIndex %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                           childNum, this->number, this->leIndex, _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
                           child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
                           child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3),
                           child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6) );
 #endif
#endif

            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = irregular_nodes.at(iedge);
            _nodes.at(3) = nodes.at(knode);
            _nodes.at(4) = nodes.at(nnode);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Tetra(childNum, mesh, this->number, _nodes);            // number, parent, coords
            mesh->addElement(child);
            children.at(2) = childNum;

            // set neigbour info
            if ( irregulars2 ) {
                child->setNeighbor( 1, this->giveNeighbor(nside) );
                child->setNeighbor( 3, this->giveNeighbor(iside) );
                child->setNeighbor( 4, this->giveNeighbor(jside) );

                child->setIrregular( 3, irregular_nodes.at(jedge) );
                child->setIrregular( 6, irregular_nodes.at(kkedge) );
                child->setIrregular( 4, irregular_nodes.at(iiedge) );
            } else {
                child->setNeighbor( 1, -this->giveNeighbor(nside) );
                child->setNeighbor( 3, -this->giveNeighbor(iside) );
                child->setNeighbor( 4, -this->giveNeighbor(jside) );
            }

            // neihgbor2 of child2 is changed to negative during subdivision (if any) of child1
            child->setNeighbor(2, childNum - 1);
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || jshared || iishared || kkshared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 2, shared_edges.at(iedge) );
                }

                if ( jshared ) {
                    child->setSharedEdge( 3, shared_edges.at(jedge) );
                }

                if ( iishared ) {
                    child->setSharedEdge( 4, shared_edges.at(iiedge) );
                }

                if ( kkshared ) {
                    child->setSharedEdge( 6, shared_edges.at(kkedge) );
                }

                child->setSharedEdge(1, i_shared_id);
                child->setSharedEdge(5, n_shared_id);
            }

#endif

#ifdef __PARALLEL_MODE
            // set connected elements
            if ( mesh->giveNode( nodes.at(inode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(inode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(nnode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(nnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(nnode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(nnode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(jnode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(jnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(knode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(knode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum);
            }

            if ( ishared ) {
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 1);
            }

#endif

#ifdef DEBUG_INFO
 #ifdef __PARALLEL_MODE
            OOFEM_LOG_INFO("[%d] Child %d generated on parent %d [%d] (leIndex %d, nds %d %d %d %d [%d %d %d %d], ngbs %d %d %d %d, irr %d %d %d %d %d %d [%d %d %d %d %d %d])\n",
                           mesh->giveSubdivision()->giveRank(), childNum,
                           this->number, this->giveGlobalNumber(), this->leIndex,
                           _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
                           mesh->giveNode( _nodes.at(1) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(2) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(3) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(4) )->giveGlobalNumber(),
                           child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
                           child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3),
                           child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6),
                           ( child->giveIrregular(1) ) ? mesh->giveNode( child->giveIrregular(1) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(2) ) ? mesh->giveNode( child->giveIrregular(2) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(3) ) ? mesh->giveNode( child->giveIrregular(3) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(4) ) ? mesh->giveNode( child->giveIrregular(4) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(5) ) ? mesh->giveNode( child->giveIrregular(5) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(6) ) ? mesh->giveNode( child->giveIrregular(6) )->giveGlobalNumber() : 0);
 #else
            OOFEM_LOG_INFO( "Child %d generated on parent %d (leIndex %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                           childNum, this->number, this->leIndex, _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
                           child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
                           child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3),
                           child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6) );
 #endif
#endif

            // set connected elements to nonshared nodes on outer boundary
#ifdef __PARALLEL_MODE
            if ( !ishared ) {
#endif
            if ( mesh->giveNode( irregular_nodes.at(iedge) )->giveNumber() < 0 ) {                   // check for marked local irregular
                // irregular node on outer boundary
                // insert connected elements
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(iedge) )->insertConnectedElement(childNum - 1);
            }

            // update connectivity of both end nodes of iedge if not shared
            // and ONLY if there is already list of connected elements
#ifdef __PARALLEL_MODE
            if ( mesh->giveNode( nodes.at(jnode) )->giveParallelMode() != DofManager_shared ) {                // update already done
#endif
            if ( mesh->giveNode( nodes.at(jnode) )->giveConnectedElements()->giveSize() ) {
                mesh->giveNode( nodes.at(jnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum - 1);
            }

#ifdef __PARALLEL_MODE
        }

#endif
#ifdef __PARALLEL_MODE
            if ( mesh->giveNode( nodes.at(knode) )->giveParallelMode() != DofManager_shared ) {                // update already done
#endif
            if ( mesh->giveNode( nodes.at(knode) )->giveConnectedElements()->giveSize() ) {
                mesh->giveNode( nodes.at(knode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum);
            }

#ifdef __PARALLEL_MODE
        }
#endif
#ifdef __PARALLEL_MODE
        }
#endif
#ifdef __PARALLEL_MODE
            if ( !iishared ) {
#endif
            // update connectivity of both end nodes of iiedge if not shared
            // and ONLY if there is already list of connected elements
#ifdef __PARALLEL_MODE
            if ( mesh->giveNode( nodes.at(inode) )->giveParallelMode() != DofManager_shared ) {                // update already done
#endif
            if ( mesh->giveNode( nodes.at(inode) )->giveConnectedElements()->giveSize() ) {
                mesh->giveNode( nodes.at(inode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 1);
            }

#ifdef __PARALLEL_MODE
        }
#endif
#ifdef __PARALLEL_MODE
            if ( mesh->giveNode( nodes.at(nnode) )->giveParallelMode() != DofManager_shared ) {                // update already done
#endif
            if ( mesh->giveNode( nodes.at(nnode) )->giveConnectedElements()->giveSize() ) {
                mesh->giveNode( nodes.at(nnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(nnode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(nnode) )->insertConnectedElement(childNum - 1);
            }

#ifdef __PARALLEL_MODE
        }
#endif
#ifdef __PARALLEL_MODE
        }
#endif
        } else {
            // count number of irregulars on each child;
            // if there is one irregular only, the set leIndex is valid
            // otherwise it is corrected later
            if ( irregular_nodes.at(iedge) ) {
                irregulars1++;
                irregulars2++;
                leIndex1 = leIndex2 = 2;
            }

            if ( irregular_nodes.at(jedge) ) {
                irregulars1++;
                leIndex1 = 3;
            }

            if ( irregular_nodes.at(kedge) ) {
                irregulars1++;
                leIndex1 = 1;
            }

            if ( irregular_nodes.at(jjedge) ) {
                irregulars2++;
                leIndex2 = 5;
            }

            if ( irregular_nodes.at(kkedge) ) {
                irregulars2++;
                leIndex2 = 6;
            }

            if ( irregulars1 > 1 ) {
                // use parent information about the longest edge on nside
                leIndex1 = this->side_leIndex.at(nside);

                if ( leIndex1 == jedge ) {
                    leIndex1 = 3;
                } else if ( leIndex1 == kedge ) {
                    leIndex1 = 1;
                } else {
#ifdef DEBUG_CHECK
                    if ( leIndex1 != iedge ) {
                        OOFEM_ERROR("side longest edge inconsistency on %d", this->number);
                    }

#endif
                    leIndex1 = 2;
                }
            }

            if ( irregulars2 > 1 ) {
                // use parent information about the longest edge on iside
                leIndex2 = this->side_leIndex.at(iside);

                if ( leIndex2 == kkedge ) {
                    leIndex2 = 6;
                } else if ( leIndex2 == jjedge ) {
                    leIndex2 = 5;
                } else {
#ifdef DEBUG_CHECK
                    if ( leIndex2 != iedge ) {
                        OOFEM_ERROR("side longest edge inconsistency on %d", this->number);
                    }

#endif
                    leIndex2 = 2;
                }
            }

#ifdef __PARALLEL_MODE
            int j_shared_id = 0, k_shared_id = 0;

            // check whether new edges are potentially shared
            if ( iishared ) {
                if ( mesh->giveNode( nodes.at(jnode) )->giveParallelMode() == DofManager_shared ) {
                    if ( !this->giveNeighbor(kside) ) {
                        _edge = new Subdivision :: RS_SharedEdge(this->mesh);
                        _edge->setEdgeNodes( irregular_nodes.at(iiedge), nodes.at(jnode) );

                        this->mesh->addEdge(_edge);
                        eNum++;

                        // get the tentative partitions
                        // they must be confirmed via mutual communication
                        IntArray partitions;
                        if ( _edge->giveSharedPartitions(partitions) ) {
                            j_shared_id = eNum;
                            _edge->setPartitions(partitions);
                            // put edge number into queue of shared edges to resolve remote partitions
                            sharedEdgesQueue.push_back(eNum);
                        }
                    }
                }

                if ( mesh->giveNode( nodes.at(knode) )->giveParallelMode() == DofManager_shared ) {
                    if ( !this->giveNeighbor(jside) ) {
                        _edge = new Subdivision :: RS_SharedEdge(this->mesh);
                        _edge->setEdgeNodes( irregular_nodes.at(iiedge), nodes.at(knode) );

                        this->mesh->addEdge(_edge);
                        eNum++;

                        // get the tentative partitions
                        // they must be confirmed via mutual communication
                        IntArray partitions;
                        if ( _edge->giveSharedPartitions(partitions) ) {
                            k_shared_id = eNum;
                            _edge->setPartitions(partitions);
                            // put edge number into queue of shared edges to resolve remote partitions
                            sharedEdgesQueue.push_back(eNum);
                        }
                    }
                }
            }

#endif

            _nodes.at(1) = nodes.at(inode);
            _nodes.at(2) = nodes.at(jnode);
            _nodes.at(3) = nodes.at(knode);
            _nodes.at(4) = irregular_nodes.at(iiedge);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Tetra(childNum, mesh, this->number, _nodes);            // number, parent, coords
            mesh->addElement(child);
            children.at(1) = childNum;
            // set neigbour info
            if ( irregulars1 ) {
                child->setNeighbor( 1, this->giveNeighbor(nside) );
                child->setNeighbor( 2, this->giveNeighbor(kside) );
                child->setNeighbor( 4, this->giveNeighbor(jside) );

                child->setIrregular( 1, irregular_nodes.at(kedge) );
                child->setIrregular( 3, irregular_nodes.at(jedge) );
                child->setIrregular( 2, irregular_nodes.at(iedge) );
            } else {
                child->setNeighbor( 1, -this->giveNeighbor(nside) );
                child->setNeighbor( 2, -this->giveNeighbor(kside) );
                child->setNeighbor( 4, -this->giveNeighbor(jside) );
            }

            // neihgbor3 of child1 is changed to negative during subdivision (if any) of child2
            child->setNeighbor(3, childNum + 1);
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || jshared || kshared || iishared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 2, shared_edges.at(iedge) );
                }

                if ( jshared ) {
                    child->setSharedEdge( 3, shared_edges.at(jedge) );
                }

                if ( kshared ) {
                    child->setSharedEdge( 1, shared_edges.at(kedge) );
                }

                if ( iishared ) {
                    child->setSharedEdge( 4, shared_edges.at(iiedge) );
                }

                child->setSharedEdge(5, j_shared_id);
                child->setSharedEdge(6, k_shared_id);
            }

#endif

#ifdef DEBUG_INFO
 #ifdef __PARALLEL_MODE
            OOFEM_LOG_INFO("[%d] Child %d generated on parent %d [%d] (leIndex %d, nds %d %d %d %d [%d %d %d %d], ngbs %d %d %d %d, irr %d %d %d %d %d %d [%d %d %d %d %d %d])\n",
                           mesh->giveSubdivision()->giveRank(), childNum,
                           this->number, this->giveGlobalNumber(), this->leIndex,
                           _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
                           mesh->giveNode( _nodes.at(1) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(2) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(3) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(4) )->giveGlobalNumber(),
                           child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
                           child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3),
                           child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6),
                           ( child->giveIrregular(1) ) ? mesh->giveNode( child->giveIrregular(1) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(2) ) ? mesh->giveNode( child->giveIrregular(2) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(3) ) ? mesh->giveNode( child->giveIrregular(3) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(4) ) ? mesh->giveNode( child->giveIrregular(4) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(5) ) ? mesh->giveNode( child->giveIrregular(5) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(6) ) ? mesh->giveNode( child->giveIrregular(6) )->giveGlobalNumber() : 0);
 #else
            OOFEM_LOG_INFO( "Child %d generated on parent %d (leIndex %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                           childNum, this->number, this->leIndex, _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
                           child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
                           child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3),
                           child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6) );
 #endif
#endif

            _nodes.at(1) = irregular_nodes.at(iiedge);
            _nodes.at(2) = nodes.at(jnode);
            _nodes.at(3) = nodes.at(knode);
            _nodes.at(4) = nodes.at(nnode);
            childNum = mesh->giveNumberOfElements() + 1;
            child = new Subdivision :: RS_Tetra(childNum, mesh, this->number, _nodes);            // number, parent, coords
            mesh->addElement(child);
            children.at(2) = childNum;

            // set neigbour info
            if ( irregulars2 ) {
                child->setNeighbor( 2, this->giveNeighbor(kside) );
                child->setNeighbor( 3, this->giveNeighbor(iside) );
                child->setNeighbor( 4, this->giveNeighbor(jside) );

                child->setIrregular( 5, irregular_nodes.at(jjedge) );
                child->setIrregular( 6, irregular_nodes.at(kkedge) );
                child->setIrregular( 2, irregular_nodes.at(iedge) );
            } else {
                child->setNeighbor( 2, -this->giveNeighbor(kside) );
                child->setNeighbor( 3, -this->giveNeighbor(iside) );
                child->setNeighbor( 4, -this->giveNeighbor(jside) );
            }

            // neihgbor1 of child2 is changed to negative during subdivision (if any) of child1
            child->setNeighbor(1, childNum - 1);
#ifdef __PARALLEL_MODE
            // set shared info
            if ( ishared || jjshared || kkshared || iishared ) {
                child->makeSharedEdges();
                if ( ishared ) {
                    child->setSharedEdge( 2, shared_edges.at(iedge) );
                }

                if ( jjshared ) {
                    child->setSharedEdge( 5, shared_edges.at(jjedge) );
                }

                if ( kkshared ) {
                    child->setSharedEdge( 6, shared_edges.at(kkedge) );
                }

                if ( iishared ) {
                    child->setSharedEdge( 4, shared_edges.at(iiedge) );
                }

                child->setSharedEdge(1, j_shared_id);
                child->setSharedEdge(3, k_shared_id);
            }

#endif

#ifdef __PARALLEL_MODE
            // set connected elements
            if ( mesh->giveNode( nodes.at(jnode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(jnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(knode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(knode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(inode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(inode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 1);
            }

            if ( mesh->giveNode( nodes.at(nnode) )->giveParallelMode() == DofManager_shared ) {
                mesh->giveNode( nodes.at(nnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(nnode) )->insertConnectedElement(childNum);
            }

            if ( iishared ) {
                mesh->giveNode( irregular_nodes.at(iiedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(iiedge) )->insertConnectedElement(childNum - 1);
            }

#endif

#ifdef DEBUG_INFO
 #ifdef __PARALLEL_MODE
            OOFEM_LOG_INFO("[%d] Child %d generated on parent %d [%d] (leIndex %d, nds %d %d %d %d [%d %d %d %d], ngbs %d %d %d %d, irr %d %d %d %d %d %d [%d %d %d %d %d %d])\n",
                           mesh->giveSubdivision()->giveRank(), childNum,
                           this->number, this->giveGlobalNumber(), this->leIndex,
                           _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
                           mesh->giveNode( _nodes.at(1) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(2) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(3) )->giveGlobalNumber(),
                           mesh->giveNode( _nodes.at(4) )->giveGlobalNumber(),
                           child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
                           child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3),
                           child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6),
                           ( child->giveIrregular(1) ) ? mesh->giveNode( child->giveIrregular(1) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(2) ) ? mesh->giveNode( child->giveIrregular(2) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(3) ) ? mesh->giveNode( child->giveIrregular(3) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(4) ) ? mesh->giveNode( child->giveIrregular(4) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(5) ) ? mesh->giveNode( child->giveIrregular(5) )->giveGlobalNumber() : 0,
                           ( child->giveIrregular(6) ) ? mesh->giveNode( child->giveIrregular(6) )->giveGlobalNumber() : 0);
 #else
            OOFEM_LOG_INFO( "Child %d generated on parent %d (leIndex %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                           childNum, this->number, this->leIndex, _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
                           child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
                           child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3),
                           child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6) );
 #endif
#endif

            // set connected elements to nonshared nodes on outer boundary
#ifdef __PARALLEL_MODE
            if ( !iishared ) {
#endif
            if ( mesh->giveNode( irregular_nodes.at(iiedge) )->giveNumber() < 0 ) {                   // check for marked local irregular
                // irregular node on outer boundary
                // insert connected elements
                mesh->giveNode( irregular_nodes.at(iiedge) )->insertConnectedElement(childNum);
                mesh->giveNode( irregular_nodes.at(iiedge) )->insertConnectedElement(childNum - 1);
            }

            // update connectivity of both end nodes of iiedge if not shared
            // and ONLY if there is already list of connected elements
#ifdef __PARALLEL_MODE
            if ( mesh->giveNode( nodes.at(inode) )->giveParallelMode() != DofManager_shared ) {                // update already done
#endif
            if ( mesh->giveNode( nodes.at(inode) )->giveConnectedElements()->giveSize() ) {
                mesh->giveNode( nodes.at(inode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(inode) )->insertConnectedElement(childNum - 1);
            }

#ifdef __PARALLEL_MODE
        }
#endif
#ifdef __PARALLEL_MODE
            if ( mesh->giveNode( nodes.at(nnode) )->giveParallelMode() != DofManager_shared ) {                // update already done
#endif
            if ( mesh->giveNode( nodes.at(nnode) )->giveConnectedElements()->giveSize() ) {
                mesh->giveNode( nodes.at(nnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(nnode) )->insertConnectedElement(childNum);
            }

#ifdef __PARALLEL_MODE
        }
#endif
#ifdef __PARALLEL_MODE
        }
#endif
#ifdef __PARALLEL_MODE
            if ( !ishared ) {
#endif
            // update connectivity of both end nodes of iedge if not shared
            // and ONLY if there is already list of connected elements
#ifdef __PARALLEL_MODE
            if ( mesh->giveNode( nodes.at(jnode) )->giveParallelMode() != DofManager_shared ) {                // update already done
#endif
            if ( mesh->giveNode( nodes.at(jnode) )->giveConnectedElements()->giveSize() ) {
                mesh->giveNode( nodes.at(jnode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(jnode) )->insertConnectedElement(childNum - 1);
            }

#ifdef __PARALLEL_MODE
        }
#endif
#ifdef __PARALLEL_MODE
            if ( mesh->giveNode( nodes.at(knode) )->giveParallelMode() != DofManager_shared ) {                // update already done
#endif
            if ( mesh->giveNode( nodes.at(knode) )->giveConnectedElements()->giveSize() ) {
                mesh->giveNode( nodes.at(knode) )->eraseConnectedElement( this->giveNumber() );
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum);
                mesh->giveNode( nodes.at(knode) )->insertConnectedElement(childNum - 1);
            }

#ifdef __PARALLEL_MODE
        }
#endif
#ifdef __PARALLEL_MODE
        }
#endif
        }

        if ( irregulars1 ) {
#ifdef DEBUG_CHECK
            if ( !mesh->giveElement( children.at(1) )->giveIrregular(leIndex1) ) {
                OOFEM_ERROR("leIndex incosistency for child 1 of %d", this->number);
            }

#endif
            mesh->giveElement( children.at(1) )->setLeIndex(leIndex1);
            mesh->giveElement( children.at(1) )->generate(sharedEdgesQueue);
        }

        if ( irregulars2 ) {
#ifdef DEBUG_CHECK
            if ( !mesh->giveElement( children.at(2) )->giveIrregular(leIndex2) ) {
                OOFEM_ERROR("leIndex incosistency for child 2 of %d", this->number);
            }

#endif
            mesh->giveElement( children.at(2) )->setLeIndex(leIndex2);
            mesh->giveElement( children.at(2) )->generate(sharedEdgesQueue);
        }

        // if there is neighbor of "this" (local or remote)  not designated for bisection
        // its neihgbor corresponding to "this" must be made negative to enforce update_neighbors
        for ( i = 1; i <= 4; i++ ) {
            if ( this->neghbours_base_elements.at(i) ) {
#ifdef DEBUG_CHECK
                if ( this->neghbours_base_elements.at(i) < 0 ) {
                    OOFEM_ERROR("negative neighbor of %d not expected", this->number);
                }

#endif

                ngb = mesh->giveElement( this->neghbours_base_elements.at(i) );
                if ( !ngb->hasIrregulars() ) {
                    ind = ngb->giveNeighbors()->findFirstIndexOf(this->number);
                    if ( ind ) {
                        ngb->setNeighbor(ind, -this->number);
                    }
                }
            }
        }
    }
}



void
Subdivision :: RS_Triangle :: update_neighbours()
{
    bool found;
    int i, iside, parentNeighbor;
    const IntArray *ca;

    /*
     * neghbours_base_elements pre-set in generate using this rule:
     * value > 0 .. real neighbour known, value equals to real neighbor on new mesh
     * value = 0 .. boundary
     * value < 0 .. neighbour set to parent element, actual neighbour has to be determined from parent children
     *             or if no child exists then to parent new counterpart.
     */

    for ( iside = 1; iside <= 3; iside++ ) {
        if ( this->neghbours_base_elements.at(iside) < 0 ) {
            parentNeighbor = -this->neghbours_base_elements.at(iside);
            // test if parentNeighbor has been subdivided
            if ( mesh->giveElement(parentNeighbor)->hasIrregulars() ) {                  // HUHU replace by !isTerminal
                // old neighbour bisected -> we loop over its children to find an appropriate new neighbor
                ca = mesh->giveElement(parentNeighbor)->giveChildren();
                found = false;
                for ( i = 1; i <= ca->giveSize(); i++ ) {
                    if ( this->isNeighborOf( mesh->giveElement( ca->at(i) ) ) ) {
                        this->neghbours_base_elements.at(iside) = ca->at(i);
                        found = true;
                        break;
                    }
                }

                if ( !found ) {
                    OOFEM_ERROR("update_neighbours failed for element %d", this->number);
                }
            } else {
                // parent neighbour remains actual neighbour
                this->neghbours_base_elements.at(iside) = parentNeighbor;
            }
        } else if ( this->neghbours_base_elements.at(iside) > 0 ) {
            // neighbor element already set
        }
    } // end loop over element sides
}


void
Subdivision :: RS_Tetra :: update_neighbours()
{
    bool found;
    int i, j, k, iside, parentNeighbor;
    const IntArray *ca1, *ca2, *ca3;

    /*
     * neghbours_base_elements pre-set in generate using this rule:
     * value > 0 .. real neighbour known, value equals to real neighbor on new mesh
     * value = 0 .. boundary
     * value < 0 .. neighbour set to parent element, actual neighbour has to be determined from terminal children of parent
     *             or if no child exists then to parent new counterpart.
     */

    for ( iside = 1; iside <= 4; iside++ ) {
        if ( this->neghbours_base_elements.at(iside) < 0 ) {
            parentNeighbor = -this->neghbours_base_elements.at(iside);
            // test if parentNeighbor has been subdivided
            if ( mesh->giveElement(parentNeighbor)->hasIrregulars() ) {                // HUHU replace by !isTerminal
                // old neighbour bisected -> we loop over its terminal children to find an appropriate new neighbor;
                // there may be at maximum 3 levels of parent tetra subdivision
                ca1 = mesh->giveElement(parentNeighbor)->giveChildren();
                found = false;
                for ( i = 1; i <= ca1->giveSize(); i++ ) {
                    if ( mesh->giveElement( ca1->at(i) )->isTerminal() ) {
                        if ( this->isNeighborOf( mesh->giveElement( ca1->at(i) ) ) ) {
                            this->neghbours_base_elements.at(iside) = ca1->at(i);
                            found = true;
                            break;
                        }
                    } else {
                        ca2 = mesh->giveElement( ca1->at(i) )->giveChildren();
                        for ( j = 1; j <= ca2->giveSize(); j++ ) {
                            if ( mesh->giveElement( ca2->at(j) )->isTerminal() ) {
                                if ( this->isNeighborOf( mesh->giveElement( ca2->at(j) ) ) ) {
                                    this->neghbours_base_elements.at(iside) = ca2->at(j);
                                    found = true;
                                    break;
                                }
                            } else {
                                ca3 = mesh->giveElement( ca2->at(j) )->giveChildren();
                                for ( k = 1; k <= ca3->giveSize(); k++ ) {
                                    if ( mesh->giveElement( ca3->at(k) )->isTerminal() ) {
                                        if ( this->isNeighborOf( mesh->giveElement( ca3->at(k) ) ) ) {
                                            this->neghbours_base_elements.at(iside) = ca3->at(k);
                                            found = true;
                                            break;
                                        }
                                    }

#ifdef DEBUG_CHECK
                                    else {
                                        OOFEM_ERROR("too many subdivision levels for element %d", this->number);
                                    }
#endif
                                }
                            }

                            if ( found ) {
                                break;
                            }
                        }
                    }

                    if ( found ) {
                        break;
                    }
                }

                if ( !found ) {
                    OOFEM_ERROR("failed for element %d (side %d, elem %d)",
                                 this->number, iside, parentNeighbor);
                }
            } else {
                // parent neighbour remains actual neighbour
                this->neghbours_base_elements.at(iside) = parentNeighbor;
            }
        } else if ( this->neghbours_base_elements.at(iside) > 0 ) {
            // neighbor element already set
        }
    } // end loop over element side faces

#ifdef DEBUG_CHECK
    // check updated neighbors
    IntArray snodes1, snodes2;
    RS_Element *ngb;
    for ( i = 1; i <= 4; i++ ) {
        if ( this->neghbours_base_elements.at(i) ) {
            if ( this->neghbours_base_elements.at(i) > 0 ) {
                this->giveSideNodes(i, snodes1);
                ngb = mesh->giveElement( this->neghbours_base_elements.at(i) );
                if ( ngb->giveNumber() < 0 ) {
                    OOFEM_ERROR("neighbour %d of %d is temporary",
                                 this->neghbours_base_elements.at(i), this->number);
                }

                if ( !ngb->isTerminal() ) {
                    OOFEM_ERROR("neighbour %d of %d not terminal",
                                 this->neghbours_base_elements.at(i), this->number);
                }

                if ( ngb->giveNumber() < this->number ) {
                    // ngb has been already processed
                    j = ngb->giveNeighbors()->findFirstIndexOf(this->number);
                    if ( j ) {
                        ngb->giveSideNodes(j, snodes2);
                        for ( k = 1; k <= 3; k++ ) {
                            if ( snodes2.findFirstIndexOf( snodes1.at(k) ) ) {
                                continue;
                            }

                            OOFEM_ERROR( "element %d is not neighbor (%d) of element %d",
                                         this->number, i, this->neghbours_base_elements.at(i) );
                        }
                    } else {
                        OOFEM_ERROR( "element %d is not neighbor (%d) of element %d",
                                     this->number, i, this->neghbours_base_elements.at(i) );
                    }
                } else {
                    // ngb has not been yet processed ==> I cannot check particular side of ngb
                    for ( k = 1; k <= 3; k++ ) {
                        if ( ngb->giveNodes()->findFirstIndexOf( snodes1.at(k) ) ) {
                            continue;
                        }

                        OOFEM_ERROR( "element %d is not neighbor of element %d",
                                     this->number, this->neghbours_base_elements.at(i) );
                    }
                }
            } else {
                OOFEM_ERROR("negative neighbour %d of %d not expected",
                             this->neghbours_base_elements.at(i), this->number);
            }
        }
    }

#endif
}


bool
Subdivision :: RS_Triangle :: isNeighborOf(Subdivision :: RS_Element *elem)
{
    // Simplified implementation, considering only one element type - triangles
    int i, _c = 0;
    for ( i = 1; i <= 3; i++ ) {
        _c +=  elem->containsNode( this->nodes.at(i) );
    }

    return ( _c == 2 );
}


bool
Subdivision :: RS_Tetra :: isNeighborOf(Subdivision :: RS_Element *elem)
{
    // Simplified implementation, considering only one element type - tetras
    int i, _c = 0;
    for ( i = 1; i <= 4; i++ ) {
        _c +=  elem->containsNode( this->nodes.at(i) );
    }

    return ( _c == 3 );
}


void
Subdivision :: RS_Triangle :: giveSideNodes(int iside, IntArray &snodes)
{
    int inode, jnode;

    inode = iside;
    jnode = ( iside < 3 ) ? iside + 1 : 1;

    snodes.resize(2);
    snodes.at(1) = nodes.at(inode);
    snodes.at(2) = nodes.at(jnode);
}


void
Subdivision :: RS_Tetra :: giveSideNodes(int iside, IntArray &snodes)
{
    int inode, jnode, knode;

    if ( iside == 1 ) {
        inode = 1;
        jnode = 2;
        knode = 3;
    } else {
        inode = iside - 1;
        jnode = ( iside < 4 ) ? iside : 1;
        knode = 4;
    }

    snodes.resize(3);
    snodes.at(1) = nodes.at(inode);
    snodes.at(2) = nodes.at(jnode);
    snodes.at(3) = nodes.at(knode);
}



double
Subdivision :: RS_Triangle :: giveDensity()
{
    double x1, x2, x3, y1, y2, y3;
    x1 = mesh->giveNode( nodes.at(1) )->giveCoordinate(1);
    x2 = mesh->giveNode( nodes.at(2) )->giveCoordinate(1);
    x3 = mesh->giveNode( nodes.at(3) )->giveCoordinate(1);

    y1 = mesh->giveNode( nodes.at(1) )->giveCoordinate(2);
    y2 = mesh->giveNode( nodes.at(2) )->giveCoordinate(2);
    y3 = mesh->giveNode( nodes.at(3) )->giveCoordinate(2);

    return sqrt( fabs( 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) ) );
}


double
Subdivision :: RS_Tetra :: giveDensity()
{
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    double dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3;
    double vol, d;
    x1 = mesh->giveNode( nodes.at(1) )->giveCoordinate(1);
    x2 = mesh->giveNode( nodes.at(2) )->giveCoordinate(1);
    x3 = mesh->giveNode( nodes.at(3) )->giveCoordinate(1);
    x4 = mesh->giveNode( nodes.at(4) )->giveCoordinate(1);

    y1 = mesh->giveNode( nodes.at(1) )->giveCoordinate(2);
    y2 = mesh->giveNode( nodes.at(2) )->giveCoordinate(2);
    y3 = mesh->giveNode( nodes.at(3) )->giveCoordinate(2);
    y4 = mesh->giveNode( nodes.at(4) )->giveCoordinate(2);

    z1 = mesh->giveNode( nodes.at(1) )->giveCoordinate(3);
    z2 = mesh->giveNode( nodes.at(2) )->giveCoordinate(3);
    z3 = mesh->giveNode( nodes.at(3) )->giveCoordinate(3);
    z4 = mesh->giveNode( nodes.at(4) )->giveCoordinate(3);

    dx1 = x2 - x1;
    dy1 = y2 - y1;
    dz1 = z2 - z1;
    dx2 = x3 - x1;
    dy2 = y3 - y1;
    dz2 = z3 - z1;
    dx3 = x4 - x1;
    dy3 = y4 - y1;
    dz3 = z4 - z1;

    vol = ( dx3 * ( dy1 * dz2 - dz1 * dy2 ) + dy3 * ( dz1 * dx2 - dx1 * dz2 ) + dz3 * ( dx1 * dy2 - dy1 * dx2 ) ) / 6.0;
    d = exp(log(vol) / 3.0);

    return d;
}


void
Subdivision :: RS_Triangle :: importConnectivity(ConnectivityTable *ct)
{
    IntArray me(1), conn;
    int iNode, jNode, el, i;

    neghbours_base_elements.resize(3);
    neghbours_base_elements.zero();           // initialized to have no neighbour

#if 0   // old version
    me.at(1) = this->number;
    ct->giveElementNeighbourList(conn, me);

    int iside;
    for ( iside = 1; iside <= 3; iside++ ) {
        iNode = nodes.at(iside);
        jNode = nodes.at( ( iside == 3 ) ? 1 : iside + 1 );

        // select right neighbour
        for ( i = 1; i <= conn.giveSize(); i++ ) {
            el = conn.at(i);
            if ( el == this->number ) {
                continue;
            }

 #ifdef __PARALLEL_MODE
            if ( mesh->giveElement(el)->giveParallelMode() != Element_local ) {
                continue;
            }

 #endif
            if ( mesh->giveElement(el)->containsNode(iNode) &&
                mesh->giveElement(el)->containsNode(jNode) ) {
                neghbours_base_elements.at(iside) = el;
                break;
            }
        }
    }

#endif

    me.at(1) = nodes.at(3);
    ct->giveNodeNeighbourList(conn, me);
    iNode = nodes.at(1);
    jNode = nodes.at(2);

    // select right neighbour
    for ( i = 1; i <= conn.giveSize(); i++ ) {
        el = conn.at(i);
        if ( el == this->number ) {
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( mesh->giveElement(el)->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        if ( mesh->giveElement(el)->containsNode(iNode) ) {
            neghbours_base_elements.at(3) = el;
            break;
        }
    }

    for ( i = 1; i <= conn.giveSize(); i++ ) {
        el = conn.at(i);
        if ( el == this->number ) {
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( mesh->giveElement(el)->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        if ( mesh->giveElement(el)->containsNode(jNode) ) {
            neghbours_base_elements.at(2) = el;
            break;
        }
    }

    me.at(1) = iNode;
    ct->giveNodeNeighbourList(conn, me);

    // select right neighbour
    for ( i = 1; i <= conn.giveSize(); i++ ) {
        el = conn.at(i);
        if ( el == this->number ) {
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( mesh->giveElement(el)->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        if ( mesh->giveElement(el)->containsNode(jNode) ) {
            neghbours_base_elements.at(1) = el;
            break;
        }
    }
}


void
Subdivision :: RS_Tetra :: importConnectivity(ConnectivityTable *ct)
{
    IntArray me(1), conn;
    int iNode, jNode, kNode, el, i;

    neghbours_base_elements.resize(4);
    neghbours_base_elements.zero();           // initialized to have no neighbour

#if 0   // old version
    me.at(1) = this->number;
    ct->giveElementNeighbourList(conn, me);

    iNode = nodes.at(1);
    jNode = nodes.at(2);
    kNode = nodes.at(3);

    // select right base neighbour
    for ( i = 1; i <= conn.giveSize(); i++ ) {
        el = conn.at(i);
        if ( el == this->number ) {
            continue;
        }

 #ifdef __PARALLEL_MODE
        if ( mesh->giveElement(el)->giveParallelMode() != Element_local ) {
            continue;
        }

 #endif
        if ( mesh->giveElement(el)->containsNode(iNode) &&
            mesh->giveElement(el)->containsNode(jNode) &&
            mesh->giveElement(el)->containsNode(kNode) ) {
            neghbours_base_elements.at(1) = el;
            break;
        }
    }

    kNode = nodes.at(4);
    int iside;
    for ( iside = 1; iside <= 3; iside++ ) {
        iNode = nodes.at(iside);
        jNode = nodes.at( ( iside == 3 ) ? 1 : iside + 1 );

        // select right side neighbour
        for ( i = 1; i <= conn.giveSize(); i++ ) {
            el = conn.at(i);
            if ( el == this->number ) {
                continue;
            }

 #ifdef __PARALLEL_MODE
            if ( mesh->giveElement(el)->giveParallelMode() != Element_local ) {
                continue;
            }

 #endif
            if ( mesh->giveElement(el)->containsNode(iNode) &&
                mesh->giveElement(el)->containsNode(jNode) &&
                mesh->giveElement(el)->containsNode(kNode) ) {
                neghbours_base_elements.at(iside + 1) = el;
                break;
            }
        }
    }

#endif

    bool has_i, has_j, has_k;

    me.at(1) = nodes.at(4);
    ct->giveNodeNeighbourList(conn, me);
    iNode = nodes.at(1);
    jNode = nodes.at(2);
    kNode = nodes.at(3);

    // select right neighbour
    for ( i = 1; i <= conn.giveSize(); i++ ) {
        el = conn.at(i);
        if ( el == this->number ) {
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( mesh->giveElement(el)->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        has_i = mesh->giveElement(el)->containsNode(iNode);
        has_j = mesh->giveElement(el)->containsNode(jNode);
        has_k = mesh->giveElement(el)->containsNode(kNode);

        if ( has_i && has_j ) {
            neghbours_base_elements.at(2) = el;
        }

        if ( has_j && has_k ) {
            neghbours_base_elements.at(3) = el;
        }

        if ( has_k && has_i ) {
            neghbours_base_elements.at(4) = el;
        }
    }

    me.at(1) = iNode;
    ct->giveNodeNeighbourList(conn, me);

    // select right neighbour
    for ( i = 1; i <= conn.giveSize(); i++ ) {
        el = conn.at(i);
        if ( el == this->number ) {
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( mesh->giveElement(el)->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        if ( mesh->giveElement(el)->containsNode(jNode) &&
            mesh->giveElement(el)->containsNode(kNode) ) {
            neghbours_base_elements.at(1) = el;
            break;
        }
    }

    /*
     * #ifdef DEBUG_INFO
     * OOFEM_LOG_INFO("Elem %d connectivity imported (nds %d %d %d %d ngbs %d %d %d %d)\n", this->number,
     *                                                       nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4),
     *                                                       neghbours_base_elements.at(1), neghbours_base_elements.at(2),
     *                                                       neghbours_base_elements.at(3), neghbours_base_elements.at(4));
     ******#endif
     */
}


#ifdef __OOFEG
void
Subdivision :: RS_Node :: drawGeometry()
//
// draws graphics representation of receiver
//
{
    GraphicObj *go;
    EPixel color;
    BOOLEAN suc;
    const char *colors[] = {
        "orange", "black"
    };

    WCRec p [ 1 ]; /* point */
    p [ 0 ].x = ( FPNum ) this->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveCoordinate(3);

    EASValsSetMType(FILLED_CIRCLE_MARKER);
    color = ColorGetPixelFromString(const_cast< char * >(colors [ 0 ]), & suc);
    EASValsSetColor(color);
    //EASValsSetColor( gc.getNodeColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetMSize(8);
    go = CreateMarker3D(p);
    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | MTYPE_MASK | MSIZE_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);

    char num [ 6 ];
    color = ColorGetPixelFromString(const_cast< char * >(colors [ 1 ]), & suc);
    EASValsSetColor(color);
    //EASValsSetColor( gc.getNodeColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveCoordinate(3);
    sprintf(num, "%d", this->number);
    go = CreateAnnText3D(p, num);
    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void
Subdivision :: RS_Triangle :: drawGeometry()
{
    WCRec p [ 3 ];
    GraphicObj *go;
    EPixel color;
    BOOLEAN suc;
    const char *colors[] = {
        "DodgerBlue", "black"
    };

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    color = ColorGetPixelFromString(const_cast< char * >(colors [ 0 ]), & suc);
    EASValsSetColor(color);
    color = ColorGetPixelFromString(const_cast< char * >(colors [ 1 ]), & suc);
    EASValsSetEdgeColor(color);
    //EASValsSetColor( gc.getElementColor() );
    //EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    //    EASValsSetShrink(0.8);
    p [ 0 ].x = ( FPNum ) mesh->giveNode( nodes.at(1) )->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) mesh->giveNode( nodes.at(1) )->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) mesh->giveNode( nodes.at(2) )->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) mesh->giveNode( nodes.at(2) )->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) mesh->giveNode( nodes.at(3) )->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) mesh->giveNode( nodes.at(3) )->giveCoordinate(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    //EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | SHRINK_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void
Subdivision :: RS_Tetra :: drawGeometry()
{
    WCRec p [ 4 ];
    GraphicObj *go;
    EPixel color;
    BOOLEAN suc;
    const char *colors[] = {
        "DodgerBlue", "black"
    };

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    color = ColorGetPixelFromString(const_cast< char * >(colors [ 0 ]), & suc);
    EASValsSetColor(color);
    color = ColorGetPixelFromString(const_cast< char * >(colors [ 1 ]), & suc);
    EASValsSetEdgeColor(color);
    //EASValsSetColor( gc.getElementColor() );
    //EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    //EASValsSetShrink(0.8);
    p [ 0 ].x = ( FPNum ) mesh->giveNode( nodes.at(1) )->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) mesh->giveNode( nodes.at(1) )->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) mesh->giveNode( nodes.at(1) )->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) mesh->giveNode( nodes.at(2) )->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) mesh->giveNode( nodes.at(2) )->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) mesh->giveNode( nodes.at(2) )->giveCoordinate(3);
    p [ 2 ].x = ( FPNum ) mesh->giveNode( nodes.at(3) )->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) mesh->giveNode( nodes.at(3) )->giveCoordinate(2);
    p [ 2 ].z = ( FPNum ) mesh->giveNode( nodes.at(3) )->giveCoordinate(3);
    p [ 3 ].x = ( FPNum ) mesh->giveNode( nodes.at(4) )->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) mesh->giveNode( nodes.at(4) )->giveCoordinate(2);
    p [ 3 ].z = ( FPNum ) mesh->giveNode( nodes.at(4) )->giveCoordinate(3);

    go =  CreateTetra(p);
    //EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | SHRINK_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK | FILL_MASK, go);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK | FILL_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif


MesherInterface :: returnCode
Subdivision :: createMesh(TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew)
{
    // import data from old domain
    int parent, nnodes = domain->giveNumberOfDofManagers(), nelems = domain->giveNumberOfElements();
    IntArray enodes;
    Subdivision :: RS_Node *_node;
    Subdivision :: RS_Element *_element;
    Timer timer;

    timer.startTimer();

    if ( this->mesh ) {
        delete mesh;
    }

    mesh = new Subdivision :: RS_Mesh(this);

    // import nodes
    // all nodes (including remote which are not needed) are imported to ensure consistency between
    // node number (in the mesh) and its parent number (in the domain) because the domain is used to import connectivities
    for ( int i = 1; i <= nnodes; i++ ) {
        _node = new Subdivision :: RS_Node( i, mesh, i, * ( domain->giveNode ( i )->giveCoordinates() ),
                                           domain->giveErrorEstimator ( )->giveRemeshingCrit ( )->giveRequiredDofManDensity ( i, tStep ),
                                           domain->giveNode ( i )->isBoundary() );
        _node->setGlobalNumber( domain->giveNode(i)->giveGlobalNumber() );
#ifdef __PARALLEL_MODE
        _node->setParallelMode( domain->giveNode(i)->giveParallelMode() );
        _node->setPartitions( * domain->giveNode(i)->givePartitionList() );
#endif
        this->mesh->addNode(_node);
    }

    // import elements
    // all elements (including remote which are not needed) are imported to ensure consistency between
    // element number (in the mesh) and its parent number (in the domain) because the domain is used to import connectivities
    for ( int i = 1; i <= nelems; i++ ) {
        if ( domain->giveElement(i)->giveGeometryType() == EGT_triangle_1 ) {
            enodes.resize(3);
            for ( int j = 1; j <= 3; j++ ) {
                enodes.at(j) = domain->giveElement(i)->giveDofManagerNumber(j);
            }

            _element = new Subdivision :: RS_Triangle(i, mesh, i, enodes);
        } else if ( domain->giveElement(i)->giveGeometryType() == EGT_tetra_1 ) {
            enodes.resize(4);
            for ( int j = 1; j <= 4; j++ ) {
                enodes.at(j) = domain->giveElement(i)->giveDofManagerNumber(j);
            }

            _element = new Subdivision :: RS_Tetra(i, mesh, i, enodes);
        } else {
            OOFEM_ERROR("Unsupported element geometry (element %d)", i);
            _element = NULL;
        }
        _element->setGlobalNumber( domain->giveElement(i)->giveGlobalNumber() );
#ifdef __PARALLEL_MODE
        _element->setParallelMode( domain->giveElement(i)->giveParallelMode() );
#endif
        this->mesh->addElement(_element);

    }

    // import connectivities for local elements only
    for ( int i = 1; i <= nelems; i++ ) {
#ifdef __PARALLEL_MODE
        if ( this->mesh->giveElement(i)->giveParallelMode() != Element_local ) {
            continue;
        }

#endif
        this->mesh->giveElement(i)->importConnectivity( domain->giveConnectivityTable() );
    }

#ifdef __PARALLEL_MODE
    // import connectivities for shared nodes only
    // build global shared node map (for top level only)
    this->mesh->initGlobalSharedNodeMap();
    for ( int i = 1; i <= nnodes; i++ ) {
        _node = this->mesh->giveNode(i);
        if ( _node->giveParallelMode() != DofManager_shared ) {
            continue;
        }

        _node->importConnectivity( domain->giveConnectivityTable() );
        this->mesh->insertGlobalSharedNodeMap(_node);
    }

    // build array of shared edges (for top level only)
    // this can be done after connectivity is available for all shared nodes
    for ( int i = 1; i <= nnodes; i++ ) {
        _node = this->mesh->giveNode(i);
        if ( _node->giveParallelMode() != DofManager_shared ) {
            continue;
        }

        _node->numberSharedEdges();
    }

    // initiate queue for shared edge exchange
    for ( int i = 1; i <= mesh->giveNumberOfEdges(); i++ ) {
        if ( mesh->giveEdge(i)->givePartitions()->giveSize() ) {
            sharedEdgesQueue.push_back(i);
        }
    }

    // exchange shared edges (on top level)
    this->exchangeSharedEdges();
#endif

#ifdef __OOFEG
 #ifdef DRAW_MESH_BEFORE_BISECTION
    nelems = mesh->giveNumberOfElements();
    for ( int i = 1; i <= nelems; i++ ) {
  #ifdef __PARALLEL_MODE
        if ( this->mesh->giveElement(i)->giveParallelMode() != Element_local ) {
            continue;
        }

  #endif
        mesh->giveElement(i)->drawGeometry();
    }

    ESIEventLoop(YES, "Before bisection; Press Ctrl-p to continue");
 #endif
#endif

    // bisect mesh
    this->bisectMesh();
    // smooth mesh
    if ( smoothingFlag ) {
        // smooth only if there are new irregulars;
        // this is reasonable only if the smoothing is done locally !!!
        if ( nnodes != mesh->giveNumberOfNodes() ) {
            this->smoothMesh();
        }
    }

#ifdef __OOFEG
 #ifdef DRAW_MESH_AFTER_BISECTION
    nelems = mesh->giveNumberOfElements();
    for ( int i = 1; i <= nelems; i++ ) {
        if ( !mesh->giveElement(i)->isTerminal() ) {
            continue;
        }

  #ifdef __PARALLEL_MODE
        if ( this->mesh->giveElement(i)->giveParallelMode() != Element_local ) {
            continue;
        }

  #endif
        mesh->giveElement(i)->drawGeometry();
    }

    ESIEventLoop(YES, "After bisection; Press Ctrl-p to continue");
 #endif
#endif

    Dof *dof;
    DofManager *parentNodePtr, *node;
    Element *elem;
    CrossSection *crossSection;
    Material *mat;
    NonlocalBarrier *barrier;
    GeneralBoundaryCondition *bc;
    InitialCondition *ic;
    Function *func;
    Set *set;
    std :: string name;

    // create new mesh (missing param for new mesh!)
    nnodes = mesh->giveNumberOfNodes();
    ( * dNew ) = new Domain( 2, domain->giveSerialNumber() + 1, domain->giveEngngModel() );
    ( * dNew )->setDomainType( domain->giveDomainType() );

    // copy dof managers
    ( * dNew )->resizeDofManagers(nnodes);
    const IntArray dofIDArrayPtr = domain->giveDefaultNodeDofIDArry();
    for ( int inode = 1; inode <= nnodes; inode++ ) {
        parent = mesh->giveNode(inode)->giveParent();
        if ( parent ) {
            parentNodePtr = domain->giveNode(parent);
            // inherit all data from parent (bc, ic, load, etc.)
            node = classFactory.createDofManager(parentNodePtr->giveInputRecordName(), inode, * dNew);
            node->setNumberOfDofs(0);
            node->setLoadArray( * parentNodePtr->giveLoadArray() );
            // create individual DOFs
            for ( Dof *idofPtr: *parentNodePtr ) {
                if ( idofPtr->isPrimaryDof() ) {
                    dof = new MasterDof( node, idofPtr->giveBcId(), idofPtr->giveIcId(), idofPtr->giveDofID() );
                } else {
                    SimpleSlaveDof *simpleSlaveDofPtr = dynamic_cast< SimpleSlaveDof * >(idofPtr);
                    if ( simpleSlaveDofPtr ) {
                        dof = new SimpleSlaveDof( node, simpleSlaveDofPtr->giveMasterDofManagerNum(), idofPtr->giveDofID() );
                    } else {
                        OOFEM_ERROR("unsupported DOF type");
                        dof = NULL;
                    }
                }
                node->appendDof(dof);
            }

            node->setGlobalNumber( parentNodePtr->giveGlobalNumber() );
#ifdef __PARALLEL_MODE
            node->setParallelMode( parentNodePtr->giveParallelMode() );
            node->setPartitionList( parentNodePtr->givePartitionList() );
#endif
        } else {
            // newly created node (irregular)
            node = new Node(inode, *dNew);
            //create new node with default DOFs
            int ndofs = dofIDArrayPtr.giveSize();
            node->setNumberOfDofs(0);

            // create individual DOFs
            for ( int idof = 1; idof <= ndofs; idof++ ) {
#ifdef NM
                dof = NULL;
                FloatArray *coords = mesh->giveNode(inode)->giveCoordinates();
                if ( !dof ) {
                    if ( fabs(coords->at(1) - 200.0) < 0.000001 ) {
                        if ( coords->at(2) > -0.000001 && coords->at(2) < 97.500001 ) {
                            dof = new MasterDof( node, 1, 0, dofIDArrayPtr.at ( idof ) );
                        }
                    }
                }

                if ( !dof ) {
                    if ( fabs( coords->at(2) ) < 0.000001 ) {
                        dof = new MasterDof( node, 1, 0, dofIDArrayPtr.at ( idof ) );
                    }
                }

                if ( !dof ) {
                    if ( fabs( coords->at(1) ) < 0.000001 ) {
                        if ( coords->at(2) > 102.4999999 && coords->at(2) < 199.9999999 ) {
                            dof = new SimpleSlaveDof( node, 5, ( DofID ) dofIDArrayPtr.at ( idof ) );                                // HUHU
                        }
                    }
                }

                if ( !dof ) {
                    if ( fabs(coords->at(2) - 200.0) < 0.000001 ) {
                        if ( coords->at(1) > 0.000001 && coords->at(1) < 200.000001 ) {
                            dof = new SimpleSlaveDof( node, 6, ( DofID ) dofIDArrayPtr.at ( idof ) );                                  // HUHU
                        }
                    }
                }

                if ( !dof ) {
                    dof = new MasterDof( node, 0, 0, dofIDArrayPtr.at ( idof ) );
                }

#else
 #ifdef BRAZIL_2D
                dof = NULL;
                FloatArray *coords = mesh->giveNode(inode)->giveCoordinates();
                if ( !dof ) {
                    if ( fabs( coords->at(1) ) < 0.000001 ) {
                        if ( coords->at(2) > 0.000001 ) {
                            if ( idof == 1 ) {
                                dof = new MasterDof( node, 1, 0, dofIDArrayPtr.at ( idof ) );
                            }
                        }
                    }
                }

                if ( !dof ) {
                    if ( fabs( coords->at(2) ) < 0.000001 ) {
                        if ( coords->at(1) > 0.000001 ) {
                            if ( idof == 2 ) {
                                dof = new MasterDof( node, 1, 0, dofIDArrayPtr.at ( idof ) );
                            }
                        }
                    }
                }

                if ( !dof ) {
                    dof = new MasterDof( node, 0, 0, dofIDArrayPtr.at ( idof ) );
                }

 #else
  #ifdef THREEPBT_3D
                dof = NULL;
                FloatArray *coords = mesh->giveNode(inode)->giveCoordinates();
                if ( !dof ) {
                    if ( fabs( coords->at(1) ) < 0.000001 ) {
                        if ( coords->at(2) > 0.000001 ) {
                            if ( coords->at(3) > 499.999999 ) {
                                if ( idof == 1 || idof == 3 ) {
                                    dof = new MasterDof( node, 1, 0, dofIDArrayPtr.at ( idof ) );
                                }
                            }
                        }
                    }
                }

                if ( !dof ) {
                    if ( coords->at(1) > 1999.999999 ) {
                        if ( coords->at(2) > 0.000001 ) {
                            if ( coords->at(3) > 499.999999 ) {
                                if ( idof == 3 ) {
                                    dof = new MasterDof( node, 1, 0, dofIDArrayPtr.at ( idof ) );
                                }
                            }
                        }
                    }
                }

                if ( !dof ) {
                    dof = new MasterDof( node, 0, 0, dofIDArrayPtr.at ( idof ) );
                }

  #else
   #ifdef HEADEDSTUD
                double dist, rad;
                FloatArray *coords = mesh->giveNode(inode)->giveCoordinates();

                dof = NULL;
                if ( !dof ) {
                    dist = coords->at(1) * coords->at(1) + coords->at(3) * coords->at(3);
                    if ( coords->at(2) < 88.000001 ) {
                        if ( coords->at(2) > 70.000001 ) {
                            if ( fabs(dist - 7.0 * 7.0) < 0.01 ) {                            // be very tolerant (geometry is not precise)
                                if ( idof == 1 || idof == 3 ) {
                                    dof = new MasterDof( node, 1, 0, ( DofIDItem ) dofIDArrayPtr.at ( idof ) );
                                }
                            }
                        } else if ( coords->at(2) > 67.9999999 ) {
                            rad = 18.0 - 11.0 / 5.5 * ( coords->at(2) - 64.5 );
                            if ( fabs(dist - rad * rad) < 0.01 ) {                            // be very tolerant (geometry is not precise)
                                if ( idof == 1 || idof == 3 ) {
                                    dof = new MasterDof( node, 1, 0, ( DofIDItem ) dofIDArrayPtr.at ( idof ) );
                                }

                                if ( idof == 2 ) {
                                    dof = new MasterDof( node, 2, 0, ( DofIDItem ) dofIDArrayPtr.at ( idof ) );
                                    OOFEM_LOG_INFO("Subdivision: Conrolled Node %d", inode);
                                }
                            }
                        }
                    }
                }

                if ( !dof ) {
                    if ( coords->at(1) > 299.999999 || coords->at(3) > 299.999999 ) {
                        dof = new MasterDof( node, 1, 0, ( DofIDItem ) dofIDArrayPtr.at ( idof ) );
                    }
                }

                if ( !dof ) {
                    if ( coords->at(1) < 0.00000001 ) {
                        if ( idof == 1 ) {
                            dof = new MasterDof( node, 1, 0, ( DofIDItem ) dofIDArrayPtr.at ( idof ) );
                        }
                    }
                }

                if ( !dof ) {
                    dof = new MasterDof( node, 0, 0, ( DofIDItem ) dofIDArrayPtr.at ( idof ) );
                }

   #else
                dof = new MasterDof( node, 0, 0, dofIDArrayPtr.at ( idof ) );
   #endif
  #endif
 #endif
#endif
                node->appendDof(dof);
            }

            node->setGlobalNumber( mesh->giveNode(inode)->giveGlobalNumber() );
#ifdef __PARALLEL_MODE
            node->setParallelMode( mesh->giveNode(inode)->giveParallelMode() );
            node->setPartitionList( mesh->giveNode(inode)->givePartitions() );
#endif
        }

        // set node coordinates
        static_cast< Node * >(node)->setCoordinates( * mesh->giveNode(inode)->giveCoordinates() );
        node->setBoundaryFlag( mesh->giveNode(inode)->isBoundary() );
        ( * dNew )->setDofManager(inode, node);
    } // end creating dof managers

    // create elements
    // count number of local terminal elements first
    int nterminals = 0;
    nelems = mesh->giveNumberOfElements();
    for ( int ielem = 1; ielem <= nelems; ielem++ ) {
#ifdef __PARALLEL_MODE
        if ( mesh->giveElement(ielem)->giveParallelMode() != Element_local ) {
            continue;
        }

#endif
        if ( mesh->giveElement(ielem)->isTerminal() ) {
            nterminals++;
        }
    }

#ifdef __PARALLEL_MODE
    IntArray parentElemMap(nterminals);
#endif
    ( * dNew )->resizeElements(nterminals);
    int eNum = 0;
    for ( int ielem = 1; ielem <= nelems; ielem++ ) {
#ifdef __PARALLEL_MODE
        if ( mesh->giveElement(ielem)->giveParallelMode() != Element_local ) {
            continue;
        }

#endif
        if ( !mesh->giveElement(ielem)->isTerminal() ) {
            continue;
        }

        eNum++;
        parent = mesh->giveElement(ielem)->giveTopParent();
#ifdef __PARALLEL_MODE
        parentElemMap.at(eNum) = parent;
#endif
        if ( parent ) {
            // Copy most of the existing parent element:
            DynamicInputRecord ir( *domain->giveElement ( parent ) );
            ir.setField(* mesh->giveElement(ielem)->giveNodes(), _IFT_Element_nodes);
            ir.giveRecordKeywordField(name);
            elem = classFactory.createElement(name.c_str(), eNum, * dNew);
            elem->initializeFrom(& ir);
            elem->setGlobalNumber( mesh->giveElement(ielem)->giveGlobalNumber() );
#ifdef __PARALLEL_MODE
            //ir.setRecordKeywordNumber( mesh->giveElement(ielem)->giveGlobalNumber() );
            // not subdivided elements inherit globNum, subdivided give -1
            // local elements have array partitions empty !
#endif
            ( * dNew )->setElement(eNum, elem);
        } else {
            OOFEM_ERROR("parent element missing");
        }
    } // end loop over elements

    // create the rest of the model description (BCs, CrossSections, Materials, etc)
    // cross sections
    int ncrosssect = domain->giveNumberOfCrossSectionModels();
    ( * dNew )->resizeCrossSectionModels(ncrosssect);
    for ( int i = 1; i <= ncrosssect; i++ ) {
        DynamicInputRecord ir( *domain->giveCrossSection ( i ) );
        ir.giveRecordKeywordField(name);

        crossSection = classFactory.createCrossSection(name.c_str(), i, * dNew);
        crossSection->initializeFrom(& ir);
        ( * dNew )->setCrossSection(i, crossSection);
    }

    // materials
    int nmat = domain->giveNumberOfMaterialModels();
    ( * dNew )->resizeMaterials(nmat);
    for ( int i = 1; i <= nmat; i++ ) {
        DynamicInputRecord ir( *domain->giveMaterial ( i ) );
        ir.giveRecordKeywordField(name);

        mat = classFactory.createMaterial(name.c_str(), i, * dNew);
        mat->initializeFrom(& ir);
        ( * dNew )->setMaterial(i, mat);
    }

    // barriers
    int nbarriers = domain->giveNumberOfNonlocalBarriers();
    ( * dNew )->resizeNonlocalBarriers(nbarriers);
    for ( int i = 1; i <= nbarriers; i++ ) {
        DynamicInputRecord ir( *domain->giveNonlocalBarrier ( i ) );
        ir.giveRecordKeywordField(name);
        
        barrier = classFactory.createNonlocalBarrier(name.c_str(), i, * dNew);
        barrier->initializeFrom(& ir);
        ( * dNew )->setNonlocalBarrier(i, barrier);
    }

    // boundary conditions
    int nbc = domain->giveNumberOfBoundaryConditions();
    ( * dNew )->resizeBoundaryConditions(nbc);
    for ( int i = 1; i <= nbc; i++ ) {
        DynamicInputRecord ir( *domain->giveBc ( i ) );
        ir.giveRecordKeywordField(name);


        bc = classFactory.createBoundaryCondition(name.c_str(), i, * dNew);
        bc->initializeFrom(& ir);
        ( * dNew )->setBoundaryCondition(i, bc);
    }

    // initial conditions
    int nic = domain->giveNumberOfInitialConditions();
    ( * dNew )->resizeInitialConditions(nic);
    for ( int i = 1; i <= nic; i++ ) {
        DynamicInputRecord ir( *domain->giveIc ( i ) );

        ic = new InitialCondition(i, *dNew);
        ic->initializeFrom(& ir);
        ( * dNew )->setInitialCondition(i, ic);
    }

    // load time functions
    int nfunc = domain->giveNumberOfFunctions();
    ( * dNew )->resizeFunctions(nfunc);
    for ( int i = 1; i <= nfunc; i++ ) {
        DynamicInputRecord ir( *domain->giveFunction ( i ) );
        ir.giveRecordKeywordField(name);

        func = classFactory.createFunction(name.c_str(), i, * dNew);
        func->initializeFrom(& ir);
        ( * dNew )->setFunction(i, func);
    }

    // sets
    int nset = domain->giveNumberOfSets();
    ( * dNew )->resizeSets(nset);
    for ( int i = 1; i <= nset; i++ ) {
        DynamicInputRecord ir( *domain->giveSet ( i ) );
        ir.giveRecordKeywordField(name);

        set = new Set(i, * dNew);
        set->initializeFrom(& ir);
        ( * dNew )->setSet(i, set);
    }
    

    // post initialize components
    ( * dNew )->postInitialize();

    // copy output manager settings
    ( * dNew )->giveOutputManager()->beCopyOf( domain->giveOutputManager() );

    timer.stopTimer();
#ifdef __PARALLEL_MODE
    OOFEM_LOG_INFO( "[%d] Subdivision: created new mesh (%d nodes and %d elements) in %.2fs\n",
                   ( * dNew )->giveEngngModel()->giveRank(), nnodes, eNum, timer.getUtime() );
#else
    OOFEM_LOG_INFO( "Subdivision: created new mesh (%d nodes and %d elements) in %.2fs\n",
                   nnodes, eNum, timer.getUtime() );
#endif

 


#ifdef __PARALLEL_MODE
 #ifdef __VERBOSE_PARALLEL
    nnodes = ( * dNew )->giveNumberOfDofManagers();
    for ( int inode = 1; inode <= nnodes; inode++ ) {
        if ( ( * dNew )->giveDofManager(inode)->giveParallelMode() == DofManager_shared ) {
            //OOFEM_LOG_INFO ("[%d] Shared Node %d[%d]\n", this->giveRank(), inode, (*dNew)->giveDofManager(inode)->giveGlobalNumber());
        }
    }

 #endif

    // we need to assign global numbers to newly generated elements
    this->assignGlobalNumbersToElements(* dNew);

    bool nonloc = false;
    nmat = ( * dNew )->giveNumberOfMaterialModels();
    for ( int im = 1; im <= nmat; im++ ) {
        if ( ( * dNew )->giveMaterial(im)->giveInterface(NonlocalMaterialExtensionInterfaceType) ) {
            nonloc = true;
        }
    }

    if ( nonloc ) {
        exchangeRemoteElements(* dNew, parentElemMap);
    }

    ( * dNew )->commitTransactions( ( * dNew )->giveTransactionManager() );

    // print some statistics
    if (0) {
      for (int in=1; in<=(*dNew)->giveNumberOfDofManagers(); in++) {
        DynamicInputRecord ir;
        (*dNew)->giveDofManager(in)->giveInputRecord(ir);
        OOFEM_LOG_INFO("[%d]:[%d]:%s\n", this->giveRank(), (*dNew)->giveDofManager(in)->giveGlobalNumber(), ir.giveRecordAsString().c_str());
      }
      IntArray nodes;
      for (int in=1; in<=(*dNew)->giveNumberOfElements(); in++) {
        DynamicInputRecord ir;
        (*dNew)->giveElement(in)->giveInputRecord(ir);
        nodes = (*dNew)->giveElement(in)->giveDofManArray();
        // translate local node numbers to globnums
        for (int ii=1; ii<=nodes.giveSize(); ii++)
          nodes.at(ii) = (*dNew)->giveNode(nodes.at(ii))->giveGlobalNumber();
        ir.setField(nodes, "gnodes");
        OOFEM_LOG_INFO("[%d]:[%d]:%s\n", this->giveRank(), (*dNew)->giveElement(in)->giveGlobalNumber(), ir.giveRecordAsString().c_str());
      }
    }
    nelems = ( * dNew )->giveNumberOfElements();
    int localVals [ 2 ], globalVals [ 2 ];
    localVals [ 0 ] = 0;
    for ( int ielem = 1; ielem <= nelems; ielem++ ) {
        if ( ( * dNew )->giveElement(ielem)->giveParallelMode() == Element_local ) {
            localVals [ 0 ]++;
        }
    }

    nnodes = ( * dNew )->giveNumberOfDofManagers();
    localVals [ 1 ] = 0;
    for ( int inode = 1; inode <= nnodes; inode++ ) {
        if ( ( * dNew )->giveDofManager(inode)->isLocal() ) {
            localVals [ 1 ]++;
        }
    }

    MPI_Reduce(localVals, globalVals, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if ( this->giveRank() == 0 ) {
        OOFEM_LOG_INFO("Subdivision: new mesh info: %d nodes, %d elements in total\n", globalVals [ 1 ], globalVals [ 0 ]);
    }
#else
    // we need to assign global numbers to newly generated elements
    this->assignGlobalNumbersToElements(* dNew);
#endif

    return MI_OK;
}


void
Subdivision :: bisectMesh()
{
    int ie, nelems = mesh->giveNumberOfElements(), nelems_old = 0, terminal_local_elems = nelems;
    int nnodes = mesh->giveNumberOfNodes(), nnodes_old;
    double iedensity, rdensity;
    int repeat = 1, loop = 0, max_loop = 0;     // max_loop != 0 use only for debugging
    RS_Element *elem;
    RS_Node *node;
    //std::queue<int> subdivqueue;
#ifdef __PARALLEL_MODE
    int remote_elems = 0;
    int myrank = this->giveRank();
    int problem_size = this->giveNumberOfProcesses();
    int value;
    int *partitionsIrregulars = new int[ problem_size ];
#endif
    
    // get the max globnum on the mesh
    // determine max global number of local nodes
    int maxlocalglobal = 0, maxglobalnumber;
    for ( int in = 1; in <= nnodes; in++ ) {
        maxlocalglobal = max( maxlocalglobal, mesh->giveNode(in)->giveGlobalNumber() );
    }
#ifdef __PARALLEL_MODE
    // determine max global number on all partitions
    MPI_Allreduce(& maxlocalglobal, & maxglobalnumber, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
    maxglobalnumber=maxlocalglobal;
#endif


    // repeat bisection until no new element is created
    while ( repeat && ( loop < max_loop || max_loop == 0 ) ) {


        nnodes_old = nnodes;
#ifdef __PARALLEL_MODE
        OOFEM_LOG_INFO("[%d] Subdivision::bisectMesh: entering bisection loop %d\n", myrank, ++loop);
#else
        OOFEM_LOG_INFO("Subdivision::bisectMesh: entering bisection loop %d\n", ++loop);
#endif
        repeat = 0;
        // process only newly created elements in pass 2 and more
        for ( ie = nelems_old + 1; ie <= nelems; ie++ ) {
            elem = mesh->giveElement(ie);
#ifdef __PARALLEL_MODE
            // skip bisection of remote elements (in first pass only (nelems_old = 0));
            // in pass 2 and more there should be no remote elements because
            // only newly created elements are processed
            if ( nelems_old == 0 ) {
                if ( elem->giveParallelMode() == Element_remote ) {
                    remote_elems++;
                    terminal_local_elems--;
                    continue;
                }
            }

#endif
            if ( !elem->isTerminal() ) {
                continue;
            }

#ifdef __PARALLEL_MODE
            // bisect local elements only
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

#endif

#ifdef DEBUG_CHECK
            if ( elem->giveQueueFlag() ) {
                OOFEM_ERROR("unexpected queue flag of %d ", ie);
            }

#endif

            iedensity = elem->giveDensity();
            rdensity = elem->giveRequiredDensity();

            // first select all candidates for local bisection based on required mesh density

            if ( rdensity < iedensity ) {
                subdivqueue.push(ie);
                elem->setQueueFlag(true);



#ifndef __PARALLEL_MODE
                repeat = 1;                // force repetition in seqeuntial run
#endif
                /*
                 * #ifdef __PARALLEL_MODE
                 * OOFEM_LOG_INFO("[%d] Subdivision: scheduling element %d[%d] for bisection, dens=%lf rdens=%lf\n", myrank, ie, elem->giveGlobalNumber(), iedensity, rdensity);
                 ******#else
                 * OOFEM_LOG_INFO("Subdivision: scheduling element %d for bisection, dens=%lf rdens=%lf\n", ie, iedensity, rdensity);
                 ******#endif
                 */
            }
        }

#ifdef DEBUG_INFO
 #ifdef __PARALLEL_MODE
        OOFEM_LOG_INFO("[%d] (with %d nodes and %d elems)\n", myrank, nnodes_old, terminal_local_elems + remote_elems);
 #else
        OOFEM_LOG_INFO("(with %d nodes and %d elems)\n", nnodes_old, terminal_local_elems);
 #endif
#endif

#ifdef __PARALLEL_MODE
        for ( value = 0; value == 0; value = exchangeSharedIrregulars() ) {
#endif
        // loop over subdivision queue to bisect all local elements there
        while ( !subdivqueue.empty() ) {
            elem = mesh->giveElement( subdivqueue.front() );
#ifdef DEBUG_CHECK
 #ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                OOFEM_ERROR( "nonlocal element %d not expected for bisection", elem->giveNumber() );
            }

 #endif
#endif
            elem->evaluateLongestEdge();
            elem->bisect(subdivqueue, sharedIrregularsQueue);
            subdivqueue.pop();
        }

#ifdef __PARALLEL_MODE
        // in parallel communicate with neighbours the irregular nodes on shared bondary
        }

#endif

        int in;
        nnodes = mesh->giveNumberOfNodes();

#ifndef __PARALLEL_MODE
        int myrank = 0;
#endif
        // assign global numbers to newly introduced irregulars while
        // keeping global numbering of existing (master) nodes
        // idea: first determine the max globnum already assigned
        // and start global numbering of new nodes from this value up

        // determine number of local irregulars
        // this is needed to determine offsets on each partition to start global numbering
        // the numbering of irregulars starts on rank 0 partition by numbering subseqeuntly all its irregulars
        // then we proceed with rank 1, etc.
        // the shared irregulars receive their number from partition with the lowest rank.

        // count local irregulars that receive their global number from this partition
        int localIrregulars = 0;
        for ( in = nnodes_old+1; in <= nnodes; in++ ) {
            if ( this->isNodeLocalIrregular(mesh->giveNode(in), myrank) && ( mesh->giveNode(in)->giveGlobalNumber() == 0 ) ) {
                localIrregulars++;
            }
        }

#ifdef __PARALLEL_MODE
 #ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_INFO("[%d] Subdivision::bisectMesh: number of new local irregulars is %d\n", myrank, localIrregulars);
 #endif
        int localOffset = 0;
        int irank, gnum,  globalIrregulars = 0;
        // gather number of local irregulars from all partitions
        MPI_Allgather(& localIrregulars, 1, MPI_INT, partitionsIrregulars, 1, MPI_INT, MPI_COMM_WORLD);
        // compute local offset
        for ( irank = 0; irank < myrank; irank++ ) {
            localOffset += partitionsIrregulars [ irank ];
        }
        // start to assign global numbers to local irregulars
        int availGlobNum = maxglobalnumber + localOffset;
        for ( in = nnodes_old+1; in <= nnodes; in++ ) {
            node = mesh->giveNode(in);
            if ( this->isNodeLocalIrregular(node, myrank) && ( node->giveGlobalNumber() == 0 ) ) {
                // set negative globnum to mark newly assigned nodes with globnum to participate in shared globnum data exchange
                node->setGlobalNumber( -( ++availGlobNum ) );
 #ifdef __VERBOSE_PARALLEL
                //OOFEM_LOG_INFO ("[%d] %d --> [%d]\n", myrank, in, availGlobNum);
 #endif
            }
        }
#else
        int localOffset=0;
        // start to assign global numbers to local irregulars
        int availGlobNum = maxglobalnumber+localOffset;
        for ( in = nnodes_old+1; in <= nnodes; in++ ) {
            node = mesh->giveNode(in);
            node->setGlobalNumber( ( ++availGlobNum ) );
        }

#endif

#ifdef __PARALLEL_MODE
        // finally, communicate global numbers assigned to shared irregulars
        this->assignGlobalNumbersToSharedIrregulars();
        for ( in = nnodes_old+1; in <= nnodes; in++ ) {
            node = mesh->giveNode(in);
            gnum = node->giveGlobalNumber();
            if ( gnum < 0 ) {
                // turn all globnums to positive
                node->setGlobalNumber(-gnum);
                // update global shared node map (for refined level)
                // this must be done after globnum is made positive
                this->mesh->insertGlobalSharedNodeMap(node);
            } else if ( gnum == 0 ) {
                OOFEM_ERROR("rank [%d] zero globnum identified on node %d", myrank, in);
            }
        }

        // update max globnum
        for ( irank = 0; irank < problem_size; irank++ ) {
            globalIrregulars += partitionsIrregulars [ irank ];
        }

        if ( globalIrregulars ) {
            maxglobalnumber += globalIrregulars;
            // exchange shared edges
            // this must be done after globnums are assigned to new shared irregulars
            /* BP
            this->exchangeSharedEdges();
            repeat = 1;                                    // force repetition in parallel run
            */
        }

#else
        maxglobalnumber+=localIrregulars;
#endif

        // symbolic bisection is finished;
        // now we need to create new mesh;
        // also there is a need to update element connectivities
        for ( ie = 1; ie <= nelems; ie++ ) {
            elem = mesh->giveElement(ie);
            if ( !elem->isTerminal() ) {
                continue;
            }

#ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            elem->generate(sharedEdgesQueue);
        }

#ifdef __PARALLEL_MODE
        if ( globalIrregulars ) {
            // exchange shared edges
            // this must be done after globnums are assigned to new shared irregulars
            this->exchangeSharedEdges();
            repeat = 1;                                    // force repetition in parallel run
        }
#endif
        // unmark local unshared irregulars marked for connectivity setup
        // important this may be not done before generate !!!
        for ( in = nnodes_old+1; in <= nnodes; in++ ) {
            if ( mesh->giveNode(in)->giveNumber() < 0 ) {
                mesh->giveNode(in)->setNumber( -mesh->giveNode(in)->giveNumber() );
            }
        }

        nelems_old = nelems;
        nelems = mesh->giveNumberOfElements();
        terminal_local_elems = 0;
        for ( ie = 1; ie <= nelems; ie++ ) {
            elem = mesh->giveElement(ie);
            if ( !elem->isTerminal() ) {
                continue;
            }

#ifdef __PARALLEL_MODE
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            elem->update_neighbours();
            terminal_local_elems++;
        }

#if 0
 #ifdef __PARALLEL_MODE
        int global_repeat = 0;

        // determine whether any partition needs additional bisection pass
        MPI_Allreduce(& repeat, & global_repeat, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        repeat = global_repeat;
 #endif
#endif

#ifdef __OOFEG
 #ifdef DRAW_MESH_AFTER_EACH_BISECTION_LEVEL
        for ( ie = 1; ie <= nelems; ie++ ) {
            if ( !mesh->giveElement(ie)->isTerminal() ) {
                continue;
            }

            mesh->giveElement(ie)->drawGeometry();
        }

        ESIEventLoop( YES, const_cast< char * >("Subdivision Bisection; Press Ctrl-p to continue") );
 #endif
#endif
    } // end bisection loop
#ifdef __PARALLEL_MODE
    if (partitionsIrregulars) {
      delete[] partitionsIrregulars;
    }
#endif
}


void
Subdivision :: smoothMesh()
{
    int nnodes, nelems, i, j, in, ie;
    int pos, number, reg, nd, cycles = 6;
    IntArray snodes;
    FloatArray *coords;
    RS_Element *elem;
    //bool fixed;
    IntArray node_num_elems, node_con_elems;
    IntArray node_num_nodes, node_con_nodes;

#ifdef DEBUG_SMOOTHING
    char buffer [ 1024 ];
    int len;
#endif
#ifdef __PARALLEL_MODE
    OOFEM_LOG_INFO( "[%d] Subdivision::smoothMesh\n", this->giveRank() );
#else
    OOFEM_LOG_INFO("Subdivision::smoothMesh\n");
#endif
    nnodes = mesh->giveNumberOfNodes();
    nelems = mesh->giveNumberOfElements();

    // count number of elements incident to nodes
    node_num_elems.resize(nnodes + 1);
    node_num_elems.zero();

    for ( ie = 1; ie <= nelems; ie++ ) {
        elem = mesh->giveElement(ie);
        if ( !elem->isTerminal() ) {
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( elem->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        for ( i = 1; i <= elem->giveNodes()->giveSize(); i++ ) {
            node_num_elems.at( elem->giveNode(i) )++;
        }
    }

    // recalculate the number of elements to current addresses
    pos = 1;
    for ( in = 1; in <= nnodes; in++ ) {
        number = node_num_elems.at(in);
        node_num_elems.at(in) = pos;
        pos += number;
    }

    node_num_elems.at(nnodes + 1) = pos;
    node_con_elems.resize(pos);

    // store element numbers incident to nodes
    for ( ie = 1; ie <= nelems; ie++ ) {
        elem = mesh->giveElement(ie);
        if ( !elem->isTerminal() ) {
            continue;
        }

#ifdef __PARALLEL_MODE
        if ( elem->giveParallelMode() != Element_local ) {
            continue;
        }

#endif

        for ( i = 1; i <= elem->giveNodes()->giveSize(); i++ ) {
            node_con_elems.at(node_num_elems.at( elem->giveNode(i) )++) = ie;
        }
    }

    // recalculate the addresses to address of the first element
    pos = 1;
    for ( in = 1; in <= nnodes; in++ ) {
        number = node_num_elems.at(in) - pos;
        node_num_elems.at(in) = pos;
        pos += number;
    }

#ifdef DEBUG_SMOOTHING
 #ifdef __PARALLEL_MODE
    OOFEM_LOG_INFO( "[%d] Subdivision::smoothMesh: connectivity node_id: elem_ids\n", this->giveRank() );
 #else
    OOFEM_LOG_INFO("Subdivision::smoothMesh: connectivity node_id: elem_ids\n");
 #endif
    for ( in = 1; in <= nnodes; in++ ) {
        len = 0;
        for ( i = node_num_elems.at(in); i < node_num_elems.at(in + 1); i++ ) {
            sprintf( buffer + len, " %d", node_con_elems.at(i) );
            len = strlen(buffer);
            if ( len >= 1024 ) {
                OOFEM_ERROR("too small buffer");
            }
        }

        OOFEM_LOG_INFO("%d:%s\n", in, buffer);
    }

#endif

    // count number of nodes incident to nodes
    node_num_nodes.resize(nnodes + 1);
    node_num_nodes.zero();

    for ( in = 1; in <= nnodes; in++ ) {
        for ( i = node_num_elems.at(in); i < node_num_elems.at(in + 1); i++ ) {
            elem = mesh->giveElement( node_con_elems.at(i) );
            for ( j = 1; j <= elem->giveNodes()->giveSize(); j++ ) {
                nd = elem->giveNode(j);
                if ( nd != in && mesh->giveNode(nd)->giveNumber() > 0 ) {
                    mesh->giveNode(nd)->setNumber(-nd);                          // mark connected node
                    node_num_nodes.at(in)++;
                }
            }
        }

        // unmark connected nodes
        for ( i = node_num_elems.at(in); i < node_num_elems.at(in + 1); i++ ) {
            elem = mesh->giveElement( node_con_elems.at(i) );
            for ( j = 1; j <= elem->giveNodes()->giveSize(); j++ ) {
                nd = elem->giveNode(j);
                mesh->giveNode(nd)->setNumber(nd);
            }
        }
    }

    // recalculate the number of nodes to current addresses
    pos = 1;
    for ( in = 1; in <= nnodes; in++ ) {
        number = node_num_nodes.at(in);
        node_num_nodes.at(in) = pos;
        pos += number;
    }

    node_num_nodes.at(nnodes + 1) = pos;
    node_con_nodes.resize(pos);

    // store nodes incident to nodes
    for ( in = 1; in <= nnodes; in++ ) {
        for ( i = node_num_elems.at(in); i < node_num_elems.at(in + 1); i++ ) {
            elem = mesh->giveElement( node_con_elems.at(i) );
            for ( j = 1; j <= elem->giveNodes()->giveSize(); j++ ) {
                nd = elem->giveNode(j);
                if ( nd != in && mesh->giveNode(nd)->giveNumber() > 0 ) {
                    mesh->giveNode(nd)->setNumber(-nd);                          // mark connected node
                    node_con_nodes.at(node_num_nodes.at(in)++) = nd;
                }
            }
        }

        // unmark connected nodes
        for ( i = node_num_elems.at(in); i < node_num_elems.at(in + 1); i++ ) {
            elem = mesh->giveElement( node_con_elems.at(i) );
            for ( j = 1; j <= elem->giveNodes()->giveSize(); j++ ) {
                nd = elem->giveNode(j);
                mesh->giveNode(nd)->setNumber(nd);
            }
        }
    }

    // recalculate the addresses to address of the first node
    pos = 1;
    for ( in = 1; in <= nnodes; in++ ) {
        number = node_num_nodes.at(in) - pos;
        node_num_nodes.at(in) = pos;
        pos += number;
    }

#ifdef DEBUG_SMOOTHING
 #ifdef __PARALLEL_MODE
    OOFEM_LOG_INFO( "[%d] Subdivision::smoothMesh: connectivity node_id: node_ids\n", this->giveRank() );
 #else
    OOFEM_LOG_INFO("Subdivision::smoothMesh: connectivity node_id: node_ids\n");
 #endif
    for ( in = 1; in <= nnodes; in++ ) {
        len = 0;
        for ( i = node_num_nodes.at(in); i < node_num_nodes.at(in + 1); i++ ) {
            sprintf( buffer + len, " %d", node_con_nodes.at(i) );
            len = strlen(buffer);
            if ( len >= 1024 ) {
                OOFEM_ERROR("too small buffer");
            }
        }

        OOFEM_LOG_INFO("%d:%s\n", in, buffer);
    }

#endif

    // identify fixed nodes which should not be subjected to smoothing;
    // these are nodes on boundary (geometrical or material);
    // note that some of these node may be already marked as boundary if this flag was available on input
    // or when processed during bisection or when processed (and assigned) during previous smoothing;
    // the last case is applicable only if node's method storeYourself packs the boundary flag;
    bool fixed;
    for ( in = 1; in <= nnodes; in++ ) {
        if ( mesh->giveNode(in)->isBoundary() ) {
            continue;                                                             // skip boundary node
        }

        if ( node_num_elems.at(in) == node_num_elems.at(in + 1) ) {
            continue;                                                             // skip isolated node
        }

        fixed = false;
        elem = mesh->giveElement( node_con_elems.at( node_num_elems.at(in) ) );
        // check for geometrical boundary
        for ( j = 1; j <= elem->giveNeighbors()->giveSize(); j++ ) {
            if ( elem->giveNeighbor(j) == 0 ) {
                elem->giveSideNodes(j, snodes);
                if ( snodes.findFirstIndexOf(in) != 0 ) {
                    mesh->giveNode(in)->setNumber(-in);                                             // mark fixed node on geometrical boundary
                    fixed = true;
                    break;
                }
            }
        }

        if ( !fixed ) {
            reg = domain->giveElement( elem->giveTopParent() )->giveRegionNumber();
            for ( i = node_num_elems.at(in) + 1; i < node_num_elems.at(in + 1); i++ ) {
                elem = mesh->giveElement( node_con_elems.at(i) );
                // check for material boundary
                if ( domain->giveElement( elem->giveTopParent() )->giveRegionNumber() != reg ) {
                    mesh->giveNode(in)->setNumber(-in);                                             //mark fixed node on material boundary
                    fixed = true;
                    break;
                }

                // check for geometrical boundary
                for ( j = 1; j <= elem->giveNeighbors()->giveSize(); j++ ) {
                    if ( elem->giveNeighbor(j) == 0 ) {
                        elem->giveSideNodes(j, snodes);
                        if ( snodes.findFirstIndexOf(in) != 0 ) {
                            mesh->giveNode(in)->setNumber(-in);                                                 //mark fixed node on geometrical boundary
                            fixed = true;
                            break;
                        }
                    }
                }

                if ( fixed ) {
                    break;
                }
            }
        }
    }

#ifdef QUICK_HACK
    int jn;
    IntArray orderedNodes;
    Subdivision :: RS_CompareNodePositions cmp(mesh);
    orderedNodes.resize(nnodes);
    for ( jn = 1; jn <= nnodes; jn++ ) {
        orderedNodes.at(jn) = jn;
    }

    sort(orderedNodes, cmp);
#endif

    while ( cycles-- ) {
#ifdef QUICK_HACK
        for ( jn = 1; jn <= nnodes; jn++ ) {
            in = orderedNodes.at(jn);
#else
        for ( in = 1; in <= nnodes; in++ ) {
#endif
#ifdef __PARALLEL_MODE
            if ( ( mesh->giveNode(in)->giveParallelMode() == DofManager_shared ) ||
                ( mesh->giveNode(in)->giveParallelMode() == DofManager_null ) ) {
                continue;                                                                                 // skip shared and remote node
            }

#endif
            if ( mesh->giveNode(in)->giveNumber() < 0 ) {
                continue;                                                                           // skip fixed node
            }

            if ( mesh->giveNode(in)->isBoundary() ) {
                continue;                                                                           // skip boundary node
            }

            coords = mesh->giveNode(in)->giveCoordinates();

#ifdef DEBUG_CHECK
            if ( coords ) {
                int count = 0;
                coords->zero();
                for ( i = node_num_nodes.at(in); i < node_num_nodes.at(in + 1); i++ ) {
                    if ( mesh->giveNode( node_con_nodes.at(i) ) ) {
                        if ( mesh->giveNode( node_con_nodes.at(i) )->giveCoordinates() ) {
                          coords->add( *(mesh->giveNode( node_con_nodes.at(i))->giveCoordinates()));
                            count++;
                        } else {
                            OOFEM_ERROR("node %d without coordinates", in);
                        }
                    } else {
                        OOFEM_ERROR("undefined node %d", in);
                    }
                }

                if ( !count ) {
                    OOFEM_ERROR("node %d without connectivity", in);
                }

                coords->times( 1.0 / ( node_num_nodes.at(in + 1) - node_num_nodes.at(in) ) );
            } else {
                OOFEM_ERROR("node %d without coordinates", in);
            }

#else
            coords->zero();
            for ( i = node_num_nodes.at(in); i < node_num_nodes.at(in + 1); i++ ) {
                coords->add( * mesh->giveNode( node_con_nodes.at(i) )->giveCoordinates() );
            }

            coords->times( 1.0 / ( node_num_nodes.at(in + 1) - node_num_nodes.at(in) ) );
#endif
        }
    }

    // unmark fixed nodes and marked them as boundary
    for ( in = 1; in <= nnodes; in++ ) {
        if ( mesh->giveNode(in)->giveNumber() < 0 ) {
            mesh->giveNode(in)->setNumber(in);

#ifndef QUICK_HACK
            // it is not clear what boundary flag may cause
            // therefore it is not set in current version
            //mesh->giveNode(in)->setBoundary(true);
#endif
        }
    }
}


bool
Subdivision :: isNodeLocalIrregular(Subdivision :: RS_Node *node, int myrank)
{
#ifdef __PARALLEL_MODE
    if ( node->isIrregular() ) {
        if ( node->giveParallelMode() == DofManager_local ) {
            return true;
        } else if ( node->giveParallelMode() == DofManager_shared ) {
            int i, minpart, npart;
            const IntArray *partitions = node->givePartitions();
            npart = partitions->giveSize();
            minpart = myrank;
            for ( i = 1; i <= npart; i++ ) {
                minpart = min( minpart, partitions->at(i) );
            }

            if ( minpart == myrank ) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
#else
    return node->isIrregular();
#endif
}

void
Subdivision :: assignGlobalNumbersToElements(Domain *d)
{

#ifdef __PARALLEL_MODE

    int problem_size = this->giveNumberOfProcesses();
    int myrank = this->giveRank();
    int i, nelems, numberOfLocalElementsToNumber = 0;
    int *partitionNumberOfElements = new int[ problem_size ];
    int localMaxGlobnum = 0, globalMaxGlobnum;

    // idea: first determine the number of local elements waiting for new global id
    // and also determine max global number assigned up to now
    nelems = d->giveNumberOfElements();
    for ( i = 1; i <= nelems; i++ ) {
        localMaxGlobnum = max( localMaxGlobnum, d->giveElement(i)->giveGlobalNumber() );
 #ifdef DEBUG_CHECK
        if ( d->giveElement(i)->giveParallelMode() == Element_remote ) {
            OOFEM_ERROR("unexpected remote element %d ", i);
        }

 #endif
        if ( d->giveElement(i)->giveGlobalNumber() <= 0 ) {
            numberOfLocalElementsToNumber++;
        }
    }

    // determine number of elements across all partitions
    MPI_Allgather(& numberOfLocalElementsToNumber, 1, MPI_INT,
                  partitionNumberOfElements, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allreduce(& localMaxGlobnum, & globalMaxGlobnum, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
 #ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO("[%d] Subdivision::assignGlobalNumbersToElements: max globnum %d, new elements %d\n", myrank, globalMaxGlobnum, numberOfLocalElementsToNumber);
 #endif

    // compute local offset
    int startOffset = globalMaxGlobnum, availGlobNum;
    for ( i = 0; i < myrank; i++ ) {
        startOffset += partitionNumberOfElements [ i ];
    }

    // lets assign global numbers on each partition to local elements
    availGlobNum = startOffset;
    for ( i = 1; i <= nelems; i++ ) {
        if ( d->giveElement(i)->giveGlobalNumber() <= 0 ) {
            d->giveElement(i)->setGlobalNumber(++availGlobNum);
        }
    }

 #ifdef __VERBOSE_PARALLEL
    /*
     * for (i=1; i<=nelems; i++) {
     * OOFEM_LOG_INFO ("[%d] Element %d[%d]\n", myrank, i,d->giveElement(i)->giveGlobalNumber());
     * }
     */
 #endif

    if (partitionNumberOfElements) {
      delete[] partitionNumberOfElements;
    }

#else // local case

    int i, nelems, numberOfLocalElementsToNumber = 0;
    int localMaxGlobnum = 0;

    // idea: first determine the number of local elements waiting for new global id
    // and also determine max global number assigned up to now
    nelems = d->giveNumberOfElements();
    for ( i = 1; i <= nelems; i++ ) {
        localMaxGlobnum = max( localMaxGlobnum, d->giveElement(i)->giveGlobalNumber() );
        if ( d->giveElement(i)->giveGlobalNumber() <= 0 ) {
            numberOfLocalElementsToNumber++;
        }
    }

    // compute local offset
    int startOffset = localMaxGlobnum, availGlobNum;

    // lets assign global numbers on each partition to local elements
    availGlobNum = startOffset;
    for ( i = 1; i <= nelems; i++ ) {
        if ( d->giveElement(i)->giveGlobalNumber() <= 0 ) {
            d->giveElement(i)->setGlobalNumber(++availGlobNum);
        }
    }

#endif

}



#ifdef __PARALLEL_MODE
bool
Subdivision :: exchangeSharedIrregulars()
{
    // loop over local sharedIrregularsQueue
    int globalSharedIrregularsQueueEmpty, localSharedIrregularsQueueEmpty = this->sharedIrregularsQueue.empty();
 #ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO("[%d] Subdivision :: exchangeSharedIrregulars: localSharedIrregularsQueueEmpty %d\n",
                   this->giveRank(), localSharedIrregularsQueueEmpty);
 #endif
    MPI_Allreduce(& localSharedIrregularsQueueEmpty, & globalSharedIrregularsQueueEmpty, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
 #ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO("[%d] Subdivision :: exchangeSharedIrregulars: globalSharedIrregularsQueueEmpty %d\n",
                   this->giveRank(), globalSharedIrregularsQueueEmpty);
 #endif
    if ( globalSharedIrregularsQueueEmpty ) {
        return true;
    } else {
 #ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_INFO( "[%d] Subdivision :: exchangeSharedIrregulars: started\n", this->giveRank() );
 #endif

        // there are some shared irregulars -> data exchange
        CommunicatorBuff cb(this->giveNumberOfProcesses(), CBT_dynamic);
        Communicator com(domain->giveEngngModel(), &cb, this->giveRank(), this->giveNumberOfProcesses(), CommMode_Dynamic);

        com.packAllData(this, this, & Subdivision :: packSharedIrregulars);
        com.initExchange(SHARED_IRREGULAR_DATA_TAG);
        com.unpackAllData(this, this, & Subdivision :: unpackSharedIrregulars);
        com.finishExchange();
        this->sharedIrregularsQueue.clear();

        return false;
    }
}


int
Subdivision :: packSharedIrregulars(Subdivision *s, ProcessCommunicator &pc)
{
    int pi, iNode, jNode, rproc = pc.giveRank();
    int myrank = this->giveRank();
    const IntArray *sharedPartitions;
    std :: list< int > :: const_iterator sharedIrregQueueIter;
    IntArray edgeInfo(2);

    if ( rproc == myrank ) {
        return 1;                  // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
    for ( sharedIrregQueueIter = sharedIrregularsQueue.begin();
          sharedIrregQueueIter != sharedIrregularsQueue.end();
          sharedIrregQueueIter++ ) {
        pi = ( * sharedIrregQueueIter );
        sharedPartitions = this->mesh->giveNode(pi)->givePartitions();
        if ( sharedPartitions->contains(rproc) ) {
            // the info about new local shared irregular node needs to be sent to remote partition
            // the new irregular on remote partition is identified using two nodes (glonums) defining
            // an edge on which irregular node is introduced
            ( ( RS_IrregularNode * ) this->mesh->giveNode(pi) )->giveEdgeNodes(iNode, jNode);
            edgeInfo.at(1) = this->mesh->giveNode(iNode)->giveGlobalNumber();
            edgeInfo.at(2) = this->mesh->giveNode(jNode)->giveGlobalNumber();
            pcbuff->write(SUBDIVISION_SHARED_IRREGULAR_REC_TAG);
            edgeInfo.storeYourself(*pcbuff);
 #ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_INFO("[%d] Subdivision::packSharedIrregulars: packing shared node %d [%d %d] for %d\n",
                           myrank, pi, edgeInfo.at(1), edgeInfo.at(2), rproc);
 #endif
        }
    }

    pcbuff->write(SUBDIVISION_END_DATA);

    return 1;
}

int
Subdivision :: unpackSharedIrregulars(Subdivision *s, ProcessCommunicator &pc)
{
    int myrank = this->giveRank();
    int ie, _type, iNum, iproc = pc.giveRank();
    int iNode, jNode, elems;
    double density;
    IntArray edgeInfo(2);
    const IntArray *iElems, *jElems;
    FloatArray coords;
    Subdivision :: RS_SharedEdge *edge;
    Subdivision :: RS_Element *elem;
    Subdivision :: RS_IrregularNode *irregular;
    int eIndex;

    if ( iproc == myrank ) {
        return 1;                                 // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();

    pcbuff->read(_type);
    // unpack dofman data
    while ( _type != SUBDIVISION_END_DATA ) {
        if ( _type == SUBDIVISION_SHARED_IRREGULAR_REC_TAG ) {
            edgeInfo.restoreYourself(*pcbuff);
 #ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_INFO("[%d] Subdivision::unpackSharedIrregulars: received shared node record [%d %d] from %d ...\n",
                           myrank, edgeInfo.at(1), edgeInfo.at(2), iproc);
 #endif

            iNode = mesh->sharedNodeGlobal2Local( edgeInfo.at(1) );
            jNode = mesh->sharedNodeGlobal2Local( edgeInfo.at(2) );

            // get elements incident simultaneiusly to iNode and jNode
            iElems = mesh->giveNode(iNode)->giveConnectedElements();
            jElems = mesh->giveNode(jNode)->giveConnectedElements();

            IntArray common;
            if ( iElems->giveSize() <= jElems->giveSize() ) {
                common.preallocate( iElems->giveSize() );
            } else {
                common.preallocate( jElems->giveSize() );
            }

            // I do rely on the fact that the arrays are ordered !!!
            // I am using zero chunk because array common is large enough
            elems = iElems->findCommonValuesSorted(* jElems, common, 0);
#if 0
            if ( !elems) {
               // get type of the next record
               pcbuff->read(_type);
               continue;
            }            
#else
            if ( !elems ) {
                OOFEM_ERROR("[%d] no element found sharing nodes %d[%d] and %d[%d]",
                            myrank, iNode, edgeInfo.at(1), jNode, edgeInfo.at(2));
            }
#endif
            // check on the first element whether irregular exists
            elem = mesh->giveElement( common.at(1) );
            eIndex = elem->giveEdgeIndex(iNode, jNode);
 #ifdef DEBUG_CHECK
            if ( !elem->giveSharedEdges()->giveSize() ) {
                OOFEM_ERROR("element %d sharing nodes %d %d has no shared edges",
                             elem->giveNumber(), iNode, jNode);
            }

            if ( !elem->giveSharedEdge(eIndex) ) {
                OOFEM_ERROR("element %d sharing nodes %d and %d has no shared edge %d",
                             elem->giveNumber(), iNode, jNode, eIndex);
            }

 #endif
            if ( elem->giveIrregular(eIndex) ) {
                // irregular already exists
 #ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_INFO( "...already exists as %d on element %d\n", elem->giveIrregular(eIndex), elem->giveNumber() );
 #endif
            } else {
                // irregular does not exist:
                // compute coordinates of new irregular
                coords = * ( mesh->giveNode(iNode)->giveCoordinates() );
                coords.add( * mesh->giveNode(jNode)->giveCoordinates() );
                coords.times(0.5);
 #ifdef HEADEDSTUD
                double dist, rad, rate;
                FloatArray *c;

                c = mesh->giveNode(iNode)->giveCoordinates();
                dist = c->at(1) * c->at(1) + c->at(3) * c->at(3);
                if ( c->at(2) > 69.9999999 ) {
                    rad = 7.0;
                } else if ( c->at(2) < 64.5000001 ) {
                    rad = 18.0;
                } else {
                    rad = 18.0 - 11.0 / 5.5 * ( c->at(2) - 64.5 );
                }

                if ( fabs(dist - rad * rad) < 0.01 ) {                // be very tolerant (geometry is not precise)
                    c = mesh->giveNode(jNode)->giveCoordinates();
                    dist = c->at(1) * c->at(1) + c->at(3) * c->at(3);
                    if ( c->at(2) > 69.9999999 ) {
                        rad = 7.0;
                    } else if ( c->at(2) < 64.5000001 ) {
                        rad = 18.0;
                    } else {
                        rad = 18.0 - 11.0 / 5.5 * ( c->at(2) - 64.5 );
                    }

                    if ( fabs(dist - rad * rad) < 0.01 ) {                    // be very tolerant (geometry is not precise)
                        dist = coords.at(1) * coords.at(1) + coords.at(3) * coords.at(3);
                        if ( coords.at(2) > 69.9999999 ) {
                            rad = 7.0;
                        } else if ( coords.at(2) < 64.5000001 ) {
                            rad = 18.0;
                        } else {
                            rad = 18.0 - 11.0 / 5.5 * ( coords.at(2) - 64.5 );
                        }

                        rate = rad / sqrt(dist);
                        coords.at(1) *= rate;
                        coords.at(3) *= rate;
                    }
                }

 #endif
                // compute required density of a new node
                density = 0.5 * ( mesh->giveNode(iNode)->giveRequiredDensity() +
                                 mesh->giveNode(jNode)->giveRequiredDensity() );
                // create new irregular to receiver
                iNum = mesh->giveNumberOfNodes() + 1;
                irregular = new Subdivision :: RS_IrregularNode(iNum, mesh, 0, coords, density, true);
                irregular->setParallelMode(DofManager_shared);
                irregular->setEdgeNodes(iNode, jNode);
                mesh->addNode(irregular);

 #ifdef __OOFEG
  #ifdef DRAW_IRREGULAR_NODES
                irregular->drawGeometry();
  #endif
 #endif

                // partitions are inherited from shared edge
                edge = mesh->giveEdge( elem->giveSharedEdge(eIndex) );
                irregular->setPartitions( * ( edge->givePartitions() ) );
                // add irregular to all relevant elements
                for ( ie = 1; ie <= common.giveSize(); ie++ ) {
                    elem = mesh->giveElement( common.at(ie) );
                    eIndex = elem->giveEdgeIndex(iNode, jNode);
 #ifdef DEBUG_CHECK
                    if ( !elem->giveSharedEdges()->giveSize() ) {
                        OOFEM_ERROR("element %d sharing nodes %d %d has no shared edges",
                                     elem->giveNumber(), iNode, jNode);
                    }

                    if ( !elem->giveSharedEdge(eIndex) ) {
                        OOFEM_ERROR("element %d sharing nodes %d and %d has no shared edge %d",
                                     elem->giveNumber(), iNode, jNode, eIndex);
                    }

                    if ( elem->giveIrregular(eIndex) ) {
                        OOFEM_ERROR("element %d sharing nodes %d %d already has irregular",
                                     elem->giveNumber(), iNode, jNode);
                    }

 #endif
                    elem->setIrregular(eIndex, iNum);
                    if ( !elem->giveQueueFlag() ) {
                        // schedule elem for bisection
                        subdivqueue.push( elem->giveNumber() );
                        elem->setQueueFlag(true);
                    }

 #ifdef __VERBOSE_PARALLEL
                    OOFEM_LOG_INFO( "...added as %d on element %d\n", iNum, elem->giveNumber() );
 #endif
 #ifdef DEBUG_INFO
                    if ( domain->giveElement( elem->giveTopParent() )->giveGeometryType() == EGT_tetra_1 ) {
                        // do not print global numbers of elements because they are not available (they are assigned at once after bisection);
                        // do not print global numbers of irregulars as these may not be available yet
                        OOFEM_LOG_INFO( "[%d] Shared irregular %d added on %d (edge %d, nodes %d %d [%d %d], nds %d %d %d %d [%d %d %d %d], ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
                                       myrank, iNum,
                                       elem->giveNumber(), eIndex, iNode, jNode,
                                       mesh->giveNode(iNode)->giveGlobalNumber(), mesh->giveNode(jNode)->giveGlobalNumber(),
                                       elem->giveNode(1), elem->giveNode(2), elem->giveNode(3), elem->giveNode(4),
                                       mesh->giveNode( elem->giveNode(1) )->giveGlobalNumber(),
                                       mesh->giveNode( elem->giveNode(2) )->giveGlobalNumber(),
                                       mesh->giveNode( elem->giveNode(3) )->giveGlobalNumber(),
                                       mesh->giveNode( elem->giveNode(4) )->giveGlobalNumber(),
                                       elem->giveNeighbor(1), elem->giveNeighbor(2), elem->giveNeighbor(3), elem->giveNeighbor(4),
                                       elem->giveIrregular(1), elem->giveIrregular(2), elem->giveIrregular(3),
                                       elem->giveIrregular(4), elem->giveIrregular(5), elem->giveIrregular(6) );
                    }

 #endif
                }
            }
        } else {
            OOFEM_ERROR("unknown tag received");
        }

        // get type of the next record
        pcbuff->read(_type);
    }

    return 1;
}

void
Subdivision :: assignGlobalNumbersToSharedIrregulars()
{
    // there are some shared irregulars -> data exchange
    CommunicatorBuff cb(this->giveNumberOfProcesses(), CBT_dynamic);
    Communicator com(domain->giveEngngModel(), &cb, this->giveRank(),
                     this->giveNumberOfProcesses(), CommMode_Dynamic);

    com.packAllData(this, this, & Subdivision :: packIrregularSharedGlobnums);
    com.initExchange(SHARED_IRREGULAR_DATA_TAG);
    com.unpackAllData(this, this, & Subdivision :: unpackIrregularSharedGlobnums);
    com.finishExchange();
}

int
Subdivision :: packIrregularSharedGlobnums(Subdivision *s, ProcessCommunicator &pc)
{
    int rproc = pc.giveRank();
    int in, iNode, jNode, nnodes = mesh->giveNumberOfNodes();
    int myrank = this->giveRank();
    IntArray edgeInfo(3);
    const IntArray *sharedPartitions;

    if ( rproc == myrank ) {
        return 1;                  // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
    for ( in = 1; in <= nnodes; in++ ) {
        if ( this->isNodeLocalSharedIrregular(mesh->giveNode(in), myrank) && ( mesh->giveNode(in)->giveGlobalNumber() < 0 ) ) {
            sharedPartitions = this->mesh->giveNode(in)->givePartitions();
            if ( sharedPartitions->contains(rproc) ) {
                // the info about new local shared irregular node needs to be sent to remote partition
                // the new irregular on remote partition is identified using two nodes defining
                // an edge on which irregular node is introduced
                ( ( RS_IrregularNode * ) this->mesh->giveNode(in) )->giveEdgeNodes(iNode, jNode);
                edgeInfo.at(1) = this->mesh->giveNode(iNode)->giveGlobalNumber();
                edgeInfo.at(2) = this->mesh->giveNode(jNode)->giveGlobalNumber();
                edgeInfo.at(3) = this->mesh->giveNode(in)->giveGlobalNumber();         // keep negative
 #ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_INFO("[%d] packIrregularSharedGlobnums: sending %d [%d] - [%d %d] to %d\n",
                               myrank, in, -edgeInfo.at(3), edgeInfo.at(1), edgeInfo.at(2), rproc);
 #endif
                pcbuff->write(SUBDIVISION_SHARED_IRREGULAR_REC_TAG); // KOKO asi zbytecne
                edgeInfo.storeYourself(*pcbuff);
            }
        }
    }

    pcbuff->write(SUBDIVISION_END_DATA);

    return 1;
}


int
Subdivision :: unpackIrregularSharedGlobnums(Subdivision *s, ProcessCommunicator &pc)
{
    int myrank = this->giveRank();
    int _type, iproc = pc.giveRank();
    int iNode, jNode, iNum, elems;
    IntArray edgeInfo(3);
    const IntArray *iElems, *jElems;
    RS_Element *elem;
    RS_Node *node;
    int eIndex;

    if ( iproc == myrank ) {
        return 1;                                 // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();

    pcbuff->read(_type);
    // unpack dofman data
    while ( _type != SUBDIVISION_END_DATA ) {
        if ( _type == SUBDIVISION_SHARED_IRREGULAR_REC_TAG ) {   // KOKO asi zbytecne
            edgeInfo.restoreYourself(*pcbuff);

            iNode = mesh->sharedNodeGlobal2Local( edgeInfo.at(1) );
            jNode = mesh->sharedNodeGlobal2Local( edgeInfo.at(2) );

            // get elements incident simultaneiusly to iNode and jNode
            iElems = mesh->giveNode(iNode)->giveConnectedElements();
            jElems = mesh->giveNode(jNode)->giveConnectedElements();

            IntArray common;
            if ( iElems->giveSize() <= jElems->giveSize() ) {
                common.preallocate( iElems->giveSize() );
            } else {
                common.preallocate( jElems->giveSize() );
            }

            // I do rely on the fact that the arrays are ordered !!!
            // I am using zero chunk because array common is large enough
            elems = iElems->findCommonValuesSorted(* jElems, common, 0);
#if 0
            if ( !elems ) {
               pcbuff->read(_type);
               continue;
            }
#else
            if ( !elems ) {
                OOFEM_ERROR("no element found sharing nodes %d and %d",
                             iNode, jNode);
            }
#endif

            // assign globnum to appropriate edge on the first element
            elem = mesh->giveElement( common.at(1) );
            eIndex = elem->giveEdgeIndex(iNode, jNode);
            iNum = elem->giveIrregular(eIndex);
 #ifdef DEBUG_CHECK
            if ( !elem->giveSharedEdges()->giveSize() ) {
                OOFEM_ERROR("element %d sharing nodes %d %d has no shared edges",
                             elem->giveNumber(), iNode, jNode);
            }

            if ( !elem->giveSharedEdge(eIndex) ) {
                OOFEM_ERROR("element %d sharing nodes %d and %d has no shared edge %d",
                             elem->giveNumber(), iNode, jNode, eIndex);
            }

            if ( !iNum ) {
                OOFEM_ERROR("element %d sharing nodes %d %d does not have irregular",
                             elem->giveNumber(), iNode, jNode);
            }

 #endif
            node = mesh->giveNode(iNum);
            node->setGlobalNumber( edgeInfo.at(3) );                    // keep negative
            // update global shared node map (for refined level)
            this->mesh->insertGlobalSharedNodeMap(node);

 #ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_INFO("[%d] Subdivision::unpackIrregularSharedGlobnums: received %d [%d] - [%d %d] from %d\n",
                           myrank, iNum, -edgeInfo.at(3), edgeInfo.at(1), edgeInfo.at(2), iproc);
 #endif
        } else {
            OOFEM_ERROR("unknown tag received");
        }

        pcbuff->read(_type);
    }

    return 1;
}

bool
Subdivision :: isNodeLocalSharedIrregular(Subdivision :: RS_Node *node, int myrank)
{
    if ( node->isIrregular() ) {
        if ( node->giveParallelMode() == DofManager_shared ) {
            int i, minpart, npart;
            const IntArray *partitions = node->givePartitions();
            npart = partitions->giveSize();
            minpart = myrank;
            for ( i = 1; i <= npart; i++ ) {
                minpart = min( minpart, partitions->at(i) );
            }

            if ( minpart == myrank ) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
}





int
Subdivision :: giveRank()
{
    return this->domain->giveEngngModel()->giveRank();
}

int
Subdivision :: giveNumberOfProcesses()
{
    return this->domain->giveEngngModel()->giveNumberOfProcesses();
}


/*
 * Idea is to rebuild the remote elements from scratch after subdivision
 * Reason is that sbdivision in mot cases followed by smoothing,
 * so that it becomes difficult to reconstruct the remot elements from old remote elements
 *
 * the algorith should be summarized as follows:
 * 1. For each remote partition:
 *    For each node (INODE) shared with remote partition:
 *         Add all elements sharing INODE into the queue:
 *         For each element in the queue:
 *             Mark element as remote
 *             For each node (JNODE) of this element:
 *                Compute its distance from INODE
 *                If less than nonlocal radius:
 *                    Add all elements sharing JNODE to queue
 *
 *  This is an aproximate algrithm, that works correctly for reasonably shaped grids,
 *  for grids with badly shaped elements, some interactions can be neglected, due to the fact
 *  that real distances between intagration points are not considered directly.
 *
 *  This algorithm failes if the mirror zone propagates to partition not adjacent to current partition.
 *  Therefore the algorithm has been abandaned and replaced by that which mirrors children
 *  of all originally mirrored elements;
 *  If the smoothing is applied this is also only approximate algorithm, as the mesh movement
 *  may casuse that some new elements being not child of origanally mirrored elements are shifted
 *  into the mirrored band !!!
 */
void
Subdivision :: exchangeRemoteElements(Domain *d, IntArray &parentMap)
{
    int nproc = this->giveNumberOfProcesses(), myrank = this->giveRank();
    CommunicatorBuff cb(nproc, CBT_dynamic);
    Communicator com(d->giveEngngModel(), &cb, myrank, nproc, CommMode_Dynamic);

    // move existing dofmans and elements, that will be local on current partition,
    // into local map
    RS_packRemoteElemsStruct s;
    s.d = d;
    s.parentElemMap = & parentMap;
    com.packAllData(this, & s, & Subdivision :: packRemoteElements);
    com.initExchange(SUBDIVISION_MIGRATE_REMOTE_ELEMENTS_TAG);

    // remove existing remote elements and null nodes
    int i, nelem = d->giveNumberOfElements();
    int nnodes = d->giveNumberOfDofManagers();
    DomainTransactionManager *dtm = d->giveTransactionManager();
    for ( i = 1; i <= nnodes; i++ ) {
        if ( d->giveDofManager(i)->giveParallelMode() == DofManager_null ) {
            dtm->addDofManTransaction(DomainTransactionManager :: DTT_Remove, d->giveDofManager(i)->giveGlobalNumber(), NULL);
        }
    }

    for ( i = 1; i <= nelem; i++ ) {
        if ( d->giveElement(i)->giveParallelMode() == Element_remote ) {
            dtm->addElementTransaction(DomainTransactionManager :: DTT_Remove, d->giveElement(i)->giveGlobalNumber(), NULL);
        }
    }

    // receive remote data
    com.unpackAllData(this, d, & Subdivision :: unpackRemoteElements);
    com.finishExchange();
}

int
Subdivision :: packRemoteElements(RS_packRemoteElemsStruct *s, ProcessCommunicator &pc)
{
    Domain *d = s->d;
    int rproc = pc.giveRank();
    int i, inode;
    int nn, in;
    int myrank = this->giveRank();
    DofManager *inodePtr;
    Element *elemPtr, *relemPtr;
    std :: queue< int >elemCandidates;
    std :: set< int >remoteElements, processedElements, nodesToSend;

    if ( rproc == myrank ) {
        return 1;                  // skip local partition
    }

    EngngModel *emodel = d->giveEngngModel();
    // comMap refers to original (parent) elements
    const IntArray *comMap = emodel->giveProblemCommunicator(EngngModel :: PC_nonlocal)->giveProcessCommunicator(rproc)->giveToSendMap();
 #ifdef __OOFEG
  #ifdef DRAW_REMOTE_ELEMENTS
    oofegGraphicContext gc;
    EPixel geocolor = gc.getElementColor();
    gc.setElementColor( gc.getActiveCrackColor() );
  #endif
 #endif
    for ( i = 1; i <= d->giveNumberOfElements(); i++ ) {
        // remote parent skipped - parentElemMap has zero value for them
        if ( comMap->contains( s->parentElemMap->at(i) ) ) {
            remoteElements.insert(i);
 #ifdef __OOFEG
  #ifdef DRAW_REMOTE_ELEMENTS
            TimeStep *tStep = emodel->giveCurrentStep();
            d->giveElement(i)->drawRawGeometry(gc, tStep);
  #endif
 #endif
        }
    }

 #ifdef __OOFEG
  #ifdef DRAW_REMOTE_ELEMENTS
    ESIEventLoop(YES, "Remote element packing ; Press Ctrl-p to continue");
    gc.setElementColor(geocolor);
  #endif
 #endif

    // now the list of elements to became remote on given remote partition is in remoteElements set

    // need to pack elements definition and corresponding nodes!
    /*
     * mark shecheduled nodes:
     * loop over elements in remoteElements set and add all their nodes (except those that are shared)
     */
    for ( int elNum: remoteElements ) {
        relemPtr = d->giveElement(elNum);
        nn = relemPtr->giveNumberOfNodes();
        for ( in = 1; in <= nn; in++ ) {
            inode = relemPtr->giveDofManagerNumber(in);
            if ( ( d->giveDofManager(inode)->giveParallelMode() == DofManager_local ) ||
                ( ( d->giveDofManager(inode)->giveParallelMode() == DofManager_shared ) && ( !d->giveDofManager(inode)->givePartitionList()->contains(rproc) ) ) ) {
                // nodesToSend is set, therefore duplicity is avoided
                nodesToSend.insert(inode);
            }
        }
    }

    //-----------------end here-------------------

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
    // send nodes that define remote_elements gometry
    for ( int nodeNum: nodesToSend ) {
        inodePtr = d->giveDofManager(nodeNum);

        pcbuff->write( inodePtr->giveInputRecordName() );
        pcbuff->write( inodePtr->giveGlobalNumber() );
        inodePtr->saveContext(*pcbuff, CM_Definition);
    }

    // pack end-of-element-record
    pcbuff->write("");

    // send elements
    for ( int elNum: remoteElements ) {
        elemPtr = d->giveElement(elNum);
        // pack local element (node numbers shuld be global ones!!!)
        // pack type
        pcbuff->write( elemPtr->giveInputRecordName() );
        // nodal numbers shuld be packed as global !!
        elemPtr->saveContext(*pcbuff, CM_Definition | CM_DefinitionGlobal);
        //OOFEM_LOG_INFO ("[%d] Sending Remote elem %d[%d] to rank %d\n", myrank,*si, elemPtr->giveGlobalNumber(), rproc );
    }

    // pack end-of-element-record
    pcbuff->write("");

    return 1;
}


int
Subdivision :: unpackRemoteElements(Domain *d, ProcessCommunicator &pc)
{
    int myrank = d->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();
    int _globnum;
    std :: string _type;
    DomainTransactionManager *dtm = d->giveTransactionManager();

    if ( iproc == myrank ) {
        return 1;                                 // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();

    // unpack dofman data
    do {
        pcbuff->read(_type);
        if ( _type.size() == 0 ) {
            break;
        }

        // receiving new local dofManager
        pcbuff->read(_globnum);

        DofManager *dofman;
        bool _newentry = false;
        if ( ( dofman = dtm->giveDofManager(_globnum) ) == NULL ) {
            // data not available -> create a new one
            _newentry = true;
            dofman = classFactory.createDofManager(_type.c_str(), 0, d);
        }

        dofman->setGlobalNumber(_globnum);
        // unpack dofman state (this is the local dofman, not available on remote)
        dofman->restoreContext(*pcbuff, CM_Definition);
        dofman->setParallelMode(DofManager_null);
        // add transaction if new entry allocated; otherwise existing one has been modified via returned dofman
        if ( _newentry ) {
            dtm->addDofManTransaction(DomainTransactionManager :: DTT_ADD, _globnum, dofman);
        }
    } while ( 1 );

    IntArray elemPartitions(1);
    elemPartitions.at(1) = iproc;
    int nrecv = 0;

    do {
        pcbuff->read(_type);
        if ( _type.size() == 0 ) {
            break;
        }

        Element *elem = classFactory.createElement(_type.c_str(), 0, d);
        elem->restoreContext(*pcbuff, CM_Definition | CM_DefinitionGlobal);
        elem->setParallelMode(Element_remote);
        elem->setPartitionList(elemPartitions);
        dtm->addElementTransaction(DomainTransactionManager :: DTT_ADD, elem->giveGlobalNumber(), elem);
        nrecv++;
        //OOFEM_LOG_INFO ("[%d] Received Remote elem [%d] from rank %d\n", myrank, elem->giveGlobalNumber(), iproc );
        //recvElemList.push_back(elem);
    } while ( 1 );

    return 1;
}




/* CAUTION: the exchange can be done only for not subdivided edges !
 *       subdivided edges inherit partitions from parent edge (in ::generate) */

void
Subdivision :: exchangeSharedEdges()
{
    int i, pi, iNode, jNode, elems;
    std :: list< int > :: const_iterator sharedEdgeQueueIter;
    Subdivision :: RS_SharedEdge *edge;
    Subdivision :: RS_Element *elem;
    const IntArray *iElems, *jElems;

    // loop over local sharedEdgessQueue
    int globalSharedEdgesQueueEmpty, localSharedEdgesQueueEmpty = this->sharedEdgesQueue.empty();
 #ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO("[%d] Subdivision :: exchangeSharedEdges: localSharedEdgesQueueEmpty %d\n",
                   this->giveRank(), localSharedEdgesQueueEmpty);
 #endif
    MPI_Allreduce(& localSharedEdgesQueueEmpty, & globalSharedEdgesQueueEmpty, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
 #ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO("[%d] Subdivision :: exchangeSharedEdges: globalSharedEdgesQueueEmpty %d\n",
                   this->giveRank(), globalSharedEdgesQueueEmpty);
 #endif
    if ( !globalSharedEdgesQueueEmpty ) {
 #ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_INFO( "[%d] Subdivision :: exchangeSharedEdges: started\n", this->giveRank() );
 #endif

        // there are some shared edges -> data exchange

        CommunicatorBuff cb(this->giveNumberOfProcesses(), CBT_dynamic);
        Communicator com(domain->giveEngngModel(), &cb, this->giveRank(), this->giveNumberOfProcesses(), CommMode_Dynamic);

        com.packAllData(this, this, & Subdivision :: packSharedEdges);

        // remove all tentative partitions on queued shared edges after packing relevant data
        for ( sharedEdgeQueueIter = sharedEdgesQueue.begin();
              sharedEdgeQueueIter != sharedEdgesQueue.end();
              sharedEdgeQueueIter++ ) {
            pi = ( * sharedEdgeQueueIter );

            edge = mesh->giveEdge(pi);
            edge->removePartitions();
        }

        com.initExchange(SHARED_EDGE_DATA_TAG);
        com.unpackAllData(this, this, & Subdivision :: unpackSharedEdges);
        com.finishExchange();

        // unmark unshared edges from elements after unpacking data and before clearing the queue
        for ( sharedEdgeQueueIter = sharedEdgesQueue.begin();
              sharedEdgeQueueIter != sharedEdgesQueue.end();
              sharedEdgeQueueIter++ ) {
            pi = ( * sharedEdgeQueueIter );

            edge = mesh->giveEdge(pi);
            if ( edge->givePartitions()->giveSize() ) {
                continue;
            }

            edge->giveEdgeNodes(iNode, jNode);
 #ifdef __VERBOSE_PARALLEL
            //OOFEM_LOG_INFO("edge %d %d [%d %d] not shared\n", iNode, jNode,
            //mesh->giveNode(iNode)->giveGlobalNumber(), mesh->giveNode(jNode)->giveGlobalNumber());
 #endif
            iElems = mesh->giveNode(iNode)->giveConnectedElements();
            jElems = mesh->giveNode(jNode)->giveConnectedElements();

            IntArray common;
            if ( iElems->giveSize() <= jElems->giveSize() ) {
                common.preallocate( iElems->giveSize() );
            } else {
                common.preallocate( jElems->giveSize() );
            }

            // I do rely on the fact that the arrays are ordered !!!
            // I am using zero chunk because array common is large enough
            elems = iElems->findCommonValuesSorted(* jElems, common, 0);
 #ifdef DEBUG_CHECK
            if ( !elems ) {
                OOFEM_ERROR("no element found sharing nodes %d and %d corresponding to edge %d",
                             iNode, jNode, pi);
            }

 #endif
            for ( i = 1; i <= elems; i++ ) {
                elem = mesh->giveElement( common.at(i) );
 #ifdef DEBUG_CHECK
                if ( !elem->giveSharedEdges()->giveSize() ) {
                    OOFEM_ERROR("element %d sharing nodes %d %d has no shared edges",
                                 elem->giveNumber(), iNode, jNode);
                }

 #endif
                elem->setSharedEdge(elem->giveEdgeIndex(iNode, jNode), 0);
            }
        }

        this->sharedEdgesQueue.clear();
    }
}



int
Subdivision :: packSharedEdges(Subdivision *s, ProcessCommunicator &pc)
{
    int pi, iNode, jNode, rproc = pc.giveRank();
    int myrank = this->giveRank();
    const IntArray *sharedPartitions;
    std :: list< int > :: const_iterator sharedEdgeQueueIter;
    IntArray edgeInfo(2);

    if ( rproc == myrank ) {
        return 1;                  // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
    for ( sharedEdgeQueueIter = sharedEdgesQueue.begin();
          sharedEdgeQueueIter != sharedEdgesQueue.end();
          sharedEdgeQueueIter++ ) {
        pi = ( * sharedEdgeQueueIter );

        sharedPartitions = this->mesh->giveEdge(pi)->givePartitions();
        if ( sharedPartitions->contains(rproc) ) {
            this->mesh->giveEdge(pi)->giveEdgeNodes(iNode, jNode);
            edgeInfo.at(1) = this->mesh->giveNode(iNode)->giveGlobalNumber();
            edgeInfo.at(2) = this->mesh->giveNode(jNode)->giveGlobalNumber();
 #ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_INFO("[%d] Subdivision::packSharedEdges: sending [%d %d] to %d\n", myrank, edgeInfo.at(1), edgeInfo.at(2), rproc);
 #endif
            pcbuff->write(SUBDIVISION_SHARED_EDGE_REC_TAG);                // KOKO zbytecne
            edgeInfo.storeYourself(*pcbuff);
        }
    }

    pcbuff->write(SUBDIVISION_END_DATA);

    return 1;
}


int
Subdivision :: unpackSharedEdges(Subdivision *s, ProcessCommunicator &pc)
{
    int myrank = this->giveRank();
    int _type, iproc = pc.giveRank();
    int iNode, jNode, elems;
    IntArray edgeInfo(2);
    const IntArray *iElems, *jElems;
    std :: list< int > :: const_iterator sharedEdgeQueueIter;
    Subdivision :: RS_SharedEdge *edge;
    Subdivision :: RS_Element *elem;
    int eIndex;

    if ( iproc == myrank ) {
        return 1;                                 // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();

    pcbuff->read(_type);
    // unpack dofman data
    while ( _type != SUBDIVISION_END_DATA ) {
        if ( _type == SUBDIVISION_SHARED_EDGE_REC_TAG ) {   // KOKO zbytecne
            edgeInfo.restoreYourself(*pcbuff);

 #ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_INFO("[%d] Subdivision::unpackSharedEdges: receiving [%d %d] from %d\n",
                           myrank, edgeInfo.at(1), edgeInfo.at(2), iproc);
 #endif

            iNode = mesh->sharedNodeGlobal2Local( edgeInfo.at(1) );
            jNode = mesh->sharedNodeGlobal2Local( edgeInfo.at(2) );

            if (!(iNode && jNode)) {
               pcbuff->read(_type);
               continue;
            }

            // get elements incident simultaneiusly to iNode and jNode
            iElems = mesh->giveNode(iNode)->giveConnectedElements();
            jElems = mesh->giveNode(jNode)->giveConnectedElements();

            IntArray common;
            if ( iElems->giveSize() <= jElems->giveSize() ) {
                common.preallocate( iElems->giveSize() );
            } else {
                common.preallocate( jElems->giveSize() );
            }

            // I do rely on the fact that the arrays are ordered !!!
            // I am using zero chunk because array common is large enough
            elems = iElems->findCommonValuesSorted(* jElems, common, 0);
            if ( elems ) {
                // find the relevant edge on the first element
                elem = mesh->giveElement( common.at(1) );
                eIndex = elem->giveEdgeIndex(iNode, jNode);
 #ifdef DEBUG_CHECK
                if ( !elem->giveSharedEdges()->giveSize() ) {
                    OOFEM_ERROR("element %d sharing nodes %d %d has no shared edges",
                                 elem->giveNumber(), iNode, jNode);
                }

                if ( !elem->giveSharedEdge(eIndex) ) {
                    OOFEM_ERROR("element %d sharing nodes %d and %d has no shared edge %d",
                                 elem->giveNumber(), iNode, jNode, eIndex);
                }

 #endif
                edge = mesh->giveEdge( elem->giveSharedEdge(eIndex) );
                // use zero chunk; the array is large enough from tentative set of partitions
                edge->addPartition(iproc, 0);

 #ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_INFO( "[%d] Partition %d added to shared edge %d (%d %d) [%d %d]\n",
                               myrank, iproc, elem->giveSharedEdge(eIndex), iNode, jNode, edgeInfo.at(1), edgeInfo.at(2) );
 #endif
            }
        } else {
            OOFEM_ERROR("unknown tag received");
        }

        pcbuff->read(_type);
    }

    return 1;
}



void
Subdivision :: RS_Mesh :: insertGlobalSharedNodeMap(Subdivision :: RS_Node *node)
{
    int key;
    key = node->giveGlobalNumber();
    sharedNodeMap [ key ] = node;
}


int
Subdivision :: RS_Mesh :: sharedNodeGlobal2Local(int _globnum)
{
    if ( sharedNodeMap.find(_globnum) != sharedNodeMap.end() ) {
        // node is already available -> update only
        return ( sharedNodeMap [ _globnum ]->giveNumber() );
    } else {
        return 0;
    }
}

#endif



int
Subdivision :: RS_CompareNodePositions :: operator() (int i, int j)
{
    double icoord, jcoord;

    icoord = m->giveNode(i)->giveCoordinates()->at(1);
    jcoord = m->giveNode(j)->giveCoordinates()->at(1);

    if ( icoord < jcoord ) {
        return -1;
    }

    if ( icoord > jcoord ) {
        return 1;
    }

    icoord = m->giveNode(i)->giveCoordinates()->at(2);
    jcoord = m->giveNode(j)->giveCoordinates()->at(2);

    if ( icoord < jcoord ) {
        return -1;
    }

    if ( icoord > jcoord ) {
        return 1;
    }

    icoord = m->giveNode(i)->giveCoordinates()->at(3);
    jcoord = m->giveNode(j)->giveCoordinates()->at(3);

    if ( icoord < jcoord ) {
        return -1;
    }

    if ( icoord > jcoord ) {
        return 1;
    }

    return 0;
}
} // end namespace oofem
