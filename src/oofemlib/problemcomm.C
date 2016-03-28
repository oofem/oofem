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

#include "problemcomm.h"
#include "intarray.h"
#include "error.h"
#include "engngm.h"
#include "element.h"
#include "dofmanager.h"

#ifdef __USE_MPI
 #include <mpi.h>
#endif

#define  __VERBOSE_PARALLEL

namespace oofem {
ProblemCommunicator :: ProblemCommunicator(EngngModel *emodel, CommunicatorBuff *b, int rank, int size) :
    Communicator(emodel, b, rank, size)
{
    this->initialized = false;
}


ProblemCommunicator :: ~ProblemCommunicator()
{ }


NodeCommunicator :: NodeCommunicator(EngngModel *emodel, CommunicatorBuff *b, int rank, int size) : 
    ProblemCommunicator(emodel, b, rank, size)
{ }

ElementCommunicator :: ElementCommunicator(EngngModel *emodel, CommunicatorBuff *b, int rank, int size) : 
    ProblemCommunicator(emodel, b, rank, size)
{ }

void
NodeCommunicator :: setUpCommunicationMaps(EngngModel *pm, bool excludeSelfCommFlag, bool forceReinit)
{
#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("NodeCommunicator :: setUpCommunicationMaps", "Setting up communication maps", rank);
#endif

    if ( !forceReinit && initialized ) {
        return;
    }

    Domain *domain = pm->giveDomain(1);
    int nnodes = domain->giveNumberOfDofManagers();

    //
    // receive and send maps are same and are assembled locally
    // using DofManager's partition lists.
    //

    IntArray domainNodeSendCount(size);

    for ( int i = 1; i <= nnodes; i++ ) {
        DofManager *dman = domain->giveDofManager(i);
        const IntArray *partitionList = dman->givePartitionList();
        if ( dman->giveParallelMode() == DofManager_shared ) {
            for ( int j = 1; j <= partitionList->giveSize(); j++ ) {
                if ( !( excludeSelfCommFlag && ( this->rank == partitionList->at(j) ) ) ) {
                    domainNodeSendCount.at(partitionList->at(j) + 1)++;
                }
            }
        }
    }

    // build maps simultaneously
    IntArray pos(size);
    std :: vector< IntArray >maps( size );
    for ( int i = 0; i < size; i++ ) {
        maps [ i ].resize( domainNodeSendCount.at ( i + 1 ) );
    }


    for ( int i = 1; i <= nnodes; i++ ) {
        DofManager *dman = domain->giveDofManager(i);
        // if combination node & element cut can occur, test for shared DofMan mode
        const IntArray *partitionList = dman->givePartitionList();
        if ( dman->giveParallelMode() == DofManager_shared ) {
            for ( int j = 1; j <= partitionList->giveSize(); j++ ) {
                int partition = partitionList->at(j);
                if ( !( excludeSelfCommFlag && ( this->rank == partition ) ) ) {
                    maps [ partition ].at( ++pos.at(partition + 1) ) = i;
                }
            }
        }
    }

    // set up domain communicators maps
    for ( int i = 0; i < size; i++ ) {
        this->setProcessCommunicatorToSendArry(this->giveProcessCommunicator(i), maps [ i ]);
        this->setProcessCommunicatorToRecvArry(this->giveProcessCommunicator(i), maps [ i ]);
        //this->giveDomainCommunicator(i)->setToSendArry (this->engngModel, maps[i]);
        //this->giveDomainCommunicator(i)->setToRecvArry (this->engngModel, maps[i]);
    }

    initialized = true;
}


void
ElementCommunicator :: setUpCommunicationMaps(EngngModel *pm,  bool excludeSelfCommFlag, bool forceReinit)
{
#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("ElementCommunicator :: setUpCommunicationMaps", "Setting up communication maps", rank);
#endif

    if ( !forceReinit && initialized ) {
        return;
    }

    OOFEM_LOG_RELEVANT("[%d] ElementCommunicator :: Setting up communication maps\n", rank);

    Domain *domain = pm->giveDomain(1);

    /*
     * Initially, each partition knows for which nodes a receive
     * is needed (and can therefore compute easily the recv map),
     * but does not know for which nodes it should send data to which
     * partition. Hence, the communication setup is performed by
     * broadcasting "send request" lists of nodes for which
     * a partition expects to receive data (ie. of those nodes
     * which the partition uses, but does not own) to all
     * collaborating processes. The "send request" list are
     * converted into send maps.
     */

    // receive maps can be build locally,
    // but send maps should be assembled from broadcasted lists (containing
    // expected receive nodes) of remote partitions.

    // first build local receive map
    IntArray domainNodeRecvCount(size);
    int domainRecvListSize = 0, domainRecvListPos = 0;
    int nelems;
    int result = 1;

    nelems = domain->giveNumberOfElements();
    for ( int i = 1; i <= nelems; i++ ) {
        Element *element = domain->giveElement(i);
        const IntArray *partitionList = element->givePartitionList();
        if ( element->giveParallelMode() == Element_remote ) {
            // size of partitionList should be 1 <== only ine master
            for ( int j = 1; j <= partitionList->giveSize(); j++ ) {
                if ( !( excludeSelfCommFlag && ( this->rank == partitionList->at(j) ) ) ) {
                    domainRecvListSize++;
                    domainNodeRecvCount.at(partitionList->at(j) + 1)++;
                }
            }
        }
    }

    // build maps simultaneously
    IntArray pos(size);
    std :: vector< IntArray >maps( size );
    for ( int i = 0; i < size; i++ ) {
        maps [ i ].resize( domainNodeRecvCount.at ( i + 1 ) );
    }

    // allocate also domain receive list to be broadcasted
    IntArray domainRecvList(domainRecvListSize);

    if ( domainRecvListSize ) {
        for ( int i = 1; i <= nelems; i++ ) {
            // test if element is remote one
            Element *element = domain->giveElement(i);
            if ( element->giveParallelMode() == Element_remote ) {
                domainRecvList.at(++domainRecvListPos) = element->giveGlobalNumber();

                const IntArray *partitionList = element->givePartitionList();
                // size of partitionList should be 1 <== only ine master
                for ( int j = 1; j <= partitionList->giveSize(); j++ ) {
                    if ( !( excludeSelfCommFlag && ( this->rank == partitionList->at(j) ) ) ) {
                        int partition = partitionList->at(j);
                        maps [ partition ].at( ++pos.at(partition + 1) ) = i;
                    }
                }
            }
        }
    }

    // set up domains recv communicator maps
    for ( int i = 0; i < size; i++ ) {
        this->setProcessCommunicatorToRecvArry(this->giveProcessCommunicator(i), maps [ i ]);
        //this->giveDomainCommunicator(i)->setToRecvArry(this->engngModel, maps [ i ]);
    }


#ifdef __VERBOSE_PARALLEL
    for (int i=0; i<size; i++) {
      fprintf (stderr, "domain %d-%d: domainCommRecvsize is %d\n",rank,i,this->giveProcessCommunicator(i)->giveToRecvMap()->giveSize() );
      printf ("domain %d-%d: reecv map:",rank,i);
      this->giveProcessCommunicator(i)->giveToRecvMap()->printYourself();
    }
#endif
        

    // to assemble send maps, we must analyze broadcasted remote domain send lists
    // and we must also broadcast our send list.

#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Remote Element-cut broadcasting started", rank);
#endif


    StaticCommunicationBuffer commBuff(MPI_COMM_WORLD);
    IntArray remoteDomainRecvList;
    IntArray toSendMap;
    int localExpectedSize, globalRecvSize;
    int sendMapPos, sendMapSize, globalDofManNum;

    // determine the size of receive buffer using AllReduce operation
#ifndef IBM_MPI_IMPLEMENTATION
    localExpectedSize = domainRecvList.givePackSize(commBuff);
#else
    localExpectedSize = domainRecvList.givePackSize(commBuff) + 1;
#endif


#ifdef __USE_MPI
    result = MPI_Allreduce(& localExpectedSize, & globalRecvSize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if ( result != MPI_SUCCESS ) {
        OOFEM_ERROR("MPI_Allreduce failed");
    }

#else
WARNING: NOT SUPPORTED MESSAGE PARSING LIBRARY
#endif

#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Finished reducing receiveBufferSize", rank);
#endif


    // resize to fit largest received message
    commBuff.resize(globalRecvSize);

    // resize toSend map to max possible size
    toSendMap.resize(globalRecvSize);

    for ( int i = 0; i < size; i++ ) { // loop over domains
        commBuff.init();
        if ( i == rank ) {
            //current domain has to send its receive list to all domains
            // broadcast domainRecvList

#ifdef __VERBOSE_PARALLEL
            VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Broadcasting own send list", rank);
#endif

            domainRecvList.storeYourself(commBuff);
            result = commBuff.bcast(i);
            if ( result != MPI_SUCCESS ) {
                OOFEM_ERROR("commBuff broadcast failed");
            }

#ifdef __VERBOSE_PARALLEL
            VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Broadcasting own send list finished", rank);
#endif
        } else {
#ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d\n",
                            rank, "ProblemCommunicator :: unpackAllData", i);
#endif
            // receive broadcasted lists
            result = commBuff.bcast(i);
            if ( result != MPI_SUCCESS ) {
                OOFEM_ERROR("commBuff broadcast failed");
            }

#ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d finished\n",
                            rank, "ProblemCommunicator :: unpackAllData", i);
#endif


            // unpack remote receive list
            if ( remoteDomainRecvList.restoreYourself(commBuff) != CIO_OK ) {
                OOFEM_ERROR("unpack remote receive list failed");
            }

            // find if remote elements are in local partition
            // if yes add them into send map for correcponding i-th partition
            sendMapPos = 0;
            sendMapSize = 0;
            // determine sendMap size
            for ( int j = 1; j <= nelems; j++ ) { // loop over local elements
                Element *element = domain->giveElement(j);
                if ( element->giveParallelMode() == Element_local ) {
                    globalDofManNum = element->giveGlobalNumber();
                    // test id globalDofManNum is in remoteDomainRecvList
                    if ( remoteDomainRecvList.findFirstIndexOf(globalDofManNum) ) {
                        sendMapSize++;
                    }
                }
            }

            toSendMap.resize(sendMapSize);

            for ( int j = 1; j <= nelems; j++ ) { // loop over local elements
                Element *element = domain->giveElement(j);
                if ( element->giveParallelMode() == Element_local ) {
                    globalDofManNum = element->giveGlobalNumber();
                    // test id globalDofManNum is in remoteDomainRecvList
                    if ( remoteDomainRecvList.findFirstIndexOf(globalDofManNum) ) {
                        // add this local DofManager number to sed map for active partition
                        toSendMap.at(++sendMapPos) = j;
                    }
                }
            } // end loop over local DofManagers

            // set send map to i-th process communicator
            this->setProcessCommunicatorToSendArry(this->giveProcessCommunicator(i), toSendMap);

#ifdef __VERBOSE_PARALLEL
            fprintf (stderr, "domain %d-%d: domainCommSendsize is %d\n",rank,i,this->giveProcessCommunicator(i)->giveToSendMap()->giveSize() );
            printf ("domain %d-%d: send map:",rank,i);
            this->giveProcessCommunicator(i)->giveToSendMap()->printYourself();
            
#endif
                

            //this->giveDomainCommunicator(i)->setToSendArry (this->engngModel, toSendMap);
        } // end receiving broadcasted lists

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Receiving broadcasted send maps finished", rank);
#endif
    } // end loop over domains

    initialized = true;
}

int
NodeCommunicator :: setProcessCommunicatorToSendArry(ProcessCommunicator *processComm, IntArray &map)
{
    sortCommMap(map, & ProblemCommunicator :: DofManCmp);
    processComm->setToSendArry(engngModel, map, 0); ///@todo We could call separate functions for node and elements instead of passing an argument
    return 1;
}

int
NodeCommunicator :: setProcessCommunicatorToRecvArry(ProcessCommunicator *processComm, IntArray &map)
{
    sortCommMap(map, & ProblemCommunicator :: DofManCmp);
    processComm->setToRecvArry(engngModel, map, 0);
    return 1;
}

int
ElementCommunicator :: setProcessCommunicatorToSendArry(ProcessCommunicator *processComm, IntArray &map)
{
    sortCommMap(map, & ProblemCommunicator :: ElemCmp);
    processComm->setToSendArry(engngModel, map, 1);
    return 1;
}

int
ElementCommunicator :: setProcessCommunicatorToRecvArry(ProcessCommunicator *processComm, IntArray &map)
{
    sortCommMap(map, & ProblemCommunicator :: ElemCmp);
    processComm->setToRecvArry(engngModel, map, 1);
    return 1;
}


void
ProblemCommunicator :: sortCommMap( IntArray &map, int ( ProblemCommunicator :: *cmp )( int, int ) )
{
    this->quickSortCommMap(map, 1, map.giveSize(), cmp);
}


void
ProblemCommunicator :: quickSortCommMap( IntArray &map, int l, int r, int ( ProblemCommunicator :: *cmp )( int, int ) )
{
    if ( r <= l ) {
        return;
    }

    int i = quickSortPartition(map, l, r, cmp);
    quickSortCommMap(map, l, i - 1, cmp);
    quickSortCommMap(map, i + 1, r, cmp);
}


int
ProblemCommunicator :: quickSortPartition( IntArray &map, int l, int r, int ( ProblemCommunicator :: *cmp )( int, int ) )
{
    int i = l - 1, j = r;
    int v = map.at(r);
    int swap;

    for ( ; ; ) {
        while ( ( ( this->*cmp )(map.at(++i), v) ) < 0 ) {
            ;
        }

        while ( ( ( this->*cmp )( v, map.at(--j) ) ) < 0 ) {
            if ( j == l ) {
                break;
            }
        }

        if ( i >= j ) {
            break;
        }

        swap = map.at(i);
        map.at(i) = map.at(j);
        map.at(j) = swap;
    }

    swap = map.at(i);
    map.at(i) = map.at(r);
    map.at(r) = swap;
    return i;
}


int
ProblemCommunicator :: DofManCmp(int i, int j)
{
    return ( engngModel->giveDomain(1)->giveDofManager(i)->giveGlobalNumber() -
            engngModel->giveDomain(1)->giveDofManager(j)->giveGlobalNumber() );
}
int
ProblemCommunicator :: ElemCmp(int i, int j)
{
    return ( engngModel->giveDomain(1)->giveElement(i)->giveGlobalNumber() -
            engngModel->giveDomain(1)->giveElement(j)->giveGlobalNumber() );
}
} // end namespace oofem

