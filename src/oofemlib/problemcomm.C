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

#ifdef __PARALLEL_MODE

#include "problemcomm.h"
#include "intarray.h"
#include "error.h"
#include "element.h"

#ifdef __USE_MPI
 #include <mpi.h>
#endif

namespace oofem {
ProblemCommunicator :: ProblemCommunicator(EngngModel *emodel, CommunicatorBuff *b, int rank, int size,
                                           ProblemCommunicatorMode mode) :
    Communicator(emodel, b, rank, size)
{
    this->mode = mode;
    this->initialized = false;
}


ProblemCommunicator :: ~ProblemCommunicator()
{ }


void
ProblemCommunicator :: setUpCommunicationMapsForNodeCut(EngngModel *pm, bool excludeSelfCommFlag)
{
    Domain *domain = pm->giveDomain(1);
    int nnodes = domain->giveNumberOfDofManagers();
    int i, j, partition;

    if ( this->mode != ProblemCommMode__NODE_CUT ) {
        _error("setUpCommunicationMapsForNodeCut: invalid mode");
    }

    //
    // receive and send maps are same and are assembled locally
    // using DofManager's partition lists.
    //

    IntArray domainNodeSendCount(size);
    const IntArray *partitionList;

    for ( i = 1; i <= nnodes; i++ ) {
        partitionList = domain->giveDofManager(i)->givePartitionList();
        if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
            for ( j = 1; j <= partitionList->giveSize(); j++ ) {
                if ( !( excludeSelfCommFlag && ( this->rank == partitionList->at(j) ) ) ) {
                    domainNodeSendCount.at(partitionList->at(j) + 1)++;
                }
            }
        }
    }

    // build maps simultaneously
    IntArray pos(size);
    IntArray **maps = new IntArray * [ size ];
    for ( i = 0; i < size; i++ ) {
        maps [ i ] = new IntArray( domainNodeSendCount.at(i + 1) );
    }


    for ( i = 1; i <= nnodes; i++ ) {
        // if combination node & element cut can occur, test for shared DofMan mode
        partitionList = domain->giveDofManager(i)->givePartitionList();
        if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
            for ( j = 1; j <= partitionList->giveSize(); j++ ) {
                partition = partitionList->at(j);
                if ( !( excludeSelfCommFlag && ( this->rank == partition ) ) ) {
                    maps [ partition ]->at( ++pos.at(partition + 1) ) = i;
                }
            }
        }
    }

    // set up domain communicators maps
    for ( i = 0; i < size; i++ ) {
        this->setProcessCommunicatorToSendArry(this->giveProcessCommunicator(i), * maps [ i ]);
        this->setProcessCommunicatorToRecvArry(this->giveProcessCommunicator(i), * maps [ i ]);
        //this->giveDomainCommunicator(i)->setToSendArry (this->engngModel, *maps[i]);
        //this->giveDomainCommunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
    }

    // delete local maps
    for ( i = 0; i < size; i++ ) {
        delete maps [ i ];
    }

    delete[] maps;
}


void
ProblemCommunicator :: setUpCommunicationMapsForElementCut(EngngModel *pm,
                                                           bool excludeSelfCommFlag)
{
    Domain *domain = pm->giveDomain(1);
    int nnodes = domain->giveNumberOfDofManagers();
    int i, j, partition;

    if ( this->mode == ProblemCommMode__ELEMENT_CUT ) {
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
        const IntArray *partitionList;
        DofManager *dofMan;
        //Element    *element;
        int domainRecvListSize = 0, domainRecvListPos = 0;
        //int nelems;
        int result = 1;

        for ( i = 1; i <= nnodes; i++ ) {
            partitionList = domain->giveDofManager(i)->givePartitionList();
            if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_remote ) {
                // size of partitionList should be 1 <== only ine master
                for ( j = 1; j <= partitionList->giveSize(); j++ ) {
                    if ( !( excludeSelfCommFlag && ( this->rank == partitionList->at(j) ) ) ) {
                        domainRecvListSize++;
                        domainNodeRecvCount.at(partitionList->at(j) + 1)++;
                    }
                }
            }
        }

        // build maps simultaneously
        IntArray pos(size);
        IntArray **maps = new IntArray * [ size ];
        for ( i = 0; i < size; i++ ) {
            maps [ i ] = new IntArray( domainNodeRecvCount.at(i + 1) );
        }

        // allocate also domain receive list to be broadcasted
        IntArray domainRecvList(domainRecvListSize);

        if ( domainRecvListSize ) {
            for ( i = 1; i <= nnodes; i++ ) {
                // test if node is remote DofMan
                dofMan = domain->giveDofManager(i);
                if ( dofMan->giveParallelMode() == DofManager_remote ) {
                    domainRecvList.at(++domainRecvListPos) = dofMan->giveGlobalNumber();

                    partitionList = domain->giveDofManager(i)->givePartitionList();
                    // size of partitionList should be 1 <== only ine master
                    for ( j = 1; j <= partitionList->giveSize(); j++ ) {
                        if ( !( excludeSelfCommFlag && ( this->rank == partitionList->at(j) ) ) ) {
                            partition = partitionList->at(j);
                            maps [ partition ]->at( ++pos.at(partition + 1) ) = i;
                        }
                    }
                }
            }
        }

        // set up process recv communicator maps
        for ( i = 0; i < size; i++ ) {
            this->setProcessCommunicatorToRecvArry(this->giveProcessCommunicator(i), * maps [ i ]);
            //this->giveDomainCommunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
        }

        // delete local maps
        for ( i = 0; i < size; i++ ) {
            delete maps [ i ];
        }

        delete[] maps;

        // to assemble send maps, we must analyze broadcasted remote domain send lists
        // and we must also broadcast our send list.

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Element-cut broadcasting started", rank);
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
            _error("setUpCommunicationMaps: MPI_Allreduce failed");
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

        for ( i = 0; i < size; i++ ) { // loop over domains
            commBuff.init();
            if ( i == rank ) {
                //current domain has to send its receive list to all domains
                // broadcast domainRecvList

#ifdef __VERBOSE_PARALLEL
                VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Broadcasting own send list", rank);
#endif

                commBuff.packIntArray(domainRecvList);
                result = commBuff.bcast(i);
                if ( result != MPI_SUCCESS ) {
                    _error("setUpCommunicationMaps: commBuff broadcast failed");
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
                    _error("setUpCommunicationMaps: commBuff broadcast failed");
                }

#ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d finished\n",
                                rank, "ProblemCommunicator :: unpackAllData", i);
#endif


                // unpack remote receive list
                if ( !commBuff.unpackIntArray(remoteDomainRecvList) ) {
                    _error("ProblemCommunicator::setUpCommunicationMaps: unpack remote receive list failed");
                }

                // find if remote nodes are in local partition
                // if yes add them into send map for correcponding i-th partition
                sendMapPos = 0;
                sendMapSize = 0;
                // determine sendMap size
                for ( j = 1; j <= nnodes; j++ ) { // loop over local DofManagers
                    dofMan = domain->giveDofManager(j);
                    globalDofManNum = dofMan->giveGlobalNumber();
                    // test id globalDofManNum is in remoteDomainRecvList
                    if ( remoteDomainRecvList.findFirstIndexOf(globalDofManNum) ) {
                        sendMapSize++;
                    }
                }

                toSendMap.resize(sendMapSize);

                for ( j = 1; j <= nnodes; j++ ) { // loop over local DofManagers
                    dofMan = domain->giveDofManager(j);
                    globalDofManNum = dofMan->giveGlobalNumber();
                    // test id globalDofManNum is in remoteDomainRecvList
                    if ( remoteDomainRecvList.findFirstIndexOf(globalDofManNum) ) {
                        // add this local DofManager number to sed map for active partition
                        toSendMap.at(++sendMapPos) = j;
                    }
                } // end loop over local DofManagers

                // set send map to i-th process communicator
                this->setProcessCommunicatorToSendArry(this->giveProcessCommunicator(i), toSendMap);
                //this->giveDomainCommunicator(i)->setToSendArry (this->engngModel, toSendMap);
            } // end receiving broadcasted lists

#ifdef __VERBOSE_PARALLEL
            VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Receiving broadcasted send maps finished", rank);
#endif
        } // end loop over domains

    } else {
        _error("setUpCommunicationMapsForElementCut: unknown mode");
    }
}


void
ProblemCommunicator :: setUpCommunicationMapsForRemoteElementMode(EngngModel *pm,
                                                                  bool excludeSelfCommFlag)
{
    //int nnodes = domain->giveNumberOfDofManagers();
    Domain *domain = pm->giveDomain(1);
    int i, j, partition;

    if ( this->mode == ProblemCommMode__REMOTE_ELEMENT_MODE ) {
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
        const IntArray *partitionList;
        //DofManager *dofMan;
        Element *element;
        int domainRecvListSize = 0, domainRecvListPos = 0;
        int nelems;
        int result = 1;

        nelems = domain->giveNumberOfElements();
        for ( i = 1; i <= nelems; i++ ) {
            partitionList = domain->giveElement(i)->givePartitionList();
            if ( domain->giveElement(i)->giveParallelMode() == Element_remote ) {
                // size of partitionList should be 1 <== only ine master
                for ( j = 1; j <= partitionList->giveSize(); j++ ) {
                    if ( !( excludeSelfCommFlag && ( this->rank == partitionList->at(j) ) ) ) {
                        domainRecvListSize++;
                        domainNodeRecvCount.at(partitionList->at(j) + 1)++;
                    }
                }
            }
        }

        // build maps simultaneously
        IntArray pos(size);
        IntArray **maps = new IntArray * [ size ];
        for ( i = 0; i < size; i++ ) {
            maps [ i ] = new IntArray( domainNodeRecvCount.at(i + 1) );
        }

        // allocate also domain receive list to be broadcasted
        IntArray domainRecvList(domainRecvListSize);

        if ( domainRecvListSize ) {
            for ( i = 1; i <= nelems; i++ ) {
                // test if element is remote one
                element = domain->giveElement(i);
                if ( element->giveParallelMode() == Element_remote ) {
                    domainRecvList.at(++domainRecvListPos) = element->giveGlobalNumber();

                    partitionList = domain->giveElement(i)->givePartitionList();
                    // size of partitionList should be 1 <== only ine master
                    for ( j = 1; j <= partitionList->giveSize(); j++ ) {
                        if ( !( excludeSelfCommFlag && ( this->rank == partitionList->at(j) ) ) ) {
                            partition = partitionList->at(j);
                            maps [ partition ]->at( ++pos.at(partition + 1) ) = i;
                        }
                    }
                }
            }
        }

        // set up domains recv communicator maps
        for ( i = 0; i < size; i++ ) {
            this->setProcessCommunicatorToRecvArry(this->giveProcessCommunicator(i), * maps [ i ]);
            //this->giveDomainCommunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
        }

        /*
         * #ifdef __VERBOSE_PARALLEL
         * for (i=0; i<size; i++) {
         * fprintf (stderr, "domain %d-%d: domainCommRecvsize is %d\n",rank,i,this->giveDomainCommunicator(i)->giveRecvBuff()->giveSize() );
         * printf ("domain %d-%d: reecv map:",rank,i);
         * this->giveDomainCommunicator(i)->giveToRecvMap()->printYourself();
         * }
         *#endif
         */

        // delete local maps
        for ( i = 0; i < size; i++ ) {
            delete maps [ i ];
        }

        delete [] maps;

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
            _error("setUpCommunicationMaps: MPI_Allreduce failed");
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

        for ( i = 0; i < size; i++ ) { // loop over domains
            commBuff.init();
            if ( i == rank ) {
                //current domain has to send its receive list to all domains
                // broadcast domainRecvList

#ifdef __VERBOSE_PARALLEL
                VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Broadcasting own send list", rank);
#endif

                commBuff.packIntArray(domainRecvList);
                result = commBuff.bcast(i);
                if ( result != MPI_SUCCESS ) {
                    _error("setUpCommunicationMaps: commBuff broadcast failed");
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
                    _error("setUpCommunicationMaps: commBuff broadcast failed");
                }

#ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d finished\n",
                                rank, "ProblemCommunicator :: unpackAllData", i);
#endif


                // unpack remote receive list
                if ( !commBuff.unpackIntArray(remoteDomainRecvList) ) {
                    _error("ProblemCommunicator::setUpCommunicationMaps: unpack remote receive list failed");
                }

                // find if remote elements are in local partition
                // if yes add them into send map for correcponding i-th partition
                sendMapPos = 0;
                sendMapSize = 0;
                // determine sendMap size
                for ( j = 1; j <= nelems; j++ ) { // loop over local elements
                    element = domain->giveElement(j);
                    if ( element->giveParallelMode() == Element_local ) {
                        globalDofManNum = element->giveGlobalNumber();
                        // test id globalDofManNum is in remoteDomainRecvList
                        if ( remoteDomainRecvList.findFirstIndexOf(globalDofManNum) ) {
                            sendMapSize++;
                        }
                    }
                }

                toSendMap.resize(sendMapSize);

                for ( j = 1; j <= nelems; j++ ) { // loop over local elements
                    element = domain->giveElement(j);
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

                /*
                 * #ifdef __VERBOSE_PARALLEL
                 *  fprintf (stderr, "domain %d-%d: domainCommSendsize is %d\n",rank,i,this->giveDomainCommunicator(i)->giveSendBuff()->giveSize() );
                 *  printf ("domain %d-%d: send map:",rank,i);
                 *  this->giveDomainCommunicator(i)->giveToSendMap()->printYourself();
                 *
                 *#endif
                 */

                //this->giveDomainCommunicator(i)->setToSendArry (this->engngModel, toSendMap);
            } // end receiving broadcasted lists

#ifdef __VERBOSE_PARALLEL
            VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Receiving broadcasted send maps finished", rank);
#endif
        } // end loop over domains

    } else {
        _error("setUpCommunicationMapsForRemoteElementMode: unknown mode");
    }
}

void
ProblemCommunicator :: setUpCommunicationMaps(EngngModel *pm, bool excludeSelfCommFlag,
                                              bool forceReinit)
{
#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("ProblemCommunicator :: setUpCommunicationMaps", "Setting up communication maps", rank);
#endif

    if ( !forceReinit && initialized ) {
        return;
    }

    if ( this->mode == ProblemCommMode__NODE_CUT ) {
        setUpCommunicationMapsForNodeCut(pm, excludeSelfCommFlag);
    } else if ( this->mode == ProblemCommMode__ELEMENT_CUT ) {
        setUpCommunicationMapsForElementCut(pm, excludeSelfCommFlag);
    } else if ( this->mode == ProblemCommMode__REMOTE_ELEMENT_MODE ) {
        setUpCommunicationMapsForRemoteElementMode(pm, excludeSelfCommFlag);
    } else {
        _error("setUpCommunicationMaps: unknown mode");
    }

    initialized = true;
}

/*
 *
 * //
 * // Old compact version
 * //
 * void
 * ProblemCommunicator :: setUpCommunicationMaps (Domain* domain)
 * {
 * int nnodes = domain->giveNumberOfDofManagers();
 * int i, j, partition;
 *
 *#ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("ProblemCommunicator :: setUpCommunicationMaps", "Setting up communication maps", rank);
 *#endif
 *
 *
 * if (this->mode == PC__NODE_CUT) {
 *
 * //
 * // receive and send maps are same and are assembled locally
 * // using DofManager's partition lists.
 * //
 *
 * IntArray domainNodeSendCount(size);
 * const IntArray* partitionList;
 *
 * for (i=1; i<= nnodes; i++) {
 * partitionList = domain->giveDofManager (i) -> givePartitionList();
 * if (domain->giveDofManager (i)->giveParallelMode () == DofManager_shared) {
 *  for (j=1; j<=partitionList->giveSize(); j++) {
 *   domainNodeSendCount.at(partitionList->at(j)+1) ++;
 *  }
 * }
 * }
 *
 * // build maps simultaneosly
 * IntArray pos(size);
 * IntArray **maps = new IntArray* [size];
 * for (i=0; i<size; i++) maps[i]=new IntArray(domainNodeSendCount.at(i+1));
 *
 *
 * for (i=1; i<= nnodes; i++) {
 * // if combination node & element cut can occur, test for shared DofMan mode
 * partitionList = domain->giveDofManager (i) -> givePartitionList();
 * if (domain->giveDofManager (i)->giveParallelMode () == DofManager_shared) {
 *  for (j=1; j<=partitionList->giveSize(); j++) {
 *   partition = partitionList->at(j);
 *   maps[partition]->at(++pos.at(partition+1)) = i;
 *  }
 * }
 * }
 * // set up domain communicators maps
 * for (i=0; i<size; i++) {
 * this->setDomainCommunicatorToSendArry (this->giveDomainCommunicator(i), *maps[i]);
 * this->setDomainCommunicatorToRecvArry (this->giveDomainCommunicator(i), *maps[i]);
 * //this->giveDomainCommunicator(i)->setToSendArry (this->engngModel, *maps[i]);
 * //this->giveDomainCommunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
 * }
 *
 * // delete local maps
 * for (i=0; i<size; i++) delete maps[i];
 * delete maps;
 *
 * } else if ((this->mode == PC__ELEMENT_CUT) ||
 *     (this->mode == PC__REMOTE_ELEMENT_MODE)) {
 *
 *
 * //    Initaially, each partition knows for which nodes a receive
 * //    is needed (and can therefore compute easily the recv map),
 * //    but does not know for which nodes it should send data to which
 * //    partition. Hence, the communication setup is performed by
 * //    broadcasting "send request" lists of nodes for which
 * //    a partition expects to receive data (ie. of those nodes
 * //    which the partition uses, but does not own) to all
 * //    collaborating processes. The "send request" list are
 * //    converted into send maps.
 *
 * // receive maps can be build localy,
 * // but send maps should be assembled from broadcasted lists (containing
 * // expected receive nodes) of remote partitions.
 *
 * // first build local receive map
 * IntArray domainNodeRecvCount(size);
 * const IntArray* partitionList;
 * DofManager *dofMan;
 * Element    *element;
 * int domainRecvListSize = 0, domainRecvListPos = 0;
 * int nelems;
 * int result = 1;
 *
 * if (this->mode == PC__ELEMENT_CUT) {
 * for (i=1; i<= nnodes; i++) {
 *  partitionList = domain->giveDofManager (i) -> givePartitionList();
 *  if (domain->giveDofManager (i)->giveParallelMode () == DofManager_remote) {
 *   // size of partitionList should be 1 <== only ine master
 *   for (j=1; j<=partitionList->giveSize(); j++) {
 *    domainRecvListSize++;
 *    domainNodeRecvCount.at(partitionList->at(j)+1) ++;
 *   }
 *  }
 * }
 * } else {            // PC__REMOTE_ELEMENT_MODE
 * nelems = domain->giveNumberOfElements();
 * for (i=1; i<= nelems; i++) {
 *  partitionList = domain->giveElement (i) -> givePartitionList();
 *  if (domain->giveElement (i)->giveParallelMode () == Element_remote) {
 *   // size of partitionList should be 1 <== only ine master
 *   for (j=1; j<=partitionList->giveSize(); j++) {
 *    domainRecvListSize++;
 *    domainNodeRecvCount.at(partitionList->at(j)+1) ++;
 *   }
 *  }
 * }
 * }
 *
 * // build maps simultaneosly
 * IntArray pos(size);
 * IntArray **maps = new IntArray* [size];
 * for (i=0; i<size; i++) maps[i]=new IntArray(domainNodeRecvCount.at(i+1));
 *
 * // allocate also domain receive list to be broadcasted
 * IntArray domainRecvList (domainRecvListSize);
 *
 * if (domainRecvListSize) {
 *
 * if (this->mode == PC__ELEMENT_CUT) {
 *  for (i=1; i<= nnodes; i++) {
 *   // test if node is remote DofMan
 *   dofMan = domain->giveDofManager (i);
 *   if (dofMan->giveParallelMode() == DofManager_remote) {
 *
 *    domainRecvList.at(++domainRecvListPos) = dofMan->giveGlobalNumber  ();
 *
 *    partitionList = domain->giveDofManager (i) -> givePartitionList();
 *    // size of partitionList should be 1 <== only ine master
 *    for (j=1; j<=partitionList->giveSize(); j++) {
 *     partition = partitionList->at(j);
 *     maps[partition]->at(++pos.at(partition+1)) = i;
 *    }
 *   }
 *  }
 * } else {          // PC__REMOTE_ELEMENT_MODE
 *  for (i=1; i<= nelems; i++) {
 *   // test if element is remote one
 *   element = domain->giveElement (i);
 *   if (element->giveParallelMode() == Element_remote) {
 *
 *    domainRecvList.at(++domainRecvListPos) = element->giveGlobalNumber  ();
 *
 *    partitionList = domain->giveElement (i) -> givePartitionList();
 *    // size of partitionList should be 1 <== only ine master
 *    for (j=1; j<=partitionList->giveSize(); j++) {
 *     partition = partitionList->at(j);
 *     maps[partition]->at(++pos.at(partition+1)) = i;
 *    }
 *   }
 *  }
 * }
 *
 * }
 * // set up domains recv communicator maps
 * for (i=0; i<size; i++) {
 * this->setDomainCommunicatorToRecvArry (this->giveDomainCommunicator(i), *maps[i]);
 * //this->giveDomainCommunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
 * }
 *
 * // delete local maps
 * for (i=0; i<size; i++) delete maps[i];
 * delete maps;
 *
 * // to assemble send maps, we must analyze broadcasted remote domain send lists
 * // and we must also broadcast our send list.
 *
 *#ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Element-cut broadcasting started", rank);
 *#endif
 *
 *
 * CommunicationBuffer commBuff (MPI_COMM_WORLD);
 * IntArray remoteDomainRecvList;
 * IntArray toSendMap;
 * int localExpectedSize, globalRecvSize ;
 * int sendMapPos, sendMapSize, globalDofManNum;
 *
 * // determine the size of receive buffer using AllReduce operation
 * localExpectedSize = domainRecvList.givePackSize(commBuff);
 *
 *#ifdef __USE_MPI
 * result = MPI_Allreduce (&localExpectedSize, &globalRecvSize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
 * if (result != MPI_SUCCESS) _error ("setUpCommunicationMaps: MPI_Allreduce failed");
 *#else
 * WARNING: NOT SUPPORTED MESSAGE PARSING LIBRARY
 *#endif
 *
 *#ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Finished reducing receiveBufferSize", rank);
 *#endif
 *
 *
 * // resize to fit largest received message
 * commBuff.resize (globalRecvSize);
 *
 * // resize toSend map to max possible size
 * toSendMap.resize (globalRecvSize);
 *
 * for (i=0; i<size; i++) { // loop over domains
 * commBuff.init();
 * if (i == rank) {
 *  //current domain has to send its receive list to all domains
 *  // broadcast domainRecvList
 *
 *#ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Broadcasting own send list", rank);
 *#endif
 *
 *  commBuff.packIntArray (domainRecvList);
 *  result = commBuff.bcast (i);
 *  if (result!= MPI_SUCCESS) _error ("setUpCommunicationMaps: commBuff broadcast failed");
 *
 *#ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Broadcasting own send list finished", rank);
 *#endif
 *
 * } else {
 *
 *#ifdef __VERBOSE_PARALLEL
 *  fprintf(stderr, "\n[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d",
 *      rank,"ProblemCommunicator :: unpackAllData", i);
 *#endif
 *  // receive broadcasted lists
 *  result = commBuff.bcast (i);
 *  if (result != MPI_SUCCESS) _error ("setUpCommunicationMaps: commBuff broadcast failed");
 *
 *#ifdef __VERBOSE_PARALLEL
 *  fprintf(stderr, "\n[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d finished",
 *      rank,"ProblemCommunicator :: unpackAllData", i);
 *#endif
 *
 *
 *  // unpack remote receive list
 *  if (!commBuff.unpackIntArray(remoteDomainRecvList))
 *   _error ("ProblemCommunicator::setUpCommunicationMaps: unpack remote receive list failed");
 *
 *  if (this->mode == PC__ELEMENT_CUT) {
 *   // find if remote nodes are in local partition
 *   // if yes add them into send map for correcponding i-th partition
 *   sendMapPos = 0;
 *   sendMapSize = 0;
 *   // determine sendMap size
 *   for (j=1; j<= nnodes; j++) { // loop over local DofManagers
 *    dofMan = domain->giveDofManager (j);
 *    globalDofManNum = dofMan->giveGlobalNumber();
 *    // test id globalDofManNum is in remoteDomainRecvList
 *    if (remoteDomainRecvList.findFirstIndexOf (globalDofManNum)) sendMapSize++;
 *   }
 *   toSendMap.resize (sendMapSize);
 *
 *   for (j=1; j<= nnodes; j++) { // loop over local DofManagers
 *    dofMan = domain->giveDofManager (j);
 *    globalDofManNum = dofMan->giveGlobalNumber();
 *    // test id globalDofManNum is in remoteDomainRecvList
 *    if (remoteDomainRecvList.findFirstIndexOf (globalDofManNum)) {
 *     // add this local DofManager number to sed map for active partition
 *     toSendMap.at(++sendMapPos) = j;
 *    }
 *   } // end loop over local DofManagers
 *  } else {       // PC__REMOTE_ELEMENT_MODE
 *   // find if remote elements are in local partition
 *   // if yes add them into send map for correcponding i-th partition
 *   sendMapPos = 0;
 *   sendMapSize = 0;
 *   // determine sendMap size
 *   for (j=1; j<= nelems; j++) { // loop over local elements
 *    element = domain->giveElement (j);
 *    globalDofManNum = element->giveGlobalNumber();
 *    // test id globalDofManNum is in remoteDomainRecvList
 *    if (remoteDomainRecvList.findFirstIndexOf (globalDofManNum)) sendMapSize++;
 *   }
 *   toSendMap.resize (sendMapSize);
 *
 *   for (j=1; j<= nelems; j++) { // loop over local elements
 *     element = domain->giveElement (j);
 *    globalDofManNum = element->giveGlobalNumber();
 *    // test id globalDofManNum is in remoteDomainRecvList
 *    if (remoteDomainRecvList.findFirstIndexOf (globalDofManNum)) {
 *     // add this local DofManager number to sed map for active partition
 *     toSendMap.at(++sendMapPos) = j;
 *    }
 *   } // end loop over local DofManagers
 *  }
 *
 *  // set send map to i-th domain communicator
 *  this->setDomainCommunicatorToSendArry (this->giveDomainCommunicator(i), toSendMap);
 *  //this->giveDomainCommunicator(i)->setToSendArry (this->engngModel, toSendMap);
 * } // end receiving broadcasted lists
 *#ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("ProblemCommunicator::setUpCommunicationMaps", "Receiving broadcasted send maps finished", rank);
 *#endif
 *
 * } // end loop over domains
 * } else _error ("setUpCommunicationMaps: unknown mode");
 * }
 */


int
ProblemCommunicator :: setProcessCommunicatorToSendArry(ProcessCommunicator *processComm, IntArray &map)
{
    if ( this->mode == ProblemCommMode__NODE_CUT ) {
        sortCommMap(map, & ProblemCommunicator :: DofManCmp);
    } else if ( this->mode == ProblemCommMode__REMOTE_ELEMENT_MODE ||
               this->mode == ProblemCommMode__ELEMENT_CUT ) {
        sortCommMap(map, & ProblemCommunicator :: ElemCmp);
    } else {
        _error("setDomainCommunicatorToSendArry: unknown mode");
    }

    processComm->setToSendArry(engngModel, map, this->mode);
    return 1;
}

int
ProblemCommunicator :: setProcessCommunicatorToRecvArry(ProcessCommunicator *processComm, IntArray &map)
{
    if ( this->mode == ProblemCommMode__NODE_CUT ) {
        sortCommMap(map, & ProblemCommunicator :: DofManCmp);
    } else if ( this->mode == ProblemCommMode__REMOTE_ELEMENT_MODE ||
               this->mode == ProblemCommMode__ELEMENT_CUT ) {
        sortCommMap(map, & ProblemCommunicator :: ElemCmp);
    } else {
        _error("setDomainCommunicatorToRecvArry: unknown mode");
    }

    processComm->setToRecvArry(engngModel, map, this->mode);
    return 1;
}



void
ProblemCommunicator :: sortCommMap( IntArray &map, int ( ProblemCommunicator :: *cmp ) (int, int) )
{
    this->quickSortCommMap(map, 1, map.giveSize(), cmp);
}


void
ProblemCommunicator :: quickSortCommMap( IntArray &map, int l, int r, int ( ProblemCommunicator :: *cmp ) (int, int) )
{
    if ( r <= l ) {
        return;
    }

    int i = quickSortPartition(map, l, r, cmp);
    quickSortCommMap(map, l, i - 1, cmp);
    quickSortCommMap(map, i + 1, r, cmp);
}

#define GLOBNUM(locnum) ( engngModel->giveDomain(1)->giveDofManager(locnum)->giveGlobalNumber() )


int
ProblemCommunicator :: quickSortPartition( IntArray &map, int l, int r, int ( ProblemCommunicator :: *cmp ) (int, int) )
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

#undef GLOBNUM

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
#endif
