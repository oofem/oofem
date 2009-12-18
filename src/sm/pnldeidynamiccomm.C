/* $Header: /home/cvs/bp/oofem/sm/src/pnldeidynamiccomm.C,v 1.2.4.1 2004/04/05 15:19:47 bp Exp $ */
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

#ifdef __PARALLEL_MODE

#include "pnldeidynamiccomm.h"
#include "pnldeidynamic.h"
#include "intarray.h"

#ifdef __USE_MPI
#ifndef __MAKEDEPEND
#include "mpi.h"
#endif
#endif

namespace oofem {

PNlDEIDynamicComunicator :: PNlDEIDynamicComunicator(PNlDEIDynamic *emodel, int rank, int size,
                                                     PNlDEIDynamicComunicatorMode mode) :
    Communicator< PNlDEIDynamic >(emodel, rank, size)
{
    this->mode = mode;
}


PNlDEIDynamicComunicator :: ~PNlDEIDynamicComunicator()
{ }


void
PNlDEIDynamicComunicator :: setUpCommunicationMapsForNodeCut(EngngModel *pm)
{
    Domain *domain = pm->giveDomain(1);
    int nnodes = domain->giveNumberOfDofManagers();
    int i, j, partition;

    if ( this->mode != PNlDEIDynamicComunicator__NODE_CUT ) {
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
                domainNodeSendCount.at(partitionList->at(j) + 1)++;
            }
        }
    }

    // build maps simultaneosly
    IntArray pos(size);
    IntArray **maps = new IntArray * [ size ];
    for ( i = 0; i < size; i++ ) {
        maps [ i ] = new IntArray( domainNodeSendCount.at(i + 1) );
    }


    for ( i = 1; i <= nnodes; i++ ) {
        // if combination node & eleemnt cut can occur, test for shared DofMan mode
        partitionList = domain->giveDofManager(i)->givePartitionList();
        if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
            for ( j = 1; j <= partitionList->giveSize(); j++ ) {
                partition = partitionList->at(j);
                maps [ partition ]->at( ++pos.at(partition + 1) ) = i;
            }
        }
    }

    // set up domain communicators maps
    for ( i = 0; i < size; i++ ) {
        this->setProblemComunicatorToSendArry(this->giveProblemComunicator(i), * maps [ i ]);
        this->setProblemComunicatorToRecvArry(this->giveProblemComunicator(i), * maps [ i ]);
        //this->giveDomainComunicator(i)->setToSendArry (this->engngModel, *maps[i]);
        //this->giveDomainComunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
    }

    // delete local maps
    for ( i = 0; i < size; i++ ) {
        delete maps [ i ];
    }

    delete maps;
}


void
PNlDEIDynamicComunicator :: setUpCommunicationMapsForElementCut(EngngModel *pm)
{
    Domain *domain = pm->giveDomain(1);
    int nnodes = domain->giveNumberOfDofManagers();
    int i, j, partition;

    if ( this->mode == PNlDEIDynamicComunicator__ELEMENT_CUT ) {
        /*
         * Initaially, each partition knows for which nodes a receive
         * is needed (and can therefore compute easily the recv map),
         * but does not know for which nodes it should send data to which
         * partition. Hence, the communication setup is performed by
         * broadcasting "send request" lists of nodes for which
         * a partition expects to receive data (ie. of those nodes
         * which the partition uses, but does not own) to all
         * collaborating processes. The "send request" list are
         * converted into send maps.
         */

        // receive maps can be build localy,
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
                    domainRecvListSize++;
                    domainNodeRecvCount.at(partitionList->at(j) + 1)++;
                }
            }
        }

        // build maps simultaneosly
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
                        partition = partitionList->at(j);
                        maps [ partition ]->at( ++pos.at(partition + 1) ) = i;
                    }
                }
            }
        }

        // set up problem recv communicator maps
        for ( i = 0; i < size; i++ ) {
            this->setProblemComunicatorToRecvArry(this->giveProblemComunicator(i), * maps [ i ]);
            //this->giveDomainComunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
        }

        // delete local maps
        for ( i = 0; i < size; i++ ) {
            delete maps [ i ];
        }

        delete maps;

        // to assemble send maps, we must analyze broadcasted remote domain send lists
        // and we must also broadcast our send list.

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Element-cut broadcasting started", rank);
#endif


        CommunicationBuffer commBuff(MPI_COMM_WORLD);
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
        VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Finished reducing receiveBufferSize", rank);
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
                VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Broadcasting own send list", rank);
#endif

                commBuff.packIntArray(domainRecvList);
                result = commBuff.bcast(i);
                if ( result != MPI_SUCCESS ) {
                    _error("setUpCommunicationMaps: commBuff broadcast failed");
                }

#ifdef __VERBOSE_PARALLEL
                VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Broadcasting own send list finished", rank);
#endif
            } else {
#ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_INFO("[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d\n",
                               rank, "PNlDEIDynamicComunicator :: unpackAllData", i);
#endif
                // receive broadcasted lists
                result = commBuff.bcast(i);
                if ( result != MPI_SUCCESS ) {
                    _error("setUpCommunicationMaps: commBuff broadcast failed");
                }

#ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_INFO("[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d finished\n",
                               rank, "PNlDEIDynamicComunicator :: unpackAllData", i);
#endif


                // unpack remote receive list
                if ( !commBuff.unpackIntArray(remoteDomainRecvList) ) {
                    _error("PNlDEIDynamicComunicator::setUpCommunicationMaps: unpack remote receive list failed");
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

                // set send map to i-th problem communicator
                this->setProblemComunicatorToSendArry(this->giveProblemComunicator(i), toSendMap);
                //this->giveDomainComunicator(i)->setToSendArry (this->engngModel, toSendMap);
            } // end receiving broadcasted lists

#ifdef __VERBOSE_PARALLEL
            VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Receiving broadcasted send maps finished", rank);
#endif
        } // end loop over domains

    } else {
        _error("setUpCommunicationMapsForElementCut: unknown mode");
    }
}


void
PNlDEIDynamicComunicator :: setUpCommunicationMapsForRemoteElementMode(EngngModel *pm)
{
    //int nnodes = domain->giveNumberOfDofManagers();
    Domain *domain = pm->giveDomain(1);
    int i, j, partition;

    if ( this->mode == PNlDEIDynamicComunicator__REMOTE_ELEMENT_MODE ) {
        /*
         * Initaially, each partition knows for which nodes a receive
         * is needed (and can therefore compute easily the recv map),
         * but does not know for which nodes it should send data to which
         * partition. Hence, the communication setup is performed by
         * broadcasting "send request" lists of nodes for which
         * a partition expects to receive data (ie. of those nodes
         * which the partition uses, but does not own) to all
         * collaborating processes. The "send request" list are
         * converted into send maps.
         */

        // receive maps can be build localy,
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
                    domainRecvListSize++;
                    domainNodeRecvCount.at(partitionList->at(j) + 1)++;
                }
            }
        }

        // build maps simultaneosly
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
                        partition = partitionList->at(j);
                        maps [ partition ]->at( ++pos.at(partition + 1) ) = i;
                    }
                }
            }
        }

        // set up domains recv communicator maps
        for ( i = 0; i < size; i++ ) {
            this->setProblemComunicatorToRecvArry(this->giveProblemComunicator(i), * maps [ i ]);
            //this->giveDomainComunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
        }

        /*
         #ifdef __VERBOSE_PARALLEL
         * for (i=0; i<size; i++) {
         * fprintf (stderr, "domain %d-%d: domainCommRecvsize is %d\n",rank,i,this->giveDomainComunicator(i)->giveRecvBuff()->giveSize() );
         * printf ("domain %d-%d: reecv map:",rank,i);
         * this->giveDomainComunicator(i)->giveToRecvMap()->printYourself();
         * }
         #endif
         */

        // delete local maps
        for ( i = 0; i < size; i++ ) {
            delete maps [ i ];
        }

        delete maps;

        // to assemble send maps, we must analyze broadcasted remote domain send lists
        // and we must also broadcast our send list.

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Remote Element-cut broadcasting started", rank);
#endif


        CommunicationBuffer commBuff(MPI_COMM_WORLD);
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
        VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Finished reducing receiveBufferSize", rank);
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
                VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Broadcasting own send list", rank);
#endif

                commBuff.packIntArray(domainRecvList);
                result = commBuff.bcast(i);
                if ( result != MPI_SUCCESS ) {
                    _error("setUpCommunicationMaps: commBuff broadcast failed");
                }

#ifdef __VERBOSE_PARALLEL
                VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Broadcasting own send list finished", rank);
#endif
            } else {
#ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_INFO("[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d\n",
                               rank, "PNlDEIDynamicComunicator :: unpackAllData", i);
#endif
                // receive broadcasted lists
                result = commBuff.bcast(i);
                if ( result != MPI_SUCCESS ) {
                    _error("setUpCommunicationMaps: commBuff broadcast failed");
                }

#ifdef __VERBOSE_PARALLEL
                OOFEM_LOG_INFO("[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d finished\n",
                               rank, "PNlDEIDynamicComunicator :: unpackAllData", i);
#endif


                // unpack remote receive list
                if ( !commBuff.unpackIntArray(remoteDomainRecvList) ) {
                    _error("PNlDEIDynamicComunicator::setUpCommunicationMaps: unpack remote receive list failed");
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

                // set send map to i-th problem communicator
                this->setProblemComunicatorToSendArry(this->giveProblemComunicator(i), toSendMap);

                /*
                 #ifdef __VERBOSE_PARALLEL
                 *  fprintf (stderr, "domain %d-%d: domainCommSendsize is %d\n",rank,i,this->giveDomainComunicator(i)->giveSendBuff()->giveSize() );
                 *  printf ("domain %d-%d: send map:",rank,i);
                 *  this->giveDomainComunicator(i)->giveToSendMap()->printYourself();
                 *
                 #endif
                 */

                //this->giveDomainComunicator(i)->setToSendArry (this->engngModel, toSendMap);
            } // end receiving broadcasted lists

#ifdef __VERBOSE_PARALLEL
            VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Receiving broadcasted send maps finished", rank);
#endif
        } // end loop over domains

    } else {
        _error("setUpCommunicationMapsForRemoteElementMode: unknown mode");
    }
}

void
PNlDEIDynamicComunicator :: setUpCommunicationMaps(EngngModel *pm)
{
#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator :: setUpCommunicationMaps", "Setting up communication maps", rank);
#endif


    if ( this->mode == PNlDEIDynamicComunicator__NODE_CUT ) {
        setUpCommunicationMapsForNodeCut(pm);
    } else if ( this->mode == PNlDEIDynamicComunicator__ELEMENT_CUT ) {
        setUpCommunicationMapsForElementCut(pm);
    } else if ( this->mode == PNlDEIDynamicComunicator__REMOTE_ELEMENT_MODE ) {
        setUpCommunicationMapsForRemoteElementMode(pm);
    } else {
        _error("setUpCommunicationMaps: unknown mode");
    }
}

/*
 *
 * //
 * // Old compact version
 * //
 * void
 * PNlDEIDynamicComunicator :: setUpCommunicationMaps (Domain* domain)
 * {
 * int nnodes = domain->giveNumberOfDofManagers();
 * int i, j, partition;
 *
 #ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator :: setUpCommunicationMaps", "Setting up communication maps", rank);
 #endif
 *
 *
 * if (this->mode == PNlDEIDynamicComunicator__NODE_CUT) {
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
 * // if combination node & eleemnt cut can occur, test for shared DofMan mode
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
 * this->setDomainComunicatorToSendArry (this->giveDomainComunicator(i), *maps[i]);
 * this->setDomainComunicatorToRecvArry (this->giveDomainComunicator(i), *maps[i]);
 * //this->giveDomainComunicator(i)->setToSendArry (this->engngModel, *maps[i]);
 * //this->giveDomainComunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
 * }
 *
 * // delete local maps
 * for (i=0; i<size; i++) delete maps[i];
 * delete maps;
 *
 * } else if ((this->mode == PNlDEIDynamicComunicator__ELEMENT_CUT) ||
 *     (this->mode == PNlDEIDynamicComunicator__REMOTE_ELEMENT_MODE)) {
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
 * if (this->mode == PNlDEIDynamicComunicator__ELEMENT_CUT) {
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
 * } else {            // PNlDEIDynamicComunicator__REMOTE_ELEMENT_MODE
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
 * if (this->mode == PNlDEIDynamicComunicator__ELEMENT_CUT) {
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
 * } else {          // PNlDEIDynamicComunicator__REMOTE_ELEMENT_MODE
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
 * this->setDomainComunicatorToRecvArry (this->giveDomainComunicator(i), *maps[i]);
 * //this->giveDomainComunicator(i)->setToRecvArry (this->engngModel, *maps[i]);
 * }
 *
 * // delete local maps
 * for (i=0; i<size; i++) delete maps[i];
 * delete maps;
 *
 * // to assemble send maps, we must analyze broadcasted remote domain send lists
 * // and we must also broadcast our send list.
 *
 #ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Element-cut broadcasting started", rank);
 #endif
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
 #ifdef __USE_MPI
 * result = MPI_Allreduce (&localExpectedSize, &globalRecvSize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
 * if (result != MPI_SUCCESS) _error ("setUpCommunicationMaps: MPI_Allreduce failed");
 #else
 * WARNING: NOT SUPPORTED MESSAGE PARSING LIBRARY
 #endif
 *
 #ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Finished reducing receiveBufferSize", rank);
 #endif
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
 #ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Broadcasting own send list", rank);
 #endif
 *
 *  commBuff.packIntArray (domainRecvList);
 *  result = commBuff.bcast (i);
 *  if (result!= MPI_SUCCESS) _error ("setUpCommunicationMaps: commBuff broadcast failed");
 *
 #ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Broadcasting own send list finished", rank);
 #endif
 *
 * } else {
 *
 #ifdef __VERBOSE_PARALLEL
 *  fprintf(stderr, "\n[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d",
 *      rank,"PNlDEIDynamicComunicator :: unpackAllData", i);
 #endif
 *  // receive broadcasted lists
 *  result = commBuff.bcast (i);
 *  if (result != MPI_SUCCESS) _error ("setUpCommunicationMaps: commBuff broadcast failed");
 *
 #ifdef __VERBOSE_PARALLEL
 *  fprintf(stderr, "\n[process rank %3d]: %-30s: Receiving broadcasted send map from partition %3d finished",
 *      rank,"PNlDEIDynamicComunicator :: unpackAllData", i);
 #endif
 *
 *
 *  // unpack remote receive list
 *  if (!commBuff.unpackIntArray(remoteDomainRecvList))
 *   _error ("PNlDEIDynamicComunicator::setUpCommunicationMaps: unpack remote receive list failed");
 *
 *  if (this->mode == PNlDEIDynamicComunicator__ELEMENT_CUT) {
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
 *  } else {       // PNlDEIDynamicComunicator__REMOTE_ELEMENT_MODE
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
 *  this->setDomainComunicatorToSendArry (this->giveDomainComunicator(i), toSendMap);
 *  //this->giveDomainComunicator(i)->setToSendArry (this->engngModel, toSendMap);
 * } // end receiving broadcasted lists
 #ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("PNlDEIDynamicComunicator::setUpCommunicationMaps", "Receiving broadcasted send maps finished", rank);
 #endif
 *
 * } // end loop over domains
 * } else _error ("setUpCommunicationMaps: unknown mode");
 * }
 */


int
PNlDEIDynamicComunicator :: setProblemComunicatorToSendArry(ProblemComunicator< PNlDEIDynamic > *problemComm, IntArray &map)
{
    if ( ( this->mode == PNlDEIDynamicComunicator__NODE_CUT ) ) {
        sortCommMap(map, & PNlDEIDynamicComunicator :: DofManCmp);
    } else if ( ( this->mode == PNlDEIDynamicComunicator__REMOTE_ELEMENT_MODE ) ||
               ( this->mode == PNlDEIDynamicComunicator__ELEMENT_CUT ) ) {
        sortCommMap(map, & PNlDEIDynamicComunicator :: ElemCmp);
    } else {
        _error("setDomainComunicatorToSendArry: unknown mode");
    }

    problemComm->setToSendArry(engngModel, map, this->mode);
    return 1;
}

int
PNlDEIDynamicComunicator :: setProblemComunicatorToRecvArry(ProblemComunicator< PNlDEIDynamic > *problemComm, IntArray &map)
{
    if ( ( this->mode == PNlDEIDynamicComunicator__NODE_CUT ) ) {
        sortCommMap(map, & PNlDEIDynamicComunicator :: DofManCmp);
    } else if ( ( this->mode == PNlDEIDynamicComunicator__REMOTE_ELEMENT_MODE ) ||
               ( this->mode == PNlDEIDynamicComunicator__ELEMENT_CUT ) ) {
        sortCommMap(map, & PNlDEIDynamicComunicator :: ElemCmp);
    } else {
        _error("setDomainComunicatorToRecvArry: unknown mode");
    }

    problemComm->setToRecvArry(engngModel, map, this->mode);
    return 1;
}



void
PNlDEIDynamicComunicator :: sortCommMap( IntArray &map, int ( PNlDEIDynamicComunicator :: *cmp )( int, int ) )
{
    this->quickSortCommMap(map, 1, map.giveSize(), cmp);
}


void
PNlDEIDynamicComunicator :: quickSortCommMap( IntArray &map, int l, int r, int ( PNlDEIDynamicComunicator :: *cmp )( int, int ) )
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
PNlDEIDynamicComunicator :: quickSortPartition( IntArray &map, int l, int r, int ( PNlDEIDynamicComunicator :: *cmp )( int, int ) )
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
PNlDEIDynamicComunicator :: DofManCmp(int i, int j)
{
    return ( engngModel->giveDomain(1)->giveDofManager(i)->giveGlobalNumber() -
            engngModel->giveDomain(1)->giveDofManager(j)->giveGlobalNumber() );
}
int
PNlDEIDynamicComunicator :: ElemCmp(int i, int j)
{
    return ( engngModel->giveDomain(1)->giveElement(i)->giveGlobalNumber() -
            engngModel->giveDomain(1)->giveElement(j)->giveGlobalNumber() );
}

} // end namespace oofem
#endif
