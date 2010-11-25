/* $Header: /home/cvs/bp/oofem/sm/src/feticommunicator.C,v 1.2.4.1 2004/04/05 15:19:46 bp Exp $ */
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

#include "feticommunicator.h"
#include "engngm.h"
#include "intarray.h"
#include "dofmanager.h"

#ifdef __USE_MPI
 #ifndef __MAKEDEPEND
  #include "mpi.h"
 #endif
#endif

namespace oofem {
FETICommunicator :: FETICommunicator(EngngModel *emodel, CommunicatorBuff *b, int rank, int size) :
    Communicator(emodel, b, rank, size)
{
    if ( rank != 0 ) {
        _error("FETICommunicator: bad rank number, expected rank 0 for master");
    }
}


FETICommunicator :: ~FETICommunicator()
{ }


void
FETICommunicator :: setUpCommunicationMaps(EngngModel *pm)
{
    int i, j, l, maxRec;
    int globaldofmannum, localNumber, ndofs, npart;
    int numberOfBoundaryDofMans;
    int source, tag;
    IntArray numberOfPartitionBoundaryDofMans(size);
    StaticCommunicationBuffer commBuff(MPI_COMM_WORLD);
    EModelDefaultEquationNumbering dn;
    // FETIBoundaryDofManager *dofmanrec;
    // Map containing boundary dof managers records, the key is corresponding global number
    // value is corresponding local master dof manager number
    map< int, int, less< int > >BoundaryDofManagerMap;
    // communication maps of slaves
    IntArray **commMaps = new IntArray * [ size ];
    // location array
    IntArray locNum;
    Domain *domain = pm->giveDomain(1);

    // check if receiver is master
    if ( this->rank != 0 ) {
        _error("FETICommunicator::setUpCommunicationMaps : rank 0 (master) expected as receiver");
    }

    // resize receive buffer
    commBuff.resize( commBuff.givePackSize(MPI_INT, 1) );

    //
    // receive data
    //
    for ( i = 1; i < size; i++ ) {
        commBuff.iRecv(MPI_ANY_SOURCE, FETICommunicator :: NumberOfBoundaryDofManagersMsg);
        while ( !commBuff.testCompletion(source, tag) ) {
            ;
        }

        // unpack data
        commBuff.unpackInt(j);
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received data from partition %3d (received %d)\n",
                        rank, "FETICommunicator :: setUpCommunicationMaps : received number of boundary dofMans", source, j);
#endif
        numberOfPartitionBoundaryDofMans.at(source + 1) = j;
        commBuff.init();
    }

    MPI_Barrier(MPI_COMM_WORLD);


    // determine the total number of boundary dof managers at master
    int nnodes = domain->giveNumberOfDofManagers();
    j = 0;
    for ( i = 1; i <= nnodes; i++ ) {
        if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
            j++;
        }
    }

    numberOfPartitionBoundaryDofMans.at(1) = j;

    //
    // receive list of bounadry dof managers with corresponding number of dofs from each partition
    //

    // resize the receive buffer to fit all messages
    maxRec = 0;
    for ( i = 0; i < size; i++ ) {
        if ( numberOfPartitionBoundaryDofMans.at(i + 1) > maxRec ) {
            maxRec = numberOfPartitionBoundaryDofMans.at(i + 1);
        }
    }

    commBuff.resize( 2 * maxRec * commBuff.givePackSize(MPI_INT, 1) );
    // resize communication maps acordingly
    for ( i = 0; i < size; i++ ) {
        j = numberOfPartitionBoundaryDofMans.at(i + 1);
        commMaps [ i ] = new IntArray(j);
    }


    // add local master contribution first
    // loop over all dofmanager data received
    i = 0;
    for ( j = 1; j <= numberOfPartitionBoundaryDofMans.at(1); j++ ) {
        // fing next shared dofman
        while ( !( domain->giveDofManager(++i)->giveParallelMode() == DofManager_shared ) ) {
            ;
        }

        globaldofmannum = domain->giveDofManager(i)->giveGlobalNumber();
        domain->giveDofManager(i)->giveCompleteLocationArray(locNum, dn);
        ndofs = 0;
        for ( l = 1; l <= locNum.giveSize(); l++ ) {
            if ( locNum.at(l) ) {
                ndofs++;
            }
        }

        // add corresponding entry to master map of boundary dof managers
        if ( ( localNumber = BoundaryDofManagerMap [ globaldofmannum ] ) == 0 ) { // no local counterpart exist
            // create it
            boundaryDofManList.push_back( FETIBoundaryDofManager(globaldofmannum, 0, ndofs) );
            // remember the local number; actual position in vector is localNumber-1
            localNumber = BoundaryDofManagerMap [ globaldofmannum ] = ( boundaryDofManList.size() );
            boundaryDofManList.back().addPartition(0);
        } else { // update the corresponding record
            boundaryDofManList [ localNumber - 1 ].addPartition(0);
            if ( boundaryDofManList [ localNumber - 1 ].giveNumberOfDofs() != ndofs ) {
                _error("FETICommunicator :: setUpCommunicationMaps : ndofs size mismatch");
            }
        }

        // remember communication map for particular partition
        commMaps [ 0 ]->at(j) = localNumber;
    }

    //
    // receive data from slave partitions
    //

    for ( i = 1; i < size; i++ ) {
        commBuff.iRecv(MPI_ANY_SOURCE, FETICommunicator :: BoundaryDofManagersRecMsg);
        while ( !commBuff.testCompletion(source, tag) ) {
            ;
        }

        // unpack data
#ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received data from partition %3d\n",
                        rank, "FETICommunicator :: setUpCommunicationMaps : received boundary dofMans records", source);
#endif

        // loop over all dofmanager data received
        for ( j = 1; j <= numberOfPartitionBoundaryDofMans.at(source + 1); j++ ) {
            commBuff.unpackInt(globaldofmannum);
            commBuff.unpackInt(ndofs);

            // add corresponding entry to master map of boundary dof managers
            if ( ( localNumber = BoundaryDofManagerMap [ globaldofmannum ] ) == 0 ) { // no local counterpart exist
                // create it
                boundaryDofManList.push_back( FETIBoundaryDofManager(globaldofmannum, 0, ndofs) );
                // remember the local number; actual position in vector is localNumber-1
                localNumber = BoundaryDofManagerMap [ globaldofmannum ] = ( boundaryDofManList.size() );
                boundaryDofManList.back().addPartition(source);
            } else { // update the corresponding record
                boundaryDofManList [ localNumber - 1 ].addPartition(source);
                if ( boundaryDofManList [ localNumber - 1 ].giveNumberOfDofs() != ndofs ) {
                    _error("FETICommunicator :: setUpCommunicationMaps : ndofs size mismatch");
                }
            }

            // remember communication map for particular partition
            commMaps [ source ]->at(j) = localNumber;
        }

        commBuff.init();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //
    // assign code numbers to boundary dofs
    //
    numberOfEquations = 0;
    numberOfBoundaryDofMans = boundaryDofManList.size();
    for ( i = 1; i <= numberOfBoundaryDofMans; i++ ) {
        npart = boundaryDofManList [ i - 1 ].giveNumberOfSharedPartitions();
        ndofs = boundaryDofManList [ i - 1 ].giveNumberOfDofs();
        //neqs  = npart*ndofs;
        boundaryDofManList [ i - 1 ].setCodeNumbers(numberOfEquations); // updates numberOfEquations
    }

    // store the commMaps
    for ( i = 0; i < size; i++ ) {
        if ( i != 0 ) {
            this->giveProcessCommunicator(i)->setToSendArry(engngModel, * commMaps [ i ], 0);
            this->giveProcessCommunicator(i)->setToRecvArry(engngModel, * commMaps [ i ], 0);
        } else {
            masterCommMap = * commMaps [ i ];
        }

        delete commMaps [ i ];
    }

    delete commMaps;

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("FETICommunicator::setUpCommunicationMaps", "communication maps setup finished", rank);
#endif
}



// ==================================================================================================== //

/*
 * int
 * PNlDEIDynamicCommunicator :: setDomainCommunicatorToSendArry (DomainCommunicator<PNlDEIDynamic>* domainComm, IntArray& map)
 * {
 * if ((this->mode == PNlDEIDynamicCommunicator__NODE_CUT))
 * sortCommMap (map, &PNlDEIDynamicCommunicator::DofManCmp);
 * else if ((this->mode == PNlDEIDynamicCommunicator__REMOTE_ELEMENT_MODE)
 || (this->mode == PNlDEIDynamicCommunicator__ELEMENT_CUT))
 * sortCommMap (map, &PNlDEIDynamicCommunicator::ElemCmp);
 * else _error ("setDomainCommunicatorToSendArry: unknown mode");
 *
 * domainComm->setToSendArry (engngModel, map, this->mode);
 * return 1;
 * }
 *
 * int
 * PNlDEIDynamicCommunicator :: setDomainCommunicatorToRecvArry (DomainCommunicator<PNlDEIDynamic>* domainComm, IntArray& map)
 * {
 * if ((this->mode == PNlDEIDynamicCommunicator__NODE_CUT))
 * sortCommMap (map, &PNlDEIDynamicCommunicator::DofManCmp);
 * else if ((this->mode == PNlDEIDynamicCommunicator__REMOTE_ELEMENT_MODE) ||
 *    (this->mode == PNlDEIDynamicCommunicator__ELEMENT_CUT))
 * sortCommMap (map, &PNlDEIDynamicCommunicator::ElemCmp);
 * else _error ("setDomainCommunicatorToRecvArry: unknown mode");
 *
 * domainComm->setToRecvArry (engngModel, map, this->mode);
 * return 1;
 * }
 *
 *
 *
 * void
 * PNlDEIDynamicCommunicator :: sortCommMap (IntArray& map, int (PNlDEIDynamicCommunicator::*cmp) (int,int))
 * {
 * this->quickSortCommMap (map, 1, map.giveSize(), cmp);
 * }
 *
 *
 * void
 * PNlDEIDynamicCommunicator :: quickSortCommMap (IntArray& map, int l, int r, int (PNlDEIDynamicCommunicator::*cmp) (int,int))
 * {
 * if (r<=l) return;
 * int i = quickSortPartition (map, l, r, cmp);
 * quickSortCommMap (map, l, i-1, cmp);
 * quickSortCommMap (map, i+1, r, cmp);
 * }
 *
 *#define GLOBNUM(locnum) (engngModel->giveDomain()->giveDofManager(locnum)->giveGlobalNumber())
 *
 *
 * int
 * PNlDEIDynamicCommunicator :: quickSortPartition (IntArray& map, int l, int r, int (PNlDEIDynamicCommunicator::*cmp) (int,int))
 * {
 * int i=l-1, j=r;
 * int v = map.at(r);
 * int swap;
 *
 * for (;;) {
 * while (((this->*cmp) (map.at(++i), v)) < 0 );
 * while (((this->*cmp) (v, map.at(--j))) < 0 ) if (j==l) break;
 * if (i >= j) break;
 * swap = map.at(i); map.at(i) = map.at(j); map.at(j) = swap;
 * }
 * swap = map.at(i); map.at(i) = map.at(r); map.at(r) = swap;
 * return i;
 * }
 *
 *#undef GLOBNUM
 *
 * int
 * PNlDEIDynamicCommunicator :: DofManCmp (int i, int j)
 * {
 * return (engngModel->giveDomain()->giveDofManager(i)->giveGlobalNumber() -
 *   engngModel->giveDomain()->giveDofManager(j)->giveGlobalNumber());
 * }
 * int
 * PNlDEIDynamicCommunicator :: ElemCmp (int i, int j)
 * {
 * return (engngModel->giveDomain()->giveElement(i)->giveGlobalNumber() -
 *   engngModel->giveDomain()->giveElement(j)->giveGlobalNumber());
 * }
 */
} // end namespace oofem
#endif
