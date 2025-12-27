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

#include "sm/FETISolver/feticommunicator.h"
#include "engngm.h"
#include "intarray.h"
#include "dofmanager.h"
#include "unknownnumberingscheme.h"
#include "domain.h"

#ifdef __USE_MPI
 #include <mpi.h>
#endif

namespace oofem {
FETICommunicator :: FETICommunicator(EngngModel *emodel, CommunicatorBuff *b, int rank, int size) :
    Communicator(emodel, b, rank, size)
{
    if ( rank != 0 ) {
        OOFEM_ERROR("bad rank number, expected rank 0 for master");
    }
}


FETICommunicator :: ~FETICommunicator()
{ }


void
FETICommunicator :: setUpCommunicationMaps(EngngModel *pm)
{
    int i, j, maxRec;
    int globaldofmannum, localNumber, ndofs;
    int numberOfBoundaryDofMans;
    int source, tag;
    int size = this->processComms.size();
    IntArray numberOfPartitionBoundaryDofMans(size);
    StaticCommunicationBuffer commBuff(MPI_COMM_WORLD);
    EModelDefaultEquationNumbering dn;
    // FETIBoundaryDofManager *dofmanrec;
    // Map containing boundary dof managers records, the key is corresponding global number
    // value is corresponding local master dof manager number
    map< int, int, less< int > >BoundaryDofManagerMap;
    // communication maps of slaves
    std::vector<IntArray> commMaps(size);
    // location array
    IntArray locNum;
    Domain *domain = pm->giveDomain(1);

    // check if receiver is master
    if ( this->rank != 0 ) {
        OOFEM_ERROR("rank 0 (master) expected as receiver");
    }

    // resize receive buffer
    commBuff.resize( commBuff.givePackSizeOfInt(1) );

    //
    // receive data
    //
    for ( int i = 1; i < size; i++ ) {
        commBuff.iRecv(MPI_ANY_SOURCE, FETICommunicator :: NumberOfBoundaryDofManagersMsg);
        while ( !commBuff.testCompletion(source, tag) ) {
            ;
        }

        // unpack data
        commBuff.read(j);
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
    for ( int i = 1; i <= nnodes; i++ ) {
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
    for ( int i = 0; i < size; i++ ) {
        if ( numberOfPartitionBoundaryDofMans.at(i + 1) > maxRec ) {
            maxRec = numberOfPartitionBoundaryDofMans.at(i + 1);
        }
    }

    commBuff.resize( 2 * maxRec * commBuff.givePackSizeOfInt(1) );
    // resize communication maps acordingly
    for ( int i = 0; i < size; i++ ) {
        j = numberOfPartitionBoundaryDofMans.at(i + 1);
        commMaps [ i ].resize(j);
    }


    // add local master contribution first
    // loop over all dofmanager data received
    i = 0;
    for ( int j = 1; j <= numberOfPartitionBoundaryDofMans.at(1); j++ ) {
        // fing next shared dofman
        while ( !( domain->giveDofManager(++i)->giveParallelMode() == DofManager_shared ) ) {
            ;
        }

        globaldofmannum = domain->giveDofManager(i)->giveGlobalNumber();
        domain->giveDofManager(i)->giveCompleteLocationArray(locNum, dn);
        ndofs = 0;
        for ( int l = 1; l <= locNum.giveSize(); l++ ) {
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
                OOFEM_ERROR("ndofs size mismatch");
            }
        }

        // remember communication map for particular partition
        commMaps [ 0 ].at(j) = localNumber;
    }

    //
    // receive data from slave partitions
    //

    for ( int i = 1; i < size; i++ ) {
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
        for ( int j = 1; j <= numberOfPartitionBoundaryDofMans.at(source + 1); j++ ) {
            commBuff.read(globaldofmannum);
            commBuff.read(ndofs);

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
                    OOFEM_ERROR("ndofs size mismatch");
                }
            }

            // remember communication map for particular partition
            commMaps [ source ].at(j) = localNumber;
        }

        commBuff.init();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    //
    // assign code numbers to boundary dofs
    //
    numberOfEquations = 0;
    numberOfBoundaryDofMans = boundaryDofManList.size();
    for ( int i = 1; i <= numberOfBoundaryDofMans; i++ ) {
        boundaryDofManList [ i - 1 ].setCodeNumbers(numberOfEquations); // updates numberOfEquations
    }

    // store the commMaps
    for ( int i = 0; i < size; i++ ) {
        if ( i != 0 ) {
            this->giveProcessCommunicator(i)->setToSendArry(engngModel, commMaps [ i ], 0);
            this->giveProcessCommunicator(i)->setToRecvArry(engngModel, commMaps [ i ], 0);
        } else {
            masterCommMap = commMaps [ i ];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("FETICommunicator::setUpCommunicationMaps", "communication maps setup finished", rank);
#endif
}
} // end namespace oofem
