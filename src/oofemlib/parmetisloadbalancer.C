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

#include "parmetisloadbalancer.h"
#include "domain.h"
#include "engngm.h"
#include "element.h"
#include "dofmanager.h"
#include "connectivitytable.h"
#include "error.h"
#include "parallel.h"
#include "processcomm.h"
#include "communicator.h"
#include "classfactory.h"

#include <set>
#include <stdlib.h>

namespace oofem {
//#define ParmetisLoadBalancer_DEBUG_PRINT

REGISTER_LoadBalancer(ParmetisLoadBalancer);

ParmetisLoadBalancer :: ParmetisLoadBalancer(Domain *d) : LoadBalancer(d)
{
    elmdist = NULL;
    tpwgts  = NULL;
}

ParmetisLoadBalancer :: ~ParmetisLoadBalancer()
{
    if ( elmdist ) {
        delete[] elmdist;
    }

    if ( tpwgts ) {
        delete[] tpwgts;
    }
}


void
ParmetisLoadBalancer :: calculateLoadTransfer()
{
    idx_t *eind, *eptr, *xadj, *adjncy, *vwgt, *vsize;
    idx_t *part;
    int i, nlocalelems, eind_size, nelem = domain->giveNumberOfElements();
    int ndofman, idofman, numflag, ncommonnodes, options [ 4 ], ie, nproc;
    int edgecut, wgtflag, ncon;
    real_t ubvec [ 1 ], itr;
    Element *ielem;
    MPI_Comm communicator = MPI_COMM_WORLD;
    LoadBalancerMonitor *lbm = domain->giveEngngModel()->giveLoadBalancerMonitor();

    nproc = domain->giveEngngModel()->giveNumberOfProcesses();
    // init parmetis element numbering
    this->initGlobalParmetisElementNumbering();
    // prepare data structures for ParMETIS_V3_Mesh2Dual
    // count the size of eind array
    eind_size = 0;
    nlocalelems = 0;
    for ( i = 1; i <= nelem; i++ ) {
        ielem = domain->giveElement(i);
        if ( ielem->giveParallelMode() == Element_local ) {
            nlocalelems++;
            eind_size += ielem->giveNumberOfDofManagers();
        }
    }

    // allocate eind and eptr arrays
    eind = new idx_t [ eind_size ];
    eptr = new idx_t [ nlocalelems + 1 ];
    if ( ( eind == NULL ) || ( eptr == NULL ) ) {
        OOFEM_ERROR("failed to allocate eind and eptr arrays");
    }

    // fill in the eind and eptr (mesh graph)
    int eind_pos = 0, eptr_pos = 0;
    for ( i = 1; i <= nelem; i++ ) {
        ielem = domain->giveElement(i);
        if ( ielem->giveParallelMode() == Element_local ) {
            eptr [ eptr_pos ] = eind_pos;
            ndofman = ielem->giveNumberOfDofManagers();
            for ( idofman = 1; idofman <= ndofman; idofman++ ) {
                eind [ eind_pos++ ] = ielem->giveDofManager(idofman)->giveGlobalNumber() - 1;
            }

            eptr_pos++;
        }
    }

    // last rec
    eptr [ nlocalelems ] = eind_pos;

    // call ParMETIS_V3_Mesh2Dual to construct dual graph (in parallel)
    // dual graph: elements are vertices; element edges are graph edges
    // this is necessary, since cut runs through graph edges
    numflag = 0;
    ncommonnodes = 2;
    ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, & numflag, & ncommonnodes, & xadj, & adjncy, & communicator);

 #ifdef ParmetisLoadBalancer_DEBUG_PRINT
    int myrank = domain->giveEngngModel()->giveRank();
    // DEBUG PRINT
    fprintf(stderr, "[%d] xadj:", myrank);
    for ( i = 0; i <= nlocalelems; i++ ) {
        fprintf(stderr, " %d", xadj [ i ]);
    }

    fprintf(stderr, "\n[%d] adjncy:", myrank);
    for ( i = 0; i < xadj [ nlocalelems ]; i++ ) {
        fprintf(stderr, " %d", adjncy [ i ]);
    }

    fprintf(stderr, "\n");
 #endif


    // setup imbalance tolerance for each vertex weight - ubvec param
    ubvec [ 0 ] = 1.05;
    // setup options array
    options [ 0 ] = 1; // set to zero for default
    options [ 1 ] = 1; // get timings
    options [ 2 ] = 15; // random seed
    options [ 3 ] = 1; // sub-domains and processors are coupled
    // set ratio of inter-proc communication compared to data redistribution time
    itr = 1000.0;
    // set partition weights by quering load balance monitor
    const FloatArray &_procweights = lbm->giveProcessorWeights();
    if ( tpwgts == NULL ) {
        if ( ( tpwgts = new real_t [ nproc ] ) == NULL ) {
            OOFEM_ERROR("failed to allocate tpwgts");
        }
    }

    for ( i = 0; i < nproc; i++ ) {
        tpwgts [ i ] = _procweights(i);
    }

    /*
     * // log processor weights
     * OOFEM_LOG_RELEVANT ("[%d] ParmetisLoadBalancer: proc weights: ", myrank);
     * for (i=0; i<nproc; i++) OOFEM_LOG_RELEVANT ("%4.3f ",tpwgts[i]);
     * OOFEM_LOG_RELEVANT ("\n");
     */

    // obtain vertices weights (element weights) representing relative computational cost
    if ( ( vwgt = new idx_t [ nlocalelems ] ) == NULL ) {
        OOFEM_ERROR("failed to allocate vwgt");
    }

    if ( ( vsize = new idx_t [ nlocalelems ] ) == NULL ) {
        OOFEM_ERROR("failed to allocate vsize");
    }

    for ( ie = 0, i = 0; i < nelem; i++ ) {
        ielem = domain->giveElement(i + 1);
        if ( ielem->giveParallelMode() == Element_local ) {
            vwgt [ ie ]    = ( int ) ( ielem->predictRelativeComputationalCost() * 100.0 );
            vsize [ ie++ ] = 1; //ielem->predictRelativeRedistributionCost();
        }
    }

    wgtflag = 2;
    numflag = 0;
    ncon = 1;
    if ( ( part = new idx_t [ nlocalelems ] ) == NULL ) {
        OOFEM_ERROR("failed to allocate part");
    }

    // call ParMETIS balancing routineParMETIS_V3_AdaptiveRepart
    ParMETIS_V3_AdaptiveRepart(elmdist, xadj, adjncy, vwgt, vsize, NULL, & wgtflag, & numflag, & ncon, & nproc,
                               tpwgts, ubvec, & itr, options, & edgecut, part, & communicator);

    // part contains partition vector for local elements on receiver
    // we need to map it to domain elements (this is not the same, since
    // domain may contain not only its local elements but remote elements as well)
    int loc_num = 0;
    this->elementPart.resize(nelem);
    for ( i = 1; i <= nelem; i++ ) {
        ielem = domain->giveElement(i);
        if ( ielem->giveParallelMode() == Element_local ) {
            this->elementPart.at(i) = part [ loc_num++ ];
        } else {
            // we can not say anything about remote elements; this information is available on partition
            // that has its local counterpart
            this->elementPart.at(i) = -1;
        }
    }

    if ( part ) {
        delete[] part;
    }

 #ifdef ParmetisLoadBalancer_DEBUG_PRINT
    // debug
    fprintf(stderr, "[%d] edgecut: %d elementPart:", myrank, edgecut);
    for ( i = 1; i <= nelem; i++ ) {
        fprintf( stderr, " %d", elementPart.at(i) );
    }

    fprintf(stderr, "\n");
 #endif

    // delete allocated xadj, adjncy arrays by ParMETIS
    delete[] eind;
    delete[] eptr;
    delete[] vwgt;
    delete[] vsize;
    free(xadj);
    free(adjncy);

    this->labelDofManagers();
}

void
ParmetisLoadBalancer :: initGlobalParmetisElementNumbering()
{
    int nproc = domain->giveEngngModel()->giveNumberOfProcesses();
    int myrank = domain->giveEngngModel()->giveRank();
    IntArray procElementCounts(nproc);

    //if (procElementCounts) delete procElementCounts;
    if ( elmdist == NULL ) {
        elmdist = new idx_t [ nproc + 1 ];
        if ( elmdist == NULL ) {
            OOFEM_ERROR("failed to allocate elmdist array");
        }
    }

    // determine number of local elements for the receiver
    int i, nlocelem = 0, nelem = domain->giveNumberOfElements();
    int globnum;

    for ( i = 1; i <= nelem; i++ ) {
        if ( domain->giveElement(i)->giveParallelMode() == Element_local ) {
            nlocelem++;
        }
    }

    procElementCounts(myrank) = nlocelem;

    MPI_Allgather(& nlocelem, 1, MPI_INT, procElementCounts.givePointer(), 1, MPI_INT, MPI_COMM_WORLD);
    elmdist [ 0 ] = 0;
    for ( i = 0; i < nproc; i++ ) {
        elmdist [ i + 1 ] = elmdist [ i ] + procElementCounts(i);
    }

    // we need to number elements sequentially on each partition (and we start from rank 0)
    // compute local offset
    myGlobNumOffset = 0;
    for ( i = 0; i < myrank; i++ ) {
        myGlobNumOffset += procElementCounts(i);
    }

    /* assemble maps of local numbering
     * map is necessary since we may have remote elements that are not
     * part of local domain for load balancing purposes
     */
    globnum = myGlobNumOffset + 1;
    lToGMap.resize(nelem);
    gToLMap.resize(nelem);
    for ( i = 1; i <= nelem; i++ ) {
        if ( domain->giveElement(i)->giveParallelMode() == Element_local ) {
            lToGMap.at(i) = globnum;
            gToLMap.at(globnum - myGlobNumOffset) = i;
            globnum++;
        } else {
            lToGMap.at(i) = 0;
        }
    }
}

void
ParmetisLoadBalancer :: labelDofManagers()
{
    int idofman, ndofman = domain->giveNumberOfDofManagers();
    ConnectivityTable *ct = domain->giveConnectivityTable();
    const IntArray *dofmanconntable;
    DofManager *dofman;
    Element *ielem;
    dofManagerParallelMode dmode;
    std :: set< int, std :: less< int > >__dmanpartitions;
    int myrank = domain->giveEngngModel()->giveRank();
    int nproc = domain->giveEngngModel()->giveNumberOfProcesses();
    int ie, npart;

    // resize label array
    dofManState.resize(ndofman);
    dofManState.zero();
    // resize dof man partitions
    dofManPartitions.clear();
    dofManPartitions.resize(ndofman);

 #ifdef ParmetisLoadBalancer_DEBUG_PRINT
    int _cols = 0;
    fprintf(stderr, "[%d] DofManager labels:\n", myrank);
 #endif

    // loop over local dof managers
    for ( idofman = 1; idofman <= ndofman; idofman++ ) {
        dofman = domain->giveDofManager(idofman);
        dmode = dofman->giveParallelMode();
        if ( ( dmode == DofManager_local ) || ( dmode == DofManager_shared ) ) {
            dofmanconntable = ct->giveDofManConnectivityArray(idofman);
            __dmanpartitions.clear();
            for ( ie = 1; ie <= dofmanconntable->giveSize(); ie++ ) {
                ielem = domain->giveElement( dofmanconntable->at(ie) );
                // assemble list of partitions sharing idofman dofmanager
                // set is used to include possibly repeated partition only once
                if ( ielem->giveParallelMode() == Element_local ) {
                    __dmanpartitions.insert( giveElementPartition( dofmanconntable->at(ie) ) );
                }
            }

            npart = __dmanpartitions.size();
            dofManPartitions [ idofman - 1 ].resize( __dmanpartitions.size() );
            int i = 1;
            for ( auto &dm: __dmanpartitions ) {
                dofManPartitions [ idofman - 1 ].at(i++) = dm;
            }
        }
    }

    // handle master slave links between dofmans (master and slave required on same partition)
    this->handleMasterSlaveDofManLinks();


    /* Exchange new partitions for shared nodes */
    CommunicatorBuff cb(nproc, CBT_dynamic);
    Communicator com(domain->giveEngngModel(), &cb, myrank, nproc, CommMode_Dynamic);
    com.packAllData(this, & ParmetisLoadBalancer :: packSharedDmanPartitions);
    com.initExchange(SHARED_DOFMAN_PARTITIONS_TAG);
    com.unpackAllData(this, & ParmetisLoadBalancer :: unpackSharedDmanPartitions);
    com.finishExchange();

    /* label dof managers */
    for ( idofman = 1; idofman <= ndofman; idofman++ ) {
        dofman = domain->giveDofManager(idofman);
        dmode = dofman->giveParallelMode();
        npart = dofManPartitions [ idofman - 1 ].giveSize();
        if ( ( dmode == DofManager_local ) || ( dmode == DofManager_shared ) ) {
            // determine its state after balancing -> label
            dofManState.at(idofman) = this->determineDofManState(idofman, myrank, npart, & dofManPartitions [ idofman - 1 ]);
        } else {
            dofManState.at(idofman) = DM_NULL;
        }
    }


 #ifdef ParmetisLoadBalancer_DEBUG_PRINT
    for ( idofman = 1; idofman <= ndofman; idofman++ ) {
        fprintf(stderr, " | %d: ", idofman);
        if ( dofManState.at(idofman) == DM_NULL ) {
            fprintf(stderr, "NULL  ");
        } else if ( dofManState.at(idofman) == DM_Local ) {
            fprintf(stderr, "Local ");
        } else if ( dofManState.at(idofman) == DM_Shared ) {
            fprintf(stderr, "Shared");
        } else if ( dofManState.at(idofman) == DM_Remote ) {
            fprintf(stderr, "Remote");
        } else {
            fprintf(stderr, "Unknown");
        }

        //else if (dofManState.at(idofman) == DM_SharedExclude)fprintf (stderr, "ShdExc");
        //else if (dofManState.at(idofman) == DM_SharedNew)    fprintf (stderr, "ShdNew");
        //else if (dofManState.at(idofman) == DM_SharedUpdate) fprintf (stderr, "ShdUpd");

        if ( ( ( ++_cols % 4 ) == 0 ) || ( idofman == ndofman ) ) {
            fprintf(stderr, "\n");
        }
    }

 #endif
}

int
ParmetisLoadBalancer :: determineDofManState(int idofman, int myrank, int npart, IntArray *dofManPartitions)
{
    dofManagerParallelMode dmode = domain->giveDofManager(idofman)->giveParallelMode();
    int answer = DM_Local;

    if ( ( dmode == DofManager_local ) || ( dmode == DofManager_shared ) ) {
        if ( ( npart == 1 ) && ( dofManPartitions->at(1) == myrank ) ) {
            // local remains local
            answer = DM_Local;
        } else if ( npart == 1 ) {
            // local goes to remote partition
            answer = DM_Remote;
        } else { // npart > 1
            // local becomes newly shared
            answer = DM_Shared;
        }
    } else {
        answer = DM_NULL;
    }

    /*
     * if (dmode == DofManager_local) {
     * if ((npart == 1) && (dofManPartitions->at(1) == myrank)) {
     *  // local remains local
     *  answer = DM_Local;
     * } else if (npart == 1) {
     *  // local goes to remote partition
     *  answer = DM_Remote;
     * } else { // npart > 1
     *  // local becomes newly shared
     *  answer = DM_SharedNew;
     * }
     * } else if (dmode == DofManager_shared) {
     * // compare old and new partition list
     * int i, _same = true, containsMyRank = dofManPartitions->findFirstIndexOf (myrank);
     * const IntArray* oldpart = domain->giveDofManager(idofman)->givePartitionList();
     * for (i=1; i<=dofManPartitions->giveSize(); i++) {
     *  if ((dofManPartitions->at(i)!= myrank) &&
     *      (!oldpart->findFirstIndexOf(dofManPartitions->at(i)))) {
     *    _same=false; break;
     *  }
     * }
     * if (_same && containsMyRank) {
     *  answer = DM_Shared;
     * } else if (containsMyRank) {
     *  answer = DM_SharedUpdate;
     * } else { // !containsMyRank
     *  answer = DM_SharedExclude;
     * }
     * } else {
     * answer = DM_NULL;
     * }
     */
    return answer;
}


LoadBalancer :: DofManMode
ParmetisLoadBalancer :: giveDofManState(int idofman)
{
    return ( LoadBalancer :: DofManMode ) dofManState.at(idofman);
}


IntArray *
ParmetisLoadBalancer :: giveDofManPartitions(int idofman)
{
    return & dofManPartitions [ idofman - 1 ];
}

int
ParmetisLoadBalancer :: giveElementPartition(int ielem)
{
    return elementPart.at(ielem);
}

int
ParmetisLoadBalancer :: packSharedDmanPartitions(ProcessCommunicator &pc)
{
    int myrank = domain->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();
    int ndofman, idofman;
    DofManager *dofman;

    if ( iproc == myrank ) {
        return 1;                // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
    // loop over dofManagers and pack shared dofMan data
    ndofman = domain->giveNumberOfDofManagers();
    for ( idofman = 1; idofman <= ndofman; idofman++ ) {
        dofman = domain->giveDofManager(idofman);
        // test if iproc is in list of existing shared partitions
        if ( ( dofman->giveParallelMode() == DofManager_shared ) &&
            ( dofman->givePartitionList()->findFirstIndexOf(iproc) ) ) {
            // send new partitions to remote representation
            // fprintf (stderr, "[%d] sending shared plist of %d to [%d]\n", myrank, dofman->giveGlobalNumber(), iproc);
            pcbuff->write( dofman->giveGlobalNumber() );
            this->giveDofManPartitions(idofman)->storeYourself(*pcbuff);
        }
    }

    pcbuff->write((int)PARMETISLB_END_DATA);
    return 1;
}

int
ParmetisLoadBalancer :: unpackSharedDmanPartitions(ProcessCommunicator &pc)
{
    int myrank = domain->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();
    int _globnum, _locnum;
    IntArray _partitions;

    if ( iproc == myrank ) {
        return 1;                // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
    // init domain global2local map
    domain->initGlobalDofManMap();

    pcbuff->read(_globnum);
    // unpack dofman data
    while ( _globnum != PARMETISLB_END_DATA ) {
        _partitions.restoreYourself(*pcbuff);
        if ( ( _locnum = domain->dofmanGlobal2Local(_globnum) ) ) {
            this->addSharedDofmanPartitions(_locnum, _partitions);
        } else {
            OOFEM_ERROR("internal error, unknown global dofman %d", _globnum);
        }

        /*
         * fprintf (stderr,"[%d] Received shared plist of %d ", myrank, _globnum);
         * for (int _i=1; _i<=dofManPartitions[_locnum-1].giveSize(); _i++)
         * fprintf (stderr,"%d ", dofManPartitions[_locnum-1].at(_i));
         * fprintf (stderr,"\n");
         */
        pcbuff->read(_globnum);
    }

    return 1;
}


void ParmetisLoadBalancer :: addSharedDofmanPartitions(int _locnum, IntArray _partitions)
{
    for ( int part: _partitions ) {
        dofManPartitions [ _locnum - 1 ].insertOnce( part );
    }
}

void ParmetisLoadBalancer :: handleMasterSlaveDofManLinks()
{
    int idofman, ndofman = domain->giveNumberOfDofManagers();
    DofManager *dofman;
    //int myrank = domain->giveEngngModel()->giveRank();
    int __i, __j, __partition, _master;
    bool isSlave;
    IntArray slaveMastersDofMans;

    /*
     * We assume that in the old partitioning, the master and slave consistency was assured. This means that master is presented
     * on the same partition as slave. The master can be local (then all slaves are local) or master is shared (then slaves are on
     * partitions sharing the master).
     *
     * If master was local, then its new partitioning can be locally resolved (as all slaves were local).
     * If the master was shared, the new partitioning of master has to be communicated between old sharing partitions.
     */
    // handle master slave links between dofmans (master and slave required on same partition)

    for ( idofman = 1; idofman <= ndofman; idofman++ ) {
        dofman = domain->giveDofManager(idofman);
        isSlave = dofman->hasAnySlaveDofs();

        if ( isSlave ) {
            // ok have a look on its masters
            dofman->giveMasterDofMans(slaveMastersDofMans);
            for ( __i = 1; __i <= slaveMastersDofMans.giveSize(); __i++ ) {
                // loop over all slave masters
                _master = slaveMastersDofMans.at(__i);

                // now loop over all slave new partitions and annd then to master's partitions
                for ( __j = 1; __j <= dofManPartitions [ idofman - 1 ].giveSize(); __j++ ) {
                    __partition = dofManPartitions [ idofman - 1 ].at(__j);
                    // add slave partition to master
                    dofManPartitions [ _master - 1 ].insertOnce(__partition);
                }
            }
        }
    }
}

} // end namespace oofem
