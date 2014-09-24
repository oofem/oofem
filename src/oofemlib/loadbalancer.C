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

#include "loadbalancer.h"
#include "domain.h"
#include "engngm.h"
#include "timer.h"
#include "mathfem.h"
#include "timestep.h"
#include "floatarray.h"
#include "classfactory.h"
#include "element.h"

#include "parallel.h"
#include "processcomm.h"
#include "datastream.h"
#include "communicator.h"
#include "domaintransactionmanager.h"
#include "nonlocalmatwtp.h"

namespace oofem {
#define LoadBalancer_debug_print 0

REGISTER_LoadBalancerMonitor(WallClockLoadBalancerMonitor);


LoadBalancer :: LoadBalancer(Domain *d)  : wtpList()
{
    domain = d;
}


IRResultType
LoadBalancer :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    IntArray wtp;
    IR_GIVE_OPTIONAL_FIELD(ir, wtp, _IFT_LoadBalancer_wtp);

    this->initializeWtp(wtp);

    return IRRT_OK;
}

void
LoadBalancer :: initializeWtp(IntArray &wtp)
{
    int size = wtp.giveSize();

    if ( size ) {
        wtpList.clear();
        wtpList.reserve(size);
        for ( int iwtp: wtp ) {
            std :: unique_ptr< WorkTransferPlugin > plugin;
            if ( iwtp == 1 ) {
                plugin.reset( new NonlocalMaterialWTP(this) );
            } else {
                OOFEM_ERROR("Unknown work transfer plugin type");
            }

            wtpList.push_back(std :: move( plugin ));
        }
    }
}


void
LoadBalancer :: migrateLoad(Domain *d)
{
    // domain->migrateLoad(this);
    int nproc = d->giveEngngModel()->giveNumberOfProcesses();
    int myrank = d->giveEngngModel()->giveRank();

    OOFEM_LOG_RELEVANT("[%d] LoadBalancer: migrateLoad: migrating load\n", myrank);

    // initialize work transfer plugins before any transfer
    for ( auto &wtp: wtpList ) {
        wtp->init(d);
    }

    CommunicatorBuff cb(nproc, CBT_dynamic);
    Communicator com(d->giveEngngModel(), &cb, myrank, nproc, CommMode_Dynamic);

    // move existing dofmans and elements, that will be local on current partition,
    // into local map
    com.packAllData(this, d, & LoadBalancer :: packMigratingData);
    com.initExchange(MIGRATE_LOAD_TAG);

    // do something in between
    d->initGlobalDofManMap();
    d->initGlobalElementMap();

    this->deleteRemoteDofManagers(d);
    this->deleteRemoteElements(d);

    // receive remote data
    com.unpackAllData(this, d, & LoadBalancer :: unpackMigratingData);
    com.finishExchange();

    d->commitTransactions( d->giveTransactionManager() );


#if LoadBalancer_debug_print
    // debug print
    int nnodes = d->giveNumberOfDofManagers(), nelems = d->giveNumberOfElements();
    fprintf(stderr, "\n[%d] Nodal Table\n", myrank);
    for ( int i = 1; i <= nnodes; i++ ) {
        if ( d->giveDofManager(i)->giveParallelMode() == DofManager_local ) {
            fprintf( stderr, "[%d]: %5d[%d] local\n", myrank, i, d->giveDofManager(i)->giveGlobalNumber() );
        } else if ( d->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
            fprintf( stderr, "[%d]: %5d[%d] shared ", myrank, i, d->giveDofManager(i)->giveGlobalNumber() );
            for ( int j = 1; j <= d->giveDofManager(i)->givePartitionList()->giveSize(); j++ ) {
                fprintf( stderr, "%d ", d->giveDofManager(i)->givePartitionList()->at(j) );
            }

            fprintf(stderr, "\n");
        }
    }

    fprintf(stderr, "\n[%d] Element Table\n", myrank);
    for ( int i = 1; i <= nelems; i++ ) {
        fprintf(stderr, "%5d {", i);
        for ( int j = 1; j <= d->giveElement(i)->giveNumberOfDofManagers(); j++ ) {
            fprintf( stderr, "%d ", d->giveElement(i)->giveDofManager(j)->giveNumber() );
        }

        fprintf(stderr, "}\n");
    }

#endif

    // migrate work transfer plugin data
    for ( auto &wtp: wtpList ) {
        wtp->migrate();
    }

    // update work transfer plugin data
    for ( auto &wtp: wtpList ) {
        wtp->update();
    }

    // print some local statistics
    int nelem = domain->giveNumberOfElements();
    int nnode = domain->giveNumberOfDofManagers();
    int lnode = 0, lelem = 0;

    for ( int i = 1; i <= nnode; i++ ) {
        if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_local ) {
            lnode++;
        }
    }

    for ( int i = 1; i <= nelem; i++ ) {
        if ( domain->giveElement(i)->giveParallelMode() == Element_local ) {
            lelem++;
        }
    }

    OOFEM_LOG_RELEVANT("[%d] LB Statistics:  local elem=%d local node=%d\n", myrank, lelem, lnode);
}

int
LoadBalancer :: packMigratingData(Domain *d, ProcessCommunicator &pc)
{
    int myrank = d->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();

    //  **************************************************
    //  Pack migrating data to remote partition
    //  **************************************************

    // pack dofManagers
    if ( iproc == myrank ) {
        return 1;                // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
    ProcessCommDataStream pcDataStream(pcbuff);
    // loop over dofManagers
    int ndofman = d->giveNumberOfDofManagers();
    for ( int idofman = 1; idofman <= ndofman; idofman++ ) {
        DofManager *dofman = d->giveDofManager(idofman);
        // sync data to remote partition
        // if dofman already present on remote partition then there is no need to sync
        //if ((this->giveDofManPartitions(idofman)->findFirstIndexOf(iproc))) {
        if ( ( this->giveDofManPartitions(idofman)->findFirstIndexOf(iproc) ) &&
            ( !dofman->givePartitionList()->findFirstIndexOf(iproc) ) ) {
            pcbuff->write( dofman->giveInputRecordName() );
            pcbuff->packInt( this->giveDofManState(idofman) );
            pcbuff->packInt( dofman->giveGlobalNumber() );

            // pack dofman state (this is the local dofman, not available on remote)
            /* this is a potential performance leak, sending shared dofman to a partition,
             * in which is already shared does not require to send context (is already there)
             * here for simplicity it is always send */
            dofman->saveContext(& pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State | CM_UnknownDictState);
            // send list of new partitions
            pcbuff->packIntArray( * ( this->giveDofManPartitions(idofman) ) );
        }
    }

    // pack end-of-dofman-section record
    pcbuff->write("");

    int nelem = d->giveNumberOfElements(), nsend = 0;

    for ( int ielem = 1; ielem <= nelem; ielem++ ) { // begin loop over elements
        Element *elem = d->giveElement(ielem);
        if ( ( elem->giveParallelMode() == Element_local ) &&
            ( this->giveElementPartition(ielem) == iproc ) ) {
            // pack local element (node numbers should be global ones!!!)
            // pack type
            pcbuff->write( elem->giveInputRecordName() );
            // nodal numbers should be packed as global !!
            elem->saveContext(& pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State);
            nsend++;
        }
    } // end loop over elements

    // pack end-of-element-record
    pcbuff->write("");

    OOFEM_LOG_RELEVANT("[%d] LoadBalancer:: sending %d migrating elements to %d\n", myrank, nsend, iproc);

    return 1;
}


int
LoadBalancer :: unpackMigratingData(Domain *d, ProcessCommunicator &pc)
{
    // create temp space for dofManagers and elements
    // merging should be made by domain ?
    // maps of new dofmanagers and elements indexed by global number

    // we can put local dofManagers and elements into maps (should be done before unpacking)
    // int nproc=this->giveEngngModel()->giveNumberOfProcesses();
    int myrank = d->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();
    int _mode, _globnum;
    bool _newentry;
    std :: string _type;
    IntArray _partitions, local_partitions;
    //LoadBalancer::DofManMode dmode;
    DofManager *dofman;
    DomainTransactionManager *dtm = d->giveTransactionManager();

    //  **************************************************
    //  Unpack migrating data to remote partition
    //  **************************************************

    if ( iproc == myrank ) {
        return 1;                // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
    ProcessCommDataStream pcDataStream(pcbuff);

    // unpack dofman data
    do {
        pcbuff->read(_type);
        if ( _type.size() == 0 ) { // Empty string marks end of data
            break;
        }
        pcbuff->read(_mode);
        switch ( _mode ) {
        case LoadBalancer :: DM_Remote:
            // receiving new local dofManager
            pcbuff->read(_globnum);
            /*
             * _newentry = false;
             * if ( ( dofman = dtm->giveDofManager(_globnum) ) == NULL ) {
             *  // data not available -> create a new one
             *  _newentry = true;
             *  dofman = classFactory.createDofManager(_etype, 0, d);
             * }
             */
            _newentry = true;
            dofman = classFactory.createDofManager(_type.c_str(), 0, d);

            dofman->setGlobalNumber(_globnum);
            // unpack dofman state (this is the local dofman, not available on remote)
            dofman->restoreContext(& pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State | CM_UnknownDictState);
            // unpack list of new partitions
            pcbuff->unpackIntArray(_partitions);
            dofman->setPartitionList(& _partitions);
            dofman->setParallelMode(DofManager_local);
            // add transaction if new entry allocated; otherwise existing one has been modified via returned dofman
            if ( _newentry ) {
                dtm->addDofManTransaction(DomainTransactionManager :: DTT_ADD, _globnum, dofman);
            }

            //dmanMap[_globnum] = dofman;
            break;

        case LoadBalancer :: DM_Shared:
            // receiving new shared dofManager, that was local on sending partition
            // should be received only once (from partition where was local)
            pcbuff->read(_globnum);
            /*
             * _newentry = false;
             * if ( ( dofman = dtm->giveDofManager(_globnum) ) == NULL ) {
             *  // data not available -> mode should be SharedUpdate
             *  _newentry = true;
             *  dofman = classFactory.createDofManager(_etype, 0, d);
             * }
             */
            _newentry = true;
            dofman = classFactory.createDofManager(_type.c_str(), 0, d);


            dofman->setGlobalNumber(_globnum);
            // unpack dofman state (this is the local dofman, not available on remote)
            dofman->restoreContext(& pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State | CM_UnknownDictState);
            // unpack list of new partitions
            pcbuff->unpackIntArray(_partitions);
            dofman->setPartitionList(& _partitions);
            dofman->setParallelMode(DofManager_shared);
#ifdef __VERBOSE_PARALLEL
            fprintf(stderr, "[%d] received Shared new dofman [%d]\n", myrank, _globnum);
#endif
            // add transaction if new entry allocated; otherwise existing one has been modified via returned dofman
            if ( _newentry ) {
                dtm->addDofManTransaction(DomainTransactionManager :: DTT_ADD, _globnum, dofman);
            }

            //dmanMap[_globnum] = dofman;
            break;

        default:
            OOFEM_ERROR("unexpected dof manager mode (%d)", _mode);
        }
    } while ( 1 );

    // unpack element data
    Element *elem;
    int nrecv = 0;
    do {
        pcbuff->read(_type);
        if ( _type.size() == 0 ) {
            break;
        }

        elem = classFactory.createElement(_type.c_str(), 0, d);
        elem->restoreContext(& pcDataStream, CM_Definition | CM_State);
        elem->initForNewStep();
        dtm->addElementTransaction(DomainTransactionManager :: DTT_ADD, elem->giveGlobalNumber(), elem);
        nrecv++;
        //recvElemList.push_back(elem);
    } while ( 1 );

    OOFEM_LOG_RELEVANT("[%d] LoadBalancer:: receiving %d migrating elements from %d\n", myrank, nrecv, iproc);

    return 1;
}


/* will delete those dofmanagers, that were sent to remote partition and are locally owned here
 * so they are no longer necessary (those with state equal to DM_Remote and DM_SharedMerge)
 * This will update domain DofManager list as well as global dmanMap and physically deletes the remote dofManager
 */
void
LoadBalancer :: deleteRemoteDofManagers(Domain *d)
{
    int ndofman = d->giveNumberOfDofManagers();
    //LoadBalancer* lb = this->giveLoadBalancer();
    LoadBalancer :: DofManMode dmode;
    DofManager *dman;
    int myrank = d->giveEngngModel()->giveRank();
    DomainTransactionManager *dtm = d->giveTransactionManager();
    // loop over local nodes

    for ( int i = 1; i <= ndofman; i++ ) {
        dmode = this->giveDofManState(i);
        if ( dmode == LoadBalancer :: DM_Remote ) {
            // positive candidate found
            dtm->addDofManTransaction(DomainTransactionManager :: DTT_Remove, d->giveDofManager(i)->giveGlobalNumber(), NULL);
            // dmanMap.erase (d->giveDofManager (i)->giveGlobalNumber());
            //dman = dofManagerList->unlink (i);
            //delete dman;
        } else if ( dmode == LoadBalancer :: DM_NULL ) {
            // positive candidate found; we delete all null dof managers
            // they will be created by nonlocalmatwtp if necessary.
            // potentially, they can be reused, but this will make the code too complex
            dtm->addDofManTransaction(DomainTransactionManager :: DTT_Remove, d->giveDofManager(i)->giveGlobalNumber(), NULL);
        } else if ( dmode == LoadBalancer :: DM_Shared ) {
            dman = d->giveDofManager(i);
            dman->setPartitionList( this->giveDofManPartitions(i) );
            dman->setParallelMode(DofManager_shared);
            if ( !dman->givePartitionList()->findFirstIndexOf(myrank) ) {
                dtm->addDofManTransaction(DomainTransactionManager :: DTT_Remove, d->giveDofManager(i)->giveGlobalNumber(), NULL);
                //dmanMap.erase (this->giveDofManager (i)->giveGlobalNumber());
                //dman = dofManagerList->unlink (i);
                //delete dman;
            }
        } else if ( dmode == LoadBalancer :: DM_Local ) {
            IntArray _empty(0);
            dman = d->giveDofManager(i);
            dman->setPartitionList(& _empty);
            dman->setParallelMode(DofManager_local);
        } else {
            OOFEM_ERROR("unknown dmode encountered");
        }
    }
}

/* will delete those elements, that were sent to remote partition and are locally owned here
 * so they are no longer necessary (those with state equal to DM_Remote and DM_SharedMerge)
 * This will update domain DofManager list as well as global dmanMap and physically deletes the remote dofManager
 */
void
LoadBalancer :: deleteRemoteElements(Domain *d)
{
    int nelem =  d->giveNumberOfElements();
    int myrank = d->giveEngngModel()->giveRank();
    //LoadBalancer* lb = this->giveLoadBalancer();
    DomainTransactionManager *dtm = d->giveTransactionManager();
    //Element* elem;

    // loop over local nodes

    for ( int i = 1; i <= nelem; i++ ) {
        if ( this->giveElementPartition(i) != myrank ) {
            // positive candidate found
            // this->deleteElement (i);  // delete and set entry to NULL
            dtm->addElementTransaction(DomainTransactionManager :: DTT_Remove, d->giveElement(i)->giveGlobalNumber(), NULL);
            //elem = elementList->unlink (i);
            //dmanMap.erase (elem->giveGlobalNumber());
            //delete (elem);
        } else if ( d->giveElement(i)->giveParallelMode() != Element_local ) {
            dtm->addElementTransaction(DomainTransactionManager :: DTT_Remove, d->giveElement(i)->giveGlobalNumber(), NULL);
        }
    }
}


void
LoadBalancer :: printStatistics() const
{
    EngngModel *emodel = domain->giveEngngModel();
    int nelem, nnode;
    int lelem = 0, lnode = 0;
    int myrank = emodel->giveRank();

    nelem = domain->giveNumberOfElements();
    nnode = domain->giveNumberOfDofManagers();

    for ( int i = 1; i <= nnode; i++ ) {
        if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_local ) {
            lnode++;
        }
    }

    for ( int i = 1; i <= nelem; i++ ) {
        if ( domain->giveElement(i)->giveParallelMode() == Element_local ) {
            lelem++;
        }
    }

    double mySolutionWTime = emodel->giveTimer()->getWtime(EngngModelTimer :: EMTT_AnalysisTimer);
    double mySolutionUTime = emodel->giveTimer()->getUtime(EngngModelTimer :: EMTT_AnalysisTimer);

    OOFEM_LOG_RELEVANT("[%d] LB Statistics:  wt=%.1f ut=%.1f nelem=%d nnode=%d\n", myrank,
                       mySolutionWTime, mySolutionUTime, lelem, lnode);
}



LoadBalancerMonitor :: LoadBalancerDecisionType
WallClockLoadBalancerMonitor :: decide(TimeStep *tStep)
{
    int nproc = emodel->giveNumberOfProcesses();
    int myrank = emodel->giveRank();
    Domain *d = emodel->giveLoadBalancer()->giveDomain();
    int nelem;
    double *node_solutiontimes = new double [ nproc ];
    double *node_relcomppowers = new double [ nproc ];
    double *node_equivelements = new double [ nproc ];
    double min_st, max_st;
    double relWallClockImbalance;
    double absWallClockImbalance;
    double neqelems, sum_relcomppowers;

    if ( node_solutiontimes == NULL ) {
        OOFEM_ERROR("failed to allocate node_solutiontimes array");
    }

    if ( node_relcomppowers == NULL ) {
        OOFEM_ERROR("failed to allocate node_relcomppowers array");
    }

    if ( node_equivelements == NULL ) {
        OOFEM_ERROR("failed to allocate node_equivelements array");
    }


    // compute wall solution time of my node
    double mySolutionTime = emodel->giveTimer()->getWtime(EngngModelTimer :: EMTT_NetComputationalStepTimer);

#ifdef __LB_DEBUG
    // perturb solution time artificially if requested
    bool perturb = false;
    for ( auto perturbedStep: perturbedSteps ) {
        if ( perturbedStep.test( tStep->giveNumber() ) ) {
            perturb  = true;
            break;
        }
    }

    if ( perturb ) {
        mySolutionTime *= perturbFactor;
        OOFEM_LOG_RELEVANT("[%d] WallClockLoadBalancerMonitor: perturbed solution time by factor=%.2f\n", myrank, perturbFactor);
    }

#endif

    // collect wall clock computational time
    MPI_Allgather(& mySolutionTime, 1, MPI_DOUBLE, node_solutiontimes, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    OOFEM_LOG_RELEVANT("\nLoadBalancer:: individual processor times [sec]: (");
    for ( int i = 0; i < nproc; i++ ) {
        OOFEM_LOG_RELEVANT(" %.3f", node_solutiontimes [ i ]);
    }

    OOFEM_LOG_RELEVANT(")\n");

    // detect imbalance
    min_st = max_st = node_solutiontimes [ 0 ];
    for ( int i = 0; i < nproc; i++ ) {
        min_st = min(min_st, node_solutiontimes [ i ]);
        max_st = max(max_st, node_solutiontimes [ i ]);
    }

    absWallClockImbalance = ( max_st - min_st );
    if ( min_st ) {
        relWallClockImbalance = ( ( max_st - min_st ) / min_st );
    } else {
        relWallClockImbalance = 0.0;
    }

    // update node (processor) weights

    // compute number or equivalent elements (equavalent element has computational weight equal to 1.0)
    nelem = d->giveNumberOfElements();
    neqelems = 0.0;
    for ( int ie = 1; ie <= nelem; ie++ ) {
        if ( d->giveElement(ie)->giveParallelMode() == Element_remote ) {
            continue;
        }

        neqelems += d->giveElement(ie)->predictRelativeComputationalCost();
    }

    // exchange number or equivalent elements
    MPI_Allgather(& neqelems, 1, MPI_DOUBLE, node_equivelements, 1, MPI_DOUBLE, MPI_COMM_WORLD);


    if ( !this->staticNodeWeightFlag ) {
        // compute relative computational powers (solution_time/number_of_equivalent_elements)
        for ( int i = 0; i < nproc; i++ ) {
            node_relcomppowers [ i ] = node_equivelements [ i ] / node_solutiontimes [ i ];
        }

        // normalize computational powers
        sum_relcomppowers = 0.0;
        for ( int i = 0; i < nproc; i++ ) {
            sum_relcomppowers += node_relcomppowers [ i ];
        }

        for ( int i = 0; i < nproc; i++ ) {
            nodeWeights(i) = node_relcomppowers [ i ] / sum_relcomppowers;
        }
    }

    // log equivalent elements on nodes
    OOFEM_LOG_RELEVANT("[%d] LoadBalancer:  node equivalent elements: ", myrank);
    for ( int i = 0; i < nproc; i++ ) {
        OOFEM_LOG_RELEVANT("%6d ", ( int ) node_equivelements [ i ]);
    }

    OOFEM_LOG_RELEVANT("\n");

    // log processor weights
    OOFEM_LOG_RELEVANT("[%d] LoadBalancer: updated proc weights: ", myrank);
    for ( int i = 0; i < nproc; i++ ) {
#ifdef __LB_DEBUG
        OOFEM_LOG_RELEVANT( "%22.15e ", nodeWeights(i) );
#else
        OOFEM_LOG_RELEVANT( "%4.3f ", nodeWeights(i) );
#endif
    }

    OOFEM_LOG_RELEVANT("\n");

    delete[] node_solutiontimes;
    delete[] node_relcomppowers;
    delete[] node_equivelements;

#ifdef __LB_DEBUG
    if ( recoveredSteps.giveSize() ) {
        // recover lb if requested
        int pos;
        if ( ( pos = recoveredSteps.findFirstIndexOf( tStep->giveNumber() ) ) ) {
            double procWeight, sumWeight = 0.0, *procWeights = new double [ nproc ];

            // assign prescribed processing weight
            procWeight = processingWeights.at(pos);
            OOFEM_LOG_RELEVANT("[%d] WallClockLoadBalancerMonitor: processing weight overriden by value=%e\n", myrank, procWeight);

            // exchange processing weights
            MPI_Allgather(& procWeight, 1, MPI_DOUBLE, procWeights, 1, MPI_DOUBLE, MPI_COMM_WORLD);
            for ( int i = 0; i < nproc; i++ ) {
                nodeWeights(i) = procWeights [ i ];
                sumWeight += procWeights [ i ];
            }

            delete[] procWeights;

            if ( fabs(sumWeight - 1.0) > 1.0e-10 ) {
                OOFEM_ERROR("[%d] processing weights do not sum to 1.0 (sum = %e)\n", sumWeight);
            }

            OOFEM_LOG_RELEVANT("[%d] LoadBalancer: wall clock imbalance rel=%.2f\%,abs=%.2fs, recovering load\n", myrank, 100 * relWallClockImbalance, absWallClockImbalance);
            return LBD_RECOVER;
        } else {
            OOFEM_LOG_RELEVANT("[%d] LoadBalancer: wall clock imbalance rel=%.2f\%,abs=%.2fs, continuing\n", myrank, 100 * relWallClockImbalance, absWallClockImbalance);
            return LBD_CONTINUE;
        }
    }

#endif

    // decide
    if ( ( tStep->giveNumber() % this->lbstep == 0 ) &&
        ( ( absWallClockImbalance > this->absWallClockImbalanceTreshold ) ||
         ( ( relWallClockImbalance > this->relWallClockImbalanceTreshold ) && ( absWallClockImbalance > this->minAbsWallClockImbalanceTreshold ) ) ) ) {
        OOFEM_LOG_RELEVANT("[%d] LoadBalancer: wall clock imbalance rel=%.2f\%,abs=%.2fs, recovering load\n", myrank, 100 * relWallClockImbalance, absWallClockImbalance);
        return LBD_RECOVER;
    } else {
        OOFEM_LOG_RELEVANT("[%d] LoadBalancer: wall clock imbalance rel=%.2f\%,abs=%.2fs, continuing\n", myrank, 100 * relWallClockImbalance, absWallClockImbalance);
        return LBD_CONTINUE;
    }
}

IRResultType
LoadBalancerMonitor :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    int nproc = emodel->giveNumberOfProcesses();
    int nodeWeightMode = 0;

    nodeWeights.resize(nproc);
    for ( int i = 0; i < nproc; i++ ) {
        nodeWeights(i) = 1.0 / nproc;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, nodeWeightMode, _IFT_LoadBalancerMonitor_nodeWeightMode);
    if ( nodeWeightMode == 0 ) { // default, dynamic weights
        staticNodeWeightFlag = false;
    } else if ( nodeWeightMode == 1 ) { // equal weights for all nodes
        staticNodeWeightFlag = true;
    } else if ( nodeWeightMode == 2 ) { // user defined static weights
        IR_GIVE_OPTIONAL_FIELD(ir, nodeWeights, _IFT_LoadBalancerMonitor_initialnodeweights);
        if ( nodeWeights.giveSize() != nproc ) {
            OOFEM_ERROR("nodeWeights size not equal to number of processors");
        }

        staticNodeWeightFlag = true;
    } else {
        OOFEM_ERROR("unsupported node weight type, using default value");
        staticNodeWeightFlag = false;
    }

    return IRRT_OK;
}

IRResultType
WallClockLoadBalancerMonitor :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    result = LoadBalancerMonitor :: initializeFrom(ir);



    IR_GIVE_OPTIONAL_FIELD(ir, relWallClockImbalanceTreshold, _IFT_WallClockLoadBalancerMonitor_relwct);
    IR_GIVE_OPTIONAL_FIELD(ir, absWallClockImbalanceTreshold, _IFT_WallClockLoadBalancerMonitor_abswct);
    IR_GIVE_OPTIONAL_FIELD(ir, minAbsWallClockImbalanceTreshold, _IFT_WallClockLoadBalancerMonitor_minwct);
    IR_GIVE_OPTIONAL_FIELD(ir, lbstep, _IFT_WallClockLoadBalancerMonitor_lbstep);

#ifdef __LB_DEBUG
    perturbedSteps.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, perturbedSteps, _IFT_WallClockLoadBalancerMonitor_perturbedsteps);
    perturbFactor = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, perturbFactor, _IFT_WallClockLoadBalancerMonitor_perturbfactor);

    recoveredSteps.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, recoveredSteps, _IFT_WallClockLoadBalancerMonitor_recoveredsteps);
    processingWeights.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, processingWeights, _IFT_WallClockLoadBalancerMonitor_processingweights);
    if ( recoveredSteps.giveSize() != processingWeights.giveSize() ) {
        OOFEM_ERROR("mismatch size of lbrecoveredsteps and lbprocessingweights");
    }

#endif

    return result;
}

/*
 * void
 * LoadBalancer::migrateLoad () {}
 *
 * IRResultType
 * LoadBalancer::initializeFrom (InputRecord* ir) {
 *
 * return IRRT_OK;
 * }
 *
 * IRResultType
 * LoadBalancerMonitor::initializeFrom (InputRecord* ir) {return IRRT_OK;}
 */



LoadBalancer :: WorkTransferPlugin :: WorkTransferPlugin(LoadBalancer *_lb) {
    lb = _lb;
}
LoadBalancer :: WorkTransferPlugin :: ~WorkTransferPlugin() { }
} // end namespace oofem
