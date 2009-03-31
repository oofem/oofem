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

#include "loadbalancer.h"
#include "domain.h"
#include "engngm.h"
#include "timer.h"
#include "mathfem.h"
#include "timestep.h"
#include "usrdefsub.h"


#include "parallel.h"
#include "processcomm.h"
#include "datastream.h"
#include "communicator.h"
#include "domaintransactionmanager.h"
#include "nonlocalmatwtp.h"


#define LoadBalancer_debug_print 0

LoadBalancer :: LoadBalancer(Domain *d)  : wtpList(0)
{
    domain = d;
}


IRResultType
LoadBalancer :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    IntArray wtp;
    IR_GIVE_OPTIONAL_FIELD(ir, wtp, IFT_LoadBalancer_wtp, "wtp"); // Macro

    this->initializeWtp(wtp);

    return IRRT_OK;
}

void
LoadBalancer :: initializeWtp(IntArray &wtp) {
    int i, size = wtp.giveSize();
    WorkTransferPlugin *plugin = NULL;

    if ( size ) {
        wtpList.growTo(size);
        for ( i = 1; i <= size; i++ ) {
            if ( wtp.at(i) == 1 ) {
                plugin = new NonlocalMaterialWTP(this);
            } else {
                OOFEM_ERROR("LoadBalancer::initializeWtp: Unknown work transfer plugin type");
            }

            wtpList.put(i, plugin);
        }
    }
}


void
LoadBalancer :: migrateLoad(Domain *d)
{
    // domain->migrateLoad(this);
    int i;
    int nproc = d->giveEngngModel()->giveNumberOfProcesses();
    int myrank = d->giveEngngModel()->giveRank();

    OOFEM_LOG_RELEVANT("[%d] LoadBalancer: migrateLoad: migrating load\n", myrank);

    // initialize work transfer plugins before any transfer
    for ( i = 1; i <= wtpList.giveSize(); i++ ) {
        wtpList.at(i)->init(d);
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
    int j, nnodes = d->giveNumberOfDofManagers(), nelems = d->giveNumberOfElements();
    fprintf(stderr, "\n[%d] Nodal Table\n", myrank);
    for ( i = 1; i <= nnodes; i++ ) {
        if ( d->giveDofManager(i)->giveParallelMode() == DofManager_local ) {
            fprintf( stderr, "[%d]: %5d[%d] local\n", myrank, i, d->giveDofManager(i)->giveGlobalNumber() );
        } else if ( d->giveDofManager(i)->giveParallelMode() == DofManager_shared )  {
            fprintf( stderr, "[%d]: %5d[%d] shared ", myrank, i, d->giveDofManager(i)->giveGlobalNumber() );
            for ( j = 1; j <= d->giveDofManager(i)->givePartitionList()->giveSize(); j++ ) {
                fprintf( stderr, "%d ", d->giveDofManager(i)->givePartitionList()->at(j) );
            }

            fprintf(stderr, "\n");
        }
    }

    fprintf(stderr, "\n[%d] Element Table\n", myrank);
    for ( i = 1; i <= nelems; i++ ) {
        fprintf(stderr, "%5d {", i);
        for ( j = 1; j <= d->giveElement(i)->giveNumberOfDofManagers(); j++ ) {
            fprintf( stderr, "%d ", d->giveElement(i)->giveDofManager(j)->giveNumber() );
        }

        fprintf(stderr, "}\n");
    }

#endif

    // migrate work transfer plugin data
    for ( i = 1; i <= wtpList.giveSize(); i++ ) {
        wtpList.at(i)->migrate();
    }

    // update work transfer plugin data
    for ( i = 1; i <= wtpList.giveSize(); i++ ) {
        wtpList.at(i)->update();
    }
}

int
LoadBalancer :: packMigratingData(Domain *d, ProcessCommunicator &pc)
{
    int myrank = d->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();
    int idofman, ndofman;
    classType dtype;
    DofManager *dofman;
    LoadBalancer :: DofManMode dmode;

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
    ndofman = d->giveNumberOfDofManagers();
    for ( idofman = 1; idofman <= ndofman; idofman++ ) {
        dofman = d->giveDofManager(idofman);
        dmode = this->giveDofManState(idofman);
        dtype = dofman->giveClassID();
        // sync data to remote partition
        // if dofman already present on remote partition then there is no need to sync
        //if ((this->giveDofManPartitions(idofman)->findFirstIndexOf(iproc))) {
        if ( ( this->giveDofManPartitions(idofman)->findFirstIndexOf(iproc) ) &&
            ( !dofman->givePartitionList()->findFirstIndexOf(iproc) ) ) {
            pcbuff->packInt(dtype);
            pcbuff->packInt(dmode);
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
    pcbuff->packInt(LOADBALANCER_END_DATA);

    int ielem, nelem = d->giveNumberOfElements(), nsend = 0;

    Element *elem;

    for ( ielem = 1; ielem <= nelem; ielem++ ) { // begin loop over elements
        elem = d->giveElement(ielem);
        if ( ( elem->giveParallelMode() == Element_local ) &&
            ( this->giveElementPartition(ielem) == iproc ) ) {
            // pack local element (node numbers shuld be global ones!!!)
            // pack type
            pcbuff->packInt( elem->giveClassID() );
            // nodal numbers shuld be packed as global !!
            elem->saveContext(& pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State);
            nsend++;
        }
    } // end loop over elements

    // pack end-of-element-record
    pcbuff->packInt(LOADBALANCER_END_DATA);

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
    int _mode, _globnum, _type;
    bool _newentry;
    classType _etype;
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

    pcbuff->unpackInt(_type);
    // unpack dofman data
    while ( _type != LOADBALANCER_END_DATA ) {
        _etype = ( classType ) _type;
        pcbuff->unpackInt(_mode);
        switch ( _mode ) {
        case LoadBalancer :: DM_Remote:
            // receiving new local dofManager
            pcbuff->unpackInt(_globnum);
            /*
            _newentry = false;
            if ( ( dofman = dtm->giveDofManager(_globnum) ) == NULL ) {
                // data not available -> create a new one
                _newentry = true;
                dofman = CreateUsrDefDofManagerOfType(_etype, 0, d);
            }
            */
            _newentry = true;
            dofman = CreateUsrDefDofManagerOfType(_etype, 0, d);

            dofman->setGlobalNumber(_globnum);
            // unpack dofman state (this is the local dofman, not available on remote)
            dofman->restoreContext(& pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State | CM_UnknownDictState);
            // unpack list of new partitions
            pcbuff->unpackIntArray(_partitions);
            dofman->setPartitionList(&_partitions);
            dofman->setParallelMode(DofManager_local);
            // add transaction if new entry allocated; otherwise existing one has been modified via returned dofman
            if ( _newentry ) {
                dtm->addTransaction(DomainTransactionManager :: DTT_ADD, DomainTransactionManager :: DCT_DofManager, _globnum, dofman);
            }

            //dmanMap[_globnum] = dofman;
            break;

        case LoadBalancer :: DM_Shared:
            // receiving new shared dofManager, that was local on sending partition
            // should be received only once (from partition where was local)
            pcbuff->unpackInt(_globnum);
            /*
            _newentry = false;
            if ( ( dofman = dtm->giveDofManager(_globnum) ) == NULL ) {
                // data not available -> mode should be SharedUpdate
                _newentry = true;
                dofman = CreateUsrDefDofManagerOfType(_etype, 0, d);
            }
            */
            _newentry = true;
            dofman = CreateUsrDefDofManagerOfType(_etype, 0, d);


            dofman->setGlobalNumber(_globnum);
            // unpack dofman state (this is the local dofman, not available on remote)
            dofman->restoreContext(& pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State | CM_UnknownDictState);
            // unpack list of new partitions
            pcbuff->unpackIntArray(_partitions);
            dofman->setPartitionList(&_partitions);
            dofman->setParallelMode(DofManager_shared);
#ifdef __VERBOSE_PARALLEL
            fprintf(stderr, "[%d] received Shared new dofman [%d]\n", myrank, _globnum);
#endif
            // add transaction if new entry allocated; otherwise existing one has been modified via returned dofman
            if ( _newentry ) {
                dtm->addTransaction(DomainTransactionManager :: DTT_ADD, DomainTransactionManager :: DCT_DofManager, _globnum, dofman);
            }

            //dmanMap[_globnum] = dofman;
            break;

        default:
            OOFEM_ERROR2("LoadBalancer::unpackMigratingData: unexpected dof manager type (%d)", _type);
        }

        // get next type record
        pcbuff->unpackInt(_type);
    }

    ; // while (_type != LOADBALANCER_END_DATA);

    // unpack element data
    Element *elem;
    int nrecv = 0;
    do {
        pcbuff->unpackInt(_type);
        if ( _type == LOADBALANCER_END_DATA ) {
            break;
        }

        _etype = ( classType ) _type;
        elem = CreateUsrDefElementOfType(_etype, 0, d);
        elem->restoreContext(& pcDataStream, CM_Definition | CM_State);
	elem->initForNewStep();
        dtm->addTransaction(DomainTransactionManager :: DTT_ADD, DomainTransactionManager :: DCT_Element, elem->giveGlobalNumber(), elem);
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
    int i, ndofman =  d->giveNumberOfDofManagers();
    //LoadBalancer* lb = this->giveLoadBalancer();
    LoadBalancer :: DofManMode dmode;
    DofManager *dman;
    int myrank = d->giveEngngModel()->giveRank();
    DomainTransactionManager *dtm = d->giveTransactionManager();
    // loop over local nodes

    for ( i = 1; i <= ndofman; i++ ) {
        dmode = this->giveDofManState(i);
        if ( ( dmode == LoadBalancer :: DM_Remote ) ) {
            // positive candidate found
            dtm->addTransaction(DomainTransactionManager :: DTT_Remove, DomainTransactionManager :: DCT_DofManager, d->giveDofManager(i)->giveGlobalNumber(), NULL);
            // dmanMap.erase (d->giveDofManager (i)->giveGlobalNumber());
            //dman = dofManagerList->unlink (i);
            //delete dman;
        } else if ( ( dmode == LoadBalancer :: DM_NULL ) ) {
            // positive candidate found; we delete all null dof managers
            // they will be created by nonlocalmatwtp if necessary.
            // potentially, they can be reused, but this will make the code too complex
            dtm->addTransaction(DomainTransactionManager :: DTT_Remove, DomainTransactionManager :: DCT_DofManager, d->giveDofManager(i)->giveGlobalNumber(), NULL);
        } else if ( dmode == LoadBalancer :: DM_Shared ) {
            dman = d->giveDofManager(i);
            dman->setPartitionList(this->giveDofManPartitions(i));
            dman->setParallelMode(DofManager_shared);
            if ( !dman->givePartitionList()->findFirstIndexOf(myrank) ) {
                dtm->addTransaction(DomainTransactionManager :: DTT_Remove, DomainTransactionManager :: DCT_DofManager, d->giveDofManager(i)->giveGlobalNumber(), NULL);
                //dmanMap.erase (this->giveDofManager (i)->giveGlobalNumber());
                //dman = dofManagerList->unlink (i);
                //delete dman;
            }
        } else if ( dmode == LoadBalancer :: DM_Local ) {
            IntArray _empty(0);
            dman = d->giveDofManager(i);
            dman->setPartitionList(&_empty);
            dman->setParallelMode(DofManager_local);
        } else {
            OOFEM_ERROR("Domain::deleteRemoteDofManagers: unknown dmode encountered");
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
    int i, nelem =  d->giveNumberOfElements();
    int myrank = d->giveEngngModel()->giveRank();
    //LoadBalancer* lb = this->giveLoadBalancer();
    DomainTransactionManager *dtm = d->giveTransactionManager();
    //Element* elem;

    // loop over local nodes

    for ( i = 1; i <= nelem; i++ ) {
        if ( this->giveElementPartition(i) != myrank ) {
            // positive candidate found
            // this->deleteElement (i);  // delete and set entry to NULL
            dtm->addTransaction(DomainTransactionManager :: DTT_Remove, DomainTransactionManager :: DCT_Element, d->giveElement(i)->giveGlobalNumber(), NULL);
            //elem = elementList->unlink (i);
            //dmanMap.erase (elem->giveGlobalNumber());
            //delete (elem);
        } else if ( d->giveElement(i)->giveParallelMode() != Element_local ) {
            dtm->addTransaction(DomainTransactionManager :: DTT_Remove, DomainTransactionManager :: DCT_Element, d->giveElement(i)->giveGlobalNumber(), NULL);
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
    int i;

    nelem = domain->giveNumberOfElements();
    nnode = domain->giveNumberOfDofManagers();

    for ( i = 1; i <= nnode; i++ ) {
        if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_local ) {
            lnode++;
        }
    }

    for ( i = 1; i <= nelem; i++ ) {
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
WallClockLoadBalancerMonitor :: decide(TimeStep *atTime)
{
    int i, nproc = emodel->giveNumberOfProcesses();
    int myrank = emodel->giveRank();
    Domain *d = emodel->giveLoadBalancer()->giveDomain();
    int ie, nelem;
    double *node_solutiontimes = new double [ nproc ];
    double *node_relcomppowers = new double [ nproc ];
    double min_st, max_st;
    double relWallClockImbalance;
    double absWallClockImbalance;
    double neqelems, sum_relcomppowers;
    double myRelativeComputationalPower;
    IntArray node_equivelements(nproc);

    if ( node_solutiontimes == NULL ) {
        OOFEM_ERROR("LoadBalancer::LoadEvaluation failed to allocate node_solutiontimes array");
    }

    if ( node_relcomppowers == NULL ) {
        OOFEM_ERROR("LoadBalancer::LoadEvaluation failed to allocate node_relcomppowers array");
    }

    // compute wall solution time of my node
    double mySolutionTime = emodel->giveTimer()->getWtime(EngngModelTimer :: EMTT_NetComputationalStepTimer);

#ifdef __LB_DEBUG
    // perturb solution time artificially if requested
    bool perturb = false;
    dynaList< Range > :: iterator perturbedStepsIter;
    for ( perturbedStepsIter = perturbedSteps.begin(); perturbedStepsIter != perturbedSteps.end(); ++perturbedStepsIter ) {
        if ( ( * perturbedStepsIter ).test( atTime->giveNumber() ) ) {
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
    for ( i = 0; i < nproc; i++ ) {
        OOFEM_LOG_RELEVANT(" %.3f", node_solutiontimes [ i ]);
    }

    OOFEM_LOG_RELEVANT(")\n");



    // detect imbalance
    min_st = max_st = node_solutiontimes [ 0 ];
    for ( i = 0; i < nproc; i++ ) {
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
    for ( ie = 1; ie <= nelem; ie++ ) {
	if(d->giveElement(ie)->giveParallelMode() == Element_remote)continue;
        neqelems += d->giveElement(ie)->predictRelativeComputationalCost();
    }

    // compute relative computational power (solution_time/number_of_equivalent_elements)
    myRelativeComputationalPower = neqelems / mySolutionTime;
    // collect relative powers
    MPI_Allgather(& myRelativeComputationalPower, 1, MPI_DOUBLE, node_relcomppowers, 1, MPI_DOUBLE, MPI_COMM_WORLD);


    // get number of equivalent elements for nodes
    for ( i = 0; i < nproc; i++ ) {
        node_equivelements.at(i + 1) = ( int ) ( node_relcomppowers [ i ] * node_solutiontimes [ i ] );
    }

    // normalize computational powers
    sum_relcomppowers = 0.0;
    for ( i = 0; i < nproc; i++ ) {
        sum_relcomppowers += node_relcomppowers [ i ];
    }

    for ( i = 0; i < nproc; i++ ) {
        nodeWeights(i) = node_relcomppowers [ i ] / sum_relcomppowers;
    }

    delete[] node_solutiontimes;
    delete[] node_relcomppowers;

    // log equivalent elements on nodes
    OOFEM_LOG_RELEVANT("[%d] LoadBalancer:  node equivalent elements: ", myrank);
    for ( i = 0; i < nproc; i++ ) {
        OOFEM_LOG_RELEVANT( "%6d ", node_equivelements.at(i + 1) );
    }

    OOFEM_LOG_RELEVANT("\n");


    // log processor weights
    OOFEM_LOG_RELEVANT("[%d] LoadBalancer: updated proc weights: ", myrank);
    for ( i = 0; i < nproc; i++ ) {
        OOFEM_LOG_RELEVANT( "%4.3f ", nodeWeights(i) );
    }

    OOFEM_LOG_RELEVANT("\n");

    // decide
    if ( ( atTime->giveNumber() % this->lbstep == 0 ) && ( ( absWallClockImbalance > this->absWallClockImbalanceTreshold ) || ( relWallClockImbalance > this->relWallClockImbalanceTreshold ) ) ) {
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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    int i, nproc = emodel->giveNumberOfProcesses();

    nodeWeights.resize(nproc);
    for ( i = 0; i < nproc; i++ ) {
        nodeWeights(i) = 1.0;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, nodeWeights, IFT_LoadBalancerMonitor_initialnodeweights, "nw"); // Macro
    if ( nodeWeights.giveSize() != nproc ) {
        OOFEM_ERROR("nodeWeights size not equal to number of processors");
    }

    return IRRT_OK;
}

IRResultType
WallClockLoadBalancerMonitor :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    result = LoadBalancerMonitor :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, relWallClockImbalanceTreshold, IFT_WallClockLoadBalancerMonitor_relwct, "relwct"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, absWallClockImbalanceTreshold, IFT_WallClockLoadBalancerMonitor_abswct, "abswct"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, lbstep, IFT_WallClockLoadBalancerMonitor_lbstep, "lbstep"); // Macro

#ifdef __LB_DEBUG
    perturbedSteps.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, perturbedSteps, IFT_WallClockLoadBalancerMonitor_perturbedsteps, "lbperturbedsteps"); // Macro
    perturbFactor = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, perturbFactor, IFT_WallClockLoadBalancerMonitor_perturbfactor, "lbperturbfactor"); // Macro
#endif

    return result;
}

/*
 #else //__PARALLEL_MODE
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

#endif
