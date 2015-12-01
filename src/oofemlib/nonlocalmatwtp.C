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

#include "nonlocalmatwtp.h"
#include "nonlocalmaterialext.h"
#include "element.h"
#include "dofmanager.h"
#include "engngm.h"
#include "gausspoint.h"
#include "material.h"
#include "communicator.h"
#include "datastream.h"
#include "domaintransactionmanager.h"
#include "classfactory.h"

#include <set>

namespace oofem {
#define NonlocalMaterialWTP_DEBUG_PRINT 0

/*
 * Returns array storing nonlocal dependency
 * (in terms of global element numbers) for given element
 */
void
NonlocalMaterialWTP :: giveElementNonlocalDepArry(IntArray &answer, Domain *d, int num)
{
    std :: set< int >relems;
    std :: set< int > :: const_iterator relemsIter;
    Element *ielem = d->giveElement(num);

    if ( ielem->giveMaterial()->giveInterface(NonlocalMaterialExtensionInterfaceType) ) {
        relems.clear();
        // loop over element IRules and their IPs to retrieve remote (nonlocal) elements
        // store their global numbers in the relems set (to avoid redundancy)
        // and then keep them in nonlocTables array.
        ielem->ipEvaluator(this, & NonlocalMaterialWTP :: giveNonlocalDepArryElementPlugin, relems);

        // convert relems set into an int array
        // and store it
        answer.resize( relems.size() );
        int _i = 1;
        for ( int relem: relems ) {
            answer.at(_i++) = relem;
        }
    } else {
        answer.clear();
    }
}


void
NonlocalMaterialWTP :: giveNonlocalDepArryElementPlugin(GaussPoint *gp, std :: set< int > &s)
{
    int remoteElemNum;

    NonlocalMaterialStatusExtensionInterface *interface =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
                                                                  giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );
    if ( interface ) {
        std :: list< localIntegrationRecord > *lir = interface->giveIntegrationDomainList();

        for ( auto &intdom: *lir ) {
            remoteElemNum = ( intdom.nearGp )->giveElement()->giveGlobalNumber();
            s.insert(remoteElemNum);
        }
    }
}


/*
 * prepares the communication maps for remote elements
 * should be called immediately after load balancing,
 * before any work transfer.
 *
 */
void
NonlocalMaterialWTP :: init(Domain *domain)
{
    int ie, gie, nelem = domain->giveNumberOfElements();
    EngngModel *emodel = domain->giveEngngModel();
    Element *elem;
    int nproc = emodel->giveNumberOfProcesses();
    int myrank = emodel->giveRank();
    CommunicatorBuff cb(nproc, CBT_dynamic);
    Communicator com(emodel, &cb, myrank, nproc, CommMode_Dynamic);
    this->nonlocElementDependencyMap.clear();

    // build nonlocal element dependency array for each element
    for ( ie = 1; ie <= nelem; ie++ ) {
        elem = domain->giveElement(ie);
        if ( ( elem->giveParallelMode() == Element_local ) ) {
            gie = elem->giveGlobalNumber();
            this->giveElementNonlocalDepArry(nonlocElementDependencyMap [ gie ], domain, ie);
        }
    }

    /* send and receive nonlocElementDependencyArry of migrating elements to remote partition */
    com.packAllData(this, domain, & NonlocalMaterialWTP :: packMigratingElementDependencies);
    com.initExchange(MIGRATE_NONLOCALDEP_TAG);
    com.unpackAllData(this, domain, & NonlocalMaterialWTP :: unpackMigratingElementDependencies);
    com.finishExchange();
}



/*
 * should be called after basic local migration is finalized,
 * when all local elements are already available
 */
void
NonlocalMaterialWTP :: migrate()
{
    Domain *domain = this->lb->giveDomain();
    EngngModel *emodel = domain->giveEngngModel();
    int nproc = emodel->giveNumberOfProcesses();
    int myrank = emodel->giveRank();
    CommunicatorBuff cb(nproc, CBT_dynamic);
    Communicator com(emodel, &cb, myrank, nproc, CommMode_Dynamic);
    StaticCommunicationBuffer commBuff(MPI_COMM_WORLD);

    /*
     * build domain nonlocal element dependency list. Then exclude local elements - what remains are unsatisfied
     * remote dependencies that have to be broadcasted and received from partitions owning relevant elements
     */
    int _locsize, i, _i, ie, _size, _globnum, result, nelems = domain->giveNumberOfElements();
    int _globsize, _val;
    Element *elem;
    std :: set< int >domainElementDepSet;
    // loop over each element dep list to assemble domain list
    for ( ie = 1; ie <= nelems; ie++ ) {
        elem = domain->giveElement(ie);
        if ( ( elem->giveParallelMode() == Element_local ) ) {
            _globnum = elem->giveGlobalNumber();
            IntArray &iedep = nonlocElementDependencyMap [ _globnum ];
            _size = iedep.giveSize();
            for ( _i = 1; _i <= _size; _i++ ) {
                domainElementDepSet.insert( iedep.at(_i) );
            }

#if NonlocalMaterialWTP_DEBUG_PRINT
            fprintf(stderr, "[%d] element %d dependency:", myrank, _globnum);
            for ( _i = 1; _i <= _size; _i++ ) {
                fprintf( stderr, "%d ", iedep.at(_i) );
            }

            fprintf(stderr, "\n");
#endif
        }
    }

#if NonlocalMaterialWTP_DEBUG_PRINT
    fprintf(stderr, "[%d] nonlocal domain dependency:", myrank);
    for ( int eldep: domainElementDepSet ) {
        fprintf(stderr, "%d ", eldep);
    }

    fprintf(stderr, "\n");
#endif

    // now exclude local elements (local dependency is always satisfied)
    for ( _i = 1; _i <= nelems; _i++ ) {
        elem = domain->giveElement(_i);
        if ( elem->giveParallelMode() == Element_local ) {
            domainElementDepSet.erase( elem->giveGlobalNumber() );
        }
    }

#if NonlocalMaterialWTP_DEBUG_PRINT
    fprintf(stderr, "[%d] remote elem wish list:", myrank);
    for ( int eldep: domainElementDepSet ) {
        fprintf(stderr, "%d ", eldep);
    }

    fprintf(stderr, "\n");
#endif

    // broadcast remaining elements (unsatisfied domain nonlocal dependency) to remaining partitions
    _locsize = domainElementDepSet.size() + 1;
    result = MPI_Allreduce(& _locsize, & _globsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if ( result != MPI_SUCCESS ) {
        OOFEM_ERROR("MPI_Allreduce to determine  broadcast buffer size failed");
    }

    commBuff.resize( commBuff.givePackSizeOfInt(_globsize) );
    // remote domain wish list
    std :: set< int >remoteWishSet;

    toSendList.resize(nproc);
    for ( i = 0; i < nproc; i++ ) { // loop over partitions
        commBuff.init();
        toSendList [ i ].clear();
        if ( i == myrank ) {
            // current domain has to send its receive wish list to all domains
            commBuff.write(_locsize);
            for ( int eldep: domainElementDepSet ) {
                commBuff.write(eldep);
            }

            result = commBuff.bcast(i);
        } else {
            // unpack remote domain wish list
            remoteWishSet.clear();
            result = commBuff.bcast(i);
            // unpack size
            commBuff.read(_size);
            for ( _i = 1; _i < _size; _i++ ) {
                commBuff.read(_val);
                remoteWishSet.insert(_val);
            }

            // determine which local elements are to be sent to remotepartition
            for ( _i = 1; _i <= nelems; _i++ ) {
                elem = domain->giveElement(_i);
                if ( elem->giveParallelMode() == Element_local ) {
                    if ( remoteWishSet.find( elem->giveGlobalNumber() ) != remoteWishSet.end() ) {
                        // store local element number
                        toSendList [ i ].push_back(_i);
                    }
                }
            }
        }
    } // end loop over partitions broadcast

#if NonlocalMaterialWTP_DEBUG_PRINT
    for ( i = 0; i < nproc; i++ ) { // loop over partitions
        // print some info
        fprintf(stderr, "[%d] elements scheduled for mirroring at [%d]:",
                myrank, i);
        for ( int elnum: toSendList [ i ] ) {
            fprintf( stderr, "%d[%d] ", elnum, domain->giveElement(elnum)->giveGlobalNumber() );
        }

        fprintf(stderr, "\n");
    }

#endif





    com.packAllData(this, domain, & NonlocalMaterialWTP :: packRemoteElements);
    com.initExchange(MIGRATE_REMOTE_ELEMENTS_TAG);
    com.unpackAllData(this, domain, & NonlocalMaterialWTP :: unpackRemoteElements);
    com.finishExchange();

    domain->commitTransactions( domain->giveTransactionManager() );

#ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("NonlocalMaterialWTP::migrate", "Finished migrating remote elements", myrank);
#endif
}


void
NonlocalMaterialWTP :: update()
{
    /* Now the question is how to use nonlocElementDependencyMap, which is available for
     * each element, to fastly reinitialize nonlocal integration tables.
     *
     * if not needed, should be deleted at the end of migrate method, to free memory
     */
    this->fastRebuildNonlocalTables();
    // delete  element dep arrays
    nonlocElementDependencyMap.clear();
}


int NonlocalMaterialWTP :: packMigratingElementDependencies(Domain *d, ProcessCommunicator &pc)
{
    int myrank = d->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();

    if ( iproc == myrank ) {
        return 1;                // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();

    int ielem, nelem = d->giveNumberOfElements();
    int _globnum;
    Element *elem;

    for ( ielem = 1; ielem <= nelem; ielem++ ) { // begin loop over elements
        elem = d->giveElement(ielem);
        if ( ( elem->giveParallelMode() == Element_local ) &&
            ( lb->giveElementPartition(ielem) == iproc ) ) {
            // pack local element (node numbers shuld be global ones!!!)
            // pack type
            _globnum = elem->giveGlobalNumber();
            pcbuff->write(_globnum);
            nonlocElementDependencyMap [ _globnum ].storeYourself(*pcbuff);
        }
    } // end loop over elements

    // pack end-of-element-record
    pcbuff->write(NonlocalMaterialWTP_END_DATA);

    return 1;
}

int NonlocalMaterialWTP :: unpackMigratingElementDependencies(Domain *d, ProcessCommunicator &pc)
{
    int myrank = d->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();
    int _globnum;

    if ( iproc == myrank ) {
        return 1;                // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();

    // unpack element data
    do {
        pcbuff->read(_globnum);
        if ( _globnum == NonlocalMaterialWTP_END_DATA ) {
            break;
        }

        nonlocElementDependencyMap [ _globnum ].restoreYourself(*pcbuff);
    } while ( 1 );

    return 1;
}

int NonlocalMaterialWTP :: packRemoteElements(Domain *d, ProcessCommunicator &pc)
{
    int myrank = d->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();
    int nnodes, inode;
    DofManager *node, *dofman;
    Element *elem;

    if ( iproc == myrank ) {
        return 1;                // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();

    // here we have to pack also nodes that are shared by packed elements !!!
    // assemble set of nodes needed by those elements
    // these have to be send (except those that are shared)
    std :: set< int >nodesToSend;
    for ( int ie: toSendList [ iproc ] ) {
        //ie = d->elementGlobal2Local(gie);
        elem = d->giveElement(ie);
        nnodes = elem->giveNumberOfDofManagers();
        for ( int i = 1; i <= nnodes; i++ ) {
            node = elem->giveDofManager(i);
            if ( ( node->giveParallelMode() == DofManager_local ) ||
                ( node->isShared() && !node->givePartitionList()->contains(iproc) ) ) {
                nodesToSend.insert( node->giveGlobalNumber() );
            }
        }
    }

    // pack nodes that become null nodes on remote partition
    for ( int in: nodesToSend ) {
        inode = d->dofmanGlobal2Local(in);
        dofman = d->giveDofManager(inode);
        pcbuff->write( dofman->giveInputRecordName() );
        dofman->saveContext(*pcbuff, CM_Definition | CM_State | CM_UnknownDictState);
    }

    pcbuff->write("");

    for ( int ie: toSendList [ iproc ] ) {
        //ie = d->elementGlobal2Local(gie);
        elem = d->giveElement(ie);
        // pack local element (node numbers shuld be global ones!!!)
        // pack type
        pcbuff->write( elem->giveInputRecordName() );
        elem->saveContext(*pcbuff, CM_Definition | CM_DefinitionGlobal | CM_State);
    }

    pcbuff->write("");


    return 1;
}

int NonlocalMaterialWTP :: unpackRemoteElements(Domain *d, ProcessCommunicator &pc)
{
    int myrank = d->giveEngngModel()->giveRank();
    int iproc = pc.giveRank();
    std :: string _type;
    DofManager *dofman;
    IntArray _partitions;

    if ( iproc == myrank ) {
        return 1;                // skip local partition
    }

    // query process communicator to use
    ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();

    // unpack dofman data
    do {
        pcbuff->read(_type);
        if ( _type.size() == 0 ) {
            break;
        }
        dofman = classFactory.createDofManager(_type.c_str(), 0, d);
        dofman->restoreContext(*pcbuff, CM_Definition | CM_State | CM_UnknownDictState);
        dofman->setParallelMode(DofManager_null);
        if ( d->dofmanGlobal2Local( dofman->giveGlobalNumber() ) ) {
            // record already exist
            delete dofman;
        } else {
            d->giveTransactionManager()->addDofManTransaction(DomainTransactionManager :: DTT_ADD,
                                                              dofman->giveGlobalNumber(),
                                                              dofman);
        }
    } while ( 1 );


    // unpack element data
    Element *elem;
    _partitions.resize(1);
    _partitions.at(1) = iproc;
    do {
        pcbuff->read(_type);
        if ( _type.size() == 0 ) {
            break;
        }

        elem = classFactory.createElement(_type.c_str(), 0, d);
        elem->restoreContext(*pcbuff, CM_Definition | CM_State);
        elem->setParallelMode(Element_remote);
        elem->setPartitionList(_partitions);
        d->giveTransactionManager()->addElementTransaction(DomainTransactionManager :: DTT_ADD,
                                                           elem->giveGlobalNumber(), elem);
    } while ( 1 );

    return 1;
}

/* Now the question is how to use nonlocElementDependencyMap, which is available for
 * each element, to quickly reinitialize nonlocal integration tables.
 *
 * if not needed, should be deleted at the end of migrate method, to free memory
 *
 * first the existing data should be cleared, and new ones initialized
 * profiting from nonlocElementDependencyMap, that is available for all
 * local elements.
 */

void
NonlocalMaterialWTP :: fastRebuildNonlocalTables()
{
    Domain *d = lb->giveDomain();
    int n, i, globnum, ie, nelem = d->giveNumberOfElements();
    IntArray localElementDep;
    Element *elem;

    // build nonlocal element dependency array for each element
    for ( ie = 1; ie <= nelem; ie++ ) {
        elem = d->giveElement(ie);
        if ( ( elem->giveParallelMode() == Element_local ) ) {
            IntArray localMap;
            // translate here nonlocElementDependencyMap[_globnum] to corresponding local numbers
            globnum = elem->giveGlobalNumber();
            n = nonlocElementDependencyMap [ globnum ].giveSize();
            localElementDep.resize(n);
            for ( i = 1; i <= n; i++ ) {
                localElementDep.at(i) = d->elementGlobal2Local( nonlocElementDependencyMap [ globnum ].at(i) );
            }

            elem->ipEvaluator(this, & NonlocalMaterialWTP :: fastElementIPNonlocTableUpdater, localElementDep);
        }
    }
}



void
NonlocalMaterialWTP :: fastElementIPNonlocTableUpdater(GaussPoint *gp, IntArray &map)
{
    Element *elem = gp->giveElement();
    NonlocalMaterialExtensionInterface *iface = static_cast< NonlocalMaterialExtensionInterface * >( elem->giveMaterial()->giveInterface(NonlocalMaterialExtensionInterfaceType) );
    if ( iface ) {
        iface->rebuildNonlocalPointTable(gp, & map);
    }
}
} // end namespace oofem
