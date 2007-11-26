/* $Header: /home/cvs/bp/oofem/oofemlib/src/domain.C,v 1.31.4.2 2004/05/14 13:45:27 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/
#include "nonlocalmatwtp.h"
#include "nonlocalmaterialext.h"
#include "element.h"
#include "gausspnt.h"
#include "material.h"
#include "communicator.h"
#include "datastream.h"
#include "domaintransactionmanager.h"
#include "usrdefsub.h"
#ifndef __MAKEDEPEND
#include <set>
#endif

#define NonlocalMaterialWTP_DEBUG_PRINT 0

/*
  Returns array storing nonlocal dependency 
  (in terms of global element numbers) for given element
*/
void
NonlocalMaterialWTP::giveElementNonlocalDepArry (IntArray& answer, Domain* d, int num)
{
  int _i;
  std::set <int> relems;
  std::set<int>::const_iterator relemsIter;
  Element *ielem = d->giveElement (num);

  if (ielem->giveMaterial()->giveInterface(NonlocalMaterialExtensionInterfaceType)) {
    relems.clear();
    // loop over element IRules and their IPs to retrieve remote (nonlocal) elements
    // store their global numbers in the relems set (to avoid renundancy) 
    // and then keep them in nonlocTables array.
    ielem->ipEvaluator (this, &NonlocalMaterialWTP::giveNonlocalDepArryElementPlugin, relems);

    // convert relems set into an int array
    // and store it
    answer.resize(relems.size());
    for (relemsIter=relems.begin(), _i=1; relemsIter != relems.end(); ++relemsIter, _i++) {
      answer.at(_i) = *relemsIter;
    }
  } else {
    answer.resize(0);
  }
}


void
NonlocalMaterialWTP::giveNonlocalDepArryElementPlugin (GaussPoint* gp, std::set<int>& s)
{
  int remoteElemNum;
  
  NonlocalMaterialStatusExtensionInterface* interface = 
    (NonlocalMaterialStatusExtensionInterface*) (gp->giveElement()->giveMaterial()->
    giveStatus(gp)->giveInterface(NonlocalMaterialStatusExtensionInterfaceType));
  if (interface) {
    dynaList<localIntegrationRecord>* lir = interface->giveIntegrationDomainList();
    dynaList<localIntegrationRecord>::iterator listIter;
    
    for (listIter = lir->begin(); listIter!= lir->end(); ++listIter) {
      remoteElemNum = ((*listIter).nearGp)->giveElement()->giveGlobalNumber();
      s.insert (remoteElemNum);
    }
  }
}


/*
  prepares the communication maps for remote elements
  should be called immediately after load ballancing,
  before any work transfer.

*/
void
NonlocalMaterialWTP::init (Domain* domain)
{
  int ie, gie, nelem = domain->giveNumberOfElements();
  EngngModel* emodel = domain->giveEngngModel();
  Element* elem;
  int nproc= emodel->giveNumberOfProcesses();
  int myrank= emodel->giveRank();
  CommunicatorBuff cb (nproc, CBT_dynamic);
  Communicator com (emodel, &cb, myrank, nproc, CommMode_Dynamic);

  // build nonlocal element dependency array for each element 
  for (ie=1; ie<=nelem; ie++) {
    elem = domain->giveElement(ie);
    if ((elem->giveParallelMode() == Element_local)) {
      gie = elem->giveGlobalNumber();
      this->giveElementNonlocalDepArry (nonlocElementDependencyMap[gie], domain, ie);
    }
  }

  /* send and receive nonlocElementDependencyArry of migrating elements to remote partition */
  com.packAllData (this, domain, &NonlocalMaterialWTP::packMigratingElementDependencies);
  com.initExchange (MIGRATE_NONLOCALDEP_TAG);
  com.unpackAllData (this, domain, &NonlocalMaterialWTP::unpackMigratingElementDependencies);
  
}



/* 
   should be called after basic local migration is finalized,
   when all local elements are already available 
*/
void
NonlocalMaterialWTP::migrate ()
{
  Domain* domain = this->lb->giveDomain();
  EngngModel* emodel = domain->giveEngngModel();
  int nproc= emodel->giveNumberOfProcesses();
  int myrank= emodel->giveRank();
  CommunicatorBuff cb (nproc, CBT_dynamic);
  Communicator com (emodel, &cb, myrank, nproc, CommMode_Dynamic);
  StaticCommunicationBuffer commBuff (MPI_COMM_WORLD);

  /* 
     build domain nonlocal element dependency list. Then exclude local elements - what remains are unsatisfied 
     remote dependencies that have to be broadcasted and received from partitions owning relevant elements 
  */
  int _locsize, i, _i, ie, _size, _globnum, result, nelems=domain->giveNumberOfElements();
  int _globsize, _val;
  Element* elem;
  std::set <int> domainElementDepSet;
  std::set <int>::const_iterator sit;
  // loop over each element dep list to assemble domain list
  std::vector <IntArray>::const_iterator it;
  for (ie=1; ie <= nelems; ie++) {
    elem = domain->giveElement(ie);
    if ((elem->giveParallelMode() == Element_local)) {
      _globnum = elem->giveGlobalNumber();
      IntArray &iedep = nonlocElementDependencyMap[_globnum];
      _size=iedep.giveSize();
      for (_i=1; _i<=_size; _i++) domainElementDepSet.insert(iedep.at(_i));
#if NonlocalMaterialWTP_DEBUG_PRINT
  fprintf (stderr, "[%d] element %d dependency:", myrank, _globnum);
  for (_i=1; _i<=_size; _i++) fprintf (stderr, "%d ", iedep.at(_i));
  fprintf (stderr, "\n");
#endif
    }
  }
#if NonlocalMaterialWTP_DEBUG_PRINT
  fprintf (stderr, "[%d] nonlocal domain dependency:", myrank);
  for (sit=domainElementDepSet.begin(); 
       sit != domainElementDepSet.end(); ++sit) {
    fprintf (stderr, "%d ", *sit);
  }
  fprintf (stderr, "\n");
#endif
  
  // now exclude local elements (local dependency is always satisfied)
  for (_i=1; _i<=nelems; _i++) {
    elem = domain->giveElement(_i);
    if (elem->giveParallelMode()==Element_local)
      domainElementDepSet.erase (elem->giveGlobalNumber());
  }
#if NonlocalMaterialWTP_DEBUG_PRINT
  fprintf (stderr, "[%d] remote elem wish list:", myrank);
  for (sit=domainElementDepSet.begin(); 
       sit != domainElementDepSet.end(); ++sit) {
    fprintf (stderr, "%d ", *sit);
  }
  fprintf (stderr, "\n");
#endif

  // broadcast remaining elements (unsatisfied domain nonlocal dependency) to remaining partitions
  _locsize = domainElementDepSet.size()+1;
  result = MPI_Allreduce (&_locsize, &_globsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (result != MPI_SUCCESS) OOFEM_ERROR ("NonlocalMaterialWTP::migrate: MPI_Allreduce to determine  broadcast buffer size failed");
  commBuff.resize (commBuff.givePackSize(MPI_INT,_globsize));
  // remote domain wish list
  std::set <int> remoteWishSet;
  
  toSendList.resize(nproc);
  for (i=0; i<nproc; i++) { // loop over partitions
    commBuff.init();
    toSendList[i].clear();
    if (i == myrank) {
      // current domain has to send its receive wish list to all domains
      commBuff.packInt (_locsize);
      for (sit=domainElementDepSet.begin(); sit!=domainElementDepSet.end(); ++sit)
        commBuff.packInt(*sit);
      result = commBuff.bcast (i);
    } else {
      // unpack remote domain wish list
      remoteWishSet.clear();
      result = commBuff.bcast (i);
      // unpack size
      commBuff.unpackInt(_size);
      for (_i=1; _i<=_size; _i++) {
        commBuff.unpackInt(_val);
        remoteWishSet.insert (_val);
      }
      // determine which local elements are to be sent to remotepartition
      for (_i=1; _i<=nelems; _i++) {
        elem = domain->giveElement(_i);
        if (elem->giveParallelMode()==Element_local) {
          if (remoteWishSet.find(elem->giveGlobalNumber()) !=remoteWishSet.end()) {
            // store local element number 
            toSendList[i].push_back(_i);
          }
        }
      }
    }
  } // end loop over partitions broadcast


#if NonlocalMaterialWTP_DEBUG_PRINT
  std::list<int>::const_iterator lit;
  for (i=0; i<nproc; i++) { // loop over partitions
    // print some info
    fprintf (stderr, "[%d] elements scheduled for mirroring at [%d]:",
	     myrank,i);
    for (lit=toSendList[i].begin(); lit!=toSendList[i].end(); ++lit) {
      fprintf (stderr, "%d[%d] ", *lit, domain->giveElement(*lit)->giveGlobalNumber());
    }
    fprintf (stderr, "\n");
  }
#endif





  com.packAllData (this, domain, &NonlocalMaterialWTP::packRemoteElements);
  com.initExchange (MIGRATE_REMOTE_ELEMENTS_TAG);
  com.unpackAllData (this, domain, &NonlocalMaterialWTP::unpackRemoteElements);

  domain->commitTransactions (domain->giveTransactionManager());
#ifdef __VERBOSE_PARALLEL
   VERBOSEPARALLEL_PRINT("NonlocalMaterialWTP::migrate", "Finished migrating remote elements", myrank);
#endif
}


void
NonlocalMaterialWTP::update ()
{
  /* Now the question is how to use nonlocElementDependencyMap, which is available for
     each element, to fastly reinitialize nonlocal integration tables.

     if not needed, should be deleted at the end of migrate method, to free memory
  */
  this->fastRebuildNonlocalTables();
  // delete  element dep arrays
  nonlocElementDependencyMap.clear();
}


int NonlocalMaterialWTP::packMigratingElementDependencies (Domain* d, ProcessCommunicator& pc) 
{
  int myrank = d->giveEngngModel()->giveRank();
  int iproc = pc.giveRank();

  if (iproc == myrank) return 1; // skip local partition 
  // query process communicator to use
  ProcessCommunicatorBuff* pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream (pcbuff);

  int ielem, nelem = d->giveNumberOfElements();
  int _globnum;
  Element* elem;

  for (ielem=1; ielem<=nelem; ielem++) { // begin loop over elements
    elem = d->giveElement(ielem);
    if ((elem->giveParallelMode() == Element_local) && 
        (lb->giveElementPartition(ielem) == iproc)) {
      // pack local element (node numbers shuld be global ones!!!)
      // pack type
      _globnum = elem->giveGlobalNumber();
      pcbuff->packInt(_globnum);
      pcbuff->packIntArray (nonlocElementDependencyMap[_globnum]);
    }
  } // end loop over elements
  // pack end-of-element-record
  pcbuff->packInt (NonlocalMaterialWTP_END_DATA);
  
  return 1;

}

int NonlocalMaterialWTP::unpackMigratingElementDependencies (Domain* d, ProcessCommunicator& pc)
{
  int myrank= d->giveEngngModel()->giveRank();
  int iproc = pc.giveRank();
  int _globnum;
  
  if (iproc == myrank) return 1; // skip local partition 
  // query process communicator to use
  ProcessCommunicatorBuff* pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream (pcbuff);

  // unpack element data
  do {
    pcbuff->unpackInt (_globnum);
    if (_globnum == NonlocalMaterialWTP_END_DATA) break;
    pcbuff->unpackIntArray (nonlocElementDependencyMap[_globnum]);
  } while (1);
  
  return 1;
}

int NonlocalMaterialWTP::packRemoteElements (Domain* d, ProcessCommunicator& pc) 
{
  int myrank = d->giveEngngModel()->giveRank();
  int iproc = pc.giveRank();
  int i, ie, nnodes, inode;
  DofManager *node, *dofman;
  Element *elem;
  classType dtype;


  if (iproc == myrank) return 1; // skip local partition 
  // query process communicator to use
  ProcessCommunicatorBuff* pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream (pcbuff);
  std::list<int>::const_iterator it;

  // here we have to pack also nodes that are shared by packed elements !!!
  // assemble set of nodes needed by those elements
  // these have to be send (except those that are shared)
  std::set <int> nodesToSend;
  for (it=toSendList[iproc].begin(); it!=toSendList[iproc].end(); ++it) {
    ie = *it; //ie = d->elementGlobal2Local(*it);
    elem = d->giveElement(ie);
    nnodes = elem->giveNumberOfDofManagers();
    for (i=1; i<=nnodes; i++) {
      node = elem->giveDofManager(i);
      if ((node->giveParallelMode()==DofManager_local) ||
          (node->isShared() && !node->givePartitionList()->contains(iproc))) {
        nodesToSend.insert(elem->giveDofManager(i)->giveGlobalNumber());
      }
    }
  }
  // pack nodes that become null nodes on remote partition 
  std::set <int>::const_iterator nit;
  for (nit=nodesToSend.begin(); nit !=nodesToSend.end(); ++nit) {
    inode = d->dofmanGlobal2Local (*nit);
    dofman  = d->giveDofManager (inode);
    dtype = dofman->giveClassID();
    pcbuff->packInt (dtype);
    dofman->saveContext (&pcDataStream, CM_Definition | CM_State | CM_UnknownDictState);
  }
  pcbuff->packInt (NonlocalMaterialWTP_END_DATA);

  for (it=toSendList[iproc].begin(); it!=toSendList[iproc].end(); ++it) {
    ie=*it; //ie = d->elementGlobal2Local(*it);
    elem = d->giveElement(ie);
    // pack local element (node numbers shuld be global ones!!!)
    // pack type
    pcbuff->packInt (elem->giveClassID());
    elem->saveContext (&pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State);
  }
  pcbuff->packInt (NonlocalMaterialWTP_END_DATA);

  
  return 1;
}

int NonlocalMaterialWTP::unpackRemoteElements (Domain* d, ProcessCommunicator& pc)
{
  int myrank= d->giveEngngModel()->giveRank();
  int iproc = pc.giveRank();
  int _type;
  classType _etype;
  DofManager *dofman;
  IntArray _partitions;
  
  if (iproc == myrank) return 1; // skip local partition 
  // query process communicator to use
  ProcessCommunicatorBuff* pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream (pcbuff);

  // unpack null dofmans
  pcbuff->unpackInt (_type);
  // unpack dofman data
  while (_type != NonlocalMaterialWTP_END_DATA) {
    _etype = (classType) _type;
    dofman = ::CreateUsrDefDofManagerOfType (_etype, 0, d);
    dofman->restoreContext (&pcDataStream, CM_Definition | CM_State| CM_UnknownDictState);
    dofman->setParallelMode (DofManager_null);
    if (d->dofmanGlobal2Local(dofman->giveGlobalNumber())) {
      // record already exist
      delete dofman;
    } else {
      d->giveTransactionManager()->addTransaction (DomainTransactionManager::DTT_ADD, 
						   DomainTransactionManager::DCT_DofManager,
						   dofman->giveGlobalNumber(),
						   dofman);
    }
    
    // get next type record
    pcbuff->unpackInt (_type);
  }
  

  // unpack element data
  Element* elem;
  _partitions.resize(1); _partitions.at(1) = iproc;
  do {
    pcbuff->unpackInt (_type);
    if (_type == NonlocalMaterialWTP_END_DATA) break;
    _etype = (classType) _type;
    elem = ::CreateUsrDefElementOfType (_etype, 0, d);
    elem->restoreContext (&pcDataStream, CM_Definition | CM_State);
    elem->setParallelMode (Element_remote);
    elem->setPartitionList (_partitions);
    d->giveTransactionManager()->addTransaction (DomainTransactionManager::DTT_ADD, 
						 DomainTransactionManager::DCT_Element,
						 elem->giveGlobalNumber(),elem);
  } while (1);
  
  return 1;
}

/* Now the question is how to use nonlocElementDependencyMap, which is available for
   each element, to fastly reinitialize nonlocal integration tables.
   
   if not needed, should be deleted at the end of migrate method, to free memory

   first the existing data should be cleared, and new ones initialized
   profiting from nonlocElementDependencyMap, that is available for all
   local elements.
*/

void
NonlocalMaterialWTP::fastRebuildNonlocalTables()
{
  Domain *d = lb->giveDomain();
  int n, i, globnum, ie, nelem = d->giveNumberOfElements();
  IntArray localElementDep;
  Element *elem;

  // build nonlocal element dependency array for each element 
  for (ie=1; ie<=nelem; ie++) {
    elem = d->giveElement(ie);
    if ((elem->giveParallelMode() == Element_local)) {

      IntArray localMap;
      // translate here nonlocElementDependencyMap[_globnum] to corresponding local numbers 
      globnum = elem->giveGlobalNumber();
      n = nonlocElementDependencyMap[globnum].giveSize();
      localElementDep.resize(n);
      for (i=1; i<=n; i++) {
        localElementDep.at(i) = d->elementGlobal2Local(nonlocElementDependencyMap[globnum].at(i));
      }

      elem->ipEvaluator (this, &NonlocalMaterialWTP::fastElementIPNonlocTableUpdater,localElementDep);
    }
  }
}



void
NonlocalMaterialWTP::fastElementIPNonlocTableUpdater(GaussPoint* gp, IntArray& map)
{
  Element* elem = gp->giveElement();
  NonlocalMaterialExtensionInterface* iface = (NonlocalMaterialExtensionInterface*) elem->giveMaterial()->giveInterface(NonlocalMaterialExtensionInterfaceType);
  if (iface) iface->rebuildNonlocalPointTable (gp, &map);
}


