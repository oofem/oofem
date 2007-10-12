/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2006   Borek Patzak                                       



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

#include "loadballancer.h"
#include "domain.h"
#include "engngm.h"
#include "timer.h"
#include "mathfem.h"
#include "timestep.h"
#include "usrdefsub.h"

#ifdef __PARALLEL_MODE
#include "parallel.h"
#include "processcomm.h"
#include "datastream.h"
#include "communicator.h"
#include "domaintransactionmanager.h"
#include "nonlocalmatwtp.h"
#endif

LoadBallancer::LoadBallancer (Domain* d)  : wtpList (0)
{
  domain = d;
}

#ifdef __PARALLEL_MODE

IRResultType
LoadBallancer::initializeFrom (InputRecord* ir)
{
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                   // Required by IR_GIVE_FIELD macro

  IntArray wtp;
  IR_GIVE_OPTIONAL_FIELD (ir, wtp, IFT_LoadBallancer_wtp, "wtp"); // Macro

  this->initializeWtp (wtp);

  return IRRT_OK;
}

void
LoadBallancer::initializeWtp (IntArray& wtp) {

  int i, size = wtp.giveSize();
  WorkTransferPlugin *plugin;

  if (size) {
    wtpList.growTo (size);
    for (i=1; i<=size; i++) {
      if (wtp.at(i)==1) {
        plugin = new NonlocalMaterialWTP (this);
      } else {
        OOFEM_ERROR ("LoadBallancer::initializeWtp: Unknown work transfer plugin type");
      }
      wtpList.put(i, plugin);
    }
  }
}
                                      
                                    
void 
LoadBallancer::migrateLoad (Domain* d)
{
  // domain->migrateLoad(this);
  int i;
  int nproc= d->giveEngngModel()->giveNumberOfProcesses();
  int myrank=d->giveEngngModel()->giveRank();


  // initialize work transfer plugins before any transfer
  for (i=1; i<=wtpList.giveSize(); i++) {
    wtpList.at(i)->init (d);
  }

  CommunicatorBuff cb (nproc, CBT_dynamic);
  Communicator com (d->giveEngngModel(), &cb, myrank, nproc, CommMode_Dynamic);

  // move existing dofmans and elements, that will be local on current partition,
  // into local map
  com.packAllData (this, d, &LoadBallancer::packMigratingData);
  com.initExchange (MIGRATE_LOAD_TAG);
  
  // do something in between 
  d->initGlobalDofManMap ();
  d->initGlobalElementMap ();

  this->deleteRemoteDofManagers (d);
  this->deleteRemoteElements (d);
  
  // receive remote data
  com.unpackAllData (this, d, &LoadBallancer::unpackMigratingData);

  d->commitTransactions (d->giveTransactionManager());

#if 1
  // debug print
  int j, nnodes=d->giveNumberOfDofManagers(), nelems=d->giveNumberOfElements();
  fprintf (stderr, "\n[%d] Nodal Table\n", myrank);
  for (i=1; i<=nnodes; i++) {
    if (d->giveDofManager(i)->giveParallelMode()==DofManager_local) 
      fprintf (stderr, "[%d]: %5d[%d] local\n", myrank, i, d->giveDofManager(i)->giveGlobalNumber());
    else if (d->giveDofManager(i)->giveParallelMode()==DofManager_shared) {
      fprintf (stderr, "[%d]: %5d[%d] shared ", myrank, i, d->giveDofManager(i)->giveGlobalNumber());
      for (j=1; j<=d->giveDofManager(i)->givePartitionList()->giveSize(); j++) {
        fprintf (stderr, "%d ", d->giveDofManager(i)->givePartitionList()->at(j));
      }
      fprintf (stderr, "\n");
    }
  }
  
  fprintf (stderr, "\n[%d] Element Table\n", myrank);
  for (i=1; i<=nelems; i++) {
    fprintf (stderr, "%5d {", i);
    for (j=1; j<=d->giveElement(i)->giveNumberOfDofManagers(); j++)
      fprintf (stderr, "%d ", d->giveElement(i)->giveDofManager(j)->giveNumber());
    fprintf (stderr, "}\n");
  }
#endif

  // migrate work transfer plugin data 
  for (i=1; i<=wtpList.giveSize(); i++) {
    wtpList.at(i)->migrate ();
  }

  // update work transfer plugin data 
  for (i=1; i<=wtpList.giveSize(); i++) {
    wtpList.at(i)->update ();
  }

}

int
LoadBallancer::packMigratingData (Domain* d, ProcessCommunicator& pc) 
{
  int myrank = d->giveEngngModel()->giveRank();
  int iproc = pc.giveRank();
  int idofman, ndofman;
  classType dtype;
  DofManager* dofman;
  LoadBallancer::DofManMode dmode;

  //  **************************************************
  //  Pack migrating data to remote partition 
  //  **************************************************

  // pack dofManagers
  if (iproc == myrank) return 1; // skip local partition 
  // query process communicator to use
  ProcessCommunicatorBuff* pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream (pcbuff);
  // loop over dofManagers
  ndofman = d->giveNumberOfDofManagers();
  for (idofman=1; idofman <= ndofman; idofman++) {
    dofman = d->giveDofManager (idofman);
    dmode = this->giveDofManState(idofman);
    dtype = dofman->giveClassID();
    // sync data to remote partition 
    // if dofman already present on remote partition then there is no need to sync
    //if ((this->giveDofManPartitions(idofman)->findFirstIndexOf(iproc))) {
    if ((this->giveDofManPartitions(idofman)->findFirstIndexOf(iproc)) &&
       (!dofman->givePartitionList()->findFirstIndexOf(iproc))) { 
      pcbuff->packInt (dtype);
      pcbuff->packInt (dmode);
      pcbuff->packInt (dofman->giveGlobalNumber());

      // pack dofman state (this is the local dofman, not available on remote)
      /* this is a potential performance leak, sending shared dofman to a partition,
         in which is already shared does not require to send context (is already there)
         here for simplicity it is always send */
      dofman->saveContext (&pcDataStream, CM_Definition | CM_State | CM_UnknownDictState);
      // send list of new partitions
      pcbuff->packIntArray (*(this->giveDofManPartitions(idofman)));
    }
  }

  // pack end-of-dofman-section record
  pcbuff->packInt (LOADBALLANCER_END_DATA);
  
  int ielem, nelem = d->giveNumberOfElements();
  
  Element* elem;

  for (ielem=1; ielem<=nelem; ielem++) { // begin loop over elements
    elem = d->giveElement(ielem);
    if ((elem->giveParallelMode() == Element_local) && 
	(this->giveElementPartition(ielem) == iproc)) {
      // pack local element (node numbers shuld be global ones!!!)
      // pack type
      pcbuff->packInt (elem->giveClassID());
      // nodal numbers shuld be packed as global !!
      elem->saveContext (&pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State);
    }
  } // end loop over elements
  // pack end-of-element-record
  pcbuff->packInt (LOADBALLANCER_END_DATA);
  
  return 1;
}


int
LoadBallancer::unpackMigratingData (Domain* d, ProcessCommunicator& pc) 
{
  // create temp space for dofManagers and elements
  // merging should be made by domain ?
  // maps of new dofmanagers and elements indexed by global number

  // we can put local dofManagers and elements into maps (should be done before unpacking)
  // int nproc=this->giveEngngModel()->giveNumberOfProcesses();
  int myrank= d->giveEngngModel()->giveRank();
  int iproc = pc.giveRank();
  int _mode, _globnum, _type;
  bool _newentry;
  classType _etype;
  IntArray _partitions, local_partitions;
  //LoadBallancer::DofManMode dmode;
  DofManager* dofman;
  DomainTransactionManager* dtm = d->giveTransactionManager();

  //  **************************************************
  //  Unpack migrating data to remote partition 
  //  **************************************************

  if (iproc == myrank) return 1; // skip local partition 
  // query process communicator to use
  ProcessCommunicatorBuff* pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream (pcbuff);

  pcbuff->unpackInt (_type);
  // unpack dofman data
  while (_type != LOADBALLANCER_END_DATA) {
    _etype = (classType) _type;
    pcbuff->unpackInt (_mode);
    switch (_mode) {
    case LoadBallancer::DM_Remote: 
      // receiving new local dofManager
      pcbuff->unpackInt (_globnum);

      _newentry = false;
      if ((dofman = dtm->giveDofManager(_globnum)) == NULL) { 
        // data not available -> create a new one
        _newentry = true;
        dofman = CreateUsrDefDofManagerOfType (_etype, 0, d);
      }
      dofman->setGlobalNumber(_globnum);
      // unpack dofman state (this is the local dofman, not available on remote)
      dofman->restoreContext (&pcDataStream, CM_Definition | CM_State| CM_UnknownDictState);
      // unpack list of new partitions
      pcbuff->unpackIntArray (_partitions);
      dofman->setPartitionList (_partitions);
      dofman->setParallelMode (DofManager_local);
      // add transaction if new entry allocated; otherwise existing one has been modified via returned dofman
      if (_newentry) dtm->addTransaction (DomainTransactionManager::DTT_ADD, DomainTransactionManager::DCT_DofManager, _globnum, dofman);
      //dmanMap[_globnum] = dofman;
      break;

    case LoadBallancer::DM_Shared:      
      // receiving new shared dofManager, that was local on sending partition
      // should be received only once (from partition where was local)
      pcbuff->unpackInt (_globnum);

      _newentry = false;
      if ((dofman = dtm->giveDofManager(_globnum)) == NULL) { 
        // data not available -> mode should be SharedUpdate
        _newentry = true;
        dofman = CreateUsrDefDofManagerOfType (_etype, 0, d);
      }
      dofman->setGlobalNumber(_globnum);
      // unpack dofman state (this is the local dofman, not available on remote)
      dofman->restoreContext (&pcDataStream, CM_Definition | CM_State| CM_UnknownDictState);
      // unpack list of new partitions
      pcbuff->unpackIntArray (_partitions);
      dofman->setPartitionList (_partitions);
      dofman->setParallelMode (DofManager_shared);
#if __VERBOSE_PARALLEL
      fprintf (stderr, "[%d] received Shared new dofman [%d]\n", myrank, _globnum);
#endif
      // add transaction if new entry allocated; otherwise existing one has been modified via returned dofman
      if (_newentry) dtm->addTransaction (DomainTransactionManager::DTT_ADD, DomainTransactionManager::DCT_DofManager, _globnum, dofman);
      //dmanMap[_globnum] = dofman;
      break;
       
    default:
      OOFEM_ERROR ("LoadBallancer::unpackMigratingData: unexpected dof manager type");
    }
    // get next type record
    pcbuff->unpackInt (_type);
    
  } ; // while (_type != LOADBALLANCER_END_DATA);
  
  // unpack element data
  Element* elem;
  do {
    pcbuff->unpackInt (_type);
    if (_type == LOADBALLANCER_END_DATA) break;
    _etype = (classType) _type;
    elem = CreateUsrDefElementOfType (_etype, 0, d);
    elem->restoreContext (&pcDataStream, CM_Definition | CM_State);
    dtm->addTransaction (DomainTransactionManager::DTT_ADD, DomainTransactionManager::DCT_Element, elem->giveGlobalNumber(), elem);
    //recvElemList.push_back(elem);    
  } while (1);

  return 1;
}


/* will delete those dofmanagers, that were sent to remote partition and are locally owned here
   so they are no longer necessary (those with state equal to DM_Remote and DM_SharedMerge)
   This will update domain DofManager list as well as global dmanMap and physically deletes the remote dofManager
*/
void
LoadBallancer::deleteRemoteDofManagers (Domain* d)
{
  int i, ndofman =  d->giveNumberOfDofManagers();
  //LoadBallancer* lb = this->giveLoadBallancer();
  LoadBallancer::DofManMode dmode;
  DofManager* dman;
  int myrank= d->giveEngngModel()->giveRank();
  DomainTransactionManager* dtm = d->giveTransactionManager();
  // loop over local nodes
  
  for (i = 1; i<= ndofman; i++) {
    dmode = this->giveDofManState(i);
    if ((dmode == LoadBallancer::DM_Remote)) {
      // positive candidate found
      dtm->addTransaction (DomainTransactionManager::DTT_Remove, DomainTransactionManager::DCT_DofManager, d->giveDofManager (i)->giveGlobalNumber(), NULL);
      // dmanMap.erase (d->giveDofManager (i)->giveGlobalNumber());
      //dman = dofManagerList->unlink (i);
      //delete dman;
    } else if ((dmode == LoadBallancer::DM_NULL)) {
      // positive candidate found; we delete all null dof managers
      // they will be created by nonlocalmatwtp if necessary. 
      // potentially, they can be reused, but this will make the code too complex
      dtm->addTransaction (DomainTransactionManager::DTT_Remove, DomainTransactionManager::DCT_DofManager, d->giveDofManager (i)->giveGlobalNumber(), NULL);
    } else if (dmode == LoadBallancer::DM_Shared) {
      dman = d->giveDofManager (i);
      dman->setPartitionList (*(this->giveDofManPartitions(i)));
      dman->setParallelMode (DofManager_shared);
      if (!dman->givePartitionList()->findFirstIndexOf (myrank)) {
        dtm->addTransaction (DomainTransactionManager::DTT_Remove, DomainTransactionManager::DCT_DofManager, d->giveDofManager (i)->giveGlobalNumber(), NULL);
        //dmanMap.erase (this->giveDofManager (i)->giveGlobalNumber());
        //dman = dofManagerList->unlink (i);
        //delete dman;
      }
    } else if (dmode == LoadBallancer::DM_Local) {
      IntArray _empty(0);
      dman = d->giveDofManager (i);
      dman->setPartitionList(_empty);
      dman->setParallelMode (DofManager_local);
    } else {
      OOFEM_ERROR ("Domain::deleteRemoteDofManagers: unknown dmode encountered");
    }
  }
}

/* will delete those elements, that were sent to remote partition and are locally owned here
   so they are no longer necessary (those with state equal to DM_Remote and DM_SharedMerge)
   This will update domain DofManager list as well as global dmanMap and physically deletes the remote dofManager
*/
void
LoadBallancer::deleteRemoteElements (Domain* d)
{
  int i, nelem =  d->giveNumberOfElements();
  int myrank= d->giveEngngModel()->giveRank();
  //LoadBallancer* lb = this->giveLoadBallancer();
  DomainTransactionManager* dtm = d->giveTransactionManager();
  //Element* elem;

  // loop over local nodes
  
  for (i = 1; i<= nelem; i++) {
    if (this->giveElementPartition(i) != myrank) {
      // positive candidate found
      // this->deleteElement (i);  // delete and set entry to NULL
      dtm->addTransaction (DomainTransactionManager::DTT_Remove, DomainTransactionManager::DCT_Element, d->giveElement (i)->giveGlobalNumber(), NULL);
      //elem = elementList->unlink (i);        
      //dmanMap.erase (elem->giveGlobalNumber());
      //delete (elem);
    } else if (d->giveElement(i)->giveParallelMode() != Element_local) {
      dtm->addTransaction (DomainTransactionManager::DTT_Remove, DomainTransactionManager::DCT_Element, d->giveElement (i)->giveGlobalNumber(), NULL);
    }      
  }
}





LoadBallancerMonitor::LoadBallancerDecisionType 
WallClockLoadBallancerMonitor::decide()
{
  int i, nproc=emodel->giveNumberOfProcesses();
  int myrank=emodel->giveRank();
  long int *node_solutiontimes = new long int [nproc];
  long min_st, max_st;
  double relWallClockImbalance;
  long int absWallClockImbalance;

  if (node_solutiontimes == NULL) OOFEM_ERROR ("LoadBallancer::LoadEvaluation failed to allocate node_solutiontimes array");

  // compute wall solution time of my node
  long int mySolutionTime = emodel->giveTimer()->getWtime(EngngModelTimer::EMTT_SolutionStepTimer);
  
  // collect wall clock computational time
  MPI_Allgather (&mySolutionTime, 1, MPI_LONG, node_solutiontimes, 1, MPI_LONG, MPI_COMM_WORLD);
  // detect imbalance
  min_st = max_st = node_solutiontimes[0];
  for (i=0; i< nproc; i++) {
    min_st = min (min_st, node_solutiontimes[i]);
    max_st = max (max_st, node_solutiontimes[i]);
  }
  absWallClockImbalance = (max_st-min_st);
  if (min_st) {
    relWallClockImbalance = ((max_st-min_st)/min_st);
  } else {
    relWallClockImbalance = 0.0;
  }

  delete[] node_solutiontimes;

  // decide
  if ((absWallClockImbalance > 10) || (relWallClockImbalance > 0.1)) {
    OOFEM_LOG_RELEVANT ("[%d] LoadBallancer: wall clock imbalance rel=%.2f\%,abs=%lds, recovering load\n", myrank, relWallClockImbalance, absWallClockImbalance);
    return LBD_RECOVER;
  } else {
    OOFEM_LOG_RELEVANT ("[%d] LoadBallancer: wall clock imbalance rel=%.2f\%,abs=%lds, continuing\n", myrank, relWallClockImbalance, absWallClockImbalance);
    return LBD_CONTINUE;
  }
}



#else //__PARALLEL_MODE
void 
LoadBallancer::migrateLoad () {}

IRResultType
LoadBallancer::initializeFrom (InputRecord* ir) {return IRRT_OK;}
#endif


LoadBallancer::WorkTransferPlugin::WorkTransferPlugin (LoadBallancer* _lb) {lb=_lb;}
LoadBallancer::WorkTransferPlugin::~WorkTransferPlugin () {}
