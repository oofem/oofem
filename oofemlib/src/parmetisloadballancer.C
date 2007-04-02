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

#include "parmetisloadballancer.h"
#include "domain.h"
#include "engngm.h"
#include "element.h"
#include "dofmanager.h"
#include "conTable.h"
#include <set>

#define ParmetisLoadBallancer_DEBUG_PRINT 1
  
  
ParmetisLoadBallancer::ParmetisLoadBallancer(Domain* d) : LoadBallancer (d)
{
#ifdef __PARMETIS_MODULE
  elmdist = NULL;
#endif
}

ParmetisLoadBallancer::~ParmetisLoadBallancer ()
{
#ifdef __PARMETIS_MODULE
  if (elmdist) delete elmdist;
#endif
}


#ifdef __PARMETIS_MODULE
void
ParmetisLoadBallancer::ballanceLoad ()
{

  idxtype *eind, *eptr, *xadj, *adjncy, *part, *vwgt, *vsize;
  int i, nlocalelems, eind_size, nelem= domain->giveNumberOfElements();
  int ndofman, idofman, numflag, ncommonnodes, options[4],ie,nproc;
  int edgecut, wgtflag, ncon;
  float ubvec[1], itr;
  Element* ielem;
  MPI_Comm communicator = MPI_COMM_WORLD;

  // init parmetis element numbering 
  this->initGlobalParmetisElementNumbering();
  // prepare data structures for ParMETIS_V3_Mesh2Dual
  // count the size of eind array
  eind_size = 0; nlocalelems = 0;
  for (i=1; i<=nelem; i++) {
    ielem = domain->giveElement(i);
    if (ielem->giveParallelMode() == Element_local) {
      nlocalelems++;
      eind_size += ielem->giveNumberOfDofManagers();
    }
  }
  // allocate eind and eptr arrays
  eind = new idxtype[eind_size];
  eptr = new idxtype[nlocalelems+1];
  if ((eind==NULL) || (eptr==NULL)) OOFEM_ERROR ("ParmetisLoadBallancer::ballanceLoad: failed to allocate eind and eptr arrays");
  
  // fill in the eind and eptr (mesh graph)
  int eind_pos=0, eptr_pos=0;
  for (i=1; i<=nelem; i++) {
    ielem = domain->giveElement(i);
    if (ielem->giveParallelMode() == Element_local) {
      eptr[eptr_pos]=eind_pos;
      ndofman = ielem->giveNumberOfDofManagers();
      for (idofman = 1; idofman <= ndofman; idofman++) {
        eind[eind_pos++] = ielem->giveDofManager(idofman)->giveGlobalNumber()-1;
      }
      eptr_pos++;
    }
  }
  // last rec
  eptr[nlocalelems]=eind_pos;

  // call ParMETIS_V3_Mesh2Dual to construct dual graph (in parallel)
  numflag = 0; ncommonnodes = 2;
  ParMETIS_V3_Mesh2Dual (elmdist, eptr, eind, &numflag, &ncommonnodes, &xadj, &adjncy, &communicator);
 
#ifdef ParmetisLoadBallancer_DEBUG_PRINT
  // DEBUG PRINT
  int myrank=domain->giveEngngModel()->giveRank();
  fprintf (stderr, "[%d] xadj:", myrank);
  for (i=0; i<= nlocalelems; i++) fprintf (stderr, " %d", xadj[i]);
  fprintf (stderr, "\n[%d] adjncy:", myrank);
  for (i=0; i< xadj[nlocalelems]; i++) fprintf (stderr, " %d", adjncy[i]);
  fprintf (stderr, "\n");  
#endif


  // setup imbalance tolerance for each vertex weight - ubvec param
  ubvec[0] = 1.05;
  // setup options array
  options[0]=1; // set to zero for default
  options[1]=1; // get timings
  options[2]=15; // random seed
  options[3]=1; // sub-domains and processors are coupled
  // set ratio of inter-proc communication compared to data redistribution time
  itr = 1000.0;
  // set partition weights if not provided
  if (tpwgts == NULL) {
    if ((tpwgts = new float[nproc])==NULL) OOFEM_ERROR ("ParmetisLoadBallancer::ballanceLoad: failed to allocate tpwgts");
    for (i=0; i< nproc; i++) tpwgts[i]=1.0;
  }
  // obtain vertices weights (element weights) representing relative computational cost 
  if ((vwgt = new idxtype[nlocalelems])==NULL) OOFEM_ERROR ("ParmetisLoadBallancer::ballanceLoad: failed to allocate vwgt");
  if ((vsize= new idxtype[nlocalelems])==NULL) OOFEM_ERROR ("ParmetisLoadBallancer::ballanceLoad: failed to allocate vsize");
  for (ie = 0, i=0; i< nelem; i++) {
    ielem = domain->giveElement(i+1);
    if (ielem->giveParallelMode() == Element_local) {
      vwgt[ie]    = ielem->predictRelativeComputationalCost();
      vsize[ie++] = ielem->predictRelativeRedistributionCost();
    }
  }
  nproc=domain->giveEngngModel()->giveNumberOfProcesses();

  wgtflag = 2;
  numflag = 0;
  ncon = 1;
  if ((part = new idxtype[nlocalelems]) == NULL) OOFEM_ERROR ("ParmetisLoadBallancer::ballanceLoad: failed to allocate part");
  // call ParMETIS ballancing routineParMETIS_V3_AdaptiveRepart
  ParMETIS_V3_AdaptiveRepart (elmdist, xadj, adjncy, vwgt, vsize, NULL, &wgtflag, &numflag, &ncon, &nproc, 
                              tpwgts, ubvec, &itr, options, &edgecut, part, &communicator);
  
#ifdef ParmetisLoadBallancer_DEBUG_PRINT
  // debug
  fprintf (stderr, "[%d] edgecut: %d part:", myrank, edgecut);
  for (i=0; i< nlocalelems; i++) fprintf (stderr, " %d", part[i]);
  fprintf (stderr, "\n");
#endif
  
  // delete allocated xadj, adjncy arrays by ParMETIS
  delete xadj;
  delete adjncy;
  delete vwgt;
}

void
ParmetisLoadBallancer::initGlobalParmetisElementNumbering ()
{
  int nproc=domain->giveEngngModel()->giveNumberOfProcesses();
  int myrank=domain->giveEngngModel()->giveRank();
  IntArray procElementCounts (nproc);

  //if (procElementCounts) delete procElementCounts;
  elmdist = new idxtype[nproc+1];
  if (elmdist == NULL) OOFEM_ERROR ("ParmetisLoadBallancer::initGlobalParmetisNumbering: failed to allocate elmdist array");
  
  // determine number of local elements for the receiver
  int i, nlocelem = 0, nelem = domain->giveNumberOfElements();
  int globnum;

  for (i=1; i<=nelem; i++)
    if (domain->giveElement(i)->giveParallelMode() == Element_local) nlocelem++;
  procElementCounts(myrank) = nlocelem;

  MPI_Allgather (&nlocelem, 1, MPI_INT, procElementCounts.givePointer(), 1, MPI_INT, MPI_COMM_WORLD);
  elmdist[0]=0;
  for (i=0; i<nproc; i++) elmdist[i+1]=elmdist[i]+procElementCounts(i);

  // we need to number elements sequentially on each partition (and we start from rank 0)
  // compute local offset
  myGlobNumOffset = 0;
  for (i=0; i<myrank; i++) myGlobNumOffset+= procElementCounts(i);
  
  /* assemble maps of local numbering
     map is necessary since we may have remote elements that are not 
     part of local domain for load ballancing purposes
  */
  globnum = myGlobNumOffset+1;
  lToGMap.resize(nelem); gToLMap.resize(nelem);
  for (i=1; i<= nelem; i++) {
    if (domain->giveElement(i)->giveParallelMode() == Element_local) {
      lToGMap.at(i) = globnum;
      gToLMap.at(globnum-myGlobNumOffset) = i;
      globnum++;
    } else {
      lToGMap.at(i) = 0;
    }
  }
}

void
ParmetisLoadBallancer::labelDofManagers()
{
  int idofman, ndofman = domain->giveNumberOfDofManagers();
  ConnectivityTable* ct = domain->giveConnectivityTable();
  const IntArray* dofmanconntable;
  DofManager *dofman;
  Element* ielem;
  dofManagerParallelMode dmode;
  std::set<int, std::less<int> > __dmanpartitions;
  int myrank=domain->giveEngngModel()->giveRank();
  //int nproc=domain->giveEngngModel()->giveNumberOfProcesses();
  int ie, npart, __i;
  
  std::set<int, std::less<int> > :: iterator it;

  // resize label array
  dofManState.resize(ndofman);
  // resize dof man partitions
  dofManPartitions.resize(ndofman);

#ifdef ParmetisLoadBallancer_DEBUG_PRINT
  int _cols = 0;
  fprintf (stderr, "[%d] DofManager labels:\n", myrank);
#endif
  
  // loop over local dof managers
  for (idofman=1; idofman <= ndofman; idofman++) {
    dofman = domain->giveDofManager (idofman);
    dmode = dofman -> giveParallelMode();
    if ((dmode == DofManager_local) || (dmode == DofManager_shared)) {
      dofmanconntable = ct->giveDofManConnectivityArray (idofman);
      __dmanpartitions.clear();
      for (ie=1; ie<=dofmanconntable->giveSize(); ie++) {
        ielem = domain->giveElement(dofmanconntable->at(ie));
        // assemble list of partitions sharing idofman dofmanager
        // set is used to include possibly repeated partition only once
        if (ielem->giveParallelMode() == Element_local) {
          __dmanpartitions.insert(elmdist[ie-1]);
        }
      }
      npart = __dmanpartitions.size();
      dofManPartitions[idofman-1].resize(__dmanpartitions.size());
      for (__i=1, it=__dmanpartitions.begin(); it!=__dmanpartitions.end(); it++)
        dofManPartitions[idofman-1].at(__i++) = *it;
      // determine its state after ballancing -> label
      dofManState.at(idofman) = this->determineDofManState (idofman, myrank, npart, &dofManPartitions[idofman-1]);
    } else {
      dofManState.at(idofman) = DM_NULL;
    }
#ifdef ParmetisLoadBallancer_DEBUG_PRINT
   fprintf (stderr, " | %d: ", idofman);
   if (dofManState.at(idofman) == DM_NULL)              fprintf (stderr, "NULL ");
   else if (dofManState.at(idofman) == DM_Local)        fprintf (stderr, "Local ");
   else if (dofManState.at(idofman) == DM_Shared)       fprintf (stderr, "Shared");
   else if (dofManState.at(idofman) == DM_Remote)       fprintf (stderr, "Remote");
   else if (dofManState.at(idofman) == DM_SharedMerge)  fprintf (stderr, "ShdMer");
   else if (dofManState.at(idofman) == DM_SharedNew)    fprintf (stderr, "ShdNew");
   else if (dofManState.at(idofman) == DM_SharedUpdate) fprintf (stderr, "ShdUpd");
   if ((++_cols % 4) == 0) fprintf (stderr, "\n");
#endif   
  }
}

int
ParmetisLoadBallancer::determineDofManState (int idofman, int myrank, int npart, IntArray* dofManPartitions)
{
  dofManagerParallelMode dmode = domain->giveDofManager(idofman)->giveParallelMode();
  int answer = DM_Local;

  if (npart == 1) {
    if (dofManPartitions[idofman-1].at(1) == myrank) {
      if (dmode == DofManager_local) {
        // local remains local 
        answer = DM_Local;
      } else if (dmode == DofManager_shared) {
        // shared remains shared on local partition
        answer = DM_Shared;
      }
    } else { 
      if (dmode == DofManager_local) {
        // local node now fully on remote partition
        answer = DM_Remote;
      } else if (dmode == DofManager_shared) {
        // shared node now assigned to remote partition -> merging
        answer = DM_SharedMerge;
      }
    }
  } else { // npart > 1 : node will be newly shared, since cut runs through
    if (dmode == DofManager_local) answer = DM_SharedNew;
    else if (dmode == DofManager_shared) answer = DM_SharedUpdate;
  }
  return answer;
}


LoadBallancer::DofManMode
ParmetisLoadBallancer::giveDofManState (int idofman)
{
  return (LoadBallancer::DofManMode) dofManState.at(idofman);
}


IntArray*
ParmetisLoadBallancer::giveDofManPartitions (int idofman)
{
  return &dofManPartitions[idofman-1];
}

int
ParmetisLoadBallancer::giveElementPartition (int ielem)
{
  return elmdist[ielem-1];
}

#else
void ParmetisLoadBallancer:: ballanceLoad () {}
#endif



