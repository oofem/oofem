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

LoadBallancer::LoadBallancer (Domain* d) 
{
  domain = d;
}

#ifdef __PARALLEL_MODE
void 
LoadBallancer::migrateLoad ()
{
  domain->migrateLoad(this);
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

  // debugging -> force rebalancing in first step
  if (emodel->giveCurrentStep()->isTheFirstStep()) relWallClockImbalance = 1.0;

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
#endif
