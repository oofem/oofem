/* $Header: /home/cvs/bp/oofem/oofemlib/src/combuff.h,v 1.5 2003/04/06 14:08:23 bp Exp $ */
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
#ifndef loadballancer_h
#ifdef __PARALLEL_MODE
#include "inputrecord.h"
#include "interface.h"

class Domain;
#define MIGRATE_LOAD_TAG       9998
#define LOADBALLANCER_END_DATA 9999
/**
   Abstract base class representing general load ballancer. The task of load ballancer is to 
   recover load ballance when running in parallel. This is achieved by moving work from busy 
   nodes to other nodes to achieve an equal distribution of work. 
   In general load ballancer should repartition the problem domain, taking into account several 
   criteria:
   - It should take into account different computational requirement of different elements
   - The new partitioning should minimize the cut (to minimize the communication)
   - The new partitioning should minimize data movement (the cost of repartitioning) by 
   preserving the locality as much as possible. In other words the new and existing partitioning
   should be "similar".
   - the ballancer should decide whether the cost of rebalancing is not higher than the cost of
   continuing without rebalancing.
 */
class LoadBallancer
{
 public:
  enum LoadBallancerTimerType {LBTT_ComputationTimer, LBTT_LoadBallancingTimer, LBTT_DataTransferTimer};
  enum DofManMode { DM_NULL, DM_Local, DM_Shared, DM_SharedMerge, DM_SharedNew, DM_SharedUpdate, DM_Remote};
 protected:
  Domain* domain;

 public:

  LoadBallancer (Domain* d);
  virtual ~LoadBallancer () {}
  
  

 /**@name Profiling routines */
 //@{
  void startTimer (LoadBallancerTimerType t) {}
  void stopTimer  (LoadBallancerTimerType t) {}
 //@}

 /**@name Load evaluation and imbalance detection methods*/
 //@{
 //@}

 /**@name Work transfer calculation methods  */
 //@{
  virtual void ballanceLoad () = 0;
 //@}
  
 /**@name Work migration methods */
 //@{
  virtual void migrateLoad() = 0;
 //@}

 /**@name Query methods after Load Ballancing */
 //@{
  /// Returns the label of dofmanager after load ballancing
  virtual DofManMode giveDofManState (int idofman) = 0;
  
  /// Returns the partition list of given dofmanager after load balancing
  virtual IntArray* giveDofManPartitions (int idofman) = 0;
  
  /// Returns the new partition number assigned to local element after LB
  virtual int giveElementPartition(int ielem) = 0;

 //@}
  ///Initializes receiver acording to object description stored in input record.
  virtual IRResultType initializeFrom (InputRecord* ir) {return IRRT_OK;}

};


class LoadBallancerElementInterface : public Interface
{
 public:
  LoadBallancerElementInterface () {}
  
  virtual double predictRelativeComputationalCost ();
};

#endif
#define loadballancer_h
#endif
