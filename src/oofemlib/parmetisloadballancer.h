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
#ifndef parmetisloadbalancer_h
#ifdef __PARALLEL_MODE
#include "loadbalancer.h"

#ifdef __PARMETIS_MODULE
#ifndef __MAKEDEPEND
#include "parmetis.h"
#include <vector>
#endif
#endif

#define PARMETISLB_END_DATA 9999
#define SHARED_DOFMAN_PARTITIONS_TAG 9998

/**
   
 */
class ParmetisLoadBalancer : public LoadBalancer
{
 protected:

#ifdef __PARMETIS_MODULE
  // element numbering maps
  IntArray gToLMap, lToGMap;
  idxtype *elmdist;
  int myGlobNumOffset;
  // partition weights (user input)
  float* tpwgts;
  // array of DofManMode(s)
  IntArray dofManState;
  // array of dof man partitions
  std::vector <IntArray>  dofManPartitions;
  /// partition vector of the locally-stored elements
  IntArray elementPart;
#endif

 public:
  ParmetisLoadBalancer (Domain *d);
  ~ParmetisLoadBalancer ();
  
  virtual void calculateLoadTransfer ();

#if 1
  virtual DofManMode giveDofManState (int idofman);
  virtual IntArray* giveDofManPartitions (int idofman);
  virtual int giveElementPartition(int ielem);
#endif
protected:

#ifdef __PARMETIS_MODULE
  void initGlobalParmetisElementNumbering ();
  int  giveLocalElementNumber (int globnum) {return gToLMap.at(globnum-myGlobNumOffset);}
  int  giveGlobalElementNumber(int locnum) {return lToGMap.at(locnum);}

  /**
     Label local partition nodes (thode that are local or shared).
     Labeling consist of assigning corresponding id that characterize the
     status of local dof managar after ballancing the load. Labeling determines 
     which of local nodes remain local, or became local on other partition, 
     or became shared, etc.
  */
  void labelDofManagers ();
  int  determineDofManState (int idofman, int myrank, int npart, IntArray* dofManPartitions);

 int packSharedDmanPartitions (ProcessCommunicator& pc) ;
 int unpackSharedDmanPartitions (ProcessCommunicator& pc) ;
 void addSharedDofmanPartitions (int _locnum, IntArray _partitions);
#endif

};


#endif
#define parmetisloadbalancer_h
#endif
