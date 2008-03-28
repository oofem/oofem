/* $Header: /home/cvs/bp/oofem/sm/src/feticommunicator.h,v 1.2.4.1 2004/04/05 15:19:46 bp Exp $ */
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
#ifndef feticommunicator_h
#define feticommunicator_h

#ifdef __PARALLEL_MODE

#include "communicator.h"
#include "fetiboundarydofman.h"
#include "error.h"
#ifndef __MAKEDEPEND
#include <vector>
#endif

using namespace std;

class FETISolver;

/**
 Class representing communicator for FETI solver.
 It is attribute of FETI solver numerical method class running on master partition (rank equal to 0).
 This Communicator provides necessary services for communication with 
 associated slave partitions. It manages several domain communicators, each for communication with 
 particular partition. 
*/
class FETICommunicator : public Communicator
{
public:
 /// Enumeration used to define necessary communication tags, used to identify different messages send/received
 enum {FETICommunicatorZeroTag, NumberOfBounadryDofManagersMsg, BoundaryDofManagersRecMsg};

protected:
 /// Number of equations at master level (determined form boundary nodes)
 int numberOfEquations;
 /// list of bounadary dof managers records
 vector<FETIBoundaryDofManager> boundaryDofManList;
 /// mater communication map. Not stored in corressponding domain comm, but required in order to
 /// allow direc (memory) mapping instead of communication.
 IntArray masterCommMap;

public:


 /** 
  Constructor. Creates new communicator associated 
  to master  with rank 0.
  @param irank rank of associated partition
  @param size #number of colaborating processes
  @param mode communicator mode.
  */
 FETICommunicator (EngngModel* emodel, CommunicatorBuff* b, int rank, int size);
 /// Destructor
 ~FETICommunicator () ;


 int giveNumberOfEquations () {return numberOfEquations;}
 /**
  Service for setting up the communication patterns with other remote domains.
  Sets up the toSend and toRecv attributes in associated domain communicators.
  */
 void setUpCommunicationMaps (EngngModel* pm);
 /**
  Returns reference to i-th boundary dof manager.
  Available only on master partition.
  */
 FETIBoundaryDofManager* giveDofManager (int i) {
  return &boundaryDofManList[i-1];}
 /**
  Returns pointer to master comm map stored in receiver.
  */
 IntArray* giveMasterCommMapPtr () {return &masterCommMap;}
 /**
  Assigns  given ToSendMap  to given domain communicator. 
  Sorts map according to global entity (dofmanagers or element) numbers to ensure, that local and 
  corresponding remote domain have the identical map with same ordering. This will ensure proper packing/unpacking
  order. The corresponding domainCommunicator buffer takes care about resizing itself accordingly to hold all 
  outcomming/incomming data using engngModel->estimateMaxPackSize service.
  @param domainComm domain comm which send map will be set
  @param map send map
  */
 //int setDomainCommunicatorToSendArry (DomainCommunicator<PNlDEIDynamic>* domainComm, IntArray& map);
 /**
  Assigns  given ToRecvMap  to given domain communicator. 
  Sorts map according to global entity (dofmanagers or element) numbers to ensure, that local and 
  corresponding remote domain have the identical map with same ordering. This will ensure proper packing/unpacking
  order. The corresponding domainCommunicator buffer takes care about resizing itself  accordingly to hold all 
  outcomming/incomming data using engngModel->estimateMaxPackSize service.
  @param domainComm domain comm which Recv map will be set
  @param map recv map
  */
 //int setDomainCommunicatorToRecvArry (DomainCommunicator<PNlDEIDynamic>* domainComm, IntArray& map);
private:

 /**
  Sorts given communication map, containing local DofManager numbers acording to their 
  corresponding global numbers. It could not be sorted by standart techniques, because
  it is necessary to ask DofMAnager form domain and determine its global Number.
  @param cmp comparison function must return a negative value if first argument is less than the second,
  zero if the arguments are equal, and a positive number otherwise.
  */
 //void  sortCommMap (IntArray& map, int (PNlDEIDynamicCommunicator::*cmp) (int,int));
 /// Implementation of Quiksort algorithm
 //void  quickSortCommMap   (IntArray& map, int l, int r, int (PNlDEIDynamicCommunicator::*cmp) (int,int));
 /// Partitioning used in quiksort
 //int   quickSortPartition (IntArray& map, int l, int r, int (PNlDEIDynamicCommunicator::*cmp) (int,int));
 
 /// global dofManager number comparison function
 //int DofManCmp (int , int);
 /// global element comparison function
 //int ElemCmp (int , int);

 /**
  Service for setting up the receiver for node cut communication patterns with other remote domains.
  Sets up the toSend and toRecv attributes in associated domain communicators.
  */
 //void setUpCommunicationMapsForNodeCut (Domain* domain);
  /**
  Service for setting up the receiver for element cut communication patterns with other remote domains.
  Sets up the toSend and toRecv attributes in associated domain communicators.
  */
 //void setUpCommunicationMapsForElementCut (Domain* domain);
 /**
  Service for setting up the receiver for remote element communication patterns with other remote domains.
  Sets up the toSend and toRecv attributes in associated domain communicators.
  */
 //void setUpCommunicationMapsForRemoteElementMode (Domain* domain);


};

#endif
#endif // feticommunicator_h

