/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/communicator.h,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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

//
// Class Communicator
//

#ifndef communicator_h
#define communicator_h
#ifdef __PARALLEL_MODE

#include "processcomm.h"
#include "error.h"

#ifdef __USE_MPI
#ifndef __MAKEDEPEND
#include "mpi.h"
#endif
#endif


/**
  The Communicator and corresponding buffers (represented by this class) 
  are separated in order to allow share the same buffer by several communicators.
  Here sharing means reusing for different but NON-OVERLAPPING communications.
  if communications overlap, the different instances of CommunicatorBuff should be used!
  The CommunicatorBuff objects are registered in corresponding communicator, 
  then if maps are available, comBuff should be resized and used in subsequent ops.
  
  The registration is necessary, otherwise before each send op the buffers (given probably as parameter)  
  will be resized again (size have to be computed again) and this is probably quite cost operation. 
  When comBuff will be registered, resize is needed only when maps change, and this will not occur frequently 
  (its even quite rare).
*/
class CommunicatorBuff
{
 public:
 protected:
  /// number of processes
  int size;
  /// array of process communicators
  ProcessCommunicatorBuff** processCommBuffs;
  
 public:
  CommunicatorBuff (int s, CommBuffType t=CBT_static);
  ~CommunicatorBuff ();
  
  /**
     Returns i-th process communicator buff. The process comm buffs are numbered from rank 0.
     @param i process communicator buff index [0..size-1]
     @return pointer to corresponding process communicator buff, NULL otherwise.
  */
  ProcessCommunicatorBuff*
    giveProcessCommunicatorBuff (int i) {if (i < size) return processCommBuffs[i]; else return NULL;}
};


/**
 Class representing communicator.
 It is usually attribute of  Engng model.
 Problem comunicator provides all services for communication with 
 associated remote problems. It manages several process (task) comunicators.

 The communicator mode determines the communication:

 (Static) The mode can be static, meaning that each node can assemble its communication maps
 independently (or by independent communication). This implies that the size of
 communication buffers is known in advance. Also if no data are planned to send to remote node, there
 is no communication with this node (both sender and receiver know that there will be no data to send).
 
 (Dynamic) In this case the communication pattern and the ammount of data sent between nodes is 
 not known in advance. This requires to use dynamic (packeted) buffering.
 
*/
class Communicator 
{
 protected:
 /// rank of process 
 int rank; 
 /// number of processes
 int size;
 /// array of process communicators
 ProcessCommunicator** processComms;
 /// Engng model
 EngngModel* engngModel;
 /// mode
 CommunicatorMode mode;
 

public:
 /** 
  Constructor. Creates new communicator associated 
  to partition with number (rank) irank.
  @param irank rank of associated partition
  @param size #number of colaborating processes
  @param mode communicator mode.
  */
 Communicator (EngngModel* emodel, CommunicatorBuff* b, int rank, int size, CommunicatorMode m=CommMode_Static);
 /// Destructor
 virtual ~Communicator () ;

 /**
  Returns i-th problem communicator. The problems are numbered from rank 0.
  @param i problem communicator index [0..size-1]
  @return pointer to corresponding communicator, NULL otherwise.
  */
 ProcessCommunicator*
   giveProcessCommunicator (int i) {if (i < size) return processComms[i]; else return NULL;}

 /**
  Pack all problemCommuncators data to their send buffers. 
  @param packFunc function used to pack nodal data in to buffer. 
  @see NlDEIDynamic_Unpack_func 
  */
 template <class T> int packAllData (T* ptr, int (T::*packFunc) (ProcessCommunicator&));
 /**
  Pack all problemCommuncators data to their send buffers. 
  @param packFunc function used to pack nodal data in to buffer. 
  @see NlDEIDynamic_Unpack_func 
  */
 template <class T> int packAllData (T* ptr, FloatArray* src, int (T::*packFunc) (FloatArray*, ProcessCommunicator&));
 /**
  Unpack all problemCommuncators data  from recv buffers.
  Waits  untill receive completion before unpacking buffer.
  @param unpackFunc function used to unpack nodal data from buffer. 
  @see NlDEIDynamic_Unpack_func 
  */
 template <class T> int unpackAllData (T* ptr, int (T::*unpackFunc) (ProcessCommunicator&));
 /**
  Unpack all problemCommuncators data  from recv buffers.
  Waits  untill receive completion before unpacking buffer.
  @param unpackFunc function used to unpack nodal data from buffer. 
  @see NlDEIDynamic_Unpack_func 
  */
 template <class T> int unpackAllData (T* ptr, FloatArray* dest, int (T::*unpackFunc) (FloatArray*, ProcessCommunicator&));
 /**
  Initializes data exchange with all problems. 
  if send or receive pool is empty, communication is not preformed.
  @param tag message tag
  */
 int initExchange (int tag);
 /**
  Initializes data send exchange with all problems. 
  if send  pool is empty, communication is not preformed.
  @param tag message tag
  */
 int initSend (int tag);
 /**
  Initializes data receive exchange with all problems. 
  if  receive pool is empty, communication is not preformed.
  @param tag message tag
  */
 int initReceive (int tag);

 /** 
  Clears all buffer contens.
  */
 void clearBuffers ();
 /**
  Service for setting up the communication patterns with other remote processes.
  Sets up the toSend and toRecv attributes in associated problem communicators.
  */
 virtual void setUpCommunicationMaps (EngngModel* pm) {}

 /// prints error message and exits
 void error (char* file, int line, char *format, ...) const ;

private:
};

template <class T> int
Communicator :: packAllData (T* ptr, int (T::*packFunc) (ProcessCommunicator&))
{
 int i = size, result = 1;

 if (size) 
   for (i=0; i< size; i++) result &= giveProcessCommunicator(i)->packData (ptr, packFunc);
 return result;
}


template <class T> int
Communicator :: packAllData (T* ptr, FloatArray* src, int (T::*packFunc) (FloatArray*, ProcessCommunicator&))
{
 int i = size, result = 1;

 if (size) 
   for (i=0; i< size; i++) result &= giveProcessCommunicator(i)->packData (ptr, src, packFunc);
 return result;
}

template <class T> int
Communicator :: unpackAllData (T* ptr, int (T::*unpackFunc) (ProcessCommunicator&))
{
 int i, received, num_recv = 0, result = 1;
 IntArray recvFlag (size);
 //MPI_Status status;

 for  (i=0; i<size; i++) {
  if (giveProcessCommunicator(i)->giveToRecvMap()->giveSize()) {
   recvFlag.at(i+1) = 1;
   num_recv ++;
  }
 }
 
 while (num_recv--) { 
  
  // wait for any completion
  while (1) {
   received = 0;
   for  (i=0; i<size; i++)  {
    if (recvFlag.at(i+1)) {
      //if (giveProcessCommunicator(i)->giveRecvBuff()->testCompletion()) {
      if (giveProcessCommunicator(i)->receiveCompleted()) {
      
#ifdef __VERBOSE_PARALLEL
       OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received data from partition %3d\n",
                       rank,"Communicator :: unpackAllData", i);
#endif
      
      recvFlag.at(i+1) = 0;
      result &= giveProcessCommunicator(i)->unpackData (ptr, unpackFunc);
      received = 1;
      break;
     }
    }
   }
   if (received) break;
  }
 }

#ifdef __VERBOSE_PARALLEL
 VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier started",rank)
#endif
 
 MPI_Barrier (MPI_COMM_WORLD);

#ifdef __VERBOSE_PARALLEL
 VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier finished",rank)
#endif
 
 return result;
}



template <class T> int
Communicator :: unpackAllData (T* ptr, FloatArray* dest, int (T::*unpackFunc) (FloatArray*, ProcessCommunicator&))
{
 int i, received, num_recv = 0, result = 1;
 IntArray recvFlag (size);
 //MPI_Status status;

 for  (i=0; i<size; i++) {
  if (giveProcessCommunicator(i)->giveToRecvMap()->giveSize()) {
   recvFlag.at(i+1) = 1;
   num_recv ++;
  }
 }
 
 while (num_recv--) { 
  
  // wait for any completion
  while (1) {
   received = 0;
   for  (i=0; i<size; i++)  {
    if (recvFlag.at(i+1)) {
      //if (giveProcessCommunicator(i)->giveRecvBuff()->testCompletion()) {
      if (giveProcessCommunicator(i)->receiveCompleted()) {
      
      
#ifdef __VERBOSE_PARALLEL
       OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received data from partition %3d\n",
                       rank,"Communicator :: unpackAllData", i);
#endif
      
      recvFlag.at(i+1) = 0;
      result &= giveProcessCommunicator(i)->unpackData (ptr, dest, unpackFunc);
      received = 1;
      break;
     }
    }
   }
   if (received) break;
  }
 }

#ifdef __VERBOSE_PARALLEL
 VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier started",rank)
#endif
 
 MPI_Barrier (MPI_COMM_WORLD);

#ifdef __VERBOSE_PARALLEL
 VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier finished",rank)
#endif
 
 return result;
}




#endif
#endif
