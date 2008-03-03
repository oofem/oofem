/* $Header: /home/cvs/bp/oofem/oofemlib/src/combuff.C,v 1.8.4.1 2004/04/05 15:19:43 bp Exp $ */
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


#ifdef __PARALLEL_MODE

#ifndef __MAKEDEPEND
#include <list>
#endif

#include "compiler.h"
#include "dyncombuff.h"
#include "mathfem.h"
#include "error.h"


CommunicationPacket :: CommunicationPacket(MPI_Comm comm,int size, int num) : MPIBuffer (max(size,__CommunicationPacket_DEFAULT_SIZE), false)
{
  this->EOF_Flag = false;
  this->number = num;
  // reserve space for packet header
  curr_pos += givePackSize(comm, MPI_INT, 2);
}


CommunicationPacket :: CommunicationPacket (MPI_Comm comm,int num) : MPIBuffer (__CommunicationPacket_DEFAULT_SIZE, false)
{
  this->EOF_Flag = false;
  this->number = num;
  // reserve space for packet header
  curr_pos += givePackSize(comm, MPI_INT, 2);
}


CommunicationPacket :: ~CommunicationPacket ()
{}


void
CommunicationPacket :: init (MPI_Comm comm)
{
  MPIBuffer::init();
  this->EOF_Flag = false;
  // reserve space for packet header
  curr_pos += givePackSize(comm, MPI_INT, 2);
  
}

int 
CommunicationPacket::iSend (MPI_Comm communicator, int dest, int tag)
{
  this->packHeader(communicator);
  return (MPI_Isend (this->buff, this->curr_pos, MPI_PACKED, dest, tag, 
                     communicator, &this->request) == MPI_SUCCESS);
}


int 
CommunicationPacket::iRecv (MPI_Comm communicator, int source, int tag, int count)
{
  if (count) {
    if (count >= this->size)  {
      // reallocate itself
      if (this->resize (count) == 0) return 0;
    }
  }
  return (MPI_Irecv (this->buff,this->size, MPI_PACKED, source, tag, 
                     communicator, &this->request) == MPI_SUCCESS);
}


int 
CommunicationPacket::packHeader(MPI_Comm comm)
{
  int _arry[2];
  int _res, _pos  =0;
  
  _arry[0]=this->number;
  _arry[1]=this->EOF_Flag;

  _res = MPI_Pack (_arry, 2, MPI_INT, this->buff, size, &_pos, comm);
  
  return (_res == MPI_SUCCESS);
}

int
CommunicationPacket::unpackHeader(MPI_Comm comm)
{
  int _arry[2];
  int _res, _pos  =0;

  _res = MPI_Unpack (this->buff, this->size, &_pos, _arry, 2, MPI_INT, comm);
  this->number =   _arry[0];
  this->EOF_Flag = _arry[1];

  return (_res == MPI_SUCCESS);
}



/********DynamicCommunicationBuffer***************/

CommunicationPacketPool DynamicCommunicationBuffer::packetPool;

DynamicCommunicationBuffer::DynamicCommunicationBuffer (MPI_Comm comm, int size, bool dynamic) 
  : CommunicationBuffer (comm, size, dynamic), packet_list()
{
  number_of_packets = 0;
  /*
  // alocate first send/receive packet
  active_packet = this->allocateNewPacket (++number_of_packets);
  packet_list.push_back(active_packet);
  */
}


DynamicCommunicationBuffer::DynamicCommunicationBuffer (MPI_Comm comm, bool dynamic)
  : CommunicationBuffer (comm, dynamic), packet_list()
{
  number_of_packets = 0;
  /*
  // alocate first send/receive packet
  active_packet = this->allocateNewPacket (++number_of_packets);
  packet_list.push_back(active_packet);
  */
}

 /// Destructor.
DynamicCommunicationBuffer:: ~DynamicCommunicationBuffer ()
{
  this->clear();
}

void 
DynamicCommunicationBuffer::init ()
{
  this->clear();
}

void
DynamicCommunicationBuffer::initForPacking()
{
  this->clear();
  
  if (!active_packet) {
    active_packet = this->allocateNewPacket (++number_of_packets);
    packet_list.push_back(active_packet);
  }    
}

void 
DynamicCommunicationBuffer::initForUnpacking()
{
  recvIt = packet_list.begin();
  this->popNewRecvPacket();
}

/*
int
DynamicCommunicationBuffer::packArray (int* src, int n)
{
  int _result=1;
  int start_indx=0, end_indx, _size;

  do {
    _size = this->giveFitSize(MPI_INT, active_packet -> giveAvailableSpace(), n);
    end_indx = start_indx + _size;

    if (_size) _result &= active_packet -> packArray (communicator, src+start_indx,_size, MPI_INT);
    if (end_indx == n) break;
    // active packet full, allocate a new one
    active_packet = this->allocateNewPacket (++number_of_packets);
    packet_list.push_back(active_packet);
    start_indx = end_indx;
  } while (1);
  
  return _result;
}


int
DynamicCommunicationBuffer::unpackArray (int* dest, int n)
{
  int _result=1;
  int start_indx=0, end_indx, _size;

  do {
    _size = this->giveFitSize(MPI_INT, active_packet -> giveAvailableSpace(), n);
    end_indx = start_indx + _size;

    if (_size) _result &= active_packet->unpackArray (communicator,dest+start_indx,_size, MPI_INT);
    if (end_indx == n) break;
    // active packet exhausted, pop a new one
    this->popNewRecvPacket();
    start_indx = end_indx;
  } while (1);
  
  return _result;
}
*/

int 
DynamicCommunicationBuffer::iSend (int dest, int tag)
{
  int result = 1;
  
  std::list<CommunicationPacket*>::const_iterator it;
  /// set last (active) send packet as eof
  active_packet->setEOFFlag();

  active_rank = dest; active_tag = tag;
  for (it=packet_list.begin(); it != packet_list.end(); ++it) {
    result &= (*it)->iSend (communicator, dest, tag);
  }

  /*  
  int _myrank;
  MPI_Comm_rank (communicator, &_myrank);
  fprintf (stderr,"[%d] sending to [%d] %d packets for tag %d\n", _myrank, dest, number_of_packets, tag);
  (*(packet_list.begin()))->dump();
  */
  
 return result;

}



int 
DynamicCommunicationBuffer::iRecv (int source, int tag, int count)
{
  this->init();
  number_of_packets = 0;
  // receive first packet, but it is probably not the last one
  // create new first packet and init its receive
  active_packet = this->allocateNewPacket (++number_of_packets);
  active_rank = source; active_tag = tag;
  return active_packet->iRecv (communicator, source, tag);
}

int DynamicCommunicationBuffer::receiveCompleted ()
{
  /*
  int _myrank;
  MPI_Comm_rank (communicator, &_myrank);
  */

  if (active_packet->testCompletion()) {
    // active packet received, add it to the pool and unpach header info
    active_packet->unpackHeader(communicator);

    //fprintf (stderr, "[%d] received from [%d] packet no. %d\n", _myrank, active_rank, active_packet->getNumber());

    pushNewRecvPacket(active_packet);
    if (active_packet->hasEOFFlag()) {
      // last packet received; init for unpacking
      this->initForUnpacking();

      //fprintf (stderr,"[%d] received from [%d] %d packets for tag %d\n", _myrank, active_rank, number_of_packets, active_tag);
      //active_packet->dump();

      return 1;
    } else {
      // received next packet, but it is not the last one
      // create new packet and init its receive
      active_packet = this->allocateNewPacket (++number_of_packets);
      active_packet->iRecv (communicator, active_rank, active_tag);
      return 0;
    }
  } else {
    // active packet not yet received
    return 0;
  }
}

int DynamicCommunicationBuffer::sendCompleted ()
{
  int result = 1;
  
  std::list<CommunicationPacket*>::const_iterator it;
  
  for (it=packet_list.begin(); it != packet_list.end(); ++it) {
    result &= (*it)->testCompletion ();
  }
  
  return result;
}

int
DynamicCommunicationBuffer::giveFitSize(MPI_Datatype type, int availableSpace, int arrySize)
{
  int arrySpace, guessSize;
  MPI_Pack_size (arrySize, type, communicator, &arrySpace);
  if (availableSpace >= arrySpace) return arrySize;
  guessSize = (int) floor(((double)arrySize/(double)arrySpace)*availableSpace) + 1;
  do {
    guessSize--;
    MPI_Pack_size (guessSize, type, communicator, &arrySpace);
  } while ((availableSpace < arrySpace) && (guessSize > 0));

  return guessSize;
}
    
CommunicationPacket*
DynamicCommunicationBuffer::allocateNewPacket (int n)
{
  CommunicationPacket* result = packetPool.popPacket (communicator);
  result->init(communicator);
  result->setNumber (n);
  return result;
}

void
DynamicCommunicationBuffer::freePacket (CommunicationPacket* p)
{
  packetPool.pushPacket (p);
}

void
DynamicCommunicationBuffer:: clear ()
{
  std::list<CommunicationPacket*>::const_iterator it;
  for (it=packet_list.begin(); it != packet_list.end(); ++it) this->freePacket (*it);
  packet_list.clear();
  active_packet = NULL;
  number_of_packets = 0;
}


void
DynamicCommunicationBuffer::popNewRecvPacket ()
{
  active_packet = (*recvIt); ++recvIt;
  if (active_packet == NULL) OOFEM_ERROR ("DynamicCommunicationBuffer::popNewRecvPacket: no more packets received");
  active_packet -> init(communicator);

}

void 
DynamicCommunicationBuffer::pushNewRecvPacket (CommunicationPacket* p)
{
  packet_list.push_back (p);
}


int
DynamicCommunicationBuffer::bcast (int root) 
{
  OOFEM_ERROR ("DynamicCommunicationBuffer::bcast: not implemented");
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////

CommunicationPacket*
CommunicationPacketPool :: popPacket (MPI_Comm comm) 
{
  CommunicationPacket* result;
  
  if (available_packets.empty()) {
    // allocate new packet
    if ((result = new CommunicationPacket (comm, 0)) == NULL) 
      OOFEM_ERROR ("CommunicationPacketPool :: popPacket: allocation of new packed failed");
    allocatedPackets++; 
  } else {
    result = available_packets.front();
    available_packets.pop_front();
    freePackets--;
  }
#ifdef DEBUG
  // add packet into list of leased packets
  leased_packets.push_back(result);
#endif

  leasedPackets++;
  return result;
}

void 
CommunicationPacketPool::pushPacket (CommunicationPacket* p)
{
#ifdef DEBUG
  std::list<CommunicationPacket*>::iterator it = leased_packets.find(p);
  if (p!= leased_packets.end()) {
    // found previosly leased one
    leased_packets.erase(it);
    available_packets.push_bask(p);
  } else {
    OOFEM_ERROR ("CommunicationPacketPool::pushPacket: request to push strange packet (not allocated by pool)");
  }
#else
  available_packets.push_back(p);
#endif

  leasedPackets--;
  freePackets++;
}

void
CommunicationPacketPool::clear() 
{
  if (!leased_packets.empty()) {
    OOFEM_WARNING ("CommunicationPacketPool::clear: some packets still leased");
  }
  std::list<CommunicationPacket*>::iterator it ;
  for (it = available_packets.begin(); it != available_packets.end(); ++it) {
    if (*it) delete *it;
  }
  available_packets.clear();
  allocatedPackets=leasedPackets=freePackets=0;
}


void
CommunicationPacketPool::printInfo ()
{
  OOFEM_LOG_INFO("CommunicationPacketPool: allocated %d packets\n(packet size: %d, %d leased, %d free)\n", 
                 allocatedPackets,  __CommunicationPacket_DEFAULT_SIZE, leasedPackets,freePackets);

}
#endif
