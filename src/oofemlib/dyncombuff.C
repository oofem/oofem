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

#include <list>
#include <algorithm>

#include "dyncombuff.h"
#include "mathfem.h"
#include "error.h"

namespace oofem {
CommunicationPacket :: CommunicationPacket(MPI_Comm comm, int size, int num) : MPIBuffer(max(size, __CommunicationPacket_DEFAULT_SIZE), false)
{
    this->EOF_Flag = false;
    this->number = num;
    // reserve space for packet header
    curr_pos += givePackSize(comm, MPI_INT, 2);
}


CommunicationPacket :: CommunicationPacket(MPI_Comm comm, int num) : MPIBuffer(__CommunicationPacket_DEFAULT_SIZE, false)
{
    this->EOF_Flag = false;
    this->number = num;
    // reserve space for packet header
    curr_pos += givePackSize(comm, MPI_INT, 2);
}


CommunicationPacket :: ~CommunicationPacket()
{ }


void
CommunicationPacket :: init(MPI_Comm comm)
{
    MPIBuffer :: init();
    this->EOF_Flag = false;
    // reserve space for packet header
    curr_pos += givePackSize(comm, MPI_INT, 2);
}

int
CommunicationPacket :: iSend(MPI_Comm communicator, int dest, int tag)
{
    this->packHeader(communicator);
    return ( MPI_Isend(this->buff, this->curr_pos, MPI_PACKED, dest, tag,
                       communicator, & this->request) == MPI_SUCCESS );
}


int
CommunicationPacket :: iRecv(MPI_Comm communicator, int source, int tag, int count)
{
    if ( count ) {
        if ( count >= this->size ) {
            // reallocate itself
            if ( this->resize(count) == 0 ) {
                return 0;
            }
        }
    }

    return ( MPI_Irecv(this->buff, this->size, MPI_PACKED, source, tag,
                       communicator, & this->request) == MPI_SUCCESS );
}


int
CommunicationPacket :: testCompletion() {
    int flag;
    MPI_Status status;

    MPI_Test(& this->request, & flag, & status);
    return flag;
}

int
CommunicationPacket :: waitCompletion()
{
    MPI_Status status;

    return ( MPI_Wait(& this->request, & status) == MPI_SUCCESS );
}


int
CommunicationPacket :: packHeader(MPI_Comm comm)
{
    int _arry [ 2 ];
    int _res, _pos  = 0;

    _arry [ 0 ] = this->number;
    _arry [ 1 ] = this->EOF_Flag;

    _res = MPI_Pack(_arry, 2, MPI_INT, this->buff, size, & _pos, comm);

    return ( _res == MPI_SUCCESS );
}

int
CommunicationPacket :: unpackHeader(MPI_Comm comm)
{
    int _arry [ 2 ];
    int _res, _pos  = 0;

    _res = MPI_Unpack(this->buff, this->size, & _pos, _arry, 2, MPI_INT, comm);
    this->number =   _arry [ 0 ];
    this->EOF_Flag = _arry [ 1 ];

    return ( _res == MPI_SUCCESS );
}



/********DynamicCommunicationBuffer***************/

CommunicationPacketPool DynamicCommunicationBuffer :: packetPool;

DynamicCommunicationBuffer :: DynamicCommunicationBuffer(MPI_Comm comm, int size, bool dynamic) :
    CommunicationBuffer(comm, size, dynamic), packet_list()
{
    number_of_packets = 0;
    mode = DCB_null;
    completed = false;
    /*
     * // alocate first send/receive packet
     * active_packet = this->allocateNewPacket (++number_of_packets);
     * packet_list.push_back(active_packet);
     */
}


DynamicCommunicationBuffer :: DynamicCommunicationBuffer(MPI_Comm comm, bool dynamic) :
    CommunicationBuffer(comm, dynamic), packet_list()
{
    number_of_packets = 0;
    mode = DCB_null;
    completed = false;
    /*
     * // alocate first send/receive packet
     * active_packet = this->allocateNewPacket (++number_of_packets);
     * packet_list.push_back(active_packet);
     */
}

/// Destructor.
DynamicCommunicationBuffer :: ~DynamicCommunicationBuffer()
{
    this->clear();
}

void
DynamicCommunicationBuffer :: init()
{
    completed = false;
    this->clear();
}

void
DynamicCommunicationBuffer :: initForPacking()
{
    this->clear();

    if ( !active_packet ) {
        active_packet = this->allocateNewPacket(++number_of_packets);
        packet_list.push_back(active_packet);
    }
}

void
DynamicCommunicationBuffer :: initForUnpacking()
{
    recvIt = packet_list.begin();
    this->popNewRecvPacket();
}

/*
 * int
 * DynamicCommunicationBuffer::write (int* src, int n)
 * {
 * int _result=1;
 * int start_indx=0, end_indx, _size;
 *
 * do {
 *  _size = this->giveFitSize(MPI_INT, active_packet -> giveAvailableSpace(), n);
 *  end_indx = start_indx + _size;
 *
 *  if (_size) _result &= active_packet -> write (communicator, src+start_indx,_size, MPI_INT);
 *  if (end_indx == n) break;
 *  // active packet full, allocate a new one
 *  active_packet = this->allocateNewPacket (++number_of_packets);
 *  packet_list.push_back(active_packet);
 *  start_indx = end_indx;
 * } while (1);
 *
 * return _result;
 * }
 *
 *
 * int
 * DynamicCommunicationBuffer::read (int* dest, int n)
 * {
 * int _result=1;
 * int start_indx=0, end_indx, _size;
 *
 * do {
 *  _size = this->giveFitSize(MPI_INT, active_packet -> giveAvailableSpace(), n);
 *  end_indx = start_indx + _size;
 *
 *  if (_size) _result &= active_packet->read (communicator,dest+start_indx,_size, MPI_INT);
 *  if (end_indx == n) break;
 *  // active packet exhausted, pop a new one
 *  this->popNewRecvPacket();
 *  start_indx = end_indx;
 * } while (1);
 *
 * return _result;
 * }
 */

int
DynamicCommunicationBuffer :: iSend(int dest, int tag)
{
    int result = 1;

    /// set last (active) send packet as eof
    active_packet->setEOFFlag();

    active_rank = dest;
    active_tag = tag;
    for ( auto &packet: packet_list ) {
        result &= packet->iSend(communicator, dest, tag);
    }

    /*
     * int _myrank;
     * MPI_Comm_rank (communicator, &_myrank);
     * fprintf (stderr,"[%d] sending to [%d] %d packets for tag %d\n", _myrank, dest, number_of_packets, tag);
     * (*(packet_list.begin()))->dump();
     */
    mode = DCB_send;
    completed = false;
    return result;
}



int
DynamicCommunicationBuffer :: iRecv(int source, int tag, int count)
{
    this->init();
    number_of_packets = 0;
    // receive first packet, but it is probably not the last one
    // create new first packet and init its receive
    active_packet = this->allocateNewPacket(++number_of_packets);
    active_rank = source;
    active_tag = tag;
    mode = DCB_receive;
    completed = false;
    return active_packet->iRecv(communicator, source, tag);
}

int DynamicCommunicationBuffer :: receiveCompleted()
{
    /*
     * int _myrank;
     * MPI_Comm_rank (communicator, &_myrank);
     */
    if ( completed ) {
        return 1;
    }

    if ( active_packet->testCompletion() ) {
        // active packet received, add it to the pool and unpach header info
        active_packet->unpackHeader(communicator);

        //fprintf (stderr, "[%d] received from [%d] packet no. %d\n", _myrank, active_rank, active_packet->getNumber());

        pushNewRecvPacket(active_packet);
        if ( active_packet->hasEOFFlag() ) {
            // last packet received; init for unpacking
            this->initForUnpacking();

            //fprintf (stderr,"[%d] received from [%d] %d packets for tag %d\n", _myrank, active_rank, number_of_packets, active_tag);
            //active_packet->dump();
            completed = true;
            return 1;
        } else {
            // received next packet, but it is not the last one
            // create new packet and init its receive
            active_packet = this->allocateNewPacket(++number_of_packets);
            active_packet->iRecv(communicator, active_rank, active_tag);
            return 0;
        }
    } else {
        // active packet not yet received
        return 0;
    }
}

int DynamicCommunicationBuffer :: sendCompleted()
{
    int result = 1;

    if ( completed ) {
        return 1;
    }

    for ( auto &packet: packet_list ) {
        result &= packet->testCompletion();
    }

    completed = result;
    return result;
}

int
DynamicCommunicationBuffer :: testCompletion()
{
    if ( mode == DCB_send ) {
        return this->sendCompleted();
    } else if ( mode == DCB_receive ) {
        return this->receiveCompleted();
    } else {
        return 0;
    }
}

int
DynamicCommunicationBuffer :: waitCompletion()
{
    if ( mode == DCB_send ) {
        while ( !this->sendCompleted() ) { }

        ;
        return 1;
    } else if ( mode == DCB_receive ) {
        while ( !this->receiveCompleted() ) { }

        ;
        return 1;
    }

    return 0;
}


int
DynamicCommunicationBuffer :: giveFitSize(MPI_Datatype type, int availableSpace, int arrySize)
{
    int arrySpace, guessSize;
    MPI_Pack_size(arrySize, type, communicator, & arrySpace);
    if ( availableSpace >= arrySpace ) {
        return arrySize;
    }

    guessSize = ( int ) floor( ( ( double ) arrySize / ( double ) arrySpace ) * availableSpace ) + 1;
    do {
        guessSize--;
        MPI_Pack_size(guessSize, type, communicator, & arrySpace);
    } while ( ( availableSpace < arrySpace ) && ( guessSize > 0 ) );

    return guessSize;
}

CommunicationPacket *
DynamicCommunicationBuffer :: allocateNewPacket(int n)
{
    CommunicationPacket *result = packetPool.popPacket(communicator);
    result->init(communicator);
    result->setNumber(n);
    return result;
}

void
DynamicCommunicationBuffer :: freePacket(CommunicationPacket *p)
{
    packetPool.pushPacket(p);
}

void
DynamicCommunicationBuffer :: clear()
{
    for ( auto &packet: packet_list ) {
        this->freePacket(packet);
    }

    packet_list.clear();
    active_packet = NULL;
    number_of_packets = 0;
}


void
DynamicCommunicationBuffer :: popNewRecvPacket()
{
    active_packet = ( * recvIt );
    ++recvIt;
    if ( active_packet == NULL ) {
        OOFEM_ERROR("no more packets received");
    }

    //active_packet->init(communicator);
}

void
DynamicCommunicationBuffer :: pushNewRecvPacket(CommunicationPacket *p)
{
    packet_list.push_back(p);
}


int
DynamicCommunicationBuffer :: bcast(int root)
{
    OOFEM_ERROR("not implemented");
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////

CommunicationPacket *
CommunicationPacketPool :: popPacket(MPI_Comm comm)
{
    CommunicationPacket *result;

    if ( available_packets.empty() ) {
        // allocate new packet
        if ( ( result = new CommunicationPacket(comm, 0) ) == NULL ) {
            OOFEM_ERROR("allocation of new packed failed");
        }

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
CommunicationPacketPool :: pushPacket(CommunicationPacket *p)
{
#ifdef DEBUG
    std :: list< CommunicationPacket * > :: iterator it = std :: find(leased_packets.begin(), leased_packets.end(), p);
    if ( it != leased_packets.end() ) {
        // found previosly leased one
        leased_packets.erase(it);
        available_packets.push_back(p);
    } else {
        OOFEM_ERROR("request to push strange packet (not allocated by pool)");
    }

#else
    available_packets.push_back(p);
#endif

    leasedPackets--;
    freePackets++;
}

void
CommunicationPacketPool :: clear()
{
    if ( !leased_packets.empty() ) {
        OOFEM_WARNING("some packets still leased");
    }

    for ( auto &packet: available_packets ) {
        if ( packet ) {
            delete packet;
        }
    }

    available_packets.clear();
    allocatedPackets = leasedPackets = freePackets = 0;
}


void
CommunicationPacketPool :: printInfo()
{
    OOFEM_LOG_INFO("CommunicationPacketPool: allocated %d packets\n(packet size: %d, %d leased, %d free)\n",
                   allocatedPackets,  __CommunicationPacket_DEFAULT_SIZE, leasedPackets, freePackets);
}
} // end namespace oofem
