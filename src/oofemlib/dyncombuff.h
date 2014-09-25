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

#ifndef dyncombuff_h
#define dyncombuff_h

#include "parallel.h"
#include "combuff.h"
#include <list>

namespace oofem {
class IntArray;
class FloatArray;
class FloatMatrix;

 #define __CommunicationPacket_DEFAULT_SIZE 40960

/**
 * Class CommunicationPacket represent a data-packet, that is used to implement dynamic
 * communicator. Dynamic Communicator can pack messages into a dynamic message.
 * This dynamic message is split into a series
 * of data packets of fixed size (this is not necessary) that are send over network.
 *
 * A special header is put at the beginning of each packet buffer. This header keeps the message number
 * as well as the EOF flag indicating the last packet in message. This header is packed at the beginning
 * of each packet.
 */

class OOFEM_EXPORT CommunicationPacket : public MPIBuffer
{
protected:
    int number;
    bool EOF_Flag;

public:

 #ifdef __USE_MPI
    /// Constructor. Creates buffer of given size, using given communicator for packing.
    CommunicationPacket(MPI_Comm comm, int size, int num);
    /// Constructor. Creates empty buffer, using given communicator for packing.
    CommunicationPacket(MPI_Comm comm, int num);
 #endif
    /// Destructor.
    virtual ~CommunicationPacket();

    /**
     * Initializes buffer to empty state. All packed data are lost.
     */
    virtual void init(MPI_Comm comm);

    /**@name Services for buffer sending/receiving */
    //@{
 #ifdef __USE_MPI
    /**
     * Starts standard mode, nonblocking send.
     * @param dest Rank of destination.
     * @param tag Message tag.
     * @param communicator Communicator request handle.
     * @return Sends MIP_Succes if ok.
     */
    int iSend(MPI_Comm communicator, int dest, int tag);
    /**
     * Starts standard mode, nonblocking receive. The buffer must be large enough to receive all data.
     * @param source Rank of source.
     * @param tag Message tag.
     * @param count Number of elements to receive (bytes). Causes receive buffer to resize to count elements.
     * If zero (default value) buffer is not resized.
     * @param communicator Communicator request handle.
     * @return MIP_Succes if ok.
     */
    int iRecv(MPI_Comm communicator, int source, int tag, int count = 0);
    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @return True if operation complete, false otherwise.
     */
    virtual int testCompletion();
    /**
     * Waits until a completion of a nonblocking communication. The completion of a send operation indicates that the sender is
     * now free to update the locations in the send buffer, the completion of a receive operation indicates that the
     * receive buffer contains the received message, the receiver is now free to access it, and that the status object is set.
     * If the communication object associated with this request was created (nonblocking send or receive call),
     * then the object is deallocated by the call to MPI_WAIT and the request handle is set to MPI_REQUEST_NULL.
     * @return True if request is successful.
     */
    virtual int waitCompletion();
 #endif
    //@}

    void setNumber(int _num) { this->number = _num; }
    void setEOFFlag() { this->EOF_Flag = true; }
    int getNumber() { return number; }
    bool hasEOFFlag() { return EOF_Flag; }

    /// Packs packet header info at receiver beginning
    int packHeader(MPI_Comm);
    int unpackHeader(MPI_Comm);
};

class OOFEM_EXPORT CommunicationPacketPool
{
private:
    std :: list< CommunicationPacket * >available_packets;
    std :: list< CommunicationPacket * >leased_packets;

    int allocatedPackets, leasedPackets, freePackets;
public:
    CommunicationPacketPool() : available_packets(), leased_packets() {
        allocatedPackets = leasedPackets = freePackets;
    }
    ~CommunicationPacketPool() {
        this->clear();
    }

    CommunicationPacket *popPacket(MPI_Comm);
    void pushPacket(CommunicationPacket *);

    void printInfo();

private:
    void clear();
};


class OOFEM_EXPORT DynamicCommunicationBuffer : public CommunicationBuffer
{
protected:
    std :: list< CommunicationPacket * >packet_list;
    /// Iterator to iterate over received packets.
    std :: list< CommunicationPacket * > :: iterator recvIt;
    /// Active packet.
    CommunicationPacket *active_packet;
    /// Active rank and tag (send by initSend,initReceive, and initExchange).
    int active_tag, active_rank;
    int number_of_packets;

    /// Receiver mode.
    enum DCB_Mode { DCB_null, DCB_send, DCB_receive } mode;
    /// Static packet pool.
    static CommunicationPacketPool packetPool;
    /// Communication completion flag.
    bool completed;
public:
    /// Constructor. Creates buffer of given size, using given communicator for packing.
    DynamicCommunicationBuffer(MPI_Comm comm, int size, bool dynamic = 0);
    /// Constructor. Creates empty buffer, using given communicator for packing.
    DynamicCommunicationBuffer(MPI_Comm comm, bool dynamic = 0);
    /// Destructor.
    virtual ~DynamicCommunicationBuffer();

    virtual int resize(int newSize) { return 1; }
    virtual void init();

    /// Initialize for packing.
    virtual void initForPacking();
    /// Initialize for Unpacking (data already received).
    virtual void initForUnpacking();

    virtual int write(const int *src, int n)
    { return __write(src, n, MPI_INT); }
    virtual int write(const long *src, int n)
    { return __write(src, n, MPI_LONG); }
    virtual int write(const unsigned long *src, int n)
    { return __write(src, n, MPI_UNSIGNED_LONG); }
    virtual int write(const double *src, int n)
    { return __write(src, n, MPI_DOUBLE); }
    virtual int write(const char *src, int n)
    { return __write(src, n, MPI_CHAR); }

    virtual int read(int *dest, int n)
    { return __read(dest, n, MPI_INT); }
    virtual int read(long *dest, int n)
    { return __read(dest, n, MPI_LONG); }
    virtual int read(unsigned long *dest, int n)
    { return __read(dest, n, MPI_UNSIGNED_LONG); }
    virtual int read(double *dest, int n)
    { return __read(dest, n, MPI_DOUBLE); }
    virtual int read(char *dest, int n)
    { return __read(dest, n, MPI_CHAR); }


    virtual int iSend(int dest, int tag);
    virtual int iRecv(int source, int tag, int count = 0);
    virtual int bcast(int root);

    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @return True if operation complete, false otherwise.
     */
    int testCompletion();
    /**
     * Waits until a completion of a nonblocking communication. The completion of a send operation indicates that the sender is
     * now free to update the locations in the send buffer, the completion of a receive operation indicates that the
     * receive buffer contains the received message, the receiver is now free to access it, and that the status object is set.
     * If the communication object associated with this request was created (nonblocking send or receive call),
     * then the object is deallocated by the call to MPI_WAIT and the request handle is set to MPI_REQUEST_NULL.
     * @return True if request is successful.
     */
    virtual int waitCompletion();


    static void printInfo() { packetPool.printInfo(); }

protected:
    CommunicationPacket *allocateNewPacket(int);
    void freePacket(CommunicationPacket *);

    int receiveCompleted();
    int sendCompleted();

    void popNewRecvPacket();
    void pushNewRecvPacket(CommunicationPacket *);

    void clear();
    int giveFitSize(MPI_Datatype type, int availableSpace, int arrySize);

    /**
     * Templated low-level array packing method.
     * Templated version used since implementation is similar for different types
     * but type info is needed since implementation is relying on pointer arithmetic.
     */
    template< class T > int __write(T *src, int n, MPI_Datatype type) {
        int _result = 1;
        int start_indx = 0, end_indx, _size;
        int remaining_size = n;

        do {
            _size = this->giveFitSize(type, active_packet->giveAvailableSpace(), remaining_size);
            end_indx = start_indx + _size;

            if ( _size ) {
                _result &= active_packet->packArray(communicator, src + start_indx, _size, type);
            }

            if ( end_indx >= n ) {
                break;
            }

            // active packet full, allocate a new one
            active_packet = this->allocateNewPacket(++number_of_packets);
            packet_list.push_back(active_packet);
            start_indx = end_indx;
            remaining_size -= _size;
        } while ( 1 );

        return _result;
    }

    /**
     * Templated low-level array unpacking method.
     * Templated version used since implementation is similar for different types
     * but type info is needed since implementation is relying on pointer arithmetic.
     */
    template< class T > int __read(T *dest, int n, MPI_Datatype type) {
        int _result = 1;
        int start_indx = 0, end_indx, _size;
        int remaining_size = n;

        do {
            _size = this->giveFitSize(type, active_packet->giveAvailableSpace(), remaining_size);
            end_indx = start_indx + _size;

            if ( _size ) {
                _result &= active_packet->unpackArray(communicator, dest + start_indx, _size, type);
            }

            if ( end_indx >= n ) {
                break;
            }

            // active packet exhausted, pop a new one
            this->popNewRecvPacket();
            start_indx = end_indx;
            remaining_size -= _size;
        } while ( 1 );

        return _result;
    }
};
} // end namespace oofem

#endif // dyncombuff_h
