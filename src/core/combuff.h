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

#ifndef combuff_h
#define combuff_h

#include "oofemenv.h"
#include "datastream.h"
#include "parallel.h"

namespace oofem {

#include <stdint.h>
#include <limits.h>

#if SIZE_MAX == UCHAR_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif



 #define __CommunicationBuffer_ALLOC_CHUNK 1024
/**
 * Type with size equal to one byte (sizeof (ComBuff_BYTE_TYPE) should be 1).
 * Communication buffer buffer member is of this type.
 */
typedef char ComBuff_BYTE_TYPE;

class OOFEM_EXPORT MPIBuffer
{
protected:
    /// Size and current position in buffer in bytes (sizeof(char)).
    std::size_t size, curr_pos;
    /// Dynamic flag (if true, buffer can grow, but reallocation is needed).
    bool isDynamic;
    /// Buffer. Dynamically allocated.
    ComBuff_BYTE_TYPE *buff;
    /**
     * MPI request handle. This value is used by some message parsing functions.
     * EngngCommunicator also assembles array of commbuff handles and wait for some
     * completion (when receiving data for example).
     */
    MPI_Request request;

public:
    /// Constructor. Creates buffer of given size, using given communicator for packing.
    MPIBuffer(std::size_t size, bool dynamic = 0);
    /// Constructor. Creates empty buffer, using given communicator for packing.
    MPIBuffer(bool dynamic = 0);
    /// Destructor.
    virtual ~MPIBuffer();

    /**
     * Resizes buffer to given size. If buffer size is to be enlarged, then previously packed data
     * are kept in new buffer. Otherwise buffer is cleared using init service.
     * Current implementation only performs buffer growing, request for size decrease is ignored
     * to avoid reallocation if further request for growing is encountered.
     * @param newSize new buffer size in bytes.
     * @return nonzero if successful.
     */
    int resize(std::size_t newSize);
    /**
     * Initializes buffer to empty state. All packed data are lost.
     */
    virtual void init();

    /// @return Current buffer size.
    std::size_t giveSize() { return size; }
    /// @return Remaining space.
    std::size_t giveAvailableSpace() { return ( size - curr_pos ); }
    /**
     * Returns associated MPI request handle
     */
    MPI_Request giveRequest() { return this->request; }

    /**
     * Packs array of a values of given type into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and deallocation.
     * @param communicator Communicator handle.
     * @param src Address of first value in memory.
     * @param n Number of packed integers.
     * @param type Determines type of array values.
     * @return Nonzero if successful.
     */
    int packArray(MPI_Comm communicator, const void *src, std::size_t n, MPI_Datatype type);
    /**
     * Unpacks array of values of given type from buffer.
     * @param communicator Communicator handle.
     * @param dest Address of first value in memory, where to store values
     * @param n Number of unpacked integers.
     * @param type Determines type of array values.
     * @return Nonzero if successful.
     */
    int unpackArray(MPI_Comm communicator, void *dest, std::size_t n, MPI_Datatype type);


    /**@name Methods for determining pack size of datatype to pack/unpack to/from buffer */
    //@{
    /**
     * Returns pack size required to pack array of given type and size (c-style).
     * @param communicator Communicator handle.
     * @param type Type id.
     * @param size Size of array to pack.
     * @return Pack size required.
     */
    int givePackSize(MPI_Comm communicator, MPI_Datatype type, std::size_t size);
    //@}

    /**@name Services for buffer sending/receiving */
    //@{
    /**
     * Starts standard mode, nonblocking send.
     * @param communicator Communicator handle.
     * @param dest Rank of destination.
     * @param tag Message tag.
     * @return MPI_SUCCESS if ok.
     */
    virtual int iSend(MPI_Comm communicator, int dest, int tag);
    /**
     * Starts standard mode, nonblocking receive. The buffer must be large enough to receive all data.
     * @param source Rank of source.
     * @param tag Message tag.
     * @param count Number of elements to receive (bytes). Causes receive buffer to resize to count elements.
     * If zero (default value) buffer is not resized.
     * @param communicator Request communicator (handle).
     * @return MPI_SUCCESS if ok.
     */
    virtual int iRecv(MPI_Comm communicator, int source, int tag, std::size_t count = 0);
    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @return Nonzero if operation complete, zero otherwise.
     */
    int testCompletion();
    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @param source Source tag.
     * @param tag Tag of received message.
     * @return Nonzero if operation complete, zero otherwise.
     */
    int testCompletion(int &source, int &tag);
    /**
     * Waits until a completion of a nonblocking communication. The completion of a send operation indicates that the sender is
     * now free to update the locations in the send buffer, the completion of a receive operation indicates that the
     * receive buffer contains the received message, the receiver is now free to access it, and that the status object is set.
     * If the communication object associated with this request was created (nonblocking send or receive call),
     * then the object is deallocated by the call to MPI_WAIT and the request handle is set to MPI_REQUEST_NULL.
     * @return
     */
    virtual int waitCompletion();
    /**
     * Initializes broadcast over collaborating processes.
     * The whole buffer size is broadcasted. All buffers participating in broadcast
     * should have the same size.
     * @param communicator Communicator (handle).
     * @param root Rank of broadcast root.
     * @return MPI_SUCCESS if ok.
     */
    int bcast(MPI_Comm communicator, int root);

    void dump();

private:
    /// @return Current buffer position.
    std::size_t givePosition() { return curr_pos; }
};


/**
 * Class CommunicationBuffer provides abstraction for communication buffer.
 * buffer is used as input or output buffer to various communication
 * services provided by message parsing libraries.
 * It provides methods for buffer initialization and resizing, methods for packing and/or
 * unpacking data to/from buffer. Multiple messages can be packed/unpacked into/from buffer.
 * The services for packing/unpacking take care about multiple messages stored,
 * they maintain proper current buffer position for data inserting/retrieval.
 * Interface to low level message parsing function is provided, allowing to
 * send and receive buffer to selected destination.
 *
 */
class OOFEM_EXPORT CommunicationBuffer: public DataStream
{
protected:
    MPI_Comm communicator;
public:
    CommunicationBuffer(MPI_Comm comm, std::size_t size, bool dynamic = 0) : communicator(comm) { }
    /// Constructor. Creates empty buffer, using given communicator for packing
    CommunicationBuffer(MPI_Comm comm, bool dynamic = 0) : communicator(comm) { }
    /// Destructor.
    virtual ~CommunicationBuffer() { }

    /**
     * Resizes buffer to given size. If buffer size is to be enlarged, then previously packed data
     * are kept in new buffer. Otherwise buffer is cleared using init service.
     * Current implementation only performs buffer growing, request for size decrease is ignored
     * to avoid reallocation if further request for growing is encountered.
     * @param newSize New buffer size in bytes.
     * @return Nonzero if successful.
     */
    virtual int resize(std::size_t newSize) = 0;
    /**
     * Initializes buffer to empty state. All packed data are lost.
     */
    virtual void init() = 0;

    /// Initialize for packing
    virtual void initForPacking() = 0;
    /// Initialize for Unpacking (data already received)
    virtual void initForUnpacking() = 0;

    using DataStream::read;
    int read(bool &data) override;

    using DataStream::write;
    int write(bool data) override;

    int givePackSizeOfInt(std::size_t count) override;
    int givePackSizeOfDouble(std::size_t count) override;
    int givePackSizeOfChar(std::size_t count) override;
    int givePackSizeOfBool(std::size_t count) override;
    int givePackSizeOfLong(std::size_t count) override;
    int givePackSizeOfSizet(std::size_t count) override;


    /**@name Services for buffer sending/receiving */
    //@{
    /**
     * Starts standard mode, nonblocking send.
     * @param dest Rank of destination.
     * @param tag Message tag.
     * @return MPI_SUCCESS if ok.
     */
    virtual int iSend(int dest, int tag) = 0;
    /**
     * Starts standard mode, nonblocking receive. The buffer must be large enough to receive all data.
     * @param source Rank of source.
     * @param tag Message tag.
     * @param count Number of elements to receive (bytes). Causes receive buffer to resize to count elements.
     * If zero (default value) buffer is not resized.
     * @return MPI_SUCCESS if ok.
     */
    virtual int iRecv(int source, int tag, std::size_t count = 0) = 0;
    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @return Nonzero if operation complete, zero otherwise.
     */
    virtual int testCompletion() = 0;
    /**
     * Waits until a completion of a nonblocking communication. The completion of a send operation indicates that the sender is
     * now free to update the locations in the send buffer, the completion of a receive operation indicates that the
     * receive buffer contains the received message, the receiver is now free to access it, and that the status object is set.
     * If the communication object associated with this request was created (nonblocking send or receive call),
     * then the object is deallocated by the call to MPI_WAIT and the request handle is set to MPI_REQUEST_NULL.
     * @return Nonzero if ...
     */
    virtual int waitCompletion() = 0;
    /**
     * Initializes broadcast over collaborating processes.
     * The whole buffer size is broadcasted. All buffers participating in broadcast
     * should have the same size.
     * @param root Rank of broadcast root.
     * @return MPI_SUCCESS if ok.
     */
    virtual int bcast(int root) = 0;
    //@}
};


class OOFEM_EXPORT StaticCommunicationBuffer : public CommunicationBuffer, public MPIBuffer
{
public:
    StaticCommunicationBuffer(MPI_Comm comm, int size, bool dynamic = 0) : CommunicationBuffer(comm, size, dynamic),
        MPIBuffer(size, dynamic) { }
    /// Constructor. Creates empty buffer, using given communicator for packing
    StaticCommunicationBuffer(MPI_Comm comm, bool dynamic = 0) : CommunicationBuffer(comm, dynamic), MPIBuffer(dynamic) { }
    /// Destructor.
    virtual ~StaticCommunicationBuffer() { }

    int resize(std::size_t newSize) override { return MPIBuffer :: resize(newSize); }

    void init() override { return MPIBuffer :: init(); }
    void initForPacking() override { this->init(); }
    void initForUnpacking() override { this->init(); }

    using CommunicationBuffer::write;
    int write(const int *src, std::size_t n) override
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_INT); }
    int write(const long *src, std::size_t n) override
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_LONG); }
    int write(const unsigned long *src, std::size_t n) override
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_UNSIGNED_LONG); }
    int write(const double *src, std::size_t n) override
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_DOUBLE); }
    int write(const char *src, std::size_t n) override
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_CHAR); }

    using CommunicationBuffer::read;
    int read(int *dest, std::size_t n) override
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_INT); }
    int read(long *dest, std::size_t n) override
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_LONG); }
    int read(unsigned long *dest, std::size_t n) override
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_UNSIGNED_LONG); }
    int read(double *dest, std::size_t n) override
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_DOUBLE); }
    int read(char *dest, std::size_t n) override
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_CHAR); }

    #ifdef _MSC_VER
    int read(std::size_t* data, std::size_t count) override { return MPIBuffer :: unpackArray(this->communicator, data, count, my_MPI_SIZE_T); }
    int write(const std::size_t* data, std::size_t count) override { return MPIBuffer :: packArray(this->communicator, data, count, my_MPI_SIZE_T); }
    #endif

    int iSend(int dest, int tag) override { return MPIBuffer :: iSend(this->communicator, dest, tag); }

    int iRecv(int source, int tag, std::size_t count = 0) override { return MPIBuffer :: iRecv(this->communicator, source, tag, count); }

    int testCompletion() override { return MPIBuffer :: testCompletion(); }

    int testCompletion(int &source, int &tag) { return MPIBuffer :: testCompletion(source, tag); }

    int waitCompletion() override { return MPIBuffer :: waitCompletion(); }

    int bcast(int root) override { return MPIBuffer :: bcast(this->communicator, root); }
};
} // end namespace oofem

#endif // combuff_h
