/* $Header: /home/cvs/bp/oofem/oofemlib/src/combuff.h,v 1.5 2003/04/06 14:08:23 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//
// class CommunicationBuffer
//

#ifndef combuff_h
#define combuff_h

#ifdef __PARALLEL_MODE

 #include "parallel.h"

namespace oofem {
class IntArray;
class FloatArray;
class FloatMatrix;


 #define __CommunicationBuffer_ALLOC_CHUNK 1024
/**
 * Type with size equal to one byte (sizeof (ComBuff_BYTE_TYPE) should be 1).
 * Communication buffer buffer member is of this type.
 */
typedef char ComBuff_BYTE_TYPE;

class MPIBuffer
{
protected:
    /// Size and current position in buffer in bytes (sizeof(char)).
    int size, curr_pos;
    /// dynamic flag (if true, buffer can grow, but reallocation is needed)
    bool isDynamic;
    /// Buffer. Dynamically allocated.
    ComBuff_BYTE_TYPE *buff;
    /**
     * MPI request handle. This value is used by some message parsing functions.
     * engngcommunicator also assembles array of commbuff handles and wait for some
     * completion (when  receiveing edata for example).
     */
    MPI_Request request;


public:
    /// Constructor. Creeates buffer of given size, using given communicator for packing
    MPIBuffer(int size, bool dynamic = 0);
    /// Constructor. Creeates empty buffer, using given communicator for packing
    MPIBuffer(bool dynamic = 0);
    /// Destructor.
    virtual ~MPIBuffer();

    /**
     * Resizes buffer to given size. If bufeer size is to be enlarged, then previously packed data
     * are kept in new buffer. Otherwise buffer is cleared using init service.
     * Current implementation only performs buffer growing, request for size decrease is ignored
     * to avoid realocation if further request for groving is encountered.
     * @param newSize new buffer size in bytes.
     * @return nonzero if succesfull.
     */
    int resize(int newSize);
    /**
     * Initializes buffer to empty state. All packed data are lost.
     */
    virtual void init();

    /// Returns the  current buffer size.
    int giveSize() { return size; }
    /// returns remaining space
    int giveAvailableSpace() { return ( size - curr_pos ); }
    /**
     * Returns associated MPI request handle
     */
    MPI_Request giveRequest() { return this->request; }

    /**
     * Packs array of a values of given type into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @param src adress of first value in memory
     * @param n number of packed integers
     * @param type dermines type of array values
     * @return nonzero if succesfull
     */
    int packArray(MPI_Comm communicator, const void *src, int n, MPI_Datatype type);
    /**
     * unpacks array of values of given type from buffer.
     * @param dest adress of first value in memory, where to store values
     * @param n number of unpacked integers
     * @param type dermines type of array values
     * @return nonzero if succesfull
     */
    int unpackArray(MPI_Comm communicator, void *dest, int n, MPI_Datatype type);


    /**@name Methods for determining pack size of datatype to pack/unpack to/from buffer */
    //@{
    /**
     * Returns pack size required to pack array of given type and size (c-style).
     * @param array size
     * @param type type id
     * @return  pack size required
     */
    int givePackSize(MPI_Comm communicator, MPI_Datatype type, int size);
    //@}

    /**@name Services for buffer sending/receiving */
    //@{
    /**
     * Starts standard mode, nonblocking send.
     * @param dest rank of destination
     * @param tag message tag
     * @param communicator (handle)
     * @return sends MIP_Succes if ok
     */
    virtual int iSend(MPI_Comm communicator, int dest, int tag);
    /**
     * Starts standard mode, nonblocking receive. The buffer must be large enough to receive all data.
     * @param source rank of source
     * @param tag message tag
     * @param count number of elements to receive (bytes). Causes receive buffer to resize to count elements.
     * If zero (default value) buffer is not resized.
     * @param reguest communicator request (handle)
     * @return MIP_Succes if ok
     */
    virtual int iRecv(MPI_Comm communicator, int source, int tag, int count = 0);
    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @return true if operation complete, false otherwise.
     */
    int testCompletion();
    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @param source contain the source tag
     * @param tag contain the tag of received message
     * @return true if operation complete, false otherwise.
     */
    int testCompletion(int &source, int &tag);
    /**
     * Waits until a completion of a nonblocking communication. The completion of a send operation indicates that the sender is
     * now free to update the locations in the send buffer, the completion of a receive operation indicates that the
     * receive buffer contains the received message, the receiver is now free to access it, and that the status object is set.
     * If the communication object associated with this request was created (nonblocking send or receive call),
     * then the object is deallocated by the call to MPI_WAIT and the request handle is set to MPI_REQUEST_NULL.
     *
     */
    virtual int waitCompletion();
    /**
     * Initalizes broadcast over colaborating processes.
     * The whole buffer size is broadcasted. All buffers participating in broadcast
     * should have the same size.
     * @param root rank of broadcast root
     * @return MIP_Succes if ok
     */
    int bcast(MPI_Comm communicator, int root);


    void dump();
private:
    /// Returns current buffer position.
    int givePosition() { return curr_pos; }
};


/**
 * Class CommunicationBuffer provides absraction for comunication buffer.
 * buffer is used as input or output buffer to various comunication
 * services provided by message parsing libraries.
 * It provides methods for buffer initialization and resizing, methods for packing and/or
 * unpacking data to/from buffer. Multiple messages can be packed/unpacked into/from buffer.
 * The services for packing/unpacking take care about multiple messages stored,
 * they maintain proper current buffer position for data inserting/retrieval.
 * Interface to low level message parsing function is provided, allowing to
 * send and receive buffer to selected destination.
 *
 */
class CommunicationBuffer
{
protected:
    MPI_Comm communicator;
public:
    CommunicationBuffer(MPI_Comm comm, int size, bool dynamic = 0) { communicator = comm; }
    /// Constructor. Creeates empty buffer, using given communicator for packing
    CommunicationBuffer(MPI_Comm comm, bool dynamic = 0) { communicator = comm; }
    /// Destructor.
    virtual ~CommunicationBuffer() { }

    /**
     * Resizes buffer to given size. If bufeer size is to be enlarged, then previously packed data
     * are kept in new buffer. Otherwise buffer is cleared using init service.
     * Current implementation only performs buffer growing, request for size decrease is ignored
     * to avoid realocation if further request for groving is encountered.
     * @param newSize new buffer size in bytes.
     * @return nonzero if succesfull.
     */
    virtual int resize(int newSize) = 0;
    /**
     * Initializes buffer to empty state. All packed data are lost.
     */
    virtual void init() = 0;

    /// Initialize for packing
    virtual void initForPacking() = 0;
    /// Initialize for Unpacking (data already received)
    virtual void initForUnpacking() = 0;

    /**@name Methods for datatype packing/unpacking to/from buffer */
    //@{
    /**
     *  Packs single integer value into buffer.
     *  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     *  @return nonzero if succesfull
     */
    int packInt(int value) { return packArray(& value, 1); }
    /**
     *  Packs single double value into buffer.
     *  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     *  @return nonzero if succesfull
     */
    int packDouble(double value) { return packArray(& value, 1); }
    /**
     * Packs given IntArray  value into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @return nonzero if succesfull
     */
    int packIntArray(const IntArray &arry);
    /**
     * Packs given FloatArray  value into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @return nonzero if succesfull
     */
    int packFloatArray(const FloatArray &arry);
    /**
     * Packs given FloatMatrix  value into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @return nonzero if succesfull
     */
    int packFloatMatrix(const FloatMatrix &mtrx);

    /**
     * Packs array of values of given type into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @param src adress of first value in memory
     * @param n number of packed integers
     * @return nonzero if succesfull
     */
    //@{
    virtual int packArray(const int *src, int n) = 0;
    virtual int packArray(const long *src, int n) = 0;
    virtual int packArray(const unsigned long *src, int n) = 0;
    virtual int packArray(const double *src, int n) = 0;
    virtual int packArray(const char *src, int n) = 0;
    //@}

    /**
     *  Unpacks single integer value from buffer.
     *  @return nonzero if succesfull
     */
    int unpackInt(int &value) { return unpackArray(& value, 1); }
    /**
     *  Unpacks single double value from buffer.
     *  @return nonzero if succesfull
     */
    int unpackDouble(double &value) { return unpackArray(& value, 1); }
    /**
     * Unpacks given IntArray  value from buffer.
     * @return nonzero if succesfull
     */
    virtual int unpackIntArray(IntArray &arry);
    /**
     * Unpacks given FloatArray  value from buffer.
     * @return nonzero if succesfull
     */
    int unpackFloatArray(FloatArray &arry);
    /**
     * Unpacks given FloatMatrix  value from buffer.
     * @return nonzero if succesfull
     */
    int unpackFloatMatrix(FloatMatrix &mtrx);
    /**
     * unpacks array of values of given type from buffer.
     * @param dest adress of first value in memory, where to store values
     * @param n number of unpacked integers
     * @return nonzero if succesfull
     */
    //@{
    virtual int unpackArray(int *dest, int n) = 0;
    virtual int unpackArray(long *dest, int n) = 0;
    virtual int unpackArray(unsigned long *dest, int n) = 0;
    virtual int unpackArray(double *dest, int n) = 0;
    virtual int unpackArray(char *dest, int n) = 0;
    //@}


    /**@name Methods for determining pack size of datatype to pack/unpack to/from buffer */
    //@{
    /**
     *  Returns pack size required to pack an array of given type and size(c-style).
     *  @param array size
     *  @param type type id
     *  @return  pack size required
     */
    virtual int givePackSize(MPI_Datatype type, int size) {
        int requredSpace;
        MPI_Pack_size(size, type, communicator, & requredSpace);
        return requredSpace;
    }
    //@}

    /**@name Services for buffer sending/receiving */
    //@{
    /**
     * Starts standard mode, nonblocking send.
     * @param dest rank of destination
     * @param tag message tag
     * @param communicator (handle)
     * @return sends MIP_Succes if ok
     */
    virtual int iSend(int dest, int tag) = 0;
    /**
     * Starts standard mode, nonblocking receive. The buffer must be large enough to receive all data.
     * @param source rank of source
     * @param tag message tag
     * @param count number of elements to receive (bytes). Causes receive buffer to resize to count elements.
     * If zero (default value) buffer is not resized.
     * @param reguest communicator request (handle)
     * @return MIP_Succes if ok
     */
    virtual int iRecv(int source, int tag, int count = 0) = 0;
    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @return true if operation complete, false otherwise.
     */
    virtual int testCompletion() = 0;
    /**
     * Waits until a completion of a nonblocking communication. The completion of a send operation indicates that the sender is
     * now free to update the locations in the send buffer, the completion of a receive operation indicates that the
     * receive buffer contains the received message, the receiver is now free to access it, and that the status object is set.
     * If the communication object associated with this request was created (nonblocking send or receive call),
     * then the object is deallocated by the call to MPI_WAIT and the request handle is set to MPI_REQUEST_NULL.
     *
     */
    virtual int waitCompletion() = 0;
    /**
     * Initalizes broadcast over colaborating processes.
     * The whole buffer size is broadcasted. All buffers participating in broadcast
     * should have the same size.
     * @param root rank of broadcast root
     * @return MIP_Succes if ok
     */
    virtual int bcast(int root) = 0;
    //@}
};





class StaticCommunicationBuffer : public CommunicationBuffer, public MPIBuffer
{
public:
    StaticCommunicationBuffer(MPI_Comm comm, int size, bool dynamic = 0) : CommunicationBuffer(comm, size, dynamic),
        MPIBuffer(size, dynamic) { }
    /// Constructor. Creeates empty buffer, using given communicator for packing
    StaticCommunicationBuffer(MPI_Comm comm, bool dynamic = 0) : CommunicationBuffer(comm, dynamic), MPIBuffer(dynamic) { }
    /// Destructor.
    virtual ~StaticCommunicationBuffer() { }

    /**
     * Resizes buffer to given size. If bufeer size is to be enlarged, then previously packed data
     * are kept in new buffer. Otherwise buffer is cleared using init service.
     * Current implementation only performs buffer growing, request for size decrease is ignored
     * to avoid realocation if further request for groving is encountered.
     * @param newSize new buffer size in bytes.
     * @return nonzero if succesfull.
     */
    virtual int resize(int newSize) { return MPIBuffer :: resize(newSize); }
    /**
     * Initializes buffer to empty state. All packed data are lost.
     */
    virtual void init() { return MPIBuffer :: init(); }
    /// Initialize for packing
    virtual void initForPacking() { this->init(); }
    /// Initialize for Unpacking (data already received)
    virtual void initForUnpacking() { this->init(); }


    /**@name Methods for datatype packing/unpacking to/from buffer
     * Packs array of values if given type into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @param src adress of first value in memory
     * @param n number of packed integers
     * @return nonzero if succesfull
     */
    //@{
    virtual int packArray(const int *src, int n)
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_INT); }
    virtual int packArray(const long *src, int n)
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_LONG); }
    virtual int packArray(const unsigned long *src, int n)
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_UNSIGNED_LONG); }
    virtual int packArray(const double *src, int n)
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_DOUBLE); }
    virtual int packArray(const char *src, int n)
    { return MPIBuffer :: packArray(this->communicator, src, n, MPI_CHAR); }

    virtual int unpackArray(int *dest, int n)
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_INT); }
    virtual int unpackArray(long *dest, int n)
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_LONG); }
    virtual int unpackArray(unsigned long *dest, int n)
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_UNSIGNED_LONG); }
    virtual int unpackArray(double *dest, int n)
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_DOUBLE); }
    virtual int unpackArray(char *dest, int n)
    { return MPIBuffer :: unpackArray(this->communicator, dest, n, MPI_CHAR); }
    //@}



    /**@name Methods for determining pack size of datatype to pack/unpack to/from buffer */
    //@{
    /**
     *  Returns pack size required to pack integer array (c-style).
     *  @param array size
     *  @return  pack size required
     */
    virtual int givePackSize(MPI_Datatype type, int size) { return MPIBuffer :: givePackSize(this->communicator, type, size); }
    //@}

    /**@name Services for buffer sending/receiving */
    //@{
    /**
     * Starts standard mode, nonblocking send.
     * @param dest rank of destination
     * @param tag message tag
     * @param communicator (handle)
     * @return sends MIP_Succes if ok
     */
    virtual int iSend(int dest, int tag) { return MPIBuffer :: iSend(this->communicator, dest, tag); }
    /**
     * Starts standard mode, nonblocking receive. The buffer must be large enough to receive all data.
     * @param source rank of source
     * @param tag message tag
     * @param count number of elements to receive (bytes). Causes receive buffer to resize to count elements.
     * If zero (default value) buffer is not resized.
     * @param reguest communicator request (handle)
     * @return MIP_Succes if ok
     */
    virtual int iRecv(int source, int tag, int count = 0)  { return MPIBuffer :: iRecv(this->communicator, source, tag, count); }
    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @return true if operation complete, false otherwise.
     */
    virtual int testCompletion() { return MPIBuffer :: testCompletion(); }
    /**
     * Tests if the operation identified by this->request is complete.
     * In such case, true is returned and
     * if communication was initiated by nonblocking send/receive, then request handle
     * is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
     * @param source contain the source tag
     * @param tag contain the tag of received message
     * @return true if operation complete, false otherwise.
     */
    int testCompletion(int &source, int &tag) { return MPIBuffer :: testCompletion(source, tag); }

    /**
     * Waits until a completion of a nonblocking communication. The completion of a send operation indicates that the sender is
     * now free to update the locations in the send buffer, the completion of a receive operation indicates that the
     * receive buffer contains the received message, the receiver is now free to access it, and that the status object is set.
     * If the communication object associated with this request was created (nonblocking send or receive call),
     * then the object is deallocated by the call to MPI_WAIT and the request handle is set to MPI_REQUEST_NULL.
     *
     */
    virtual int waitCompletion() { return MPIBuffer :: waitCompletion(); };
    /**
     * Initalizes broadcast over colaborating processes.
     * The whole buffer size is broadcasted. All buffers participating in broadcast
     * should have the same size.
     * @param root rank of broadcast root
     * @return MIP_Succes if ok
     */
    virtual int bcast(int root) { return MPIBuffer :: bcast(this->communicator, root); }
    //@}
};
} // end namespace oofem
#endif
#endif // combuff_h


