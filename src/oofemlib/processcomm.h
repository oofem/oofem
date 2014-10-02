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

#ifndef processcomm_h
#define processcomm_h

#include "oofemcfg.h"
#include "combuff.h"
#include "commbufftype.h"
#include "communicatormode.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "error.h"

#include <string>

namespace oofem {
/**
 * The ProcessCommunicator and corresponding buffers (represented by this class)
 * are separated in order to allow share the same buffer by several communicators.
 * Here sharing means reusing for different but NON-OVERLAPPING communications.
 * if communications overlap, the different instances of ProcessCommunicatorBuff should be used!
 * The ProcessCommunicatorBuff objects are registered in corresponding communicator,
 * then if maps are available, comBuff should be resized and used in subsequent operations.
 *
 * The registration is necessary, otherwise before each send operation the buffers (given probably as parameter)
 * will be resized again (size have to be computed again) and this is probably quite cost operation.
 * When comBuff will be registered, resize is needed only when maps change, and this will not occur frequently
 * (its even quite rare).
 */
class OOFEM_EXPORT ProcessCommunicatorBuff
{
protected:
    /// Send buffer.
    CommunicationBuffer *send_buff;
    /// Receive buffer.
    CommunicationBuffer *recv_buff;
public:
    /// Constructor, creates empty send and receive com buffs in MPI_COMM_WORLD.
    ProcessCommunicatorBuff(CommBuffType t);
    ~ProcessCommunicatorBuff() {
        delete send_buff;
        delete recv_buff;
    }

    /**@name Methods for datatype packing/unpacking to/from buffer */
    //@{
    /**
     *  Packs single integer value into buffer.
     *  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and deallocation.
     *  @return Nonzero if successful.
     */
    int packInt(int value) { return write(& value, 1); }
    /**
     *  Packs single double value into buffer.
     *  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and deallocation.
     *  @return Nonzero if successful.
     */
    int packDouble(double value)  { return write(& value, 1); }
    /**
     * Packs array of values of given type into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and deallocation.
     * @param src Address of first value in memory.
     * @param n Number of packed integers.
     * @return Nonzero if successful.
     */
    //@{
    int write(const int *src, int n) { return send_buff->write(src, n); }
    int write(const long *src, int n) { return send_buff->write(src, n); }
    int write(const unsigned long *src, int n) { return send_buff->write(src, n); }
    int write(const double *src, int n) { return send_buff->write(src, n); }
    int write(const char *src, int n) { return send_buff->write(src, n); }
    //@}

    /**
     * Packs given IntArray  value into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and deallocation.
     * @return Nonzero if successful.
     */
    int packIntArray(const IntArray &arry) { return arry.storeYourself(send_buff); }
    /**
     * Packs given FloatArray  value into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and deallocation.
     * @return Nonzero if successful.
     */
    int packFloatArray(const FloatArray &arry) { return arry.storeYourself(send_buff); }
    /**
     * Packs given FloatMatrix  value into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and deallocation.
     * @return Nonzero if successful.
     */
    int packFloatMatrix(const FloatMatrix &mtrx) { return mtrx.storeYourself(send_buff); }
    /**
     * Packs given string into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and deallocation.
     * @return Nonzero if successful.
     */
    int write(const std :: string &str) { return this->send_buff->write(str); }
    /**
     *  Unpacks single integer value from buffer.
     *  @return Nonzero if successful.
     */
    int read(int &value)  { return read(& value, 1); }
    /**
     *  Unpacks single double value from buffer.
     *  @return Nonzero if successful.
     */
    int read(double &value) { return read(& value, 1); }
    /**
     * unpacks array of value of given type from buffer.
     * @param dest Address of first value in memory, where to store values.
     * @param n Number of unpacked integers.
     * @return Nonzero if successful.
     */
    //@{
    int read(int *dest, int n) { return static_cast<DataStream*>(recv_buff)->read(dest, n); }
    int read(long *dest, int n) { return static_cast<DataStream*>(recv_buff)->read(dest, n); }
    int read(unsigned long *dest, int n) { return static_cast<DataStream*>(recv_buff)->read(dest, n); }
    int read(double *dest, int n) { return static_cast<DataStream*>(recv_buff)->read(dest, n); }
    int read(char *dest, int n) { return static_cast<DataStream*>(recv_buff)->read(dest, n); }
    //@}
    /**
     * Unpacks given IntArray  value from buffer.
     * @return Nonzero if successful.
     */
    int unpackIntArray(IntArray &arry) { return arry.restoreYourself(recv_buff); }
    /**
     * Unpacks given FloatArray  value from buffer.
     * @return Nonzero if successful.
     */
    int unpackFloatArray(FloatArray &arry) { return arry.restoreYourself(recv_buff); }
    /**
     * Unpacks given FloatMatrix  value from buffer.
     * @return Nonzero if successful.
     */
    int unpackFloatMatrix(FloatMatrix &mtrx) { return mtrx.restoreYourself(recv_buff); }
    /**
     * Unpacks given string value from buffer.
     * @return Nonzero if successful.
     */
    int read(std :: string &str) { return recv_buff->read(str); }
    //@}

    /// Initializes send buffer to empty state. All packed data are lost.
    void initSendBuff() { send_buff->init(); }
    /// Initializes send buffer to empty state. All packed data are lost.
    void initRecvBuff() { recv_buff->init(); }
    /// Initializes receiver buffers.
    void init() {
        initSendBuff();
        initRecvBuff();
    }

    /// Initialize for packing.
    void initForPacking() { send_buff->initForPacking(); }
    /// Initialize for Unpacking (data already received).
    void initForUnpacking() { recv_buff->initForUnpacking(); }

    /**
     * Initializes data exchange with associated problem.
     * if send or receive pool is empty, the send or receive communication is not performed.
     * @param rank Partition number.
     * @param tag Message tag.
     * @return Nonzero if success.
     */
    int initExchange(int rank, int tag) {
        int result = 1;
        result &= initSend(rank, tag);
        result &= initReceive(rank, tag);

        return result;
    }
    /**
     * Initialize the send data exchange with associate problem.
     * if send pool is empty, the send communication is not performed.
     * @param rank Partition number.
     * @param tag Message tag.
     * @return Nonzero if success.
     */
    int initSend(int rank, int tag) { return send_buff->iSend(rank, tag); }
    /**
     * Initialize the receive data exchange with associate problem.
     * if receive pool is empty, the receive communication is not performed.
     * @param rank Partition number.
     * @param tag Message tag.
     * @return Nonzero if success.
     */
    int initReceive(int rank, int tag) { return recv_buff->iRecv(rank, tag); }

    /**
     *  Clears all buffer contents.
     */
    void clearBuffers() { }

    void resizeSendBuffer(int size) { send_buff->resize(size); }
    void resizeReceiveBuffer(int size) { recv_buff->resize(size); }

    /**@name Methods for manipulating/testing receiver state */
    //@{
    int sendCompleted() { return send_buff->testCompletion(); }
    int receiveCompleted() { return recv_buff->testCompletion(); }
    int testCompletion() { return ( send_buff->testCompletion() && recv_buff->testCompletion() ); }
    int waitCompletion() { return ( send_buff->waitCompletion() && recv_buff->waitCompletion() ); }
    //@}

    /**
     * Returns send buffer of receiver.
     */
    CommunicationBuffer *giveSendBuff() { return send_buff; }
    /**
     * Returns receive buffer of receiver.
     */
    CommunicationBuffer *giveRecvBuff() { return recv_buff; }
};


/**
 * Class representing process communicator for engineering model.
 * Process communicator provides all services for communication with
 * associated remote process (problem or task).
 */
class OOFEM_EXPORT ProcessCommunicator
{
protected:
    /// Associated partition (problem) number (rank)
    int rank; // remote problem id = rank

    /// Communicator buffers representation.
    ProcessCommunicatorBuff *pcBuffer;

    /// Nodes to send.
    IntArray toSend;
    /// Nodes to receive.
    IntArray toReceive;
    /// Mode.
    CommunicatorMode mode;

public:
    /**
     * Constructor. Creates new problem (partition) communicator associated
     * to partition with number (rank) irank.
     * @param b ProcessCommunicatorBuff to use.
     * @param irank Rank of associated partition.
     * @param m Mode of communicator.
     */
    ProcessCommunicator(ProcessCommunicatorBuff * b, int irank, CommunicatorMode m = CommMode_Static);
    /// Destructor
    ~ProcessCommunicator() { }

    /**
     * Returns corresponding rank of associated partition
     */
    int giveRank() { return rank; }

    ///  Returns communication buffer.
    ProcessCommunicatorBuff *giveProcessCommunicatorBuff() {
        if ( pcBuffer ) {
            return pcBuffer;
        }

        OOFEM_ERROR("ProcessCommunicatorBuff undefined");
        return NULL;
    }


    /**
     * Returns receiver to send map.
     */
    const IntArray *giveToSendMap() { return & toSend; }
    /**
     * Returns receiver to receive map.
     */
    const IntArray *giveToRecvMap() { return & toReceive; }


    /**
     * Sets receiver toSend array to src. The method assumes that toSend array is sorted according
     * to global number associated to corresponding components given in src array. THis is necessary to
     * ensure proper pack/unpack order on local and remote problem. Send buff is resized to  hold all
     * necessary data.
     * @param emodel Engng model.
     * @param src Source of toSend array.
     * @param packUnpackType Determines the type of packed quantity, used by emodel estimateMaxPackSize
     * service to estimate the size of pack/unpack buffer accordingly.
     */
    template< class T > void setToSendArry(T *emodel, const IntArray &src, int packUnpackType);
    /**
     * Sets receiver toRecv array to src. The method assumes that toRecv array is sorted according
     * to global number associated to corresponding components given in src array. THis is necessary to
     * ensure proper pack/unpack order on local and remote problem. Recv buff is resized to
     * hold all necessary data.
     * @param emodel Engng model.
     * @param src Source of toRecv array.
     * @param packUnpackType Determines the type of packed quantity, used by emodel estimateMaxPackSize
     * service to estimate the size of pack/unpack buffer accordingly.
     */
    template< class T > void setToRecvArry(T *emodel, const IntArray &src, int packUnpackType);
    /**
     * Pack nodal data to send buff.
     * @param emodel Engineering model to pack.
     * @param packFunc Function used to pack nodal data in to buffer. It uses toSend array
     * to loop over required nodes.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T > int packData( T *emodel, int ( T :: *packFunc )( ProcessCommunicator & ) )
    {
        if ( !toSend.isEmpty() || ( this->mode == CommMode_Dynamic ) ) {
            giveProcessCommunicatorBuff()->initForPacking();
            return ( emodel->*packFunc )(* this);
        } else {
            return 1;
        }
    }
    /**
     * Pack nodal data to send buff.
     * @param emodel Engineering model to pack.
     * @param src Source to pack from.
     * @param packFunc Function used to pack nodal data in to buffer. It uses toSend array
     * to loop over required nodes.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T, class P > int packData( T *emodel, P *src, int ( T :: *packFunc )( P *, ProcessCommunicator & ) )
    {
        if ( !toSend.isEmpty() || ( this->mode == CommMode_Dynamic ) ) {
            giveProcessCommunicatorBuff()->initForPacking();
            return ( emodel->*packFunc )(src, * this);
        } else {
            return 1;
        }
    }
    /**
     * Unpack nodal data from recv buff.
     * @param emodel Engineering model to unpack from.
     * @param unpackFunc Function used to unpack nodal data from buffer. It uses toRecv array
     * to loop over required nodes.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T > int unpackData( T *emodel,  int ( T :: *unpackFunc )( ProcessCommunicator & ) )
    {
        if ( !toReceive.isEmpty() || ( this->mode == CommMode_Dynamic ) ) {
            giveProcessCommunicatorBuff()->initForUnpacking();
            return ( emodel->*unpackFunc )(* this);
        } else {
            return 1;
        }
    }
    /**
     * Unpack nodal data from recv buff.
     * @param emodel Engineering model to unpack from.
     * @param dest Destination.
     * @param unpackFunc Function used to unpack nodal data from buffer. It uses toRecv array
     * to loop over required nodes.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T, class P > int unpackData( T *emodel,  P *dest, int ( T :: *unpackFunc )( P *, ProcessCommunicator & ) )
    {
        if ( !toReceive.isEmpty() || ( this->mode == CommMode_Dynamic ) ) {
            giveProcessCommunicatorBuff()->initForUnpacking();
            return ( emodel->*unpackFunc )(dest, * this);
        } else {
            return 1;
        }
    }
    /**
     * Initializes data exchange with associated problem.
     * if send or receive pool is empty, the send or receive communication is not performed.
     * @param tag Message tag.
     * @return Nonzero if success.
     */
    int initExchange(int tag);
    /**
     * Initialize the send data exchange with associate problem.
     * if send pool is empty, the send communication is not performed.
     * @param tag Message tag.
     * @return Nonzero if success.
     */
    int initSend(int tag);
    /**
     * Initialize the receive data exchange with associate problem.
     * if receive pool is empty, the receive communication is not performed.
     * @param tag Message tag.
     * @return Nonzero if success.
     */
    int initReceive(int tag);
    /// Finishes the exchange. After this call all communication buffers can be reused.
    int finishExchange();


    /**@name Methods for manipulating/testing receiver state */
    //@{
    int sendCompleted();
    int receiveCompleted();
    int testCompletion();
    int waitCompletion();

    //@}
    /// Clears all buffer contents.
    void clearBuffers();

private:
    /**
     * Resizes send buffer to needs according to toSend  array.
     * Current implementation uses EngngModel::estimateMaxPackSize function,
     * sending toSend map as parameter.
     * @param emodel Current engngModel.
     * @param packUnpackType Determines the type of packed quantity, used by emodel estimateMaxPackSize
     * service to estimate the size of pack/unpack buffer accordingly.
     * @return Nonzero if success.
     */
    template< class T > int resizeSendBuff(T *emodel, int packUnpackType);
    /**
     * Resizes receive buffer to needs according to toRecv  array.
     * Current implementation uses EngngModel::estimateMaxPackSize function,
     * sending toSend map as parameter.
     * @param emodel Current engngModel.
     * @param packUnpackType Determines the type of packed quantity, used by emodel estimateMaxPackSize
     * service to estimate the size of pack/unpack buffer accordingly.
     * @return Nonzero if success.
     */
    template< class T > int resizeRecvBuff(T *emodel, int packUnpackType);
};

template< class T > void
ProcessCommunicator :: setToSendArry(T *emodel, const IntArray &src, int packUnpackType)
{
    toSend = src;
    // toSend.sort (&NlDEIDynamicDomainComunicator::nodeSortFunct);
    //sortCommMap (toSend);
    resizeSendBuff(emodel, packUnpackType);
}


template< class T > void
ProcessCommunicator :: setToRecvArry(T *emodel, const IntArray &src, int packUnpackType)
{
    toReceive = src;
    //toReceive.sort (&NlDEIDynamicDomainComunicator::nodeSortFunct);
    //sortCommMap (toReceive);
    resizeRecvBuff(emodel, packUnpackType);
}


template< class T > int
ProcessCommunicator :: resizeSendBuff(T *emodel, int packUnpackType)
{
    int size;
    // determine space for send buffer
    size = emodel->estimateMaxPackSize(toSend, * giveProcessCommunicatorBuff()->giveSendBuff(), packUnpackType);
    giveProcessCommunicatorBuff()->resizeSendBuffer(size);
    //giveSendBuff()->resize (size);
    return 1;
}


template< class T > int
ProcessCommunicator :: resizeRecvBuff(T *emodel, int packUnpackType)
{
    int size;

    // determine space for recv buffer
    size = emodel->estimateMaxPackSize(toReceive, * giveProcessCommunicatorBuff()->giveRecvBuff(), packUnpackType);
    giveProcessCommunicatorBuff()->resizeReceiveBuffer(size);
    //giveRecvBuff()->resize (size);

    return 1;
}
} // end namespace oofem
#endif // processcomm_h
