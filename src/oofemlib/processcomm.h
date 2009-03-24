/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/processcomm.h,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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

#ifndef processcomm_h
#define processcomm_h

#ifdef __PARALLEL_MODE

#include "combuff.h"
#include "engngm.h"
#include "commbufftype.h"
#include "communicatormode.h"
#include "flotmtrx.h"


#ifndef __MAKEDEPEND
#include "mpi.h"
#endif


/**
 * The ProcessCommunicator and corresponding buffers (represented by this class)
 * are separated in order to allow share the same buffer by several communicators.
 * Here sharing means reusing for different but NON-OVERLAPPING communications.
 * if communications overlap, the different instances of ProcessCommunicatorBuff should be used!
 * The ProcessCommunicatorBuff objects are registered in corresponding communicator,
 * then if maps are available, comBuff should be resized and used in subsequent ops.
 *
 * The registration is necessary, otherwise before each send op the buffers (given probably as parameter)
 * will be resized again (size have to be computed again) and this is probably quite cost operation.
 * When comBuff will be registered, resize is needed only when maps change, and this will not occur frequently
 * (its even quite rare).
 */
class ProcessCommunicatorBuff
{
protected:
    /// Send buffer
    CommunicationBuffer *send_buff;
    /// Receive buffer
    CommunicationBuffer *recv_buff;
public:
    /// Constructor, creates empty send and receive com buffs in MPI_COMM_WORLD
    ProcessCommunicatorBuff(CommBuffType t);
    ~ProcessCommunicatorBuff() {
        if ( send_buff ) { delete send_buff; }

        if ( recv_buff ) { delete recv_buff; }
    }
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
    int packDouble(double value)  { return packArray(& value, 1); }
    /**
     * Packs array of values of given type into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @param src adress of first value in memory
     * @param n number of packed integers
     * @return nonzero if succesfull
     */
    //@{
    int packArray(const int *src, int n) { return send_buff->packArray(src, n); }
    int packArray(const long *src, int n) { return send_buff->packArray(src, n); }
    int packArray(const unsigned long *src, int n) { return send_buff->packArray(src, n); }
    int packArray(const double *src, int n) { return send_buff->packArray(src, n); }
    int packArray(const char *src, int n) { return send_buff->packArray(src, n); }
    //@}

    /**
     * Packs given IntArray  value into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @return nonzero if succesfull
     */
    int packIntArray(const IntArray &arry) { return arry.packToCommBuffer(* send_buff); }
    /**
     * Packs given FloatArray  value into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @return nonzero if succesfull
     */
    int packFloatArray(const FloatArray &arry) { return arry.packToCommBuffer(* send_buff); }
    /**
     * Packs given FloatMatrix  value into buffer.
     * Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
     * @return nonzero if succesfull
     */
    int packFloatMatrix(const FloatMatrix &mtrx) { return mtrx.packToCommBuffer(* send_buff); }
    /**
     *  Unpacks single integer value from buffer.
     *  @return nonzero if succesfull
     */
    int unpackInt(int &value)  { return unpackArray(& value, 1); }
    /**
     *  Unpacks single double value from buffer.
     *  @return nonzero if succesfull
     */
    int unpackDouble(double &value) { return unpackArray(& value, 1); }
    /**
     * unpacks array of value of given type from buffer.
     * @param dest adress of first value in memory, where to store values
     * @param n number of unpacked integers
     * @return nonzero if succesfull
     */
    //@{
    int unpackArray(int *dest, int n) { return recv_buff->unpackArray(dest, n); }
    int unpackArray(long *dest, int n) { return recv_buff->unpackArray(dest, n); }
    int unpackArray(unsigned long *dest, int n) { return recv_buff->unpackArray(dest, n); }
    int unpackArray(double *dest, int n) { return recv_buff->unpackArray(dest, n); }
    int unpackArray(char *dest, int n) { return recv_buff->unpackArray(dest, n); }
    //@}
    /**
     * Unpacks given IntArray  value from buffer.
     * @return nonzero if succesfull
     */
    int unpackIntArray(IntArray &arry) { return arry.unpackFromCommBuffer(* recv_buff); }
    /**
     * Unpacks given FloatArray  value from buffer.
     * @return nonzero if succesfull
     */
    int unpackFloatArray(FloatArray &arry) { return arry.unpackFromCommBuffer(* recv_buff); }
    /**
     * Unpacks given FloatMatrix  value from buffer.
     * @return nonzero if succesfull
     */
    int unpackFloatMatrix(FloatMatrix &mtrx) { return mtrx.unpackFromCommBuffer(* recv_buff); }
    //@}


    /**@name Methods for determining pack size of datatype to pack/unpack to/from buffer */
    //@{
    /**
     *  Returns pack size required to pack an array (c-style).
     *  @param array size
     *  @param type type id
     *  @return  pack size required
     */
    int givePackSize(MPI_Datatype type, int size) {
        int requredSpace;
        MPI_Pack_size(size, type, MPI_COMM_WORLD, & requredSpace);
        return requredSpace;
    }
    //@}


    /// Initializes send buffer to empty state. All packed data are lost.
    void initSendBuff() { send_buff->init(); }
    /// Initializes send buffer to empty state. All packed data are lost.
    void initRecvBuff() { recv_buff->init(); }
    /// initializes receiver buffers
    void init() { initSendBuff();
                  initRecvBuff(); }

    /// Initialize for packing
    void initForPacking() { send_buff->initForPacking(); }
    /// Initialize for Unpacking (data already received)
    void initForUnpacking() { recv_buff->initForUnpacking(); }

    /**
     * Initializes data exchange with associated problem.
     * if send or receive pool is empty, the send or receive communication is not preformed.
     * @param tag message tag
     * @return nonzero if success
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
     * @param tag message tag
     * @return nonzero if success
     */
    int initSend(int rank, int tag) { return send_buff->iSend(rank, tag); }
    /**
     * Initialize the receive data exchange with associate problem.
     * if receive pool is empty, the receive communication is not performed.
     * @param tag message tag
     * @return nonzero if success
     */
    int initReceive(int rank, int tag) { return recv_buff->iRecv(rank, tag); }

    /**
     *  Clears all buffer contens.
     */
    void clearBuffers() { }

    void resizeSendBuffer(int size) { send_buff->resize(size); }
    void resizeReceiveBuffer(int size) { recv_buff->resize(size); }

    /**@name Methods for manipulating/testing receiver state */
    //@{
    int sendCompleted() { return send_buff->testCompletion(); }
    int receiveCompleted() { return recv_buff->testCompletion(); }
    int testCompletion () {return (send_buff->testCompletion() && recv_buff->testCompletion());}
    int waitCompletion () {return (send_buff->waitCompletion() && recv_buff->waitCompletion());}
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
 * Class representing process communicator for NlDEIDynamic engng model.
 * Processs comunicator provides all services for communication with
 * associated remote process (problem or task).
 */
class ProcessCommunicator
{
protected:
    /// associated partition (problem) number (rank)
    int rank; // remote problem id = rank
    ///
    /// local problem
    EngngModel *localProblem;

    /// Comm Buffers representation
    ProcessCommunicatorBuff *pcBuffer;

    /// nodes to send
    IntArray toSend;
    /// nodes to receive
    IntArray toReceive;
    /// mode
    CommunicatorMode mode;

public:

    /**
     * Constructor. Creates new problem (partition) communicator associated
     * to partition with number (rank) irank.
     * @param d local problem pointer
     * @param b ProcessCommunicatorBuff to use
     * @param irank rank of associated partition
     */
    ProcessCommunicator(EngngModel *d, ProcessCommunicatorBuff *b, int irank, CommunicatorMode m = CommMode_Static);
    /// Destructor
    ~ProcessCommunicator() { }

    /**
     * Returns corresponding rank of associated partition
     */
    int giveRank() { return rank; }


    /*
     * ///  Returns send buffer of receiver.
     * CommunicationBuffer* giveSendBuff () {
     * if (pcBuffer) return pcBuffer->giveSendBuff();
     * OOFEM_ERROR ("ProcessCommunicator::giveSendBuff : ProcessCommunicatorBuff undefined");
     * return NULL;
     * }
     * /// Returns receive buffer of receiver.
     * CommunicationBuffer* giveRecvBuff () {
     * if (pcBuffer) return pcBuffer->giveRecvBuff();
     * OOFEM_ERROR ("ProcessCommunicator::giveRecvBuff : ProcessCommunicatorBuff undefined");
     * return NULL;
     * }
     */

    ///  Returns CommunicationBuffer
    ProcessCommunicatorBuff *giveProcessCommunicatorBuff() {
        if ( pcBuffer ) { return pcBuffer; }

        OOFEM_ERROR("ProcessCommunicator::giveRecvBuff : ProcessCommunicatorBuff undefined");
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
     * @param emode Engng model
     * @param src source of toSend array
     * @param packUnpackType determines the type of packed quantity, used by emodel estimateMaxPackSize
     * service to estimate the size of pack/unpack buffer accordingly.
     */
    template< class T > void setToSendArry(T *emodel, const IntArray &src, int packUnpackType);
    /**
     * Sets receiver toRecv array to src. The method assumes that toRecv array is sorted according
     * to global number associated to corresponding components given in src array. THis is necessary to
     * ensure proper pack/unpack order on local and remote problem. Recv buff is resized to
     * hold all necessary data.
     * @param emodel Engng model
     * @param src source of toRecv array
     * @param packUnpackType determines the type of packed quantity, used by emodel estimateMaxPackSize
     * service to estimate the size of pack/unpack buffer accordingly.
     */
    template< class T > void setToRecvArry(T *emodel, const IntArray &src, int packUnpackType);
    /**
     * Pack nodal data to send buff.
     * @param packFunc function used to pack nodal data in to buffer. It uses toSend array
     * to loop over required nodes.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T > int packData( T *emodel, int ( T :: *packFunc )( ProcessCommunicator & ) )
    { if ( !toSend.isEmpty() || ( this->mode == CommMode_Dynamic ) ) { giveProcessCommunicatorBuff()->initForPacking();
                                                                       return ( emodel->*packFunc )(* this); } else { return 1; } }
    /**
     * Pack nodal data to send buff.
     * @param packFunc function used to pack nodal data in to buffer. It uses toSend array
     * to loop over required nodes.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T, class P > int packData( T *emodel, P *src, int ( T :: *packFunc )( P *, ProcessCommunicator & ) )
    { if ( !toSend.isEmpty() || ( this->mode == CommMode_Dynamic ) ) { giveProcessCommunicatorBuff()->initForPacking();
                                                                       return ( emodel->*packFunc )(src, * this); } else { return 1; } }
    /**
     * Unpack nodal data from recv buff.
     * @param unpackFunc function used to unpack nodal data from buffer. It uses toRecv array
     * to loop over required nodes.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T > int unpackData( T *emodel,  int ( T :: *unpackFunc )( ProcessCommunicator & ) )
    { if ( !toReceive.isEmpty() || ( this->mode == CommMode_Dynamic ) ) { giveProcessCommunicatorBuff()->initForUnpacking();
                                                                          return ( emodel->*unpackFunc )(* this); } else { return 1; } }
    /**
     * Unpack nodal data from recv buff.
     * @param unpackFunc function used to unpack nodal data from buffer. It uses toRecv array
     * to loop over required nodes.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T, class P > int unpackData( T *emodel,  P *dest, int ( T :: *unpackFunc )( P *, ProcessCommunicator & ) )
    { if ( !toReceive.isEmpty() || ( this->mode == CommMode_Dynamic ) ) { giveProcessCommunicatorBuff()->initForUnpacking();
                                                                          return ( emodel->*unpackFunc )(dest, * this); } else { return 1; } }
    /**
     * Initializes data exchange with associated problem.
     * if send or receive pool is empty, the send or receive communication is not preformed.
     * @param tag message tag
     * @return nonzero if success
     */
    int initExchange(int tag);
    /**
     * Initialize the send data exchange with associate problem.
     * if send pool is empty, the send communication is not performed.
     * @param tag message tag
     * @return nonzero if success
     */
    int initSend(int tag);
    /**
     * Initialize the receive data exchange with associate problem.
     * if receive pool is empty, the receive communication is not performed.
     * @param tag message tag
     * @return nonzero if success
     */
    int initReceive(int tag);
    /*
     * Finishes the exchange. After this call all communication buffers can be reused.
     *
     */
    int finishExchange ();
    
    
    /**@name Methods for manipulating/testing receiver state */
    //@{
    int sendCompleted() { return giveProcessCommunicatorBuff()->sendCompleted(); }
    int receiveCompleted() { return giveProcessCommunicatorBuff()->receiveCompleted(); }
    int testCompletion() {return giveProcessCommunicatorBuff()->testCompletion();}
    int waitCompletion() {return giveProcessCommunicatorBuff()->waitCompletion();}

    //@}
    /**
     * Clears all buffer contens.
     */
    void clearBuffers();

private:
    /**
     * Resizes send buffer to needs according to toSend  array.
     * Current implementation uses EngngModel::estimateMaxPackSize function,
     * sending toSend map as parameter.
     * @param emodel current engngModel.
     * @param packUnpackType determines the type of packed quantity, used by emodel estimateMaxPackSize
     * service to estimate the size of pack/unpack buffer accordingly.
     * @return nonzero if success.
     */
    template< class T > int resizeSendBuff(T *emodel, int packUnpackType);
    /**
     * Resizes receive buffer to needs according to toRecv  array.
     * Current implementation uses EngngModel::estimateMaxPackSize function,
     * sending toSend map as parameter.
     * @param emodel current engngModel.
     * @param packUnpackType determines the type of packed quantity, used by emodel estimateMaxPackSize
     * service to estimate the size of pack/unpack buffer accordingly.
     * @return nonzero if success.
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

#endif
#endif // processcomm_h


