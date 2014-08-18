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

#ifndef communicator_h
#define communicator_h

#include "oofemcfg.h"
#include "processcomm.h"
#include "commbufftype.h"
#include "communicatormode.h"
#include "error.h"

namespace oofem {
class EngngModel;

/**
 * The Communicator and corresponding buffers (represented by this class)
 * are separated in order to allow share the same buffer by several communicators.
 * Here sharing means reusing for different but <em>non-overlapping</em> communications.
 * if communications overlap, the different instances of CommunicatorBuff should be used!
 * The CommunicatorBuff objects are registered in corresponding communicator,
 * then if maps are available, comBuff should be resized and used in subsequent ops.
 *
 * The registration is necessary, otherwise before each send op the buffers (given probably as parameter)
 * will be resized again (size have to be computed again) and this is probably quite cost operation.
 * When comBuff will be registered, resize is needed only when maps change, and this will not occur frequently
 * (its even quite rare).
 */
class OOFEM_EXPORT CommunicatorBuff
{
protected:
    /// Number of processes.
    int size;
    /// Array of process communicators.
    ProcessCommunicatorBuff **processCommBuffs;

public:
    CommunicatorBuff(int s, CommBuffType t = CBT_static);
    ~CommunicatorBuff();

    /**
     * Returns i-th process communicator buff. The process comm buffs are numbered from rank 0.
     * @param i Process communicator buff index [0..size-1].
     * @return Pointer to corresponding process communicator buff, NULL otherwise.
     */
    ProcessCommunicatorBuff *
    giveProcessCommunicatorBuff(int i) {
        if ( i < size ) {
            return processCommBuffs [ i ];
        } else {
            return NULL;
        }
    }
};


/**
 * Class representing communicator.
 * It is usually attribute of an engineering model.
 * Problem communicator provides all services for communication with
 * associated remote problems. It manages several process (task) communicators.
 *
 * The communicator mode determines the communication:
 *
 * (Static) The mode can be static, meaning that each node can assemble its communication maps
 * independently (or by independent communication). This implies that the size of
 * communication buffers is known in advance. Also if no data are planned to send to remote node, there
 * is no communication with this node (both sender and receiver know that there will be no data to send).
 *
 * (Dynamic) In this case the communication pattern and the amount of data sent between nodes is
 * not known in advance. This requires to use dynamic (packeted) buffering.
 *
 */
class OOFEM_EXPORT Communicator
{
protected:
    /// Rank of process.
    int rank;
    /// Number of processes.
    int size;
    /// Array of process communicators.
    ProcessCommunicator **processComms;
    /// Engineering model.
    EngngModel *engngModel;
    /// Mode.
    CommunicatorMode mode;

public:
    /**
     * Constructor. Creates new communicator associated
     * to partition with number (rank) irank.
     * @param emodel Associated engineering model.
     * @param buff Communicator buffer.
     * @param rank Rank of associated partition.
     * @param size Number of collaborating processes.
     * @param mode Communicator mode.
     */
    Communicator(EngngModel * emodel, CommunicatorBuff * buff, int rank, int size, CommunicatorMode mode = CommMode_Static);
    /// Destructor
    virtual ~Communicator();

    /**
     * Returns i-th problem communicator. The problems are numbered from rank 0.
     * @param i Problem communicator index [0..size-1].
     * @return Pointer to corresponding communicator, NULL otherwise.
     */
    ProcessCommunicator *
    giveProcessCommunicator(int i) {
        if ( i < size ) {
            return processComms [ i ];
        } else {
            return NULL;
        }
    }

    /**
     * Pack all problemCommunicators data to their send buffers.
     * @param ptr Pointer problem communicator.
     * @param packFunc Function used to pack nodal data in to buffer.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T > int packAllData( T *ptr, int ( T :: *packFunc )( ProcessCommunicator & ) );
    /**
     * Pack all problemCommuncators data to their send buffers.
     * @param ptr Pointer problem communicator.
     * @param src Pointer to source.
     * @param packFunc Function used to pack nodal data in to buffer.
     * @see NlDEIDynamic_Unpack_func
     */
    //template <class T> int packAllData (T* ptr, FloatArray* src, int (T::*packFunc) (FloatArray*, ProcessCommunicator&));
    template< class T, class P > int packAllData( T *ptr, P *src, int ( T :: *packFunc )( P *, ProcessCommunicator & ) );
    /**
     * Unpack all problemCommuncators data from recv buffers.
     * Waits until receive completion before unpacking buffer.
     * @param ptr Pointer problem communicator.
     * @param unpackFunc Function used to unpack nodal data from buffer.
     * @see NlDEIDynamic_Unpack_func
     */
    template< class T > int unpackAllData( T *ptr, int ( T :: *unpackFunc )( ProcessCommunicator & ) );
    /**
     * Unpack all problemCommuncators data from recv buffers.
     * Waits until receive completion before unpacking buffer.
     * @param ptr Pointer problem communicator.
     * @param src Pointer to source.
     * @param unpackFunc Function used to unpack nodal data from buffer.
     * @see NlDEIDynamic_Unpack_func
     */
    //template <class T> int unpackAllData (T* ptr, FloatArray* dest, int (T::*unpackFunc) (FloatArray*, ProcessCommunicator&));
    template< class T, class P > int unpackAllData( T *ptr, P *src, int ( T :: *unpackFunc )( P *, ProcessCommunicator & ) );
    /**
     * Initializes data exchange with all problems.
     * if send or receive pool is empty, communication is not performed.
     * @param tag message tag
     */
    int initExchange(int tag);
    /**
     * Initializes data send exchange with all problems.
     * if send  pool is empty, communication is not performed.
     * @param tag Message tag.
     * @return Nonzero if successful.
     */
    int initSend(int tag);
    /**
     * Initializes data receive exchange with all problems.
     * if receive pool is empty, communication is not performed.
     * @param tag Message tag.
     * @return Nonzero if successful.
     */
    int initReceive(int tag);
    /**
     * Finishes the exchange. After this call all communication buffers can be reused.
     * @return Nonzero if successful.
     */
    int finishExchange();

    /**
     * Clears all buffer content.
     */
    void clearBuffers();
    /**
     * Service for setting up the communication patterns with other remote processes.
     * Sets up the toSend and toRecv attributes in associated problem communicators.
     * @param pm Engineering model to use for setup.
     */
    virtual void setUpCommunicationMaps(EngngModel *pm) { }

    /// Returns string for prepending output (used by error reporting macros).
    std :: string errorInfo(const char *func) const;
};

template< class T > int
Communicator :: packAllData( T *ptr, int ( T :: *packFunc )( ProcessCommunicator & ) )
{
    int i = size, result = 1;

    if ( size ) {
        for ( i = 0; i < size; i++ ) {
            result &= giveProcessCommunicator(i)->packData(ptr, packFunc);
        }
    }

    return result;
}

/*
 * template <class T> int
 * Communicator :: packAllData (T* ptr, FloatArray* src, int (T::*packFunc) (FloatArray*, ProcessCommunicator&))
 * {
 * int i = size, result = 1;
 *
 * if (size)
 * for (i=0; i< size; i++) result &= giveProcessCommunicator(i)->packData (ptr, src, packFunc);
 * return result;
 * }
 */
template< class T, class P > int
Communicator :: packAllData( T *ptr, P *src, int ( T :: *packFunc )( P *, ProcessCommunicator & ) )
{
    int i = size, result = 1;

    if ( size ) {
        for ( i = 0; i < size; i++ ) {
            result &= giveProcessCommunicator(i)->packData(ptr, src, packFunc);
        }
    }

    return result;
}

template< class T > int
Communicator :: unpackAllData( T *ptr, int ( T :: *unpackFunc )( ProcessCommunicator & ) )
{
    int i, received, num_recv = 0, result = 1;
    IntArray recvFlag(size);
    //MPI_Status status;

    for  ( i = 0; i < size; i++ ) {
        // receive if receive map is not empty or mode is dynamic
        if ( ( giveProcessCommunicator(i)->giveToRecvMap()->giveSize() ) ||
            ( this->mode == CommMode_Dynamic ) ) {
            recvFlag.at(i + 1) = 1;
            num_recv++;
        }
    }

    while ( num_recv-- ) {
        // wait for any completion
        while ( 1 ) {
            received = 0;
            for  ( i = 0; i < size; i++ ) {
                if ( recvFlag.at(i + 1) ) {
                    //if (giveProcessCommunicator(i)->giveRecvBuff()->testCompletion()) {
                    if ( giveProcessCommunicator(i)->receiveCompleted() ) {
 #ifdef __VERBOSE_PARALLEL
                        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received data from partition %3d\n",
                                        rank, "Communicator :: unpackAllData", i);
 #endif

                        recvFlag.at(i + 1) = 0;
                        result &= giveProcessCommunicator(i)->unpackData(ptr, unpackFunc);
                        received = 1;
                        break;
                    }
                }
            }

            if ( received ) {
                break;
            }
        }
    }

 #ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier started", rank)
 #endif

    MPI_Barrier(MPI_COMM_WORLD);

 #ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier finished", rank)
 #endif

    return result;
}


/*
 * template <class T> int
 * Communicator :: unpackAllData (T* ptr, FloatArray* dest, int (T::*unpackFunc) (FloatArray*, ProcessCommunicator&))
 * {
 * int i, received, num_recv = 0, result = 1;
 * IntArray recvFlag (size);
 * //MPI_Status status;
 *
 * for  (i=0; i<size; i++) {
 * // receive if receive map is not empty or mode is dynamic
 * if ((giveProcessCommunicator(i)->giveToRecvMap()->giveSize()) ||
 *     (this->mode == CommMode_Dynamic)) {
 *   recvFlag.at(i+1) = 1;
 *   num_recv ++;
 * }
 * }
 *
 * while (num_recv--) {
 *
 * // wait for any completion
 * while (1) {
 * received = 0;
 * for  (i=0; i<size; i++)  {
 *  if (recvFlag.at(i+1)) {
 *    //if (giveProcessCommunicator(i)->giveRecvBuff()->testCompletion()) {
 *    if (giveProcessCommunicator(i)->receiveCompleted()) {
 *
 *
 ****#ifdef __VERBOSE_PARALLEL
 *     OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received data from partition %3d\n",
 *                     rank,"Communicator :: unpackAllData", i);
 ****#endif
 *
 *    recvFlag.at(i+1) = 0;
 *    result &= giveProcessCommunicator(i)->unpackData (ptr, dest, unpackFunc);
 *    received = 1;
 *    break;
 *   }
 *  }
 * }
 * if (received) break;
 * }
 * }
 *
 ****#ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier started",rank)
 ****#endif
 *
 * MPI_Barrier (MPI_COMM_WORLD);
 *
 ****#ifdef __VERBOSE_PARALLEL
 * VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier finished",rank)
 ****#endif
 *
 * return result;
 * }
 */

template< class T, class P > int
Communicator :: unpackAllData( T *ptr, P *dest, int ( T :: *unpackFunc )( P *, ProcessCommunicator & ) )
{
    int i, received, num_recv = 0, result = 1;
    IntArray recvFlag(size);
    //MPI_Status status;

    for  ( i = 0; i < size; i++ ) {
        // receive if receive map is not empty or mode is dynamic
        if ( ( giveProcessCommunicator(i)->giveToRecvMap()->giveSize() ) ||
            ( this->mode == CommMode_Dynamic ) ) {
            recvFlag.at(i + 1) = 1;
            num_recv++;
        }
    }

    while ( num_recv-- ) {
        // wait for any completion
        while ( 1 ) {
            received = 0;
            for  ( i = 0; i < size; i++ ) {
                if ( recvFlag.at(i + 1) ) {
                    //if (giveProcessCommunicator(i)->giveRecvBuff()->testCompletion()) {
                    if ( giveProcessCommunicator(i)->receiveCompleted() ) {
 #ifdef __VERBOSE_PARALLEL
                        OOFEM_LOG_DEBUG("[process rank %3d]: %-30s: Received data from partition %3d\n",
                                        rank, "Communicator :: unpackAllData", i);
 #endif

                        recvFlag.at(i + 1) = 0;
                        result &= giveProcessCommunicator(i)->unpackData(ptr, dest, unpackFunc);
                        received = 1;
                        break;
                    }
                }
            }

            if ( received ) {
                break;
            }
        }
    }

 #ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier started", rank)
 #endif

    MPI_Barrier(MPI_COMM_WORLD);

 #ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT("Communicator :: unpackAllData", "Synchronize barrier finished", rank)
 #endif

    return result;
}
} // end namespace oofem
#endif // communicator_h
