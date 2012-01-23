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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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


#ifndef pnldeidynamiccomm_h
#define pnldeidynamiccomm_h

#ifdef __PARALLEL_MODE

#include "communicator.h"
#include "pnldeidynamic.h"
#include "error.h"

namespace oofem {
/**
 * Class representing communicator for PNlDEIDynamic engng model.
 * It is attribute of pnlDeiDynamic Engng model.
 * Domain comunicator provides all services for communication with
 * associated remote domains. It manages several problem comunicators.
 */
class PNlDEIDynamicComunicator : public Communicator< PNlDEIDynamic >
{
protected:
    // setup mode
    PNlDEIDynamicComunicatorMode mode;

public:
    /**
     * Constructor. Creates new communicator associated
     * to partition with number (rank) irank.
     * @param irank rank of associated partition
     * @param size #number of colaborating processes
     * @param mode comunicator mode.
     */
    PNlDEIDynamicComunicator(PNlDEIDynamic *emodel, int rank, int size, PNlDEIDynamicComunicatorMode mode);
    /// Destructor
    virtual ~PNlDEIDynamicComunicator();

    /**
     * Service for setting up the communication patterns with other remote problems.
     * Sets up the toSend and toRecv attributes in associated problem communicators.
     */
    void setUpCommunicationMaps(EngngModel *pm);
    /**
     * Assigns  given ToSendMap  to given problem communicator.
     * Sorts map according to global entity (dofmanagers or element) numbers to ensure, that local and
     * corresponding remote problem have the identical map with same ordering. This will ensure proper packing/unpacking
     * order. The corresponding problemCommunicator buffer takes care about resizing itself accordingly to hold all
     * outcomming/incomming data using engngModel->estimateMaxPackSize service.
     * @param problemComm domain comm which send map will be set
     * @param map send map
     */
    int setProblemComunicatorToSendArry(ProblemComunicator< PNlDEIDynamic > *problemComm, IntArray &map);
    /**
     * Assigns  given ToRecvMap  to given problem communicator.
     * Sorts map according to global entity (dofmanagers or element) numbers to ensure, that local and
     * corresponding remote problem have the identical map with same ordering. This will ensure proper packing/unpacking
     * order. The corresponding problemCommunicator buffer takes care about resizing itself  accordingly to hold all
     * outcomming/incomming data using engngModel->estimateMaxPackSize service.
     * @param domainComm problem comm which Recv map will be set
     * @param map recv map
     */
    int setProblemComunicatorToRecvArry(ProblemComunicator< PNlDEIDynamic > *problemComm, IntArray &map);

private:
    /**
     * Sorts given communication map, containing local DofManager numbers acording to their
     * corresponding global numbers. It could not be sorted by standard techniques, because
     * it is necessary to ask DofMAnager form domain and determine its global Number.
     * @param cmp comparison function must return a negative value if first argument is less than the second,
     * zero if the arguments are equal, and a positive number otherwise.
     */
    void sortCommMap( IntArray & map, int ( PNlDEIDynamicComunicator :: *cmp )( int, int ) );
    /// Implementation of Quiksort algorithm
    void quickSortCommMap( IntArray & map, int l, int r, int ( PNlDEIDynamicComunicator :: *cmp )( int, int ) );
    /// Partitioning used in quiksort
    int quickSortPartition( IntArray & map, int l, int r, int ( PNlDEIDynamicComunicator :: *cmp )( int, int ) );

    /// Global dofManager number comparison function
    int DofManCmp(int, int);
    /// Global element comparison function
    int ElemCmp(int, int);

    /**
     * Service for setting up the receiver for node cut communication patterns with other remote problems.
     * Sets up the toSend and toRecv attributes in associated problem communicators.
     */
    void setUpCommunicationMapsForNodeCut(EngngModel *domain);
    /**
     * Service for setting up the receiver for element cut communication patterns with other remote problems.
     * Sets up the toSend and toRecv attributes in associated problem communicators.
     */
    void setUpCommunicationMapsForElementCut(EngngModel *domain);
    /**
     * Service for setting up the receiver for remote element communication patterns with other remote problems.
     * Sets up the toSend and toRecv attributes in associated problem communicators.
     */
    void setUpCommunicationMapsForRemoteElementMode(EngngModel *domain);
};
} // end namespace oofem
#endif
#endif // pnldeidynamiccomm_h
