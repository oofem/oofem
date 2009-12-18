/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/problemcomm.h,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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
// Class ProblemComunicator
//

#ifndef problemcomm_h
#define problemcomm_h

#ifdef __PARALLEL_MODE

#include "problemcommunicatormode.h"
#include "communicator.h"
#include "error.h"

namespace oofem {

/**
 * Class representing communicator for engng model.
 * It is assumed to be an attribute of  Engng model.
 * Domain comunicator provides all services for communication with
 * associated remote domains. It manages several process comunicators.
 */
class ProblemCommunicator : public Communicator
{
protected:
    // setup mode
    ProblemCommunicatorMode mode;
    bool initialized;

public:
    /**
     * Constructor. Creates new communicator associated
     * to partition with number (rank) irank.
     * @param irank rank of associated partition
     * @param size #number of colaborating processes
     * @param mode communicator mode.
     */
    ProblemCommunicator(EngngModel *emodel, CommunicatorBuff *b, int rank, int size, ProblemCommunicatorMode mode);
    /// Destructor
    ~ProblemCommunicator();

    /**
     * Service for setting up the communication patterns with other remote processs.
     * Sets up the toSend and toRecv attributes in associated process communicators.
     * @param excludeSelfCommFlag if set to true, the communication map of receiver with itself will be forced to be empty, otherwise it will be assembled.
     * @param forceReinit forces reinitilaization
     */
    void setUpCommunicationMaps(EngngModel *pm, bool excludeSelfCommFlag, bool forceReinit = false);
    /**
     * Assigns  given ToSendMap  to given process communicator.
     * Sorts map according to global entity (dofmanagers or element) numbers to ensure, that local and
     * corresponding remote process have the identical map with same ordering. This will ensure proper packing/unpacking
     * order. The corresponding processCommunicator buffer takes care about resizing itself accordingly to hold all
     * outcomming/incomming data using engngModel->estimateMaxPackSize service.
     * @param processComm domain comm which send map will be set
     * @param map send map
     */
    int setProcessCommunicatorToSendArry(ProcessCommunicator *processComm, IntArray &map);
    /**
     * Assigns  given ToRecvMap  to given process communicator.
     * Sorts map according to global entity (dofmanagers or element) numbers to ensure, that local and
     * corresponding remote process have the identical map with same ordering. This will ensure proper packing/unpacking
     * order. The corresponding processCommunicator buffer takes care about resizing itself  accordingly to hold all
     * outcomming/incomming data using engngModel->estimateMaxPackSize service.
     * @param domainComm process comm which Recv map will be set
     * @param map recv map
     */
    int setProcessCommunicatorToRecvArry(ProcessCommunicator *processComm, IntArray &map);
private:

    /**
     * Sorts given communication map, containing local DofManager numbers according to their
     * corresponding global numbers. It could not be sorted by standard techniques, because
     * it is necessary to ask DofMAnager form domain and determine its global Number.
     * @param cmp comparison function must return a negative value if first argument is less than the second,
     * zero if the arguments are equal, and a positive number otherwise.
     */
    void sortCommMap( IntArray & map, int ( ProblemCommunicator :: * cmp )( int, int ) );
    /// Implementation of Quiksort algorithm
    void quickSortCommMap( IntArray & map, int l, int r, int ( ProblemCommunicator :: * cmp )( int, int ) );
    /// Partitioning used in quiksort
    int quickSortPartition( IntArray & map, int l, int r, int ( ProblemCommunicator :: * cmp )( int, int ) );

    /// global dofManager number comparison function
    int DofManCmp(int, int);
    /// global element comparison function
    int ElemCmp(int, int);

    /**
     * Service for setting up the receiver for node cut communication patterns with other remote processs.
     * Sets up the toSend and toRecv attributes in associated process communicators.
     * @param excludeSelfCommFlag if set to true, the communication map of receiver with itself will be forced to be empty, otherwise it will be assembled.
     */
    void setUpCommunicationMapsForNodeCut(EngngModel *domain, bool excludeSelfCommFlag);
    /**
     * Service for setting up the receiver for element cut communication patterns with other remote processs.
     * Sets up the toSend and toRecv attributes in associated process communicators.
     * @param excludeSelfCommFlag if set to true, the communication map of receiver with itself will be forced to be empty, otherwise it will be assembled.
     */
    void setUpCommunicationMapsForElementCut(EngngModel *domain, bool excludeSelfCommFlag);
    /**
     * Service for setting up the receiver for remote element communication patterns with other remote processs.
     * Sets up the toSend and toRecv attributes in associated process communicators.
     * @param excludeSelfCommFlag if set to true, the communication map of receiver with itself will be forced to be empty, otherwise it will be assembled.
     */
    void setUpCommunicationMapsForRemoteElementMode(EngngModel *domain, bool excludeSelfCommFlag);
};

} // end namespace oofem
#endif
#endif // problemcomm_h
