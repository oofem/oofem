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

#ifndef feticommunicator_h
#define feticommunicator_h

#include "../sm/FETISolver/fetiboundarydofman.h"
#include "communicator.h"
#include "error.h"

#include <vector>

using namespace std;

namespace oofem {
class FETISolver;

/**
 * Class representing communicator for FETI solver.
 * It is attribute of FETI solver numerical method class running on master partition (rank equal to 0).
 * This Communicator provides necessary services for communication with
 * associated slave partitions. It manages several domain communicators, each for communication with
 * particular partition.
 */
class FETICommunicator : public Communicator
{
public:
    /// Enumeration used to define necessary communication tags, used to identify different messages send/received.
    enum { FETICommunicatorZeroTag, NumberOfBoundaryDofManagersMsg, BoundaryDofManagersRecMsg };

protected:
    /// Number of equations at master level (determined form boundary nodes).
    int numberOfEquations;
    /// List of boundary dof managers records.
    vector< FETIBoundaryDofManager >boundaryDofManList;
    /**
     * Master communication map. Not stored in corresponding domain comm, but required in order to
     * allow direct (memory) mapping instead of communication.
     */
    IntArray masterCommMap;

public:
    /**
     * Creates new communicator.
     * @param emodel Engineering model for communication.
     * @param b Communicator buffer.
     * @param rank Rank of associated partition.
     * @param size Number of collaborating processes.
     */
    FETICommunicator(EngngModel * emodel, CommunicatorBuff * b, int rank, int size);
    /// Destructor
    virtual ~FETICommunicator();

    int giveNumberOfDomainEquations() { return numberOfEquations; }

    void setUpCommunicationMaps(EngngModel *pm);
    /**
     * Returns reference to i-th boundary dof manager.
     * Available only on master partition.
     */
    FETIBoundaryDofManager *giveDofManager(int i) {
        return & boundaryDofManList [ i - 1 ];
    }
    /**
     * Returns pointer to master comm map stored in receiver.
     */
    IntArray *giveMasterCommMapPtr() { return & masterCommMap; }
};
} // end namespace oofem

#endif // feticommunicator_h
