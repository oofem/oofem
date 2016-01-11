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

#include "nodalaveragingrecoverymodel.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "engngm.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "processcomm.h"
 #include "communicator.h"
#endif

namespace oofem {
NodalAveragingRecoveryModel :: NodalAveragingRecoveryModel(Domain *d) : NodalRecoveryModel(d)
{ }

NodalAveragingRecoveryModel :: ~NodalAveragingRecoveryModel()
{ }

int
NodalAveragingRecoveryModel :: recoverValues(Set elementSet, InternalStateType type, TimeStep *tStep)
{
    int nnodes = domain->giveNumberOfDofManagers();
    IntArray regionNodalNumbers(nnodes);
    IntArray regionDofMansConnectivity;
    FloatArray lhs, val;


    if ( ( this->valType == type ) && ( this->stateCounter == tStep->giveSolutionStateCounter() ) ) {
        return 1;
    }

#ifdef __PARALLEL_MODE
    bool parallel = this->domain->giveEngngModel()->isParallel();
    if ( parallel ) {
        this->initCommMaps();
    }
#endif

    // clear nodal table
    this->clear();

    int regionValSize = 0;
    int regionDofMans;

    // loop over elements and determine local region node numbering and determine and check nodal values size
    if ( this->initRegionNodeNumbering(regionNodalNumbers, regionDofMans, elementSet) == 0 ) {
        return 0;
    }

    regionDofMansConnectivity.resize(regionDofMans);
    regionDofMansConnectivity.zero();

    IntArray elements = elementSet.giveElementList();
    // assemble element contributions
    for ( int i = 1; i <= elements.giveSize(); i++ ) {
        int ielem = elements.at(i);
        NodalAveragingRecoveryModelInterface *interface;
        Element *element = domain->giveElement(ielem);

        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }

        // If an element doesn't implement the interface, it is ignored.
        if ( ( interface = static_cast< NodalAveragingRecoveryModelInterface * >
                           ( element->giveInterface(NodalAveragingRecoveryModelInterfaceType) ) ) == NULL ) {
            //abort();
            continue;
        }

        int elemNodes = element->giveNumberOfDofManagers();
        // ask element contributions
        for ( int elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
            int node = element->giveDofManager(elementNode)->giveNumber();
            interface->NodalAveragingRecoveryMI_computeNodalValue(val, elementNode, type, tStep);
            // if the element cannot evaluate this variable, it is ignored
            if ( val.giveSize() == 0 ) {
                continue;
            } else if ( regionValSize == 0 ) {
                regionValSize = val.giveSize();
                lhs.resize(regionDofMans * regionValSize);
                lhs.zero();
            } else if ( val.giveSize() != regionValSize ) {
                OOFEM_LOG_RELEVANT("NodalAveragingRecoveryModel :: size mismatch for InternalStateType %s, ignoring all elements that doesn't use the size %d\n", __InternalStateTypeToString(type), regionValSize);
                continue;
            }
            int eq = ( regionNodalNumbers.at(node) - 1 ) * regionValSize;
            for ( int j = 1; j <= regionValSize; j++ ) {
                lhs.at(eq + j) += val.at(j);
            }

            regionDofMansConnectivity.at( regionNodalNumbers.at(node) )++;
        }
    } // end assemble element contributions

#ifdef __PARALLEL_MODE
    if ( parallel ) {
        this->exchangeDofManValues(lhs, regionDofMansConnectivity, regionNodalNumbers, regionValSize);
    }
#endif

    // solve for recovered values of active region
    for ( int inode = 1; inode <= nnodes; inode++ ) {
        if ( regionNodalNumbers.at(inode) ) {
            int eq = ( regionNodalNumbers.at(inode) - 1 ) * regionValSize;
            for ( int i = 1; i <= regionValSize; i++ ) {
                if ( regionDofMansConnectivity.at( regionNodalNumbers.at(inode) ) > 0 ) {
                    lhs.at(eq + i) /= regionDofMansConnectivity.at( regionNodalNumbers.at(inode) );
                } else {
                    OOFEM_WARNING("values of dofmanager %d undetermined", inode);
                    lhs.at(eq + i) = 0.0;
                }
            }
        }
    }

    // update recovered values
    this->updateRegionRecoveredValues(regionNodalNumbers, regionValSize, lhs);

    this->valType = type;
    this->stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}

#ifdef __PARALLEL_MODE

void
NodalAveragingRecoveryModel :: initCommMaps()
{
    if ( initCommMap ) {
        EngngModel *emodel = domain->giveEngngModel();
        commBuff = new CommunicatorBuff(emodel->giveNumberOfProcesses(), CBT_dynamic);
        communicator = new NodeCommunicator(emodel, commBuff, emodel->giveRank(),
                                            emodel->giveNumberOfProcesses());
        communicator->setUpCommunicationMaps(domain->giveEngngModel(), true, true);
        OOFEM_LOG_INFO("NodalAveragingRecoveryModel :: initCommMaps: initialized comm maps\n");
        initCommMap = false;
    }
}

void
NodalAveragingRecoveryModel :: exchangeDofManValues(FloatArray &lhs, IntArray &regionDofMansConnectivity,
                                                    IntArray &regionNodalNumbers, int regionValSize)
{
    parallelStruct ls( &lhs, &regionDofMansConnectivity, &regionNodalNumbers, regionValSize);

    // exchange data for shared nodes
    communicator->packAllData(this, & ls, & NodalAveragingRecoveryModel :: packSharedDofManData);
    communicator->initExchange(789);
    communicator->unpackAllData(this, & ls, & NodalAveragingRecoveryModel :: unpackSharedDofManData);
    communicator->finishExchange();
}

int
NodalAveragingRecoveryModel :: packSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1, size;
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    IntArray const *toSendMap = processComm.giveToSendMap();

    size = toSendMap->giveSize();
    for ( int i = 1; i <= size; i++ ) {
        // toSendMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node value is available for given region
        int indx = s->regionNodalNumbers->at( toSendMap->at(i) );
        if ( indx ) {
            // pack "1" to indicate that for given shared node this is a valid contribution
            result &= pcbuff->write(1);
            result &= pcbuff->write( s->regionDofMansConnectivity->at(indx) );
            int eq = ( indx - 1 ) * s->regionValSize;
            for ( int j = 1; j <= s->regionValSize; j++ ) {
                result &= pcbuff->write( s->lhs->at(eq + j) );
            }
        } else {
            // ok shared node is not in active region (determined by s->regionNodalNumbers)
            result &= pcbuff->write(0);
        }
    }

    return result;
}

int
NodalAveragingRecoveryModel :: unpackSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm)
{
    int result = 1;
    int size, flag, intValue;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    double value;

    size = toRecvMap->giveSize();
    for ( int i = 1; i <= size; i++ ) {
        int indx = s->regionNodalNumbers->at( toRecvMap->at(i) );
        // toRecvMap contains all shared dofmans with remote partition
        // one has to check, if particular shared node received contribution is available for given region
        result &= pcbuff->read(flag);
        if ( flag ) {
            // "1" to indicates that for given shared node this is a valid contribution
            result &= pcbuff->read(intValue);
            // now check if we have a valid number
            if ( indx ) {
                s->regionDofMansConnectivity->at(indx) += intValue;
            }

            int eq = ( indx - 1 ) * s->regionValSize;
            for ( int j = 1; j <= s->regionValSize; j++ ) {
                result &= pcbuff->read(value);
                if ( indx ) {
                    s->lhs->at(eq + j) += value;
                }
            }
        }
    }

    return result;
}

#endif
} // end namespace oofem
