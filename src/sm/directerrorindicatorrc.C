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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "errorestimator.h"
#include "directerrorindicatorrc.h"
#include "domain.h"
#include "element.h"
#include "conTable.h"
#include "mathfem.h"
#include "timestep.h"

namespace oofem {
DirectErrorIndicatorRC :: DirectErrorIndicatorRC(int n, ErrorEstimator *e) : RemeshingCriteria(n, e)
{
    stateCounter = -1;
#ifdef __PARALLEL_MODE
    dofManDensityExchangeFlag = true;
#endif
}

DirectErrorIndicatorRC :: ~DirectErrorIndicatorRC()
{
#ifdef __PARALLEL_MODE
#endif
}

void
DirectErrorIndicatorRC :: giveNodeChar(int inode, TimeStep *tStep, double &indicatorVal, double &currDensity)
{
    currDensity = this->giveDofManDensity(inode);
    indicatorVal =  this->giveDofManIndicator(inode, tStep);
}


double
DirectErrorIndicatorRC :: giveDofManDensity(int num)
{
#ifdef __PARALLEL_MODE
    Domain *d = this->giveDomain();
    if ( d->giveDofManager(num)->isShared() ) {
        return this->sharedDofManDensities [ num ];
    } else {
        return this->giveLocalDofManDensity(num);
    }

#else
    return this->giveLocalDofManDensity(num);

#endif
}

double
DirectErrorIndicatorRC :: giveLocalDofManDensity(int num)
{
    int isize, i;
    double currDensity = 0.0;
    const IntArray *con;
    Domain *d = this->giveDomain();
    ConnectivityTable *ct = d->giveConnectivityTable();
    DirectErrorIndicatorRCInterface *interface;

    con = ct->giveDofManConnectivityArray(num);
    isize = con->giveSize();

    for ( i = 1; i <= isize; i++ ) {
        // ask for indicator variable value and determine current mesh density
        interface = ( DirectErrorIndicatorRCInterface * )
                    d->giveElement( con->at(i) )->giveInterface(DirectErrorIndicatorRCInterfaceType);
        if ( !interface ) {
            OOFEM_WARNING2( "DirectErrorIndicatorRC::giveRequiredDofManDensity: elem %d does not support DirectErrorIndicatorRCInterface", con->at(i) );
        }

        if ( i == 1 ) {
            currDensity = interface->DirectErrorIndicatorRCI_giveCharacteristicSize();
        } else {
            currDensity =  max( currDensity, interface->DirectErrorIndicatorRCI_giveCharacteristicSize() );
        }
    }

    return currDensity;
}


double
DirectErrorIndicatorRC :: giveDofManIndicator(int num, TimeStep *tStep)
{
#ifdef __PARALLEL_MODE
    Domain *d = this->giveDomain();
    if ( d->giveDofManager(num)->isShared() ) {
        return this->sharedDofManIndicatorVals [ num ];
    } else {
        return this->giveLocalDofManIndicator(num, tStep);
    }

#else
    return this->giveLocalDofManIndicator(num, tStep);

#endif
}

double
DirectErrorIndicatorRC :: giveLocalDofManIndicator(int inode, TimeStep *tStep)
{
    int isize, i;
    const IntArray *con;
    Domain *d = this->giveDomain();
    ConnectivityTable *ct = d->giveConnectivityTable();
    DirectErrorIndicatorRCInterface *interface;
    double indicatorVal = 0.0;

    con = ct->giveDofManConnectivityArray(inode);
    isize = con->giveSize();

    for ( i = 1; i <= isize; i++ ) {
        // ask for indicator variable value and determine current mesh density
        interface = ( DirectErrorIndicatorRCInterface * )
                    d->giveElement( con->at(i) )->giveInterface(DirectErrorIndicatorRCInterfaceType);
        if ( !interface ) {
            OOFEM_WARNING2( "DirectErrorIndicatorRC::giveRequiredDofManDensity: element %d does not support DirectErrorIndicatorRCInterface", con->at(i) );
        }

        if ( i == 1 ) {
            indicatorVal = ee->giveElementError(indicatorET, d->giveElement( con->at(i) ), tStep);
        } else {
            indicatorVal = max( indicatorVal, ee->giveElementError(indicatorET, d->giveElement( con->at(i) ), tStep) );
        }
    }

    return indicatorVal;
}



int
DirectErrorIndicatorRC :: estimateMeshDensities(TimeStep *tStep)
{
    Domain *d = this->giveDomain();
    int inode, nnodes = d->giveNumberOfDofManagers();
    double indicatorVal, currDensity, proposedDensity;

    if ( stateCounter == tStep->giveSolutionStateCounter() ) {
        return 1;
    }

#ifdef __PARALLEL_MODE
    if ( initCommMap ) {
        communicator->setUpCommunicationMaps(d->giveEngngModel(), true, true);
        OOFEM_LOG_INFO("DirectErrorIndicatorRC :: estimateMeshDensities: initialized comm maps\n");
        initCommMap = false;
    }


    this->exchangeDofManDensities();
    this->exchangeDofManIndicatorVals(tStep);
#endif

    this->currStrategy  = NoRemeshing_RS;
    this->nodalDensities.resize(nnodes);

    for ( inode = 1; inode <= nnodes; inode++ ) {
        this->giveNodeChar(inode, tStep, indicatorVal, currDensity);

        if ( indicatorVal < minIndicatorLimit ) {
            this->nodalDensities.at(inode) = zeroIndicatorDensity; //currDensity;
        } else if ( indicatorVal >= maxIndicatorLimit ) {
            this->nodalDensities.at(inode) = proposedDensity = maxIndicatorDensity;

            if ( proposedDensity < ( currDensity * this->remeshingDensityRatioToggle ) ) {
                this->currStrategy  = RemeshingFromCurrentState_RS;
            }
        } else {
            // evaluate the required size
            proposedDensity = minIndicatorDensity +
                              ( indicatorVal - minIndicatorLimit ) * ( maxIndicatorDensity - minIndicatorDensity ) / ( maxIndicatorLimit - minIndicatorLimit );
            //proposedDensity = min (currDensity, proposedDensity);
            this->nodalDensities.at(inode) = proposedDensity;

            if ( proposedDensity < ( currDensity * this->remeshingDensityRatioToggle ) ) {
                this->currStrategy  = RemeshingFromCurrentState_RS;
            }
        }
    }

#ifdef __PARALLEL_MODE
    // exchange strategies between nodes to ensure consistency
    int myStrategy = this->currStrategy, globalStrategy;
    MPI_Allreduce(& myStrategy, & globalStrategy, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    this->currStrategy = ( RemeshingStrategy ) globalStrategy;
#endif

    // remember time stamp
    stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}


double
DirectErrorIndicatorRC :: giveRequiredDofManDensity(int num, TimeStep *tStep, int relative)
{
    this->estimateMeshDensities(tStep);
    if ( relative ) {
        double currDens;
        currDens = this->giveDofManDensity(num);
        return this->nodalDensities.at(num) / currDens;
    } else {
        return this->nodalDensities.at(num);
    }
}


IRResultType
DirectErrorIndicatorRC :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, minIndicatorLimit, IFT_DirectErrorIndicatorRC_minlim, "minlim"); // Macro
    IR_GIVE_FIELD(ir, maxIndicatorLimit, IFT_DirectErrorIndicatorRC_maxlim, "maxlim"); // Macro
    IR_GIVE_FIELD(ir, minIndicatorDensity, IFT_DirectErrorIndicatorRC_mindens, "mindens"); // Macro
    IR_GIVE_FIELD(ir, maxIndicatorDensity, IFT_DirectErrorIndicatorRC_maxdens, "maxdens"); // Macro
    IR_GIVE_FIELD(ir, zeroIndicatorDensity, IFT_DirectErrorIndicatorRC_defdens, "defdens"); // Macro

    remeshingDensityRatioToggle = 0.80;
    IR_GIVE_OPTIONAL_FIELD(ir, remeshingDensityRatioToggle, IFT_DirectErrorIndicatorRC_remeshingdensityratio, "remeshingdensityratio"); // Macro

#ifdef __PARALLEL_MODE
    EngngModel *emodel = domain->giveEngngModel();
    ProblemCommunicatorMode commMode = emodel->giveProblemCommMode();
    if ( commMode == ProblemCommMode__NODE_CUT ) {
        commBuff = new CommunicatorBuff(emodel->giveNumberOfProcesses(), CBT_dynamic);
        communicator = new ProblemCommunicator(emodel, commBuff, emodel->giveRank(),
                                               emodel->giveNumberOfProcesses(),
                                               commMode);
    }

#endif
    return IRRT_OK;
}

RemeshingStrategy
DirectErrorIndicatorRC :: giveRemeshingStrategy(TimeStep *tStep)
{
    this->estimateMeshDensities(tStep);
    return this->currStrategy;
}

void
DirectErrorIndicatorRC :: reinitialize()
{
    stateCounter = -1;
#ifdef __PARALLEL_MODE
    dofManDensityExchangeFlag  = true;
    initCommMap = true;
#endif
}


void
DirectErrorIndicatorRC :: setDomain(Domain *d)
{
    RemeshingCriteria :: setDomain(d);
#ifdef __PARALLEL_MODE
    dofManDensityExchangeFlag = true;
    initCommMap = true;
#endif
}


#ifdef __PARALLEL_MODE
void
DirectErrorIndicatorRC :: exchangeDofManDensities()
{
    Domain *domain = this->giveDomain();
    EngngModel *emodel = domain->giveEngngModel();
    ProblemCommunicatorMode commMode = emodel->giveProblemCommMode();

    if ( this->dofManDensityExchangeFlag ) {
        if ( commMode == ProblemCommMode__NODE_CUT ) {
            ///@todo Compute local shared dofman densities
            sharedDofManDensities.clear();
            int i, size = domain->giveNumberOfDofManagers();
            for ( i = 1; i <= size; i++ ) {
                if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
                    sharedDofManDensities [ i ] = this->giveLocalDofManDensity(i);
                }
            }

            // exchange them
            communicator->packAllData(this, & DirectErrorIndicatorRC :: packSharedDofManLocalDensities);
            communicator->initExchange(999);
            communicator->unpackAllData(this, & DirectErrorIndicatorRC :: unpackSharedDofManLocalDensities);
            communicator->finishExchange();
        }

        this->dofManDensityExchangeFlag = false;
    } // if (this->dofManDensityExchangeFlag)

}

int
DirectErrorIndicatorRC :: packSharedDofManLocalDensities(ProcessCommunicator &processComm)
{
    int result = 1, i, size;
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    IntArray const *toSendMap = processComm.giveToSendMap();

    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        result &= pcbuff->packDouble(this->sharedDofManDensities [ toSendMap->at(i) ]);
    }

    return result;
}

int
DirectErrorIndicatorRC :: unpackSharedDofManLocalDensities(ProcessCommunicator &processComm)
{
    int result = 1;
    int i, size;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    double value;

    size = toRecvMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        result &= pcbuff->unpackDouble(value);
        this->sharedDofManDensities [ toRecvMap->at(i) ] = max(value, this->sharedDofManDensities [ toRecvMap->at(i) ]);
 #ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_INFO("unpackSharedDofManLocalDensities: node %d[%d], value %f\n", toRecvMap->at(i), domain->giveDofManager( toRecvMap->at(i) )->giveGlobalNumber(), this->sharedDofManDensities [ toRecvMap->at(i) ]);
 #endif
    }

    return result;
}




void
DirectErrorIndicatorRC :: exchangeDofManIndicatorVals(TimeStep *tStep)
{
    Domain *domain = this->giveDomain();
    EngngModel *emodel = domain->giveEngngModel();
    ProblemCommunicatorMode commMode = emodel->giveProblemCommMode();

    if ( commMode == ProblemCommMode__NODE_CUT ) {
        ///@todo Compute local shared dofman indicator values
        sharedDofManIndicatorVals.clear();
        int i, size = domain->giveNumberOfDofManagers();
        for ( i = 1; i <= size; i++ ) {
            if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
                sharedDofManIndicatorVals [ i ] = this->giveLocalDofManIndicator(i, tStep);
            }
        }

        // exchange them
        communicator->packAllData(this, & DirectErrorIndicatorRC :: packSharedDofManLocalIndicatorVals);
        communicator->initExchange(999);
        communicator->unpackAllData(this, & DirectErrorIndicatorRC :: unpackSharedDofManLocalIndicatorVals);
        communicator->finishExchange();
    }
}

int
DirectErrorIndicatorRC :: packSharedDofManLocalIndicatorVals(ProcessCommunicator &processComm)
{
    int result = 1, i, size;
    //Domain *d = this->giveDomain();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    IntArray const *toSendMap = processComm.giveToSendMap();

    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        result &= pcbuff->packDouble(this->sharedDofManIndicatorVals [ toSendMap->at(i) ]);
    }

    return result;
}

int
DirectErrorIndicatorRC :: unpackSharedDofManLocalIndicatorVals(ProcessCommunicator &processComm)
{
    int result = 1;
    int i, size;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    double value;

    size = toRecvMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        result &= pcbuff->unpackDouble(value);
        this->sharedDofManIndicatorVals [ toRecvMap->at(i) ] = max(value, this->sharedDofManIndicatorVals [ toRecvMap->at(i) ]);
 #ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_INFO("unpackSharedDofManLocalIndicatorVals: node %d[%d], value %f\n", toRecvMap->at(i), domain->giveDofManager( toRecvMap->at(i) )->giveGlobalNumber(), this->sharedDofManIndicatorVals [ toRecvMap->at(i) ]);
 #endif
    }

    return result;
}


#endif
} // end namespace oofem
