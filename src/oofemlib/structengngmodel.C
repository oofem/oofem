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

#include "structengngmodel.h"
#include "dofmanager.h"
#include "dof.h"
#include "element.h"
#include "timestep.h"
#include "outputmanager.h"

namespace oofem {
StructuralEngngModel :: ~StructuralEngngModel()
{
#ifdef __PARALLEL_MODE
    delete communicator;
    delete nonlocCommunicator;
    delete commBuff;
#endif
}


void
StructuralEngngModel :: printReactionForces(TimeStep *tStep, int di)
//
// computes and prints reaction forces in all supported or restrained dofs
//
{
    //
    // Sum all equivalent forces for all connected elements
    //
    int i, numRestrDofs = 0;

    IntArray ielemDofMask;
    FloatArray reactions;
    IntArray dofManMap, dofMap, eqnMap;

    Domain *domain = this->giveDomain(di);

    // test if solution step output is active
    if ( !domain->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return;
    }

    FILE *outputStream = this->giveOutputStream();
    // map contains corresponding dofmanager and dofs numbers corresponding to prescribed equations
    // sorted according to dofmanger number and as a minor crit. according to dof number
    // this is necessary for extractor, since the sorted output is expected
    this->buildReactionTable(dofManMap, dofMap, eqnMap, tStep, di);

    //
    // print header
    //
    fprintf(outputStream, "\n\n\tR E A C T I O N S  O U T P U T:\n\t_______________________________\n\n\n");

    numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance);
    reactions.resize(numRestrDofs);
    // compute reaction forces
    this->computeReaction(reactions, tStep, di);
    //
    // loop over reactions and print them
    //
    for ( i = 1; i <= numRestrDofs; i++ ) {
        if ( domain->giveOutputManager()->testDofManOutput(dofManMap.at(i), tStep) ) {
            //fprintf(outputStream,"\tNode %8d iDof %2d reaction % .4e\n",dofManMap.at(i), dofMap.at(i),reactions.at(eqnMap.at(i)));
#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
            fprintf( outputStream, "\tNode %8d iDof %2d reaction % .4e    [bc-id: %d]\n",
                    domain->giveDofManager( dofManMap.at(i) )->giveLabel(),
                    dofMap.at(i), reactions.at( eqnMap.at(i) ),
                    domain->giveDofManager( dofManMap.at(i) )->giveDof( dofMap.at(i) )->giveBcId() );
#else
            fprintf( outputStream, "\tNode %8d iDof %2d reaction % .4e    [bc-id: %d]\n",
                    dofManMap.at(i), dofMap.at(i), reactions.at( eqnMap.at(i) ),
                    domain->giveDofManager( dofManMap.at(i) )->giveDof( dofMap.at(i) )->giveBcId() );
#endif
        }
    }
}


void
StructuralEngngModel :: computeReaction(FloatArray &answer, TimeStep *tStep, int di)
{
    int numRestrDofs = 0;
    FloatArray contribution;
    FloatArray EquivForces;

    numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance);
    answer.resize(numRestrDofs);
    answer.zero();

    // Internal forces contribution
    this->computeInternalForceReactionContribution(contribution, tStep, di);
    answer.add(contribution);
    // External loading contribution
    this->computeExternalLoadReactionContribution(contribution, tStep, di);
    answer.subtract(contribution);

#ifdef __PARALLEL_MODE
    this->updateSharedPrescribedDofManagers( answer, ReactionExchangeTag  );
#endif
}


void
StructuralEngngModel :: computeInternalForceReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    reactions.resize( this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance) );
    reactions.zero();
    this->assembleVector( reactions, tStep, EID_MomentumBalance, LastEquilibratedInternalForcesVector, VM_Total,
                          EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di) );
}


void
StructuralEngngModel :: computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    int numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance);
    reactions.resize(numRestrDofs);
    reactions.zero();
    FloatArray contribution(numRestrDofs);

    reactions.resize( this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance) );
    reactions.zero();
    this->assembleVector( reactions, tStep, EID_MomentumBalance, ExternalForcesVector, VM_Total,
                          EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di) );
}


void
StructuralEngngModel :: giveInternalForces(FloatArray &answer, bool normFlag, int di, TimeStep *stepN)
{
    // Simply assembles contributions from each element in domain
    Element *element;
    IntArray loc;
    FloatArray charVec;
    FloatMatrix R;
    int nelems;
    EModelDefaultEquationNumbering dn;

    answer.resize( this->giveNumberOfEquations( EID_MomentumBalance ) );
    answer.zero();
    if ( normFlag ) internalForcesEBENorm = 0.0;

    // Update solution state counter
    stepN->incrementStateCounter();

#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        exchangeRemoteElementData( RemoteElementExchangeTag  );
    }
#endif

    nelems = this->giveDomain(di)->giveNumberOfElements();
    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

    for ( int i = 1; i <= nelems; i++ ) {
        element = this->giveDomain(di)->giveElement(i);
#ifdef __PARALLEL_MODE
        // Skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. Their introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) continue;

#endif
        element->giveLocationArray( loc, EID_MomentumBalance, dn );
        element->giveCharacteristicVector( charVec, InternalForcesVector, VM_Total, stepN );
        if ( element->giveRotationMatrix(R, EID_MomentumBalance) ) {
            charVec.rotatedWith(R, 't');
        }
        answer.assemble(charVec, loc);

        // Compute element norm contribution.
        if ( normFlag ) internalForcesEBENorm += charVec.computeSquaredNorm();
    }

    this->timer.pauseTimer( EngngModelTimer :: EMTT_NetComputationalStepTimer );

#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers(answer, InternalForcesExchangeTag);

    if ( normFlag ) {
        // Exchange norm contributions.
        double localEBENorm = internalForcesEBENorm;
        internalForcesEBENorm = 0.0;
        MPI_Allreduce(& localEBENorm, & internalForcesEBENorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif

    // Remember last internal vars update time stamp.
    internalVarUpdateStamp = stepN->giveSolutionStateCounter();
}

void
StructuralEngngModel :: updateInternalState(TimeStep *stepN)
{
    int j, nnodes;
    Domain *domain;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);

        nnodes = domain->giveNumberOfDofManagers();
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( j = 1; j <= nnodes; j++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(j), stepN);
            }
        }


        if ( internalVarUpdateStamp != stepN->giveSolutionStateCounter() ) {
            int nelem = domain->giveNumberOfElements();
            for ( j = 1; j <= nelem; j++ ) {
                domain->giveElement(j)->updateInternalState(stepN);
            }

            internalVarUpdateStamp = stepN->giveSolutionStateCounter();
        }
    }
}


void
StructuralEngngModel :: buildReactionTable(IntArray &restrDofMans, IntArray &restrDofs,
                                           IntArray &eqn, TimeStep *tStep, int di)
{
    // determine number of restrained dofs
    Domain *domain = this->giveDomain(di);
    int numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance);
    int ndofMan = domain->giveNumberOfDofManagers();
    int i, j, indofs, rindex, count = 0;
    DofManager *inode;
    Dof *jdof;

    // initialize corresponding dofManagers and dofs for each restrained dof
    restrDofMans.resize(numRestrDofs);
    restrDofs.resize(numRestrDofs);
    eqn.resize(numRestrDofs);

    for ( i = 1; i <= ndofMan; i++ ) {
        inode = domain->giveDofManager(i);
        indofs = inode->giveNumberOfDofs();
        for ( j = 1; j <= indofs; j++ ) {
            jdof = inode->giveDof(j);
            if ( ( jdof->giveClassID() != SimpleSlaveDofClass ) && ( jdof->hasBc(tStep) ) ) { // skip slave dofs
                rindex = jdof->__givePrescribedEquationNumber();
                if ( rindex ) {
                    count++;
                    restrDofMans.at(count) = i;
                    restrDofs.at(count) = j;
                    eqn.at(count) = rindex;
                } else {
                    // NullDof has no equation number and no prescribed equation number
                    //_error ("No prescribed equation number assigned to supported DOF");
                }
            }
        }
    }
}


#ifdef __PARALLEL_MODE
int
StructuralEngngModel :: updateSharedDofManagers(FloatArray &answer, int ExchangeTag)
{
    int result = 1;


    if ( isParallel() ) {
#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: updateSharedDofManagers", "Packing data", this->giveRank() );
#endif

        result &= communicator->packAllData( ( StructuralEngngModel * ) this, & answer, & StructuralEngngModel :: packDofManagers );

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: updateSharedDofManagers", "Exchange started", this->giveRank() );
#endif

        result &= communicator->initExchange(ExchangeTag);

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: updateSharedDofManagers", "Receiving and unpacking", this->giveRank() );
#endif

        result &= communicator->unpackAllData( ( StructuralEngngModel * ) this, & answer, & StructuralEngngModel :: unpackDofManagers );
        result &= communicator->finishExchange();
    }

    return result;
}

int
StructuralEngngModel :: updateSharedPrescribedDofManagers(FloatArray &answer, int ExchangeTag)
{
    int result = 1;

    if ( isParallel() ) {
#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: updateSharedPrescribedDofManagers", "Packing data", this->giveRank() );
#endif

        result &= communicator->packAllData( ( StructuralEngngModel * ) this, & answer, & StructuralEngngModel :: packPrescribedDofManagers );

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: updateSharedPrescribedDofManagers", "Exchange started", this->giveRank() );
#endif

        result &= communicator->initExchange(ExchangeTag);

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: updateSharedDofManagers", "Receiving and unpacking", this->giveRank() );
#endif

        result &= communicator->unpackAllData( ( StructuralEngngModel * ) this, & answer, & StructuralEngngModel :: unpackPrescribedDofManagers );
        result &= communicator->finishExchange();
    }

    return result;
}


void
StructuralEngngModel :: initializeCommMaps(bool forceInit)
{
    // Set up communication patterns.
    communicator->setUpCommunicationMaps(this, true, forceInit);
    if ( nonlocalExt ) {
        nonlocCommunicator->setUpCommunicationMaps(this, true, forceInit);
    }
}


int
StructuralEngngModel :: exchangeRemoteElementData(int ExchangeTag)
{
    int result = 1;

    if ( isParallel() && nonlocalExt ) {
 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: exchangeRemoteElementData", "Packing remote element data", this->giveRank() );
 #endif

        result &= nonlocCommunicator->packAllData( ( StructuralEngngModel * ) this, & StructuralEngngModel :: packRemoteElementData );

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: exchangeRemoteElementData", "Remote element data exchange started", this->giveRank() );
 #endif

        result &= nonlocCommunicator->initExchange(ExchangeTag);

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: exchangeRemoteElementData", "Receiveng and Unpacking remote element data", this->giveRank() );
 #endif

        if ( !( result &= nonlocCommunicator->unpackAllData( ( StructuralEngngModel * ) this, & StructuralEngngModel :: unpackRemoteElementData ) ) ) {
            _error("StructuralEngngModel :: exchangeRemoteElementData: Receiveng and Unpacking remote element data");
        }

        result &= nonlocCommunicator->finishExchange();
    }

    return result;
}


int
StructuralEngngModel :: packRemoteElementData(ProcessCommunicator &processComm)
{
    int result = 1;
    int i, size;
    IntArray const *toSendMap = processComm.giveToSendMap();
    CommunicationBuffer *send_buff = processComm.giveProcessCommunicatorBuff()->giveSendBuff();
    Domain *domain = this->giveDomain(1);


    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        result &= domain->giveElement( toSendMap->at(i) )->packUnknowns( * send_buff, this->giveCurrentStep() );
    }

    return result;
}


int
StructuralEngngModel :: unpackRemoteElementData(ProcessCommunicator &processComm)
{
    int result = 1;
    int i, size;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    CommunicationBuffer *recv_buff = processComm.giveProcessCommunicatorBuff()->giveRecvBuff();
    Element *element;
    Domain *domain = this->giveDomain(1);


    size = toRecvMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        element = domain->giveElement( toRecvMap->at(i) );
        if ( element->giveParallelMode() == Element_remote ) {
            result &= element->unpackAndUpdateUnknowns( * recv_buff, this->giveCurrentStep() );
        } else {
            _error("unpackRemoteElementData: element is not remote");
        }
    }

    return result;
}


int
StructuralEngngModel :: packDofManagers(FloatArray *src, ProcessCommunicator &processComm)
{
    int result = 1;
    result &= this->packDofManagers(src, processComm, false);
    return result;
}


int
StructuralEngngModel :: packPrescribedDofManagers(FloatArray *src, ProcessCommunicator &processComm)
{
    int result = 1;
    result &= this->packDofManagers(src, processComm, true);
    return result;
}


int
StructuralEngngModel :: packDofManagers(FloatArray *src, ProcessCommunicator &processComm, bool prescribedEquations)
{
    int result = 1;
    int i, size;
    int j, ndofs, eqNum;
    Domain *domain = this->giveDomain(1);
    IntArray const *toSendMap = processComm.giveToSendMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    DofManager *dman;
    Dof *jdof;

    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        dman = domain->giveDofManager( toSendMap->at(i) );
        ndofs = dman->giveNumberOfDofs();
        for ( j = 1; j <= ndofs; j++ ) {
            jdof = dman->giveDof(j);
            if ( prescribedEquations ) {
                eqNum = jdof->__givePrescribedEquationNumber();
            } else {
                eqNum = jdof->__giveEquationNumber();
            }
            if ( jdof->isPrimaryDof() && eqNum ) {
                    result &= pcbuff->packDouble( src->at(eqNum) );
            }
        }
    }

    return result;
}


int
StructuralEngngModel :: unpackDofManagers(FloatArray *src, ProcessCommunicator &processComm)
{
    int result = 1;
    result &= this->unpackDofManagers(src, processComm, false);
    return result;
}


int
StructuralEngngModel :: unpackPrescribedDofManagers(FloatArray *src, ProcessCommunicator &processComm)
{
    int result = 1;
    result &= this->unpackDofManagers(src, processComm, true);
    return result;
}


int
StructuralEngngModel :: unpackDofManagers(FloatArray *dest, ProcessCommunicator &processComm, bool prescribedEquations)
{
    int result = 1;
    int i, size;
    int j, ndofs, eqNum;
    Domain *domain = this->giveDomain(1);
    dofManagerParallelMode dofmanmode;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    DofManager *dman;
    Dof *jdof;
    double value;


    size = toRecvMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        dman = domain->giveDofManager( toRecvMap->at(i) );
        ndofs = dman->giveNumberOfDofs();
        dofmanmode = dman->giveParallelMode();
        for ( j = 1; j <= ndofs; j++ ) {
            jdof = dman->giveDof(j);
            if ( prescribedEquations ) {
                eqNum = jdof->__givePrescribedEquationNumber();
            } else {
                eqNum = jdof->__giveEquationNumber();
            }
            if ( jdof->isPrimaryDof() && eqNum ) {
                result &= pcbuff->unpackDouble(value);
                if ( dofmanmode == DofManager_shared ) {
                    dest->at(eqNum) += value;
                } else if ( dofmanmode == DofManager_remote ) {
                    dest->at(eqNum)  = value;
                } else {
                    _error("unpackReactions: unknown dof namager parallel mode");
                }
            }
        }
    }

    return result;
}
#endif


#ifdef __PETSC_MODULE
void
StructuralEngngModel :: initPetscContexts()
{
    PetscContext *petscContext;

    int i;
    petscContextList->growTo(ndomains);
    for ( i = 0; i < this->ndomains; i++ ) {
        petscContext =  new PetscContext(this, EID_MomentumBalance);
        petscContextList->put(i + 1, petscContext);
    }
}
#endif


#ifdef __OOFEG
void
StructuralEngngModel :: showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime)
{
    Domain *domain = this->giveDomain(1);
    CharType ctype;
    int i;

    if ( type != 1 ) {
        return;
    }

    ctype = StiffnessMatrix;

    int nelems = domain->giveNumberOfElements();
    for ( i = 1; i <= nelems; i++ ) {
        domain->giveElement(i)->showSparseMtrxStructure(ctype, context, atTime);
    }
}
#endif

} // end namespace oofem
