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
#include "structuralelement.h"
#include "structuralelementevaluator.h"

namespace oofem {

StructuralEngngModel::StructuralEngngModel(int i, EngngModel* _master) : EngngModel(i, _master),
    internalVarUpdateStamp(0), internalForcesEBENorm()
{ }


StructuralEngngModel :: ~StructuralEngngModel()
{ }


void
StructuralEngngModel :: printReactionForces(TimeStep *tStep, int di)
//
// computes and prints reaction forces in all supported or restrained dofs
//
{
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

    // compute reaction forces
    this->computeReaction(reactions, tStep, di);

    //
    // loop over reactions and print them
    //
    for ( int i = 1; i <= dofManMap.giveSize(); i++ ) {
        if ( domain->giveOutputManager()->testDofManOutput(dofManMap.at(i), tStep) ) {
            fprintf( outputStream, "\tNode %8d iDof %2d reaction % .4e    [bc-id: %d]\n",
                    domain->giveDofManager( dofManMap.at(i) )->giveLabel(),
                    dofMap.at(i), reactions.at( eqnMap.at(i) ),
                    domain->giveDofManager( dofManMap.at(i) )->giveDof( dofMap.at(i) )->giveBcId() );
        }
    }
}


void
StructuralEngngModel :: computeReaction(FloatArray &answer, TimeStep *tStep, int di)
{
    int numRestrDofs = 0;
    FloatArray contribution;
    FloatArray EquivForces;

    numRestrDofs = this->giveNumberOfDomainEquations(di, EModelDefaultPrescribedEquationNumbering());
    answer.resize(numRestrDofs);
    answer.zero();

    // Add internal forces
    this->assembleVector( answer, tStep, EID_MomentumBalance, LastEquilibratedInternalForcesVector, VM_Total,
                        EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di) );
    // Subtract external loading
    ///@todo All engineering models should be using this (for consistency)
    //this->assembleVector( answer, tStep, EID_MomentumBalance, ExternalForcesVector, VM_Total,
    //                    EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di) );
    ///@todo This method is overloaded in some functions, it needs to be generalized.
    this->computeExternalLoadReactionContribution(contribution, tStep, di);
    answer.subtract(contribution);

#ifdef __PARALLEL_MODE
    this->updateSharedPrescribedDofManagers( answer, ReactionExchangeTag  );
#endif
}


void
StructuralEngngModel :: computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    int numRestrDofs = this->giveNumberOfDomainEquations(di, EModelDefaultPrescribedEquationNumbering());
    reactions.resize(numRestrDofs);
    reactions.zero();

    reactions.resize( this->giveNumberOfDomainEquations(di, EModelDefaultPrescribedEquationNumbering()) );
    reactions.zero();
    this->assembleVector( reactions, tStep, EID_MomentumBalance, ExternalForcesVector, VM_Total,
                          EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di) );
}


void
StructuralEngngModel :: giveInternalForces(FloatArray &answer, bool normFlag, int di, TimeStep *stepN)
{
    // Simply assembles contributions from each element in domain
    Domain *domain = this->giveDomain(di);
    // Update solution state counter
    stepN->incrementStateCounter();

#ifdef __PARALLEL_MODE
    ///@todo Move this into assembleVector
    if ( this->isParallel() ) {
        // Copies data from remote elements to make sure they have all information necessary for nonlocal averaging.
        exchangeRemoteElementData( RemoteElementExchangeTag  );
    }
#endif

answer.resize( this->giveNumberOfDomainEquations(di, EModelDefaultEquationNumbering()) );
    answer.zero();
    this->assembleVector( answer, stepN, EID_MomentumBalance, InternalForcesVector, VM_Total,
                          EModelDefaultEquationNumbering(), domain, normFlag ? &this->internalForcesEBENorm : NULL );

#ifdef __PARALLEL_MODE
    // Redistributes answer so that every process have the full values on all shared equations
    ///@todo This is basically "scatterN2L" which we should have available for all dofs (not just from dofmanagers). Should be moved into engngmodel as well
    this->updateSharedDofManagers(answer, InternalForcesExchangeTag);
#endif

    // Remember last internal vars update time stamp.
    internalVarUpdateStamp = stepN->giveSolutionStateCounter();
}


void
StructuralEngngModel :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    EngngModel :: updateYourself(stepN);
}


int
StructuralEngngModel :: checkConsistency()
{
    Domain *domain = this->giveDomain(1);
    int nelem = domain->giveNumberOfElements();
    // check for proper element type

    for ( int i = 1; i <= nelem; i++ ) {
        Element *ePtr = domain->giveElement(i);
        StructuralElement *sePtr = dynamic_cast< StructuralElement * >(ePtr);
        StructuralElementEvaluator *see = dynamic_cast< StructuralElementEvaluator * >(ePtr);

        if ( sePtr == NULL && see == NULL ) {
            _warning2("checkConsistency: element %d has no Structural support", i);
            return 0;
        }
    }

    EngngModel :: checkConsistency();

    return 1;
}


void
StructuralEngngModel :: updateInternalState(TimeStep *stepN)
{
    int nnodes;
    Domain *domain;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);

        nnodes = domain->giveNumberOfDofManagers();
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( int j = 1; j <= nnodes; j++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(j), stepN);
            }
        }

        if ( internalVarUpdateStamp != stepN->giveSolutionStateCounter() ) {
            int nelem = domain->giveNumberOfElements();
            for ( int j = 1; j <= nelem; j++ ) {
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
    int numRestrDofs = this->giveNumberOfDomainEquations(di, EModelDefaultPrescribedEquationNumbering());
    int ndofMan = domain->giveNumberOfDofManagers();
    int rindex, count = 0;

    // initialize corresponding dofManagers and dofs for each restrained dof
    restrDofMans.resize(numRestrDofs);
    restrDofs.resize(numRestrDofs);
    eqn.resize(numRestrDofs);

    for ( int i = 1; i <= ndofMan; i++ ) {
        DofManager *inode = domain->giveDofManager(i);
        int indofs = inode->giveNumberOfDofs();
        for ( int j = 1; j <= indofs; j++ ) {
            Dof *jdof = inode->giveDof(j);
            if ( jdof->isPrimaryDof() && ( jdof->hasBc(tStep) ) ) { // skip slave dofs
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
    // Trim to size.
    restrDofMans.resize(count);
    restrDofs.resize(count);
    eqn.resize(count);
}


#ifdef __PETSC_MODULE
void
StructuralEngngModel :: initPetscContexts()
{
    PetscContext *petscContext;

    petscContextList->growTo(ndomains);
    for ( int i = 0; i < this->ndomains; i++ ) {
        petscContext =  new PetscContext(this);
        petscContextList->put(i + 1, petscContext);
    }
}
#endif


#ifdef __OOFEG
void
StructuralEngngModel :: showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime)
{
    Domain *domain = this->giveDomain(1);

    if ( type != 1 ) {
        return;
    }

    int nelems = domain->giveNumberOfElements();
    for ( int i = 1; i <= nelems; i++ ) {
        domain->giveElement(i)->showSparseMtrxStructure(StiffnessMatrix, context, atTime);
    }
}
#endif

} // end namespace oofem
