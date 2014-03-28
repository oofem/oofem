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

#include "structengngmodel.h"
#include "dofmanager.h"
#include "dof.h"
#include "element.h"
#include "timestep.h"
#include "outputmanager.h"
#include "structuralelement.h"
#include "structuralelementevaluator.h"
#include "activebc.h"
#include "structuralinterfaceelement.h"
namespace oofem {
StructuralEngngModel :: StructuralEngngModel(int i, EngngModel *_master) : EngngModel(i, _master),
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
    FloatArray contribution;

    answer.resize( this->giveNumberOfDomainEquations( di, EModelDefaultPrescribedEquationNumbering() ) );
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
    this->updateSharedDofManagers(answer, EModelDefaultPrescribedEquationNumbering(), ReactionExchangeTag);
#endif
}


void
StructuralEngngModel :: computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    reactions.resize( this->giveNumberOfDomainEquations( di, EModelDefaultPrescribedEquationNumbering() ) );
    reactions.zero();
    this->assembleVector( reactions, tStep, EID_MomentumBalance, ExternalForcesVector, VM_Total,
                         EModelDefaultPrescribedEquationNumbering(), this->giveDomain(di) );
}


void
StructuralEngngModel :: giveInternalForces(FloatArray &answer, bool normFlag, int di, TimeStep *tStep)
{
    // Simply assembles contributions from each element in domain
    Domain *domain = this->giveDomain(di);
    // Update solution state counter
    tStep->incrementStateCounter();

    answer.resize( this->giveNumberOfDomainEquations( di, EModelDefaultEquationNumbering() ) );
    answer.zero();
    this->assembleVector(answer, tStep, EID_MomentumBalance, InternalForcesVector, VM_Total,
                         EModelDefaultEquationNumbering(), domain, normFlag ? & this->internalForcesEBENorm : NULL);

#ifdef __PARALLEL_MODE
    // Redistributes answer so that every process have the full values on all shared equations
    this->updateSharedDofManagers(answer, EModelDefaultEquationNumbering(), InternalForcesExchangeTag);
#endif

    // Remember last internal vars update time stamp.
    internalVarUpdateStamp = tStep->giveSolutionStateCounter();
}


void
StructuralEngngModel :: updateYourself(TimeStep *tStep)
{
    this->updateInternalState(tStep);
    EngngModel :: updateYourself(tStep);
}


int
StructuralEngngModel :: checkConsistency()
{
    Domain *domain = this->giveDomain(1);
    int nelem = domain->giveNumberOfElements();
    // check for proper element type

    for ( int i = 1; i <= nelem; i++ ) {
        BaseElement *ePtr = domain->giveElement(i);
        StructuralElement *sePtr = dynamic_cast< StructuralElement * >(ePtr);
        StructuralInterfaceElement *siePtr = dynamic_cast< StructuralInterfaceElement * >(ePtr);
        StructuralElementEvaluator *see = dynamic_cast< StructuralElementEvaluator * >(ePtr);

        if ( sePtr == NULL && see == NULL && siePtr == NULL ) {
            OOFEM_WARNING("checkConsistency: element %d has no Structural support", i);
            return 0;
        }
    }

    EngngModel :: checkConsistency();

    return 1;
}


void
StructuralEngngModel :: updateInternalState(TimeStep *tStep)
{
    int nnodes;
    Domain *domain;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);

        nnodes = domain->giveNumberOfDofManagers();
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( int j = 1; j <= nnodes; j++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(j), tStep);
            }
        }

        int nbc = domain->giveNumberOfBoundaryConditions();
        for ( int i = 1; i <= nbc; ++i ) {
            GeneralBoundaryCondition *bc = domain->giveBc(i);
            ActiveBoundaryCondition *abc;

            if ( ( abc = dynamic_cast< ActiveBoundaryCondition * >(bc) ) ) {
                int ndman = abc->giveNumberOfInternalDofManagers();
                for ( int j = 1; j <= ndman; j++ ) {
                    this->updateDofUnknownsDictionary(abc->giveInternalDofManager(j), tStep);
                }
            }
        }

        if ( internalVarUpdateStamp != tStep->giveSolutionStateCounter() ) {
            int nelem = domain->giveNumberOfElements();
            for ( int j = 1; j <= nelem; j++ ) {
                domain->giveElementGeometry(j)->updateInternalState(tStep);
            }

            internalVarUpdateStamp = tStep->giveSolutionStateCounter();
        }
    }
}


void
StructuralEngngModel :: buildReactionTable(IntArray &restrDofMans, IntArray &restrDofs,
                                           IntArray &eqn, TimeStep *tStep, int di)
{
    // determine number of restrained dofs
    Domain *domain = this->giveDomain(di);
    int numRestrDofs = this->giveNumberOfDomainEquations( di, EModelDefaultPrescribedEquationNumbering() );
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
    restrDofMans.resizeWithValues(count);
    restrDofs.resizeWithValues(count);
    eqn.resizeWithValues(count);
}


#ifdef __PARALLEL_MODE
void
StructuralEngngModel :: initParallelContexts()
{
    ParallelContext *parallelContext;

    parallelContextList->growTo(ndomains);
    for ( int i = 0; i < this->ndomains; i++ ) {
        parallelContext =  new ParallelContext(this);
        parallelContextList->put(i + 1, parallelContext);
    }
}
#endif


#ifdef __OOFEG
void
StructuralEngngModel :: showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);

    if ( type != 1 ) {
        return;
    }

    int nelems = domain->giveNumberOfElements();
    for ( int i = 1; i <= nelems; i++ ) {
        domain->giveElement(i)->showSparseMtrxStructure(StiffnessMatrix, context, tStep);
    }
}
#endif
} // end namespace oofem
