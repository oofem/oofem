/* $Header: /home/cvs/bp/oofem/oofemlib/src/structengngmodel.C,v 1.15.4.1 2004/04/05 15:19:44 bp Exp $ */
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


#include "structengngmodel.h"
#include "dofmanager.h"
#include "elementside.h"
#include "rigidarmnode.h"
#include "dof.h"
#include "slavedof.h"
#include "structuralelement.h"
#include "timestep.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

#include "verbose.h"
#include "conTable.h"
#include "outputmanager.h"

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
    FloatArray reactions, contribution;
    FloatArray EquivForces;
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

    //numRestrDofs = domain -> giveConnectivityTable () -> giveNumRestrDofs();
    numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance);
    reactions.resize(numRestrDofs);
    reactions.zero();

    // Internal forces contribution

    this->computeInternalForceReactionContribution(contribution, tStep, di);
    reactions.add(contribution);
    // External loading contribution
    this->computeExternalLoadReactionContribution(contribution, tStep, di);
    reactions.substract(contribution);

#ifdef __PARALLEL_MODE

    if ( commMode == ProblemCommMode__NODE_CUT ) {
#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: printReactionForces", "Packing reactions", this->giveRank() );
#endif

        communicator->packAllData( ( StructuralEngngModel * ) this, & reactions, & StructuralEngngModel :: packReactions );

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: printReactionForces", "Exchange of reactions started", this->giveRank() );
#endif

        communicator->initExchange(999);

#ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "StructuralEngngModel :: printReactionForces", "Receiving and unpacking of reactions started", this->giveRank() );
#endif

        communicator->unpackAllData( ( StructuralEngngModel * ) this, & reactions, & StructuralEngngModel :: unpackReactions );
        communicator->finishExchange();
    }

#endif

    //
    // loop over reactions and print them
    //
    for ( i = 1; i <= numRestrDofs; i++ ) {
        if ( domain->giveOutputManager()->testDofManOutput(dofManMap.at(i), tStep) ) {
            //fprintf(outputStream,"\tNode %8d iDof %2d reaction % .4e\n",dofManMap.at(i), dofMap.at(i),reactions.at(eqnMap.at(i)));
#if defined(__PARALLEL_MODE) || defined(__ENABLE_COMPONENT_LABELS)
	  fprintf( outputStream, "\tNode %8d iDof %2d reaction % .4e    [bc-id: %d]\n", 
		   domain->giveDofManager( dofManMap.at(i))->giveLabel(), 
		   dofMap.at(i), reactions.at( eqnMap.at(i) ), 
		   domain->giveDofManager( dofManMap.at(i) )->giveDof( dofMap.at(i) )->giveBcId() );
#else
            fprintf( outputStream, "\tNode %8d iDof %2d reaction % .4e    [bc-id: %d]\n", 
		     dofManMap.at(i), dofMap.at(i), reactions.at( eqnMap.at(i) ), 
		     domain->giveDofManager( dofManMap.at(i) )->giveDof( dofMap.at(i) )->giveBcId() );
#endif
        }
    }

    return;
}


#if 1

void
StructuralEngngModel :: computeInternalForceReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    reactions.resize( this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance) );
    reactions.zero();
    this->assemblePrescribedVectorFromElements( reactions, tStep, EID_MomentumBalance,
                                               LastEquilibratedNodalInternalForcesVector,
                                               VM_Total, this->giveDomain(di) );
}

void
StructuralEngngModel :: computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    int numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance);
    reactions.resize(numRestrDofs);
    reactions.zero();
    FloatArray contribution(numRestrDofs);

    this->computeElementLoadReactionContribution(contribution, tStep, di);
    reactions.add(contribution);
    this->computeNodalLoadReactionContribution(contribution, tStep, di);
    reactions.add(contribution);
}


void
StructuralEngngModel :: computeElementLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    reactions.resize( this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance) );
    reactions.zero();
    this->assemblePrescribedVectorFromElements( reactions, tStep, EID_MomentumBalance, ElementForceLoadVector, VM_Total, this->giveDomain(di) );
}

void
StructuralEngngModel :: computeNodalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    reactions.resize( this->giveNumberOfPrescribedDomainEquations(di, EID_MomentumBalance) );
    reactions.zero();
    assemblePrescribedVectorFromDofManagers( reactions, tStep, EID_MomentumBalance, NodalLoadVector, VM_Total, this->giveDomain(di) );
}

#else

void
StructuralEngngModel :: computeInternalForceReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    int i, j, k, numRestrDofs = 0;
    int dofNum, elementEquivForcesFlag, knodes, kdofNum, ielemDofManDofs;
    int reactionIndx;
    ValueModeType mode = VM_Total;

    DofManager *ielemDofMan, *master;
    Dof *kDof;
    Element *ielem;

    IntArray ielemDofMask;
    FloatArray EquivForces;

    Domain *domain = this->giveDomain(di);

    //numRestrDofs = domain -> giveConnectivityTable () -> giveNumRestrDofs();
    numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di);

    reactions.resize(numRestrDofs);
    reactions.zero();

    int nelems = domain->giveNumberOfElements();
    for ( i = 1; i <= nelems; i++ ) { // loop over elements
        elementEquivForcesFlag = 0; // indicates that element equiv forces has not been computed
        ielem = domain->giveElement(i);
        knodes = ielem->giveNumberOfDofManagers();

        dofNum = 0;
        for ( j = 1; j <= knodes; j++ ) { // loop over element nodes
            ielem->giveDofManDofIDMask(j, ielemDofMask);
            ielemDofMan = ielem->giveDofManager(j);
            ielemDofManDofs = ielemDofMask.giveSize();
            if ( ielemDofMan->giveClassID() == RigidArmNodeClass ) {
                IntArray slaveDofMask;

                master = ( ( RigidArmNode * ) ielemDofMan )->giveMasterDofMngr();
                ( ( RigidArmNode * ) ielemDofMan )->giveSlaveDofMask(slaveDofMask);

                // rigid arm node contributes to all master dofs
                // the number of dofs at nodal c.s (in which we are computing reactions)
                // is not determined by size of ielemDofMask, but by total number of dofs of master
                // because single local dof can contribute to several DOFS in master.
                for ( k = 1; k <= ielemDofMan->giveNumberOfDofs(); k++ ) { // loop over element dofManager dofs
                    dofNum++;
                    if ( slaveDofMask.at(k) ) {
                        reactionIndx = master->giveDof(k)->givePrescribedEquationNumber();
                    } else {
                        reactionIndx = ielemDofMan->giveDof(k)->givePrescribedEquationNumber();
                    }

                    // select only dof managers required to do output
                    if ( domain->giveOutputManager()->testDofManOutput(master->giveNumber(), tStep) ) {
                        if ( reactionIndx ) { // prescribed
                            // compute element equivForces if not computed previously
                            if ( !elementEquivForcesFlag ) {
                                ( ( StructuralElement * ) ielem )->giveInternalForcesVector(EquivForces, tStep, mode);
                                elementEquivForcesFlag = 1;
                            }

                            reactions.at(reactionIndx) +=  EquivForces.at(dofNum);
                        }
                    }
                }
            } else {
                for ( k = 1; k <= ielemDofManDofs; k++ ) { // loop over element dofManager dofs
                    dofNum++;
                    kdofNum = ielemDofMan->findDofWithDofId( ielemDofMask.at(k) );
                    if ( kdofNum ) {
                        kDof = ielemDofMan->giveDof(kdofNum);
                        /*
                         *    if (kDof->hasBc(tStep)) { // kDof is supported
                         *     // find corresponding master dofManager and dof number
                         *     if (kDof->giveClassID() == SimpleSlaveDofClass) {
                         *      restrDofManNum = ((SlaveDof*) kDof)->giveMasterDofManagerNum();
                         *      restrDofNum    = ((SlaveDof*) kDof)->giveMasterDofIndx();
                         *     } else {
                         *      restrDofManNum = ielemDofMan->giveNumber();
                         *      restrDofNum    = kDof->giveNumber();
                         *     }
                         *     // select only dof managers required to do output
                         *     if (domain->giveOutputManager()->testDofManOutput (ielemDofMan->giveNumber(), tStep)) {
                         *      // find corresponding reaction index
                         *      reactionIndx = domain->giveConnectivityTable() -> giveReactionIndx (restrDofManNum, restrDofNum);
                         */
                        {
                            if ( domain->giveOutputManager()->testDofManOutput(ielemDofMan->giveNumber(), tStep) ) {
                                reactionIndx = kDof->givePrescribedEquationNumber();
                                if ( reactionIndx ) {
                                    // compute element equivForces if not computed previously
                                    if ( !elementEquivForcesFlag ) {
                                        ( ( StructuralElement * ) ielem )->giveInternalForcesVector(EquivForces, tStep, mode);
                                        elementEquivForcesFlag = 1;
                                    }

                                    reactions.at(reactionIndx) +=  EquivForces.at(dofNum);
                                }
                            }
                        }
                    } else {
                        _error("printReactionForces: element dof not found");
                    }
                }
            }
        } // end of loop over element nodes

    } // end of loop over elements

}


void
StructuralEngngModel :: computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    int numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di);
    reactions.resize(numRestrDofs);
    reactions.zero();
    FloatArray contribution(numRestrDofs);

    this->computeElementLoadReactionContribution(contribution, tStep, di);
    reactions.add(contribution);
    this->computeNodalLoadReactionContribution(contribution, tStep, di);
    reactions.add(contribution);
}


void
StructuralEngngModel :: computeElementLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    int i, j, k, numRestrDofs;
    int dofNum, elementEquivForcesFlag, knodes, kdofNum, ielemDofManDofs;
    int reactionIndx;
    ValueModeType mode = VM_Total;

    DofManager *ielemDofMan, *master;
    Dof *kDof;
    Element *ielem;

    IntArray ielemDofMask;
    FloatArray EquivForces;

    Domain *domain = this->giveDomain(di);

    //numRestrDofs = domain -> giveConnectivityTable () -> giveNumRestrDofs();
    numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di);
    reactions.resize(numRestrDofs);
    reactions.zero();

    int nelems = domain->giveNumberOfElements();
    for ( i = 1; i <= nelems; i++ ) { // loop over elements
        elementEquivForcesFlag = 0; // indicates that element equiv forces has not been computed
        ielem = domain->giveElement(i);
        knodes = ielem->giveNumberOfDofManagers();

        dofNum = 0;
        for ( j = 1; j <= knodes; j++ ) { // loop over element nodes
            ielem->giveDofManDofIDMask(j, ielemDofMask);
            ielemDofMan = ielem->giveDofManager(j);
            ielemDofManDofs = ielemDofMask.giveSize();
            if ( ielemDofMan->giveClassID() == RigidArmNodeClass ) {
                IntArray slaveDofMask;

                master = ( ( RigidArmNode * ) ielemDofMan )->giveMasterDofMngr();
                ( ( RigidArmNode * ) ielemDofMan )->giveSlaveDofMask(slaveDofMask);

                // rigid arm node contributes to all master dofs
                // the number of dofs at nodal c.s (in which we are computing reactions)
                // is not determined by size of ielemDofMask, but by total number of dofs of master
                // because single local dof can contribute to several DOFS in master.
                for ( k = 1; k <= ielemDofMan->giveNumberOfDofs(); k++ ) { // loop over element dofManager dofs
                    dofNum++;
                    if ( slaveDofMask.at(k) ) {
                        reactionIndx = master->giveDof(k)->givePrescribedEquationNumber();
                    } else {
                        reactionIndx = ielemDofMan->giveDof(k)->givePrescribedEquationNumber();
                    }

                    // select only dof managers required to do output
                    if ( domain->giveOutputManager()->testDofManOutput(master->giveNumber(), tStep) ) {
                        if ( domain->giveOutputManager()->testDofManOutput(master->giveNumber(), tStep) ) {
                            reactionIndx = master->giveDof(k)->givePrescribedEquationNumber();
                            if ( reactionIndx ) {
                                // compute element equivForces if not computed previously
                                if ( !elementEquivForcesFlag ) {
                                    ( ( StructuralElement * ) ielem )->computeForceLoadVector(EquivForces, tStep, mode);
                                    elementEquivForcesFlag = 1;
                                }

                                if ( EquivForces.giveSize() ) {
                                    reactions.at(reactionIndx) +=  EquivForces.at(dofNum);
                                }
                            }
                        }
                    }
                }
            } else {
                for ( k = 1; k <= ielemDofManDofs; k++ ) { // loop over element dofManager dofs
                    dofNum++;
                    kdofNum = ielemDofMan->findDofWithDofId( ielemDofMask.at(k) );
                    if ( kdofNum ) {
                        kDof = ielemDofMan->giveDof(kdofNum);
                        /*
                         * if (kDof->hasBc(tStep)) { // kDof is supported
                         * // find corresponding master dofManager and dof number
                         * if (kDof->giveClassID() == SimpleSlaveDofClass) {
                         * restrDofManNum = ((SlaveDof*) kDof)->giveMasterDofManagerNum();
                         * restrDofNum    = ((SlaveDof*) kDof)->giveMasterDofIndx();
                         * } else {
                         * restrDofManNum = ielemDofMan->giveNumber();
                         * restrDofNum    = kDof->giveNumber();
                         * }
                         * // select only dof managers required to do output
                         * if (domain->giveOutputManager()->testDofManOutput (ielemDofMan->giveNumber(), tStep)) {
                         * // find corresponding reaction index
                         * reactionIndx = domain->giveConnectivityTable() -> giveReactionIndx (restrDofManNum, restrDofNum);
                         */
                        {
                            if ( domain->giveOutputManager()->testDofManOutput(ielemDofMan->giveNumber(), tStep) ) {
                                // find corresponding reaction index
                                reactionIndx = kDof->givePrescribedEquationNumber();
                                if ( reactionIndx ) {
                                    // compute element equivForces if not computed previously
                                    if ( !elementEquivForcesFlag ) {
                                        ( ( StructuralElement * ) ielem )->computeForceLoadVector(EquivForces, tStep, mode);
                                        elementEquivForcesFlag = 1;
                                    }

                                    if ( EquivForces.giveSize() ) {
                                        reactions.at(reactionIndx) +=  EquivForces.at(dofNum);
                                    }
                                }
                            }
                        }
                    } else {
                        _error("printReactionForces: element dof not found");
                    }
                }
            }
        } // end of loop over element nodes

    } // end of loop over elements

}


void
StructuralEngngModel :: computeNodalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di)
{
    int i, j, numRestrDofs;
    int reactionIndx;
    int nDofMan;
    ValueModeType mode = VM_Total;

    DofManager *master;

    IntArray ielemDofMask;
    FloatArray EquivForces;
    FloatArray nodalLoadVector, nonnodalLoadVector;

    Domain *domain = this->giveDomain(di);

    //numRestrDofs = domain -> giveConnectivityTable () -> giveNumRestrDofs();
    numRestrDofs = this->giveNumberOfPrescribedDomainEquations(di);
    reactions.resize(numRestrDofs);
    reactions.zero();


    // loop over dofManagers and find if they are directly loaded
    DofManager *iDofMan;
    FloatArray dofManLoad;
    int hasBc; // dofManNum;
    nDofMan = domain->giveNumberOfDofManagers();
    for ( i = 1; i <= nDofMan; i++ ) {
        iDofMan = domain->giveDofManager(i);
        if ( iDofMan->giveClassID() == RigidArmNodeClass ) {
            IntArray slaveDofMask;
            master = ( ( RigidArmNode * ) iDofMan )->giveMasterDofMngr();
            ( ( RigidArmNode * ) iDofMan )->giveSlaveDofMask(slaveDofMask);
            // determine if any BC
            // loop over linked DOF at master
            for ( hasBc = 0, j = 1; j <= master->giveNumberOfDofs(); j++ ) {
                if ( slaveDofMask.at(j) ) {
                    hasBc += master->giveDof(j)->hasBc(tStep);
                }
            }

            // loop over primary dof if any
            for ( j = 1; j <= iDofMan->giveNumberOfDofs(); j++ ) {
                if ( slaveDofMask.at(j) == 0 ) {
                    hasBc += iDofMan->giveDof(j)->hasBc(tStep);
                }
            }

            if ( hasBc ) {
                // rigid Arm node and master has bc
                iDofMan->computeLoadVectorAt(dofManLoad, tStep, mode);
                if ( dofManLoad.isEmpty() ) {
                    continue;
                }

                // loop over master dofs to find those with bc
                for ( j = 1; j <= iDofMan->giveNumberOfDofs(); j++ ) {
                    if ( slaveDofMask.at(j) ) { // slave linked
                        if ( master->giveDof(j)->hasBc(tStep) ) {
                            // reactionIndx = domain->giveConnectivityTable() -> giveReactionIndx (master->giveNumber(),
                            //                                       master->giveDof(j)->giveNumber());
                            reactionIndx = master->giveDof(j)->givePrescribedEquationNumber();
                        } else {
                            continue;
                        }
                    } else { // primary dof
                        if ( iDofMan->giveDof(j)->hasBc(tStep) ) {
                            //       reactionIndx = domain->giveConnectivityTable() -> giveReactionIndx (dofManNum,
                            //                                                                           iDofMan->giveDof(j)->giveNumber());
                            reactionIndx = iDofMan->giveDof(j)->givePrescribedEquationNumber();
                        } else {
                            continue;
                        }
                    }

                    if ( reactionIndx ) {
                        reactions.at(reactionIndx) = dofManLoad.at(j);
                    }
                }
            }
        } else { // end rigid arm node
            for ( hasBc = 0, j = 1; j <= iDofMan->giveNumberOfDofs(); j++ ) {
                hasBc += iDofMan->giveDof(j)->hasBc(tStep);
            }

            if ( hasBc ) {
                iDofMan->computeLoadVectorAt(dofManLoad, tStep, mode);
                if ( dofManLoad.isEmpty() ) {
                    continue;
                }

                // loop over  dofs to find those with bc
                for ( j = 1; j <= iDofMan->giveNumberOfDofs(); j++ ) {
                    if ( iDofMan->giveDof(j)->hasBc(tStep) ) {
                        /*
                         *    if (iDofMan->giveDof(j)->isPrimaryDof()) {
                         *     dofManNum = iDofMan->giveNumber();
                         *    } else if (iDofMan->giveDof(j)->giveClassID() == SimpleSlaveDofClass) {
                         *     dofManNum = ((SlaveDof*)iDofMan->giveDof(j))->giveMasterDofManagerNum();
                         *    } else {
                         *     _error ("printReactionForces: unknown non-primary dof type");
                         *    }
                         */
                        //      reactionIndx = domain->giveConnectivityTable() -> giveReactionIndx (dofManNum,
                        //                                        iDofMan->giveDof(j)->giveNumber());
                        reactionIndx = iDofMan->giveDof(j)->givePrescribedEquationNumber();
                        if ( reactionIndx ) {
                            reactions.at(reactionIndx) = dofManLoad.at(j);
                        }
                    }
                }
            }
        }
    }
}


#endif

/*
 * void
 * StructuralEngngModel:: computeElementEquivForces (FloatArray& answer, TimeStep *tStep, StructuralElement* ielem, ValueModeType mode) const
 * {
 * FloatArray nonnodalLoadVector  ;
 *
 * // compute element equivForces
 *
 * // use updated gp record
 * ((StructuralElement*)ielem)-> giveInternalForcesVector (answer, tStep,1);
 * //
 * // substract part corresponding to non-nodal loading.
 * // (temperature substracted in giveInternalForcesVector)
 * //
 * ((StructuralElement*)ielem)->
 * computeForceLoadVector(nonnodalLoadVector, tStep, mode) ;
 * if (nonnodalLoadVector.giveSize()) {
 * nonnodalLoadVector.negated();
 * answer.add (nonnodalLoadVector) ;
 * }
 *
 * }
 */

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
                rindex = jdof->givePrescribedEquationNumber();
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

    return;
}



#ifdef __PARALLEL_MODE
int
StructuralEngngModel :: packInternalForces(FloatArray *src, ProcessCommunicator &processComm)
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
            if ( jdof->isPrimaryDof() && ( eqNum = jdof->giveEquationNumber() ) ) {
                //fprintf (stderr, "[%d->%d] %d:%d -> %lf\n", rank,  processComm.giveRank(), toSendMap->at(i), j, src->at(eqNum));
                result &= pcbuff->packDouble( src->at(eqNum) );
            }
        }
    }

    return result;
}


int
StructuralEngngModel :: unpackInternalForces(FloatArray *dest, ProcessCommunicator &processComm)
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
            if ( jdof->isPrimaryDof() && ( eqNum = jdof->giveEquationNumber() ) ) {
                result &= pcbuff->unpackDouble(value);
                //fprintf (stderr, "[%d->%d] %d:%d <- %lf\n", rank,  processComm.giveRank(), toRecvMap->at(i), j, value);
                if ( dofmanmode == DofManager_shared ) {
                    dest->at(eqNum) += value;
                } else if ( dofmanmode == DofManager_remote ) {
                    dest->at(eqNum)  = value;
                } else {
                    _error("unpackInternalForces: unknown dof namager parallel mode");
                }
            }
        }
    }

    return result;
}

int
StructuralEngngModel :: packLoad(FloatArray *src, ProcessCommunicator &processComm)
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
            if ( jdof->isPrimaryDof() && ( eqNum = jdof->giveEquationNumber() ) ) {
                result &= pcbuff->packDouble( src->at(eqNum) );
            }
        }
    }

    return result;
}


int
StructuralEngngModel :: unpackLoad(FloatArray *dest, ProcessCommunicator &processComm)
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
            if ( jdof->isPrimaryDof() && ( eqNum = jdof->giveEquationNumber() ) ) {
                result &= pcbuff->unpackDouble(value);
                if ( dofmanmode == DofManager_shared ) {
                    dest->at(eqNum) += value;
                } else if ( dofmanmode == DofManager_remote ) {
                    dest->at(eqNum)  = value;
                } else {
		  _error3("unpackLoad: DofMamager %d[%d]: unknown parallel mode", dman->giveNumber(), dman->giveGlobalNumber());
                }
            }
        }
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
StructuralEngngModel :: packReactions(FloatArray *src, ProcessCommunicator &processComm)
{
    int result = 1;
    int i, size;
    int j, ndofs, eqNum;
    Domain *domain = this->giveDomain(1);
    IntArray const *toSendMap = processComm.giveToSendMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    //CommunicationBuffer* send_buff = processComm.giveSendBuff();
    DofManager *dman;
    Dof *jdof;

    size = toSendMap->giveSize();
    for ( i = 1; i <= size; i++ ) {
        dman = domain->giveDofManager( toSendMap->at(i) );
        ndofs = dman->giveNumberOfDofs();
        for ( j = 1; j <= ndofs; j++ ) {
            jdof = dman->giveDof(j);
            if ( jdof->isPrimaryDof() && ( eqNum = jdof->givePrescribedEquationNumber() ) ) {
                result &= pcbuff->packDouble( src->at(eqNum) );
            }
        }
    }

    return result;
}


int
StructuralEngngModel :: unpackReactions(FloatArray *dest, ProcessCommunicator &processComm)
{
    int result = 1;
    int i, size;
    int j, ndofs, eqNum;
    Domain *domain = this->giveDomain(1);
    dofManagerParallelMode dofmanmode;
    IntArray const *toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    //CommunicationBuffer* recv_buff = processComm.giveRecvBuff();
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
            if ( jdof->isPrimaryDof() && ( eqNum = jdof->givePrescribedEquationNumber() ) ) {
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
