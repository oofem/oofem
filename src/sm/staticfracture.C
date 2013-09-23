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

#include "staticfracture.h"
#include "fracturemanager.h"
#include "dofmanager.h"
#include "enrichmentdomain.h"
#include "classfactory.h"

#include <vector>

namespace oofem {

REGISTER_EngngModel( StaticFracture );

StaticFracture :: StaticFracture(int i, EngngModel *_master) : NonLinearStatic(i, _master)
{
    updateStructureFlag = false; // if true, then the internal structure needs to be updated
}


void
StaticFracture :: solveYourselfAt(TimeStep *tStep)
{

    // Initiates the total displacement to zero in the UnknownsDictionary at the first time step.
    // this must be done in order to support a dynamic equation system
    if ( tStep->isTheFirstStep() ) {
        printf("Initializing DofUnknownsDictionary... \n");
        this->initializeDofUnknownsDictionary(tStep);
    }
    
    // Initialization
    int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering()); // 1 stands for domain?
    if ( totalDisplacement.giveSize() != neq ) {
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();
        this->setTotalDisplacementFromUnknownsInDictionary(EID_MomentumBalance, VM_Total, tStep);
    }
    
    // Instanciate fracture manager
    // should be made in a more proper way with input and the like and moved to another part
    if ( tStep->isTheFirstStep() ) {
        IntArray criteriaList(1);
        criteriaList.at(1) = 1; // criteria 1 (only one supported)
        Domain *domain= this->giveDomain(1);
        this->fMan = new FractureManager( domain );
        
        // initialize failure criteria managers
        int numEl = domain->giveNumberOfElements();
        
        this->fMan->criteriaManagers.resize( criteriaList.giveSize() );

        for ( int i = 1; i <= criteriaList.giveSize(); i++ ) {
            // if local
            this->fMan->criteriaManagers.at(i-1) = new FailureCriteriaManager(Local, this->fMan);
            FailureCriteriaManager *cMan = this->fMan->criteriaManagers.at(i-1);
            cMan->list.resize(numEl);
            for ( int j = 1; j <= numEl; j++ ) { 
                cMan->list.at(j - 1) = new DamagedNeighborLayered(Local, this->fMan);
                cMan->list.at(j - 1)->thresholds.resize(1);
                cMan->list.at(j - 1)->thresholds.at(1) = -10.0;
                cMan->list.at(j - 1)->el = domain->giveElement(j);           
            }
        }
    }


    this->setUpdateStructureFlag(false);
    NonLinearStatic :: solveYourselfAt(tStep);

}

void 
StaticFracture :: updateYourself(TimeStep *tStep)
{
#if 1
    this->fMan->evaluateYourself(tStep);

    // Update XFEM structure based on the fracture manager
    this->fMan->updateXFEM(tStep); 
#endif    
    NonLinearStatic :: updateYourself(tStep);

    // Fracture/failure mechanics evaluation

    this->setUpdateStructureFlag( this->fMan->giveUpdateFlag() ); // if the internal structure need to be updated

    // Update the UnknownsDictionary if needed
    if ( this->needsStructureUpdate() ) {
        printf(" Updating DofUnknownsDictionary... \n");
        for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
            Domain *domain = this->giveDomain(idomain);
            int nnodes = domain->giveNumberOfDofManagers();
            for ( int inode = 1; inode <= nnodes; inode++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(inode), tStep);
            }
        }
    }


}

// remove
void
StaticFracture :: terminate(TimeStep *tStep)
{
    NonLinearStatic :: terminate(tStep);
}



void
StaticFracture :: updateLoadVectors(TimeStep *tStep)
{
    MetaStep *mstep = this->giveMetaStep( tStep->giveMetaStepNumber() );
    bool isLastMetaStep = ( tStep->giveNumber() == mstep->giveLastStepNumber() );

    if ( controlMode == nls_indirectControl ) { //todo@: not checked 
        #if 0
        //if ((tStep->giveNumber() == mstep->giveLastStepNumber()) && ir->hasField("fixload")) {
        if ( isLastMetaStep ) {
            if ( !mstep->giveAttributesRecord()->hasField(_IFT_NonLinearStatic_donotfixload) ) {
                OOFEM_LOG_INFO("Fixed load level\n");

                //update initialLoadVector
                if ( initialLoadVector.isEmpty() ) {
                    initialLoadVector.resize( incrementalLoadVector.giveSize() );
                }

                incrementalLoadVector.times(loadLevel);
                initialLoadVector.add(incrementalLoadVector);

                incrementalLoadVectorOfPrescribed.times(loadLevel);
                initialLoadVectorOfPrescribed.add(incrementalLoadVectorOfPrescribed);

                incrementalLoadVector.zero();
                incrementalLoadVectorOfPrescribed.zero();

                this->loadInitFlag = 1;
            }

            //if (!mstep->giveAttributesRecord()->hasField("keepll")) this->loadLevelInitFlag = 1;
        }
        #endif

    } else { // direct control
        //update initialLoadVector after each step of direct control
        //(here the loading is not proportional)

        OOFEM_LOG_DEBUG("Fixed load level\n");

        incrementalLoadVector.times(loadLevel);
        if ( initialLoadVector.giveSize() != incrementalLoadVector.giveSize() ) {
            initialLoadVector.resize( 0 );
        }

        incrementalLoadVectorOfPrescribed.times(loadLevel);
        
        //initialLoadVectorOfPrescribed.zero();
        initialLoadVectorOfPrescribed.add(incrementalLoadVectorOfPrescribed);


        incrementalLoadVector.zero();
        incrementalLoadVectorOfPrescribed.zero();

        this->loadInitFlag = 1;
    }


    // if (isLastMetaStep) {
    if ( isLastMetaStep && !mstep->giveAttributesRecord()->hasField(_IFT_NonLinearStatic_donotfixload) ) {
#ifdef VERBOSE
        OOFEM_LOG_INFO("Reseting load level\n");
#endif
        if ( mstepCumulateLoadLevelFlag ) {
            cumulatedLoadLevel += loadLevel;
        } else {
            cumulatedLoadLevel = 0.0;
        }

        this->loadLevel = 0.0;
    }
}



// Updating dofs and evaluating unknowns
#if 1

double 
StaticFracture ::  giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    // Returns the unknown quantity corresponding to the dof
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        int hash = this->giveUnknownDictHashIndx(mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else { // Value is not initiated in UnknownsDictionary
            return 0.0; ///@todo: how should one treat newly created dofs?
            //OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    } else {
        return NonLinearStatic ::  giveUnknownComponent(mode, tStep, d, dof);
    }
    
}


void
StaticFracture :: initializeDofUnknownsDictionary(TimeStep *tStep) 
{
    // Initializes all dof values to zero
    
    Domain *domain;
    Dof *iDof;
    DofManager *node;
    
    int nDofs;
    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        int nnodes = domain->giveNumberOfDofManagers();
        for ( int inode = 1; inode <= nnodes; inode++ ) {
            node = domain->giveDofManager(inode);
            nDofs = node->giveNumberOfDofs();
            for ( int i = 1; i <= nDofs; i++ ) {
                iDof = node->giveDof(i);
                iDof->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, 0.0);
            }
        }
    }
}



void
StaticFracture :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{

    // update DoF unknowns dictionary. 
    Dof *iDof;
    double val;
    for ( int i = 1; i <= inode->giveNumberOfDofs(); i++ ) {
        iDof = inode->giveDof(i);
        int eqNum = iDof->__giveEquationNumber();
        if ( iDof->hasBc(tStep) ) { 
            val = iDof->giveBcValue(VM_Total, tStep);
        } else {
            if ( eqNum > 0 ) {
                val = totalDisplacement.at(eqNum);
            } else { // new eq number
                val = 0.0;
            }
        }

        iDof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}


void
StaticFracture :: setTotalDisplacementFromUnknownsInDictionary(EquationID type, ValueModeType mode, TimeStep *tStep) 
{
    // Sets the values in the displacement vector based on stored values in the unknowns dictionaries.
    // Used in the beginning of each time step.
    Domain *domain;
    DofManager *inode;
    Dof *iDof;
    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        for ( int j = 1; j <= domain->giveNumberOfDofManagers(); j++ ) {
            inode = domain->giveDofManager(j);
            int eqNum;
            for ( int i = 1; i <= inode->giveNumberOfDofs(); i++ ) {
                iDof = inode->giveDof(i);
                eqNum = iDof->giveEqn();
                if ( eqNum > 0 ) {
                    double val = iDof->giveUnknown(mode, tStep);
                    totalDisplacement.at(eqNum) = val;
                }
            }
        }   
    }

}

#endif




} // end namespace oofem
