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
// Should be moved to future FractureManager class later
//#include "gausspnt.h"
#include <vector>

namespace oofem {
StaticFracture :: StaticFracture(int i, EngngModel *_master) : NonLinearStatic(i, _master)
{
    crackGrowthFlag = true; // if true, then the internal structure needs to be updated
}


void
StaticFracture :: solveYourselfAt(TimeStep *tStep)
{

    // Initiates the total displacement to zero in the UnknownsDictionary.
    if ( tStep->isTheFirstStep() ) {
        printf("Initializing DofUnknownsDictionary... \n");
        this->initializeDofUnknownsDictionary(tStep);
    }
    
    // Initialization
    int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());
    if ( totalDisplacement.giveSize() != neq ) {
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();
        this->setTotalDisplacementFromUnknownsInDictionary(EID_MomentumBalance, VM_Total, tStep);
    }
    
    // Instanciate fracture manager
    if ( tStep->isTheFirstStep() ) {
    
        this->fMan =  new FractureManager( this->giveDomain(1) );
        this->fMan->failureCriterias = new AList< FailureCriteria >(1);
        FailureCriteria *fc = new FailureCriteria(FC_MaxShearStress);
        fc->thresholds.resize(1);
        fc->thresholds.at(1) = 0.0009;
        this->fMan->failureCriterias->put(1, fc);
    }


    crackGrowthFlag = false;
    NonLinearStatic :: solveYourselfAt(tStep);

}

void 
StaticFracture :: updateYourself(TimeStep *stepN)
{
    
    NonLinearStatic :: updateYourself(stepN);

    this->fMan->evaluateFailureCriterias(stepN);


    this->evaluatePropagationLaw(stepN); 

    // Update the UnknownsDictionary if needed
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        printf(" Updating DofUnknownsDictionary... \n");
        for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
            Domain *domain = this->giveDomain(idomain);
            int nnodes = domain->giveNumberOfDofManagers();
            for ( int inode = 1; inode <= nnodes; inode++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(inode), stepN);
            }
        }
    }
   
}


void
StaticFracture :: terminate(TimeStep *tStep)
{
    NonLinearStatic :: terminate(tStep);
}



void
StaticFracture :: updateLoadVectors(TimeStep *stepN)
{
    MetaStep *mstep = this->giveMetaStep( stepN->giveMetaStepNumber() );
    bool isLastMetaStep = ( stepN->giveNumber() == mstep->giveLastStepNumber() );

    if ( controlMode == nls_indirectControl ) { //todo@: not checked 
        //if ((stepN->giveNumber() == mstep->giveLastStepNumber()) && ir->hasField("fixload")) {
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
    } else { // direct control
        //update initialLoadVector after each step of direct control
        //(here the loading is not proportional)
        
        /*if ( initialLoadVector.isEmpty() ) {
            initialLoadVector.resize( incrementalLoadVector.giveSize() );
        }
        */
        OOFEM_LOG_DEBUG("Fixed load level\n");

        incrementalLoadVector.times(loadLevel);
        if ( initialLoadVector.giveSize() != incrementalLoadVector.giveSize() ) {
            //initialLoadVector.resize( incrementalLoadVector.giveSize() );
            initialLoadVector.resize( 0 );
        }
        //initialLoadVector.add(incrementalLoadVector);

        incrementalLoadVectorOfPrescribed.times(loadLevel);
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



// Fracture mechanics evaluation 
// Should be in a searate class probably - FractureManager maybe

void 
StaticFracture :: evaluatePropagationLaw(TimeStep *tStep)
{
    // For element wise evaluation of prop. law
    Domain *domain;
    XfemManager *xMan; 
    Element *el;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);        
        xMan = domain->giveXfemManager(1);
        
        for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) {
            el = domain->giveElement(i);
            EnrichmentItem *ei;
            for ( int j = 1; j <= xMan->giveNumberOfEnrichmentItems(); j++ ) {    
                ei = xMan->giveEnrichmentItem(j);
                // Different actions depending on ei
                // Delamination 
                if ( Delamination *dei = dynamic_cast< Delamination * > (ei) )  {
                    for ( int k = 1; k <= ei->giveNumberOfEnrichmentDomains(); k++ ) {                
                        EnrichmentDomain *ed; 
                        ed = ei->giveEnrichmentDomain(k);
                        this->evaluatePropagationLawForDelamination(el, ed, tStep); 
                    }
                }
            }
        }

    }

    // Create additional dofs if any Enrichment Domain is updated
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        xMan->createEnrichedDofs();
    }
}


void 
StaticFracture :: evaluatePropagationLawForDelamination(Element *el, EnrichmentDomain *ed, TimeStep *tStep)
{
    
    int interfaceNum = 1;
    std::vector < FloatArray > interLamStresses;
    this->computeInterLaminarStressesAt(interfaceNum, el, tStep, interLamStresses);
    bool propagateFlag = false;
    this->evaluateFractureCriterion( interLamStresses, propagateFlag );

    if ( propagateFlag ) {
        crackGrowthFlag = true;
        printf( "\n -------------------------------\n");
        printf( "Crack growth in element %d \n", el->giveNumber() );
        

        
        IntArray dofManNumbers, elDofMans;
        elDofMans = el->giveDofManArray();
        for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {  
            // ugly piece of code that will skip enrichment of dofmans that have any bc's
            // which is not generally what you want
            #if 1
            bool hasBc= false;
            for ( int j = 1; j <= el->giveDofManager(i)->giveNumberOfDofs(); j++ ) {
                if ( el->giveDofManager(i)->giveDof(j)->hasBc(tStep) ) {
                    hasBc = true;
                    continue;
                }
            }
            #endif
            if ( !hasBc) {
                dofManNumbers.followedBy(elDofMans.at(i));
            }
        }
        
        //dofManNumbers.printYourself();
        if ( DofManList *ded = dynamic_cast< DofManList * > (ed) )  {
            ded->updateEnrichmentDomain(dofManNumbers);
        }
    }


}




void 
StaticFracture :: computeInterLaminarStressesAt(int interfaceNum, Element *el, TimeStep *tStep, std::vector < FloatArray > &interLamStresses)
{
    
    IntegrationRule *irLower = el->giveIntegrationRule(interfaceNum-1); // index from 0
    IntegrationRule *irUpper = el->giveIntegrationRule(interfaceNum);
    GaussPoint *gp;

    // find gp-pairs at the interface - should be simplified by storing this info during gp creation
    // assumes that the number of points are the same and the ordering is the same for each layer
    int numGP = irLower->getNumberOfIntegrationPoints();
    double xiMax = -1.0;
    double xiMin =  1.0;
    
    // find min and max xi coord
    for ( int i = 1; i <= numGP; i++ ) {
        gp = irLower->getIntegrationPoint(i-1);
        if ( gp->giveCoordinate(3) > xiMax ) {
            xiMax = gp->giveCoordinate(3); 
        }
        gp = irUpper->getIntegrationPoint(i-1);
        if ( gp->giveCoordinate(3) < xiMin ) {
            xiMin = gp->giveCoordinate(3); 
        }
    }

    const double tol = 1.0e-8;
    IntArray upper(0), lower(0);
    for ( int i = 0; i < numGP; i++ ) {
        gp = irLower->getIntegrationPoint(i);
        if ( abs( gp->giveCoordinate(3)-xiMax) < tol ) { // upper gp
            upper.followedBy(i);   
        }
        gp = irUpper->getIntegrationPoint(i);
        if ( abs( gp->giveCoordinate(3)-xiMin) < tol ) { // lower gp
            lower.followedBy(i); 
        }
    }

    // compute stresses
    int numInterfaceGP = lower.giveSize();
    interLamStresses.resize(numInterfaceGP);
    FloatArray vSLower, vSUpper;
    // compute mean interface stress for a gp-pair
    for ( int i = 1; i <= numInterfaceGP; i++ ) {
        gp = irUpper->getIntegrationPoint( upper.at(i) );
        el->giveIPValue(vSUpper, gp, IST_CauchyStressTensor, tStep);
        gp = irLower->getIntegrationPoint( lower.at(i) );
        el->giveIPValue(vSLower, gp, IST_CauchyStressTensor, tStep);

        interLamStresses.at(i-1).resize( vSUpper.giveSize() );
        interLamStresses.at(i-1) = 0.5 * ( vSUpper + vSLower );
    }

}



void 
StaticFracture :: evaluateFractureCriterion(std::vector < FloatArray > &interLamStresses, bool &propagateFlag)
{
    propagateFlag = false;
    for (int  i = 1; i <= interLamStresses.size(); i++ ) {

        if ( interLamStresses[i-1].at(3) > 0.0009 ) {
            propagateFlag = true;
            
            return;
        }
    }
}


} // end namespace oofem
