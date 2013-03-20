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
        this->initializeDofUnknownsDictionary(tStep);
    }
    
    // Initialization
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);
    if ( totalDisplacement.giveSize() != neq ) {
        totalDisplacement.resize(neq);
        totalDisplacement.zero();
        incrementOfDisplacement.resize(neq);
        incrementOfDisplacement.zero();
        this->setTotalDisplacementFromUnknownsInDictionary(EID_MomentumBalance, VM_Total, tStep);
    }
    
    crackGrowthFlag = false;
    NonLinearStatic :: solveYourselfAt(tStep);
    //crackGrowthFlag = true;
    /* 1) Compute fracture mechanics quantities
          What should be included
          - Element wise evaluation for cohesive zones
          - Evaluation of crack tip quantities: K, J, G
       2) Evaluate propagation model based on quantities above
       3) Update crack(s) according to propagation model

       Or should one turn it around. Set propagation model and then the model tries to evaluated what it needs.

       Need some class to keep track of the active propagation models and such
    */
    this->evaluatePropagationLaw(tStep);
    //this->forceEquationNumbering();

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
            if ( !mstep->giveAttributesRecord()->hasField(IFT_NonLinearStatic_donotfixload, "donotfixload") ) {
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
    if ( isLastMetaStep && !mstep->giveAttributesRecord()->hasField(IFT_NonLinearStatic_donotfixload, "donotfixload") ) {
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



double 
StaticFracture ::  giveUnknownComponent(EquationID type, ValueModeType mode,
                                               TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        int hash = this->giveUnknownDictHashIndx(type, mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            return 0.0; ///@todo: how should one treat newly created dofs?
            //OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    } else {
        return NonLinearStatic ::  giveUnknownComponent(type, mode, tStep, d, dof);
    }
    
}





void
StaticFracture :: initializeDofUnknownsDictionary(TimeStep *tStep) {
    //
    int nnodes, nDofs;
    double val;
    Domain *domain;
    Dof *iDof;
    DofManager *node;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        nnodes = domain->giveNumberOfDofManagers();
        for ( int inode = 1; inode <= nnodes; inode++ ) {
            node = domain->giveDofManager(inode);
            nDofs = node->giveNumberOfDofs();
            for ( int i = 1; i <= nDofs; i++ ) {
                iDof = node->giveDof(i);
                iDof->updateUnknownsDictionary(tStep->givePreviousStep(), EID_MomentumBalance, VM_Total, 0.0);
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
        //} else if ( tStep->isTheFirstStep() ) { // initialize to zero
        //     val = 0.0;
        } else {
            if ( eqNum > 0 ) {
                val = totalDisplacement.at(eqNum);
            } else { // new eq number
                val = 0.0;
            }
        }

        iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total, val);
    }
}


void
StaticFracture :: setTotalDisplacementFromUnknownsInDictionary(EquationID type, ValueModeType mode, TimeStep *tStep) 
{
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
                    double val = iDof->giveUnknown(type, mode, tStep);
                    totalDisplacement.at(eqNum) = val;
                }
            }
        }   
    }

}






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
    // temporary 
    crackGrowthFlag = true;  
    if ( el->giveNumber()== 1 && tStep->isTheFirstStep() ) {

        // updateEnrichmentDomain
        // prop law gives:
        IntArray dofManNumbers;
        

        // for example compute average stress and compare to criteria
        //for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
            for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {     
                // ugly piece of code that will skip enrichment of dofmans that have any bc's
                // which is not generally what you want
                bool hasBc= false;
                for ( int j = 1; j <= el->giveDofManager(i)->giveNumberOfDofs(); j++ ) {
                    if ( el->giveDofManager(i)->giveDof(j)->hasBc(tStep) ) {
                        hasBc = true;
                        continue;
                    }
                }
                if ( !hasBc) {
                dofManNumbers.followedBy(i);
                }
        }

        if ( DofManList *ded = dynamic_cast< DofManList * > (ed) )  {
            ded->updateEnrichmentDomain(dofManNumbers);
        }
    }


}

} // end namespace oofem
