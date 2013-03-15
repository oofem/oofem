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
    crackGrowthFlag = false;
}


void
StaticFracture :: solveYourselfAt(TimeStep *tStep)
{
    NonLinearStatic :: solveYourselfAt(tStep);

    /* 1) Compute fracture mechanics quantities
          What should be included
          - Element wise evaluation for cohesive zones
          - Evaluation of crack tip quantities: K, J, G
       2) Evaluate propagation model based on quantities above
       3) Update crack(s) according to propagation model

       Or should one turn it around. Set propagation model and then the model tries to evaluated what it needs.

       Need some class to keep track of the active propagation models and such
    */
    //this->evaluatePropagationLaw(tStep);
    //this->forceEquationNumbering();

}









void
StaticFracture :: terminate(TimeStep *tStep)
{
    NonLinearStatic :: terminate(tStep);


    this->evaluatePropagationLaw(tStep);
    /*
    this->forceEquationNumbering();


    // Just to make the simulation run through at the moment
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);
    totalDisplacement.resize(neq);
    totalDisplacement.zero();
    incrementOfDisplacement.resize(neq);
    incrementOfDisplacement.zero();

    initialLoadVector.resize(neq); 
    */

}



double 
StaticFracture ::  giveUnknownComponent(EquationID type, ValueModeType mode,
                                               TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
    return NonLinearStatic ::  giveUnknownComponent(type, mode, tStep, d, dof);

    
    if ( this->requiresUnknownsDictionaryUpdate() ) {
       if (mode == VM_Incremental) { //get difference between current and previous time variable
            return dof->giveUnknowns()->at(0) - dof->giveUnknowns()->at(1);
        }
        int hash = this->giveUnknownDictHashIndx(type, mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else {
            OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    }
    
}





void
StaticFracture :: createPreviousSolutionInDofUnknownsDictionary(TimeStep *tStep) {
    //Copy the last known temperature to be a previous solution
    int nnodes, nDofs;
    double val;
    Domain *domain;
    Dof *iDof;
    DofManager *node;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        nnodes = domain->giveNumberOfDofManagers();
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( int inode = 1; inode <= nnodes; inode++ ) {
                node = domain->giveDofManager(inode);
                nDofs = node->giveNumberOfDofs();
                for ( int i = 1; i <= nDofs; i++ ) {
                    iDof = node->giveDof(i);
                    val = iDof->giveUnknown(EID_ConservationEquation, VM_Total, tStep); //get number on hash=0(current)
                    iDof->updateUnknownsDictionary(tStep->givePreviousStep(), EID_MomentumBalance, VM_Total, val);
                }
            }
        }
    }
}




/*
int
NLTransientTransportProblem :: giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN) {
    if ( mode == VM_Total ) { //Nodal temperature
        if ( stepN->giveNumber() == this->giveCurrentStep()->giveNumber() ) { //current time
            return 0;
        } else if ( stepN->giveNumber() == this->giveCurrentStep()->giveNumber() - 1 ) { //previous time
            return 1;
        } else {
            _error5( "No history available at TimeStep %d = %f, called from TimeStep %d = %f", stepN->giveNumber(), stepN->giveTargetTime(), this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveTargetTime() );
        }
    } else {
        _error2( "ValueModeType %s undefined", __ValueModeTypeToString(mode) );
    }

    return 0;
}

void
NLTransientTransportProblem :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    // update DoF unknowns dictionary. Store the last and previous temperature only, see giveUnknownDictHashIndx
    int i, ndofs = inode->giveNumberOfDofs();
    int eqNum;
    Dof *iDof;
    double val;
    FloatArray *vect;

    for ( i = 1; i <= ndofs; i++ ) {
        iDof = inode->giveDof(i);
        eqNum = iDof->__giveEquationNumber();
        if ( iDof->hasBc(tStep) ) { // boundary condition
            val = iDof->giveBcValue(VM_Total, tStep);
        } else {
            vect = this->UnknownsField->giveSolutionVector(tStep);
            val = vect->at(eqNum);
        }

        //update temperature, which is present in every node
        iDof->updateUnknownsDictionary(tStep, EID_MomentumBalance, VM_Total, val);
    }
}


void
NLTransientTransportProblem :: copyUnknownsInDictionary(EquationID type, ValueModeType mode, TimeStep *fromTime, TimeStep *toTime) {
    int i, j, ndofs;
    double val;
    Domain *domain = this->giveDomain(1);
    int nnodes = domain->giveNumberOfDofManagers();
    DofManager *inode;
    Dof *iDof;

    for ( j = 1; j <= nnodes; j++ ) {
        inode = domain->giveDofManager(j);
        ndofs = inode->giveNumberOfDofs();
        for ( i = 1; i <= ndofs; i++ ) {
            iDof = inode->giveDof(i);
            val = iDof->giveUnknown(type, mode, fromTime);
            iDof->updateUnknownsDictionary(toTime, type, mode, val);
        }
    }
}


void
NLTransientTransportProblem :: updateInternalState(TimeStep *stepN)
{
    int j, nnodes;
    Domain *domain;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        nnodes = domain->giveNumberOfDofManagers();
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( j = 1; j <= nnodes; j++ ) {
                //update dictionary entry or add a new pair if the position is missing
                this->updateDofUnknownsDictionary(domain->giveDofManager(j), stepN);
            }
        }

        int nelem = domain->giveNumberOfElements();
        for ( j = 1; j <= nelem; j++ ) {
            domain->giveElement(j)->updateInternalState(stepN);
        }
    }
}



*/





void 
StaticFracture :: evaluatePropagationLaw(TimeStep *tStep)
{
    // For element wise evaluation of prop. law
    int ndomains = this->giveNumberOfDomains();
    int nnodes;
    Domain *domain;
    XfemManager *xMan; 


    for ( int idomain = 1; idomain <= ndomains; idomain++ ) {
        domain = this->giveDomain(idomain);        
        xMan = domain->giveXfemManager(1);
        
        EnrichmentItem *ei;
        for ( int j = 1; j <= xMan->giveNumberOfEnrichmentItems(); j++ ) {    
            ei = xMan->giveEnrichmentItem(j);
            // Different actions depending on ei

            // Delamination 
            if ( Delamination *dei = dynamic_cast< Delamination * > (ei) )  {
                EnrichmentDomain *ed; 
                for ( int k = 1; k <= ei->giveNumberOfEnrichmentDomains(); k++ ) {
                    ed = ei->giveEnrichmentDomain(k);
                    Element *el;
                    for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) {
                        el = domain->giveElement(i);

                        crackGrowthFlag = true;  
                        if ( i== 1 && tStep->isTheFirstStep() ) {
                            // updateEnrichmentDomain
                            IntArray dofManNumbers;
                            dofManNumbers.setValues(3, 3,5,6 );
                            if ( DofManList *ded = dynamic_cast< DofManList * > (ed) )  {
                                ded->addDofManagers(dofManNumbers);
                                IntArray dofIdArray;
                                for ( int m = 1; m <= dofManNumbers.giveSize(); m++ ) {
                                    DofManager *dMan = domain->giveDofManager( dofManNumbers.at(m) );
                                    ei->computeDofManDofIdArray(dofIdArray, dMan, k);
                                    xMan->addEnrichedDofsTo( dMan, dofIdArray );
                                }

                            }
                        }

                    }

                }

            }

        }
    }


}

} // end namespace oofem
