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
StaticFracture :: StaticFracture(int i, EngngModel *_master) : NonLinearStatic(i, _master){}


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
    this->forceEquationNumbering();
}

void 
StaticFracture :: evaluatePropagationLaw(TimeStep *tStep)
{
    // For element wise evaluation of prop. law
    int ndomains = this->giveNumberOfDomains();
    int nnodes;
    Domain *domain;
    XfemManager *xMan; 

    bool propagateCrack = false;
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
                        //if ( dei->isElementEnrichedByEnrichmentDomain(el,k) ) { 
                            propagateCrack = true;  
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
                        //}
                    }

                }

            }

        }
    }


}

} // end namespace oofem
