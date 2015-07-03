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

#include "contact/contactmanager.h"
#include "contact/contactdefinition.h"
#include "contact/contactelement.h"
#include "intarray.h"
#include "domain.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "masterdof.h"
#include "classfactory.h"

namespace oofem {
REGISTER_ContactDefinition(ContactDefinition)


ContactDefinition :: ContactDefinition(ContactManager *cMan)
{
    this->cMan = cMan;
    this->numberOfConstraintEq = 0;
}


ContactDefinition :: ~ContactDefinition()
{
}


int
ContactDefinition :: instanciateYourself(DataReader *dr)
{
    for ( ContactElement *cEl : this->masterElementList ) {
        cEl->instanciateYourself(dr);
        cEl->setupIntegrationPoints();
    }
    
    return 1;
}


void 
ContactDefinition :: createContactDofs()
{
    // Creates new dofs due associated with the contact (Lagrange multipliers) and appends them to the dof managers
    // Creates new dofs due associated with the contact (Lagrange multipliers) and appends them to the dof managers

    //TODO This is a bit ugly, find a better solution than asking the contact el
    if ( int numDofs = this->giveNumberOfConstraintEqToAdd() ) {
        
        // get an array with dof ids' to append to 
        IntArray dofIdArray(numDofs), dofMans;
        for ( int i = 1; i <= numDofs; i++ ) {
            dofIdArray.at(i) = this->cMan->giveDomain()->giveNextFreeDofID();
        }

        for ( ContactElement *cEl : this->masterElementList ) { 
            
            cEl->giveDofManagersToAppendTo(dofMans);
            if ( dofMans.giveSize() ) { // if the contact element adds extra dofs, store them. Maybe just store in cDef?
                cEl->setDofIdArray(dofIdArray);
            }
            
            for ( int i = 1; i <= dofMans.giveSize(); i++ ) {
                DofManager *dMan = this->cMan->giveDomain()->giveDofManager(dofMans.at(i));
                for ( auto &dofid: dofIdArray ) {
                    if ( !dMan->hasDofID( ( DofIDItem ) ( dofid ) ) ) {
                    
                        dMan->appendDof( new MasterDof( dMan, ( DofIDItem ) dofid ) );
                    }
                }
            }
        }  
    }
}


void
ContactDefinition :: computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
    //Loop through all the contact elements and let them return their internal forces vector
    FloatArray Fc;
    IntArray locArray;
    
    // TODO ask masters that are potentially in contact and not everyone
    for ( auto &master : this->masterElementList ) {
        
        // These acts as external forces so move them to the lhs
        master->computeContactForces(Fc, tStep, type, mode, s, domain, eNorms);
        Fc.negated();
        
        if ( Fc.giveSize() ) {
            master->giveLocationArray(locArray, s);
            answer.assemble(Fc, locArray);
        
            if ( eNorms ) {
                eNorms->assembleSquared( Fc, locArray );
            }
        }
    }
}


void
ContactDefinition :: computeContactTangent(SparseMtrx &answer, TimeStep *tStep,
                      CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
  
    FloatMatrix Kc;
    IntArray locArrayR, locArrayC;
    
    for ( auto &master : this->masterElementList ) {
        
        //if ( master->isInContact() ) { // tangent becomes singular with this
            //printf("node in contact: computeContactTangent\n\n");
            master->computeContactTangent(Kc, type, tStep);
            // do this in contact element?
            Kc.negated(); // should be negated!
           
            master->giveLocationArray(locArrayR, r_s);
            master->giveLocationArray(locArrayC, c_s);
            
            answer.assemble(locArrayR, locArrayC, Kc);
        //}
    }
}


}
