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
#include "Contact/contactdefinition.h"
#include "Contact/contactelement.h"
#include "intarray.h"
#include "domain.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"

namespace oofem {


ContactDefinition :: ContactDefinition(ContactManager *cMan)
{
    this->cMan = cMan;
}
    /// Destructor.
ContactDefinition :: ~ContactDefinition()
{
}




IRResultType
ContactDefinition :: initializeFrom(InputRecord *ir)
{
  
  // Example of node to node contact with two elements
  
    /*  4 - 3  <->  8 - 7
     *  |   |       |   |
     *  1 - 2  <->  5 - 6
     */

    Domain *domain = this->cMan->giveDomain();
    
    IntArray masterNodes = {2, 3};
    IntArray slaveNodes  = {5, 8};
    this->masterElementList.resize( masterNodes.giveSize() );
    for( int i = 1; i<= masterNodes.giveSize(); i++ ) {

        ContactElement *master = new Node2NodeContactL( domain->giveDofManager(masterNodes.at(i)),
                                                       domain->giveDofManager(slaveNodes.at(i)) );
        this->masterElementList[i-1] = std :: move(master);
    }
      
    return IRRT_OK;
}


int
ContactDefinition :: instanciateYourself(DataReader *dr)
{
  
    for ( ContactElement *cEl : this->masterElementList ) { 
        cEl->instanciateYourself(dr);
    }
    
  return 1;
}


void
ContactDefinition :: computeContactForces(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
    //Loop through all the contact elements and let them return their internal forces vector
    FloatArray Fc;
    IntArray locArray;
    
    for ( ContactElement *master : this->masterElementList ) {
          
        master->computeContactForces(Fc, tStep, type, mode, s, domain, eNorms);
        if ( master->isInContact() ) {
          master->giveLocationArray(locArray, s);
          Fc.negated();
          answer.assemble(Fc, locArray);
        
          if ( eNorms ) {
              eNorms->assembleSquared( Fc, locArray );
          }
          
        }
    }
  
  
}

void
ContactDefinition :: computeContactTangent(SparseMtrx *answer, TimeStep *tStep,
                      CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
  
    FloatMatrix Kc;
    IntArray locArrayR, locArrayC;
    
    for ( ContactElement *master : this->masterElementList ) {
        
        if ( master->isInContact() ) {
            //printf("node in contact: computeContactTangent\n\n");
            master->computeContactTangent(Kc, type, tStep);
            
           
            master->giveLocationArray(locArrayR, r_s);
            master->giveLocationArray(locArrayC, c_s);

            answer->assemble(locArrayR, locArrayC, Kc);
        }
    }
    
  
}






}
