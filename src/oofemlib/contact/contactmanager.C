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
#include "classfactory.h"
#include "numericalcmpn.h"

namespace oofem {
REGISTER_ContactManager(ContactManager)

ContactManager :: ContactManager(Domain *domain)
{
    this->domain = domain;
    numberOfContactDefinitions = -1;
}

ContactManager :: ~ContactManager()
{
    
    for ( ContactDefinition *cDef : contactDefinitionList ) {
        delete cDef;
    }
    
}



IRResultType
ContactManager :: initializeFrom(InputRecord *ir)
{
  
    // define one contact
    this->numberOfContactDefinitions = 1;
    this->contactDefinitionList.resize(this->numberOfContactDefinitions);
    

    for ( int i = 1; i <= numberOfContactDefinitions; i++ ) {
     
        ContactDefinition *cDef = new ContactDefinition(this);
        cDef->initializeFrom(ir);
        this->contactDefinitionList[i-1] = std :: move(cDef);
        
    }
    
    
    return IRRT_OK;
}


int 
ContactManager :: instanciateYourself(DataReader *dr)
{
    
    for ( ContactDefinition *cDef : this->contactDefinitionList ) { 
        cDef->instanciateYourself(dr);
    }
    
    return 1;
}



void 
ContactManager :: assembleVectorFromContacts(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
  if ( type == InternalForcesVector) {
    
      //printf("\n Add forces due to contact... \n");
      for ( ContactDefinition *cDef : contactDefinitionList ) {
          cDef->computeContactForces(answer, tStep, type, mode, s, domain, eNorms);
      }
  }
};


void 
ContactManager ::assembleTangentFromContacts(SparseMtrx *answer, TimeStep *tStep,
                      CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) 
{ 
  
  //printf("\n Add tangents due to contact... \n");
  for ( ContactDefinition *cDef : contactDefinitionList ) {
      cDef->computeContactTangent(answer, tStep, type, r_s, c_s);
  }
  
  
}








}
