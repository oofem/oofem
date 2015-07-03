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
#include "classfactory.h"
#include "numericalcmpn.h"
#include "error.h"

namespace oofem {
REGISTER_ContactManager(ContactManager)

ContactManager :: ContactManager(Domain *domain)
{
    this->domain = domain;
    numberOfContactDefinitions = -1;
}

ContactManager :: ~ContactManager()
{
}


IRResultType
ContactManager :: initializeFrom(InputRecord *ir)
{

    IRResultType result; // Required by IR_GIVE_FIELD macro
    
    this->numberOfContactDefinitions = 0;
    IR_GIVE_FIELD(ir, this->numberOfContactDefinitions, _IFT_ContactManager_NumberOfContactDefinitions);
  
    
    this->contactDefinitionList.resize(this->numberOfContactDefinitions);
#if 0
    for ( int i = 1; i <= numberOfContactDefinitions; i++ ) {
        std::string name;
        result = ir->giveRecordKeywordField(name);
        this->contactDefinitionList[i-1] = new ContactDefinition(this);
        this->contactDefinitionList[i-1] = classFactory.createContactDefinition( name.c_str(), this );
    }
#endif

    return IRRT_OK;
}


int 
ContactManager :: instanciateYourself(DataReader *dr)
{
    IRResultType result = IRRT_OK; // Required by IR_GIVE_FIELD macro
    std :: string name;

    // Create and instantiate contact definitions
    for ( int i = 1; i <= this->giveNumberOfContactDefinitions(); i++ ) {
        InputRecord *ir = dr->giveInputRecord(DataReader :: IR_contactDefRec, i);
        result = ir->giveRecordKeywordField(name);  
        this->contactDefinitionList[i-1].reset( classFactory.createContactDefinition( name.c_str(), this ) );
        if ( this->contactDefinitionList[i-1] ) {
            this->contactDefinitionList[i-1]->initializeFrom(ir);
            this->contactDefinitionList[i-1]->instanciateYourself(dr);
        } else {
            OOFEM_ERROR("Failed to create contact definition (%s)", name.c_str() );
        }
    }

    return result;
}


void 
ContactManager :: assembleVectorFromContacts(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
    if ( type == InternalForcesVector) {
        //printf("\n Add forces due to contact... \n");
        for ( auto &cDef : contactDefinitionList ) {
            cDef->computeContactForces(answer, tStep, type, mode, s, domain, eNorms);
        }
    }
}


void 
ContactManager :: assembleTangentFromContacts(SparseMtrx &answer, TimeStep *tStep,
                      CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) 
{ 
    //printf("\n Add tangents due to contact... \n");
    for ( auto &cDef : contactDefinitionList ) {
        cDef->computeContactTangent(answer, tStep, type, r_s, c_s);
    }
}


void
ContactManager :: createContactDofs()
{
    // Creates new dofs contacts and appends them to the dof managers
    for ( auto &cDef : contactDefinitionList ) {
        cDef->createContactDofs();
    }
}

}
