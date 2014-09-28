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

#include "Contact/cdefnode2node.h"
#include "Contact/celnode2node.h"
#include "domain.h"
#include "classfactory.h"

namespace oofem {
REGISTER_ContactDefinition(ContactDefinitionNode2Node)
REGISTER_ContactDefinition(ContactDefinitionNode2NodeL)


ContactDefinitionNode2Node :: ContactDefinitionNode2Node(ContactManager *cMan) : ContactDefinition(cMan){}
    

IRResultType
ContactDefinitionNode2Node :: initializeFrom(InputRecord *ir)
{
 
    IRResultType result; // Required by IR_GIVE_FIELD macro
    
    IntArray masterNodes;
    IntArray slaveNodes;
    IR_GIVE_FIELD(ir, masterNodes, _IFT_ContactDefinitionNode2Node_MasterNodes);
    IR_GIVE_FIELD(ir, slaveNodes, _IFT_ContactDefinitionNode2Node_SlaveNodes);
    //this->epsN = 1.0e6;
    //IR_GIVE_OPTIONAL_FIELD(ir, this->epsN, _IFT_ContactDefinitionNode2Node_PenaltyN);
    
    Domain *domain = this->giveContactManager()->giveDomain();
    for( int i = 1; i<= masterNodes.giveSize(); i++ ) {
        ContactElement *master = new Node2NodeContact( domain->giveDofManager(masterNodes.at(i)),
                                                       domain->giveDofManager(slaveNodes.at(i)));

        this->addContactElement(master);
    }
    
    return IRRT_OK;
}




// Same version but with Lagrange multipliers
ContactDefinitionNode2NodeL :: ContactDefinitionNode2NodeL(ContactManager *cMan) : ContactDefinitionNode2Node(cMan)
{
    this->setNumberOfConstraintEqToAdd(1);
}
    

IRResultType
ContactDefinitionNode2NodeL :: initializeFrom(InputRecord *ir)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro
    
    IntArray masterNodes;
    IntArray slaveNodes;
    IR_GIVE_FIELD(ir, masterNodes, _IFT_ContactDefinitionNode2Node_MasterNodes);
    IR_GIVE_FIELD(ir, slaveNodes, _IFT_ContactDefinitionNode2Node_SlaveNodes);
    
    Domain *domain = this->giveContactManager()->giveDomain();
    for( int i = 1; i<= masterNodes.giveSize(); i++ ) {
        ContactElement *master = new Node2NodeContactL( domain->giveDofManager(masterNodes.at(i)),
                                                        domain->giveDofManager(slaveNodes.at(i)) );

        this->addContactElement(master);
    }
    
    
      
    return IRRT_OK;
}





}
