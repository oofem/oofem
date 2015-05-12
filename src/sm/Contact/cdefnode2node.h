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

#ifndef contactdefinitionnode2node_h
#define contactdefinitionnode2node_h

#include "contact/contactdefinition.h"


///@name Input fields for _IFT_ContactDefinitionNode2Node
//@{
#define _IFT_ContactDefinitionNode2Node_Name "cdef_node2nodep"
#define _IFT_ContactDefinitionNode2Node_MasterNodes "masternodes"
#define _IFT_ContactDefinitionNode2Node_SlaveNodes "slavenodes"
#define _IFT_ContactDefinitionNode2Node_PenaltyN "penaltyn"

#define _IFT_ContactDefinitionNode2NodeL_Name "cdef_node2nodeL"
//@}

namespace oofem {
class Domain;
class ContactManager;
class ContactObject;
class ContactElement;

class ContactMaterial; // write this


/**
 * This class manages a particular contact definition. 
 * This keeps track of the discretization, how the contact constraints are enforced and so on
 *
 * @author Jim Brouzoulis
 */
class OOFEM_EXPORT ContactDefinitionNode2Node : public ContactDefinition
{
private:
    double epsN;
    double epsT; // these should be 'contactmaterial' par
    
public:

    /// Constructor.
    ContactDefinitionNode2Node(ContactManager *cMan);
    /// Destructor.
    virtual ~ContactDefinitionNode2Node(){};

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "ContactDefinitionNode2Node"; }
    virtual const char *giveInputRecordName() const { return _IFT_ContactDefinitionNode2Node_Name; }
 
};





/**
 * This class manages a none to node contact definition with enforcement using Lagrange multipliers.
 *
 * @author Jim Brouzoulis
 */
class OOFEM_EXPORT ContactDefinitionNode2NodeL : public ContactDefinitionNode2Node
{

public:

    /// Constructor.
    ContactDefinitionNode2NodeL(ContactManager *cMan);
    /// Destructor.
    virtual ~ContactDefinitionNode2NodeL(){};

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "ContactDefinitionNode2NodeL"; }
    virtual const char *giveInputRecordName() const { return _IFT_ContactDefinitionNode2NodeL_Name; }

};


} // end namespace oofem
#endif // contactdefinitionnode2node_h
