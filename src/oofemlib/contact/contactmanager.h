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

#ifndef contactmanager_h
#define contactmanager_h

#include "chartype.h"
#include "valuemodetype.h"
#include "oofemcfg.h"
#include "datareader.h"
#include "inputrecord.h"

#include <vector>
#include <memory>

///@name Input fields for _IFT_ContactManager
//@{
#define _IFT_ContactManager_Name "contactmanager"
#define _IFT_ContactManager_NumberOfContactDefinitions "numcontactdef"

//@}

namespace oofem {
class Domain;
class ContactDefinition;
class TimeStep;
class UnknownNumberingScheme;
class SparseMtrx;
/**
 * This class manages all the contacts in a domain 
 *
 * @author Jim Brouzoulis
 */
class OOFEM_EXPORT ContactManager
{
protected:
    Domain *domain;

private:
    std :: vector< std :: unique_ptr< ContactDefinition > > contactDefinitionList;

public:

    /// Constructor.
    ContactManager(Domain *domain);
    /// Destructor.
    virtual ~ContactManager();

    ContactManager(const ContactManager& src) = delete;
    ContactManager &operator = (const ContactManager &src) = delete;

    void createContactDofs();
    
    /// Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const { return "ContactManager"; }

    Domain *giveDomain() { return this->domain; }
    int numberOfContactDefinitions;
    ContactDefinition *giveContactDefinition(const int num) { return this->contactDefinitionList[num-1].get(); }
    int giveNumberOfContactDefinitions() const { return (int)contactDefinitionList.size(); }

    
    void assembleVectorFromContacts(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                    const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms = NULL);
    

    void assembleTangentFromContacts(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

};

} // end namespace oofem
#endif // contactmanager_h
