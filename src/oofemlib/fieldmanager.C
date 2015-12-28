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

#include "fieldmanager.h"

namespace oofem {
FieldManager :: ~FieldManager()
{ }

void
FieldManager :: registerField(std :: shared_ptr< Field >eField, FieldType key)
{
    this->externalFields.insert({key, eField});
}




bool
FieldManager :: isFieldRegistered(FieldType key)
{
    return ( this->externalFields.find(key) != this->externalFields.end() );
}

std :: shared_ptr< Field >
FieldManager :: giveField(FieldType key)
{
    auto i = this->externalFields.find(key);
    if ( i == this->externalFields.end() ) {
        std :: shared_ptr< Field >p; // std::shared_ptr<Field> p(nullptr);
        return p;
    }

    return i->second;
}

void
FieldManager :: unregisterField(FieldType key)
{
    auto i = this->externalFields.find(key);
    if ( i == this->externalFields.end() ) {
        return;
    }

    this->externalFields.erase(i);
}

std::vector<FieldType>
FieldManager :: giveRegisteredKeys()
{
    std::vector<FieldType> ret;
    for(const auto& keyField: this->externalFields) ret.push_back(keyField.first);
    return ret;
}

} // end namespace oofem
