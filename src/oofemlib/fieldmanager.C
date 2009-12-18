/* $Header: /home/cvs/bp/oofem/oofemlib/src/fieldmanager.C,v 1.1 2003/04/06 14:08:24 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#include "fieldmanager.h"

namespace oofem {

void
FieldManager :: registerField(Field *eField, FieldType key)
{
    std :: map< FieldType, Field * > :: iterator i;
    if ( ( i = this->externalFields.find(key) ) != this->externalFields.end() ) {
        // delete old entry, since map contains only pointers, not fields themselves
        delete ( * i ).second;
    }

    this->externalFields [ key ] = eField;
}

bool
FieldManager :: isFieldRegistered(FieldType key)
{
    return ( this->externalFields.find(key) != this->externalFields.end() );
}

Field *
FieldManager :: giveField(FieldType key)
{
    std :: map< FieldType, Field * > :: iterator i;
    if ( ( i = this->externalFields.find(key) ) == this->externalFields.end() ) {
        return NULL;
    }

    return ( * i ).second;
}

void
FieldManager :: unregisterField(FieldType key)
{
    std :: map< FieldType, Field * > :: iterator i;
    i = this->externalFields.find(key);
    if ( i == this->externalFields.end() ) {
        return;
    }

    delete ( * i ).second;
    this->externalFields.erase(i);
}

} // end namespace oofem
