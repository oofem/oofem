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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef fieldmanager_h
#define fieldmanager_h

#include "field.h"

#include <map>

namespace oofem {

class FieldManager
{
protected:

    /**
     * Internal datastructure to keep reference (pointer) to registered field and
     * flag, indicating, whether the field pointer is managed by field manager or not.
     * In case by managed field, the field manager is assumed to own the field and is
     * responsible for its deallocation.
     */
    class fieldRecord {
    protected:
        Field *field;
        bool   isManaged;
    public:
        /// Creates new field record, containing reference to given field.
        fieldRecord (Field* f, bool managed): field(f), isManaged(managed) { }
        fieldRecord (): field(NULL), isManaged(false) { }
        /// Destructor. Deletes managed field.
        ~fieldRecord () {if (isManaged) delete field;}

        /// Return reference to field.
        Field* giveField() {return field;}
    };


    /**
     * Field container. Stores only pointers to objects (not object themselves)
     * to avoid copying elements and to preserve the use of polymorphic types.
     */
    std :: map< FieldType, fieldRecord* >externalFields;


public:
    FieldManager(): externalFields() { }
    ~FieldManager();

    /**
     * Registers the given field (the receiver is not assumed to own given field).
     * The field is registered under given key. Using this key, it can be later accessed.
     * If managedFlag set to true, the receiver is assumed to own the field, so it is
     * responsible for its deallocation).
     */
    void registerField(Field *eField, FieldType key, bool managedFlag = false);
    /** Returns true if field is registered under key */
    bool isFieldRegistered(FieldType key);
    /**
     * Returns the previously registered field under given key; NULL otherwise
     */
    Field *giveField(FieldType key);
    /**
     * Unregisters (deletes) the field registered under given key.
     */
    void unregisterField(FieldType key);
};
} // end namespace oofem
#endif // fieldmanager_h
