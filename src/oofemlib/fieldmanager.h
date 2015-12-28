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

#ifndef fieldmanager_h
#define fieldmanager_h

#include "oofemcfg.h"
#include "field.h"

#include <map>
#include <vector>
#include <memory>

namespace oofem {
typedef std :: shared_ptr< Field > FM_FieldPtr;

class OOFEM_EXPORT FieldManager
{
protected:

    /**
     * Field container. Stores smart pointers to objects (not object themselves)
     * to avoid copying elements and to preserve the use of polymorphic types.
     * The use of shared_ptr is essential here, as some registered fields may be
     * ovned (and maintained) by emodel, some may be cretead on demand and thus
     * managed only by field manager.
     */
    std :: map< FieldType, std :: shared_ptr< Field > >externalFields;

public:
    FieldManager() : externalFields() { }
    ~FieldManager();

    /**
     * Registers the given field (the receiver is not assumed to own given field).
     * The field is registered under given key. Using this key, it can be later accessed.
     */
    void registerField(FM_FieldPtr eField, FieldType key);

    /**
     * Returns the previously registered field under given key; NULL otherwise
     */
    FM_FieldPtr giveField(FieldType key);

    /** Returns true if field is registered under key */
    bool isFieldRegistered(FieldType key);
    /**
     * Unregisters (deletes) the field registered under given key.
     */
    void unregisterField(FieldType key);
    /**
     * Returns list of registered field keys, which can be obtained by calling giveField.
     */
    std::vector<FieldType> giveRegisteredKeys();

};
} // end namespace oofem
#endif // fieldmanager_h
