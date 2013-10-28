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

//#define FIELDMANAGER_USE_SHARED_PTR
#include "oofemcfg.h"
#include "field.h"

#include <map>
#ifdef FIELDMANAGER_USE_SHARED_PTR
#include <tr1/memory>
#endif

namespace oofem {

#ifdef FIELDMANAGER_USE_SHARED_PTR
typedef std::tr1::shared_ptr<Field> FM_FieldPtr ;
#else
typedef Field* FM_FieldPtr ;
#endif

class OOFEM_EXPORT FieldManager
{
protected:

#ifdef FIELDMANAGER_USE_SHARED_PTR
  /**
   * Field container. Stores smart pointers to objects (not object themselves)
   * to avoid copying elements and to preserve the use of polymorphic types.
   * The use of shared_ptr is essential here, as some registered fields may be
   * ovned (and maintained) by emodel, some may be cretead on demand and thus 
   * managed only by field manager.
   */
  std :: map<FieldType, std::tr1::shared_ptr<Field> > externalFields;
#else  
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
#endif

public:
    FieldManager(): externalFields() { }
    ~FieldManager();

#ifdef FIELDMANAGER_USE_SHARED_PTR
    /**
     * Registers the given field (the receiver is not assumed to own given field).
     * The field is registered under given key. Using this key, it can be later accessed.
     */
    void registerField(FM_FieldPtr eField, FieldType key);
#else
    /**
     * Registers the given field (the receiver is not assumed to own given field).
     * The field is registered under given key. Using this key, it can be later accessed.
     * If managedFlag set to true, the receiver is assumed to own the field, so it is
     * responsible for its deallocation).
     */
    void registerField(FM_FieldPtr eField, FieldType key, bool managedFlag = false);
#endif

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
    
};
} // end namespace oofem
#endif // fieldmanager_h
