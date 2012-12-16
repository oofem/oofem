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

#ifndef alist_h
#define alist_h

#include "error.h"
#ifndef __MAKEDEPEND
 #include <cstddef>
#endif


namespace oofem {

/**
 * Class implementing generic list (or more precisely array).
 * It maintains the array of generic pointers to objects of type T using 1-based numbering.
 *
 * This class maintains only the links (pointers) to particular objects, objects themselfs are not contained within
 * this array. They have to be created outside (in memory, usually on heap) and then their pointers can be added to
 * array. This is sometimes called non-intrusive approach. When destructor is called, the linked objects
 * are <em>deleted</em>. To prevent the deletion, objects should be unlinked before deconstructor is called.
 *
 * The links to particular objects in array are stored in pointer array, therefore the access to particular
 * component is very efficient. On the other hand, the resizing of array is relative time expensive (the whole
 * existing pointer table must be transfered) and is recommended to set size of the array to the final size.
 */
template< class T >class AList
{
protected:
    /// Array or list size (number of components to store).
    int size;
    /// Real allocated size (may be larger than size, to prevent often reallocation).
    int allocatedSize;
    /// List size allocation increment.
    int sizeIncrement;
    /// Array of pointers to particular components.
    T **values;

public:
    /**
     * Creates list of size s.
     * @param s Initial size.
     * @param sizeIncrement Increases array in blocks of this size when necessary.
     */
    AList(int s = 0, int sizeIncrement = 0);
    /// Destructor
    ~AList();

    /**
     * @param i Index of object.
     * @return Object at given index.
     */
    T *at(int i) const;
    /**
     * Expands the receiver from its current size to newSize, in order to accommodate new entries.
     * @param newSize Size that receiver must fit.
     */
    void growTo(int newSize);
    /**
     * Checks if receiver includes object as given position.
     * @param i Index of object.
     * @return true if object is non-null at i-th entry, otherwise false.
     */
    bool includes(int i) const;
    /// @return Size of array.
    int giveSize() const { return size; }
    /// @return True if receiver is empty, otherwise false.
    int isEmpty() const { return ( size == 0 ); }
    /// @return True if receiver is not empty, otherwise false.
    int isNotEmpty() const { return ( size != 0 ); }
    /// Prints the receiver on screen.
    void printYourself() const;
    /**
     * Stores anObject at position i.
     * Enlarges the receiver if too small and delates the old value if it exists.
     * @param i Index to put object
     * @param anObject Object to put as position i.
     */
    void put(int i, T *anObject);
    /**
     * Clears receiver.
     * Objects are deleted if depending on the delete flag.
     * @param deleteObjectFlag Shows if objects should be deleted.
     */
    void clear(bool deleteObjectFlag = true);
    /// Deletes the object at i-th position.
    void remove(int i);
    /**
     * Unlinks the object a i-th position.
     * The object is returned, and its entry is
     * unlinked, so there will be no further reference to this object. Does not delete
     * the object, its pointer is returned.
     * @param i Index where to unlink.
     * @return Object at index i.
     */
    T *unlink(int i);
};

template< class T >AList< T > ::  AList(int s, int sizeIncrement)
// Constructor : creates a list of size s.
{
    register int i;
    T **p;

    allocatedSize = size = s;
    if ( size ) {
        values = new T * [ size ];
        p      = values;
        i      = size;
        while ( i-- ) {
            * p++ = NULL;
        }
    } // initialize 'values'
    else {
        values = NULL;
    }

    this->sizeIncrement = sizeIncrement;
}


template< class T >AList< T > :: ~AList()
{
    this->clear(true);
}

template< class T >void
AList< T > :: clear(bool deleteObjectFlag)
{
    int i = size;

    if ( size ) {
        if ( deleteObjectFlag ) {
            while ( i-- ) {
                delete ( values [ i ] );
            }
        }

        //      delete [size] values ;}
        delete[]  values;
    }

    allocatedSize = size = 0;
    values = NULL;
}

template< class T >void
AList< T > :: growTo(int newSize)
// Expands the receiver from its current size to newSize, in order to accommodate new entries.
{
    register int i;
    T **newValues, **p1, **p2;


    if ( newSize < size ) {
#ifdef DEBUG
        OOFEM_WARNING3("AList::growTo : new list size (%d) not larger than current size (%d)", newSize, size);
#endif
        // delete entities in indexes in the range (newSize, size)
        i = size;
        if ( size ) {
            while ( ( --i ) >= newSize ) {
                delete ( values [ i ] );
                values [ i ] = NULL;
            }
        }
    } else if ( newSize > allocatedSize ) {
        this->allocatedSize = newSize + this->sizeIncrement;
        newValues = new T * [ this->allocatedSize ];
        p1        = values;
        p2        = newValues;
        for ( i = 0; i < size; i++ ) {
            * p2++ = * p1++;
        }

        for ( i = size; i < this->allocatedSize; i++ ) {
            * p2++ = NULL;
        }

        if ( values ) {
            //      delete [size] values ;
            delete[]  values;
        }

        values = newValues;
        this->allocatedSize = newSize + this->sizeIncrement;
    }

    size = newSize;
}

template< class T >T*
AList< T > :: at(int i) const
{
#ifdef DEBUG
    if ( i <= 0 ) {
        OOFEM_ERROR2("AList :: at - Asking for negative or zero indices (%d)", i);
    }
#endif
    return values [ i - 1 ];
}

template< class T >bool
AList< T > :: includes(int i) const
// Returns True if the receiver has a non-null i-th entry, else returns
// False.
{
#ifdef DEBUG
    if ( i <= 0 ) {
        OOFEM_ERROR2("AList :: includes - Asking for negative or zero indices (%d)", i);
    }
#endif
    if ( i > size ) {
        return false;
    } else {
        return ( values [ i - 1 ] != NULL );
    }
}

template< class T >void
AList< T > :: printYourself() const
// Prints the receiver on screen.
{
    printf("List of components of size %d\n", size);
    for (int i = 1; i <= size; i++ ) {
        printf("%d : %p\n",i, values [ i - 1 ] );
    }
}

template< class T >void AList< T > :: put(int i, T *anObject)
// Stores anObject at position i. Enlarge the receiver if too small.
{
#ifdef DEBUG
    if ( i <= 0 ) {
        OOFEM_ERROR2("AList :: put - Trying to write to zero or negative indices (%d)", i);
    }
#endif
    if ( size < i ) {
        this->growTo(i);
    }

    // delete old value
    if ( values [ i - 1 ] ) {
        delete values [ i - 1 ];
    }

    values [ i - 1 ] = anObject;
}

template< class T >void AList< T > :: remove(int i)
{
#ifdef DEBUG
    if ( i < 0 ) {
        OOFEM_ERROR2("AList :: remove - Trying to remove at zero or negative indices (%d)", i);
    }
#endif

    if ( size < i ) {
        return;
    }

    if ( values [ i - 1 ] ) {
        delete values [ i - 1 ];
        values [ i - 1 ] = NULL;
    }
}


template< class T >T *AList< T > :: unlink(int i)
{
#ifdef DEBUG
    if ( i <= 0 ) {
        OOFEM_ERROR2("AList :: unlink - Trying to unlink at zero or negative indices (%d)", i);
    }
#endif
    if ( size < i ) {
        return NULL;
    }

    T *answer = values [ i - 1 ];
    values [ i - 1 ] = NULL;
    return answer;
}
} // end namespace oofem
#endif // alist_h



