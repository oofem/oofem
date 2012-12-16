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

#ifndef dynalist_h
#define dynalist_h

namespace oofem {
template< class T >class dynaList;
template< class T >class dynaListIterator;

template< class T >class listItem
{
public:
    T data;
    listItem< T > *next;
    listItem< T > *prev;
public:

    listItem(const T &value);
    ~listItem();
    listItem< T > *giveNext() { return next; }
    listItem< T > *givePrev() { return prev; }
    T &giveData() { return data; }
    T *giveDataPtr() { return & data; }

    friend class dynaList< T >;
};

/**
 * Container managing its elements as a doubly linked list.
 * The elements of a list may have any type T that is assignable, copyable
 * and comes with default constructor with no parameters.
 * A list does not provide random access. For example to access the fifth element, you must navigate
 * the first forth following the chain of elements. Thus accessing the arbitrary element is very slow.
 * THE INSERTING AND REMOVING ELEMENTS AT EACH POSITION IS FAST. The insertion and deletion takes always
 * constant time, no elements are moved, only internal pointers are manipulated.
 * To access elements in a list, you must use iterators. List provide efficient bidirectional iterator.
 */
template< class T >class dynaList
{
public:
    /// List iterator type.
    typedef dynaListIterator< T >iterator;
    /// Link type.
    typedef listItem< T > *linkType;

protected:
    linkType last; // last->next is head of list

public:
    /// Constructor, creates the empty list.
    dynaList() {
        T val;
        last = new listItem< T >(val);
        last->next = last;
        last->prev = last;
    }
    /// Destructor.
    ~dynaList();

    /**
     * Inserts at iterator position a copy of value and returns the position of new element.
     * @param position Where to insert value.
     * @param value Value to insert.
     */
    iterator insert(iterator position, const T &value);
    /**
     * Inserts a copy of value at the beginning
     * @param value Value to insert.
     */
    void pushFront(const T &value);
    /**
     * Inserts a copy of value at the end.
     * @param value Value to insert.
     */
    void pushBack(const T &value);
    /// Removes all elements (makes container empty)
    void clear();
    /**
     * Removes the element at iterator position and returns the position of next element.
     * @param position Position where to erase.
     * @return Iterator to next position.
     */
    iterator erase(iterator position);
    /// @return Size of receiver.
    int size();

    /// @return True is receiver is empty, otherwise false.
    bool isEmpty() { return ( last->next == last ); }
    /// @return Iterator for the first element.
    iterator begin() { return ( * last ).next; }
    iterator begin() const { return ( * last ).next; }
    /// @return Iterator for the position after the last element.
    iterator end() { return last; }
    iterator end() const { return last; }

protected:
    friend class dynaListIterator< T >;
};

/**
 * Bidirectional iterator for dynaList.
 */
template< class T >class dynaListIterator
{
protected:
    listItem< T > *node;

public:
    /// Constructor.
    dynaListIterator() { }
    /**
     * Constructor.
     * @param x Nodes that to iterate over.
     */
    dynaListIterator(listItem< T > *x) : node(x) { }
    /// Copy constructor.
    dynaListIterator(const dynaListIterator< T > &x) : node(x.node) { }

    /// Equality comparison operator.
    bool operator==(const dynaListIterator< T > &x) const { return node == x.node; }
    /// Inequality comparison operator.
    bool operator!=(const dynaListIterator< T > &x) const { return node != x.node; }
    /// Value access operator. Returns the element of the actual position.
    T & operator*() const { return ( * node ).data; }

    /// Lets the iterator step forward to the next element.
    dynaListIterator< T > &operator++() {
        node = ( * node ).next;
        return * this;
    }

    /// Lets the iterator step backward to the previous element.
    dynaListIterator< T > &operator--() {
        node = ( * node ).prev;
        return * this;
    }

    friend class dynaList< T >;
};


template< class T >listItem< T > :: listItem(const T &value) : data(value)
{
    next = 0;
}

template< class T >listItem< T > :: ~listItem()
{
    // delete data;
}

template< class T >dynaList< T > :: ~dynaList()
{
    this->clear();
    delete last;
}

template< class T >dynaListIterator< T >
dynaList< T > :: insert(iterator position, const T &value)
{
    listItem< T > *tmp = new listItem< T >(value);

    tmp->next = position.node;
    tmp->prev = position.node->prev;
    ( position.node->prev )->next = tmp;
    position.node->prev = tmp;
    return tmp;
}



template< class T >void dynaList< T > :: pushFront(const T &value)
{
    insert(begin(), value);
}

template< class T >void dynaList< T > :: pushBack(const T &value)
{
    insert(end(), value);
}


template< class T >void dynaList< T > :: clear()
{
    linkType cur = last->next;
    while ( cur != last ) {
        linkType tmp = cur;
        cur = cur->next;
        delete ( tmp );
    }

    last->next = last;
    last->prev = last;
}

template< class T >int dynaList< T > :: size()
{
    int size = 0;
    linkType cur = last->next;
    while ( cur != last ) {
        cur = cur->next;
        size++;
    }

    return size;
}



template< class T >dynaListIterator< T >
dynaList< T > :: erase(iterator position)
{
    linkType next_node = position.node->next;
    linkType prev_node = position.node->prev;
    prev_node->next = next_node;
    next_node->prev = prev_node;
    delete ( position.node );
    return iterator(next_node);
}
} // end namespace oofem
#endif // dynalist_h
