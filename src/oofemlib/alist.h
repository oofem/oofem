/* $Header: /home/cvs/bp/oofem/oofemlib/src/alist.h,v 1.2.4.1 2004/04/05 15:19:43 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/


//   ******************
//   *** CLASS LIST ***
//   ******************
 
#ifndef alist_h

#include "compiler.h"
#include "error.h"
#ifndef __MAKEDEPEND
#include <stdlib.h> // for NULL
#endif

class FEMComponent ;

/**
 Class implementing generic list (or more precisely array).
 It maintains the array of generic pointers to objects of type T.
 Since the generic pointer to FEMComponent can  point to 
 any derived class instance, this class can be used as list or array containing
 elements, nodes, materials, loads or load-time functions, which are represented by classes
 derived from base FEMComponent.
 This class maitains only the links (pointers) to particular objects, objects itselfs are not contained within 
 this array. They have to be created outside (in memory, usually on heap) and then their pointers can be added to
 array. This is sometimes called non-intrusive approach. When destructor is called, the pointed object 
 are DELETED.
 The links to particular objects in array are stored in pointer array, therefore the access to particular
 component is wery efficient. On the other hand, the resizing of array is relative time expensive (the whole
 existing pointer table must be transfered) and is recomended to set size of the array to the final size.
*/
template <class T> class AList
{
/*
   This class implements an array which contains elements, nodes, materials,
   loads or load-time functions.
 DESCRIPTION :
   The objects are stored in 'values', an array of poinerst to T of size
   'size'.
 TASKS :
   - storing (method 'put') and returning (method 'at') component ;
   - expanding itself, in order to accomodate more components.
*/


   protected:
   /// Array or list size (number of components to store)
      int             size ;
   /// Real allocated size (may be larger than size, to prevent often realocation)
   int             allocatedSize;
   /// List size allocation increment 
   int             sizeIncrement;
   /// Array of pointers to particular components
      T**  values ;

   public:
   /// Constructor - creates list (array) of  size s
    AList (int s = 0, int sizeIncrement = 0) ;                // constructor
   /// Destructor
    ~AList () ;                  // destructor

   /// Returns component at given index
      T*              at (int i)           {return values[i-1] ;}
      int             giveSize ()          {return size ;}
   /**
    Expands the receiver from its current size to newSize, in order to acco-
    modate new entries.
    */
      void            growTo (int newSize) ;
   /**
    Returns True if the receiver has a non-null i-th entry, else returns False.
    */
      int             includes (int) ;
   /// Returns True if receiver is empty
      int             isEmpty ()           {return (size==0) ;}
   /// Returns True if receiver is not empty
      int             isNotEmpty ()        {return (size!=0) ;}
   /// Prints the receiver on screen.
      void            printYourself () ;
   /**
    Stores anObject at position i. Enlarge the receiver if too small.
    Deletes the old value, if exist.
    */
      void            put (int i, T* anObject) ;
   /// forces receiver to be empty; objects are DELETED if deleteObjectflag is true (default) .
   void            clear (bool deleteObjectflag = true);
   /// Deleletes the object at i-th position.
   void            remove (int i);
   /**
    Unlinks the object a i-th position. The object is returned, and its entry is
    unlinked, so there will be no further reference to this object. Does not delete 
    the object, its pointer is returned.
    */
   T*        unlink (int i);
} ;


template <class T> AList<T> ::  AList (int s, int sizeIncrement)
   // Constructor : creates a list of size s.
{
   register int   i ;
   T** p ;

   allocatedSize = size = s ;
   if (size) {
      values = new T* [size] ;
      p      = values ;
      i      = size ;
      while (i--)
  *p++ = NULL ;}                              // initialize 'values'
   else
      values = NULL ;

  this->sizeIncrement = sizeIncrement;
}


template <class T> AList<T> :: ~AList ()
   // Destructor.
{
 this->clear (true);
}

template <class T> void
AList<T> :: clear (bool deleteObjectFlag)
{
   int i = size;

   if (size) {
     if (deleteObjectFlag) {
       while (i--)
	 delete (values[i]) ;
     }
     //      delete [size] values ;}
     delete[]  values ;
   }
  
   allocatedSize = size = 0;
   values = NULL ;
   
}

template <class T> void
AList<T> :: growTo (int newSize)
   // Expands the receiver from its current size to newSize, in order to acco-
   // modate new entries.
{
   register int i ;
   T **newValues,**p1,**p2 ;

#ifdef DEBUG
   if (newSize <= size) {
     OOFEM_WARNING3 ("AList::growTo : new list size (%d) not larger than current size (%d)", newSize,size) ;
   }
#endif

  if (newSize > allocatedSize) {

   this->allocatedSize = newSize + this->sizeIncrement;
   newValues = new T* [this->allocatedSize] ;
   p1        = values ;
   p2        = newValues ;
   for (i=0 ; i<size ; i++)
    *p2++ = *p1++ ;
   for (i=size ; i<this->allocatedSize ; i++)
    *p2++ = NULL ;

   if (values)
    //      delete [size] values ;
    delete[]  values ;

   values = newValues ;
   this->allocatedSize = newSize + this->sizeIncrement;
  } 

   size   = newSize ;
}



template <class T> int
AList<T> :: includes (int i)
   // Returns True if the receiver has a non-null i-th entry, else returns
   // False.
{
   if (i > size)
      return FALSE ;
   else
      return (values[i-1]!=NULL) ;
}



template <class T> void
AList<T> :: printYourself ()
   // Prints the receiver on screen.
{
   register int i ;

   printf ("list of components of size %d\n",size) ;
   for (i=1 ; i<=size ; i++) {
      if (values[i-1] == NULL)
  printf ("%d  Nil \n",i) ;
      else
  printf ("%d %lX\n",i,(long int)(values[i-1])) ;}
}


template <class T> void AList<T> :: put (int i, T* anObject)
   // Stores anObject at position i. Enlarge the receiver if too small.
{
   if (size < i)
      this -> growTo(i) ;
  // delete old value
  if (values[i-1]) delete values[i-1];
   values[i-1] = anObject ;
}

template <class T> void   AList<T> :: remove (int i)
{
 if (size < i) return;
 if (values[i-1]) {
  delete values[i-1];
  values[i-1] = NULL;
 }
}


template <class T> T*   AList<T> :: unlink (int i)
{
 if (size < i) return NULL;

 T* answer = values[i-1];
 values[i-1] = NULL;
 return answer;
}


#define alist_h
#endif



