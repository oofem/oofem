/* $Header: /home/cvs/bp/oofem/oofemlib/src/intarray.C,v 1.9.4.1 2004/04/05 15:19:43 bp Exp $ */
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
/*
 The original idea for this class comes from 
  Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 PhD Thesis, EPFL, Lausanne, 1992.
*/


//   file INTARRAY.CC

#include "intarray.h"
#include "freestor.h"
#include "cltypes.h"
#include "error.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

#ifdef __PARALLEL_MODE
#include "combuff.h"
#endif

IntArray :: IntArray (int n)
   // Constructor : creates an array of size n (filled with garbage).
{
   allocatedSize = size = n ;
   if (n)
      values = allocInt(size) ;
   else
      values = NULL ;
}


IntArray::IntArray (const IntArray& src)
{
 // copy constructor
 allocatedSize = size = src.size ;
 if (size) values = allocInt(size); else values = NULL;
 
 int *p2        = src.values ;
  int *p3        = values ;

 int i = size ;
 while (i--)
  *p3++ = *p2++ ;

}


IntArray& 
IntArray :: operator=  (const IntArray& src)
{
 // assignment: cleanup and copy
 if (values) freeInt(values);

 
 allocatedSize = size = src.size ;
 if (size) {
  values = allocInt(size) ;
 
  int *p2        = src.values ;
  int *p3        = values ;
 
  int i = size ;
  while (i--)
   *p3++ = *p2++ ;
 } else { values = NULL;}
 return *this;
}


void
IntArray :: zero () 
{
  int *p1        = values ;
 
 int i = size ;
 while (i--)
  *p1++ = 0;
}


void 
IntArray :: add (int value)
{
  int *p1        = values ;
  
  int i = size ;
  while (i--) {
    *p1+= value;
    p1++;
  }
} 


#ifdef DEBUG
int&  IntArray :: at (int i)
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}

int  IntArray :: at (int i) const
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}
#endif


#ifdef DEBUG
void  IntArray :: checkBounds (int i) const
   // Checks that the receiver includes an index i.
{
  if (i<0) OOFEM_ERROR2 ("IntArray::checkBounds : array error on index : %d < 0",i) ;
  if (i>size) 
    OOFEM_ERROR3 ("IntArray::checkBounds : array error on index : %d > %d",i,size) ;
}
#endif

void 
IntArray :: resize (int n, int allocChunk) 
{
 int *p1, *p2, *newValues, i;

 if (n <= allocatedSize) {size = n; return;}

 if (allocChunk < 0) allocChunk = 0;
 newValues = allocInt(n + allocChunk) ;

 p1 = values ;
 p2 = newValues ;
 i  = size ;
 while (i--)
  *p2++ = *p1++ ;
 
 if (values)
  freeInt (values) ;
 values = newValues ;
 allocatedSize = n + allocChunk;
 size   = n ;
}


void  IntArray :: followedBy (const IntArray& b, int allocChunk)
   // Appends the array 'b' the receiver. Returns the receiver.
{
 register int i ;
 int          newSize ;
 int          *newValues,*p1,*p2,*p3 ;
 
 newSize = size + b.size ;
 if (newSize == size)
  return ;

 if (allocChunk < 0) allocChunk = 0;

 if (newSize > allocatedSize) {
  newValues = allocInt(newSize+allocChunk) ;
  p3        = newValues ;

  p1        = values ;
  p2        = b.values ;
  p3        = newValues ;
 
  i = size ;
  while (i--)
   *p3++ = *p1++ ;
  
  i = b.size ;
  while (i--)
   *p3++ = *p2++ ;
  
  if (values)
   freeInt (values) ;
  values = newValues ;
  allocatedSize = newSize+allocChunk;
  size   = newSize ;
  
 } else {
  
  p1        = values+size ;
  p2        = b.values ;
  
  i = b.size ;
  while (i--)
   *p1++ = *p2++ ;
  
  size   = newSize;
 }
}



void  IntArray :: followedBy (const int b, int allocChunk)
   // Appends the array 'b' the receiver. Returns the receiver.
{
 register int i ;
 int          newSize ;
 int          *newValues,*p1,*p3 ;
 
 newSize = size + 1 ;

 if (newSize > allocatedSize) {
  newValues = allocInt(newSize+allocChunk) ;

  p1        = values ;
  p3        = newValues ;
 
  i = size ;
  while (i--)
   *p3++ = *p1++ ;
  
    *p3 = b ;
  
  if (values)
   freeInt (values) ;
  values = newValues ;
  allocatedSize = newSize + allocChunk;
  size   = newSize ;
  
 } else {
  
  *(values+size) = b ;
  size   = newSize;
 }
}
int
IntArray::containsOnlyZeroes () const 
{
  for (int i=0; i<size; i++)
    if (values[i]) return 0;
  return 1;
}

void  IntArray :: printYourself () const
   // Prints the receiver on screen.
{
   printf ("IntArray of size : %d\n",size) ;
   for (int i=1 ; i<=size ; ++i) {
      if (i>15) {
  printf ("   (other components not printed)") ;
  break ;}
      else
  printf ("%d  ",this->at(i)) ;}
   printf ("\n") ;
}


contextIOResultType  IntArray :: storeYourself (FILE* stream) const
// writes receiver's binary image into stream
// use id to distinguish some instances
// return value >0 succes
//              =0 file i/o error
{
  int type_id = IntArrayClass;
  // write class header
  if (fwrite(&type_id,sizeof(int),1,stream)!=1)  return (CIO_IOERR);
  // write size 
  if (fwrite(&size,sizeof(int),1,stream)!=1) return (CIO_IOERR);
  // write raw data
  if (fwrite(values,sizeof(int),size,stream)!= (size_t)size) return (CIO_IOERR);
  // return result back
  return CIO_OK;
}

contextIOResultType  IntArray :: restoreYourself (FILE* stream)
// reads receiver from stream
// warning - overwrites existing data!
// returns 0 if file i/o error
//        -1 if id od class id is not correct
{
  int class_id;
  // read class header
  if (fread(&class_id,sizeof(int),1,stream) != 1)  return (CIO_IOERR);
  if (class_id != IntArrayClass) return CIO_BADVERSION;
  // read size 
  if (fread(&size,sizeof(int),1,stream) != 1)  return (CIO_IOERR);
  if (values!=NULL) freeInt(values);
  if (size)
    values = allocInt(size) ;
  else
    values = NULL ;
  // write raw data
  if (fread(values,sizeof(int),size,stream) != (size_t)size) return (CIO_IOERR);
  // return result back
  return CIO_OK;
}


int
IntArray :: findFirstIndexOf (int value)   const
{
 // finds index of value in receiver
 // if such value  does not exists, returns zero index
 int i;
  for (i=0 ; i<size ; i++) {
  if (values[i] == value) return i+1;
 }
 // nothing found
 return 0;
}



#ifdef __PARALLEL_MODE
int 
IntArray :: packToCommBuffer (CommunicationBuffer& buff) const
{
 int result = 1;
 // pack size
 result &= buff.packInt (size);
 // pack data
 result &= buff.packIntArray (this->values, size);
 
 return result;
}

int 
IntArray :: unpackFromCommBuffer (CommunicationBuffer& buff)
{
 int newSize, result = 1;
 // unpack size
 result &= buff.unpackInt (newSize);
 // resize yourself
 this->resize(newSize);
 result &= buff.unpackIntArray (this->values, newSize);

 return result;
}

int 
IntArray :: givePackSize (CommunicationBuffer& buff)
{
 return buff.giveIntVecPackSize (1) + buff.giveIntVecPackSize (this->size);
}
#endif
