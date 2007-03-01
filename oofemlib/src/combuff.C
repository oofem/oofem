/* $Header: /home/cvs/bp/oofem/oofemlib/src/combuff.C,v 1.8.4.1 2004/04/05 15:19:43 bp Exp $ */
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


#ifdef __PARALLEL_MODE

#ifndef __MAKEDEPEND
#include <stdlib.h>
// include string.h for memmove
#include <string.h>  
#endif
#include "compiler.h"
#include "combuff.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "error.h"

#ifdef __USE_MPI
CommunicationBuffer :: CommunicationBuffer(MPI_Comm comm, int size, bool dynamic) 
{
 this->size = 0;
 curr_pos = 0;
 communicator = comm;
 buff = NULL;
 isDynamic = dynamic;
 request = MPI_REQUEST_NULL;

 this->resize (size);
}


CommunicationBuffer :: CommunicationBuffer (MPI_Comm comm, bool dynamic)
{
 this->size = 0;
 curr_pos = 0;
 communicator = comm;
 buff = NULL;
 isDynamic = dynamic;
 request = MPI_REQUEST_NULL;
}

#endif

CommunicationBuffer :: ~CommunicationBuffer ()
{
 if (buff) delete buff;
}


int 
CommunicationBuffer :: resize (int newSize)
{
 // do not shrink
 if (size >= newSize) return 1;


/*
 // first try realloc current buffer
 this->buff = (ComBuff_BYTE_TYPE*) 
  realloc (buff, newSize*sizeof (ComBuff_BYTE_TYPE));
 if (this->buff == NULL) {
  // realloc failed -> memory error
   OOFEM_ERROR ("CommunicationBuffer::resize : resize failed");
 }
 size = newSize;
*/

 ComBuff_BYTE_TYPE *newBuff;
 // allocate new memory
 if ((newBuff = (ComBuff_BYTE_TYPE*) 
      malloc (newSize*sizeof (ComBuff_BYTE_TYPE))) == NULL) {
   // alloc failed -> memory error
   OOFEM_ERROR ("CommunicationBuffer :: resize failed");
 }
 // copy old buffer into new one
 memmove (newBuff, this->buff, curr_pos);
 // dealocate old buffer
 free (this->buff);
 this->buff = newBuff;
 size = newSize;
 
 return 1;
}

void
CommunicationBuffer :: init ()
{
 curr_pos = 0;
 request = MPI_REQUEST_NULL;
}

#ifdef __USE_MPI

int 
CommunicationBuffer :: packIntArray (int* src, int n)
{
 int _size;
 // ask MPI for packing size for integer
 _size = this->giveIntVecPackSize (n);
 
 if ((this->curr_pos + _size > this->size))  {
   // reallocate itself
   if (isDynamic) {
     if (this->resize (this->curr_pos + _size + __CommunicationBuffer_ALLOC_CHUNK) == 0) return 0;
   } else { 
     OOFEM_WARNING ("CommunicationBuffer :: packIntArray: Resize requested in static mode");
     return 0;
   }
 }
 return (MPI_Pack (src, n, MPI_INT, this->buff, this->size, 
          &this->curr_pos, this->communicator) == MPI_SUCCESS);
}

int 
CommunicationBuffer :: packDoubleArray (double* src, int n)
{
 int _size;
 // ask MPI for packing size for double
 _size = this->giveDoubleVecPackSize (n);
 
 if ((this->curr_pos + _size > this->size))  {
   // reallocate itself
   if (isDynamic) {
     if (this->resize (this->curr_pos + _size + __CommunicationBuffer_ALLOC_CHUNK) == 0) return 0;
   } else {
     OOFEM_WARNING ("CommunicationBuffer :: packIntArray: Resize requested in static mode");
     return 0;
   }
 }
 return (MPI_Pack (src, n, MPI_DOUBLE, this->buff, this->size, 
          &this->curr_pos, this->communicator) == MPI_SUCCESS);
}

int
CommunicationBuffer :: packIntArray (const IntArray &arry)
{
 return arry.packToCommBuffer (*this);
}

int 
CommunicationBuffer::packFloatArray (const FloatArray &arry)
{
 return arry.packToCommBuffer (*this);
}

int 
CommunicationBuffer::packFloatMatrix (const FloatMatrix &mtrx)
{
 return mtrx.packToCommBuffer (*this);
}

int 
CommunicationBuffer :: unpackIntArray (int* dest, int n)
{
 return (MPI_Unpack (this->buff, this->size, &this->curr_pos, 
           dest, n, MPI_INT, this->communicator) == MPI_SUCCESS);
}

int 
CommunicationBuffer :: unpackDoubleArray (double* dest, int n)
{
 return (MPI_Unpack (this->buff, this->size, &this->curr_pos, 
           dest, n, MPI_DOUBLE, this->communicator) == MPI_SUCCESS);
}

int 
CommunicationBuffer :: unpackIntArray (IntArray &arry) 
{ return arry.unpackFromCommBuffer (*this);}

int CommunicationBuffer :: unpackFloatArray (FloatArray &arry) 
{ return arry.unpackFromCommBuffer (*this);}

int CommunicationBuffer :: unpackFloatMatrix (FloatMatrix&mtrx)
{ return mtrx.unpackFromCommBuffer (*this);}


int 
CommunicationBuffer::iSend (int dest, int tag)
{
 return (MPI_Isend (this->buff, this->curr_pos, MPI_PACKED, dest, tag, 
           this->communicator, &this->request) == MPI_SUCCESS);
}


int 
CommunicationBuffer::iRecv (int source, int tag, int count)
{
 if (count) {
  if (count >= this->size)  {
   // reallocate itself
   if (this->resize (count) == 0) return 0;
  }
 }
 return (MPI_Irecv (this->buff,this->size, MPI_PACKED, source, tag, 
           this->communicator, &this->request) == MPI_SUCCESS);
}

int
CommunicationBuffer :: testCompletion ()
{
 int flag;
 MPI_Status status;
 
 MPI_Test (&this->request, &flag, &status);
 return flag;
}

int
CommunicationBuffer :: testCompletion (int &source, int& tag)
{
 int flag;
 MPI_Status status;
 
 MPI_Test (&this->request, &flag, &status);

 source = status.MPI_SOURCE;
 tag = status.MPI_TAG;

 return flag;
}

int
CommunicationBuffer :: bcast (int root)
{
 return MPI_Bcast (this->buff, this->size, MPI_PACKED, root, this->communicator);
}


int 
CommunicationBuffer :: giveIntVecPackSize(int size) 
{
 int requredSpace;
 MPI_Pack_size (size, MPI_INT, this->communicator, &requredSpace);
 return requredSpace;
}

int 
CommunicationBuffer :: giveDoubleVecPackSize(int size) 
{
 int requredSpace;
 MPI_Pack_size (size, MPI_DOUBLE, this->communicator, &requredSpace);
 return requredSpace;
}


#endif

#endif
