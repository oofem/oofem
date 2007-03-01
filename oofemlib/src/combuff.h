/* $Header: /home/cvs/bp/oofem/oofemlib/src/combuff.h,v 1.5 2003/04/06 14:08:23 bp Exp $ */
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

//
// class CommunicationBuffer
//

#ifndef combuff_h
#ifdef __PARALLEL_MODE

#include "parallel.h"

class IntArray;
class FloatArray;
class FloatMatrix;

#define __CommunicationBuffer_ALLOC_CHUNK 1024
/**
 Type with size equal to one byte (sizeof (ComBuff_BYTE_TYPE) should be 1).
 Communication buffer buffer member is of this type.
*/
typedef char ComBuff_BYTE_TYPE;
/**
 Class CommunicationBuffer provides absraction for comunication buffer.
 buffer is used as input or output buffer to various comunication 
 services provided by message parsing libraries. 
 It provides methods for buffer initialization and resizing, methods for packing and/or
 unpacking data to/from buffer. Multiple messages can be packed/unpacked into/from buffer.
 The services for packing/unpacking take care about multiple messages stored, 
 they maintain proper current buffer position for data inserting/retrieval.
 Interface to low level message parsing function is provided, allowing to
 send and receive buffer to selected destination.
 
 */

class CommunicationBuffer
{
 protected:
  
  /// Size and current position in buffer in bytes (sizeof(char)).
  int size, curr_pos;
  /// dynamic flag (if true, buffer can grow, but reallocation is needed)
  bool isDynamic;
  /// Buffer. Dynamically allocated.
  ComBuff_BYTE_TYPE* buff;
  
#ifdef __USE_MPI
 /// MPI_Communicator used for packing
 MPI_Comm communicator;
 /**
  MPI request handle. This value is used by some message parsing functions.
  engngcommunicator also assembles array of commbuff handles and wait for some
  completion (when  receiveing edata for example).
  */
 MPI_Request request;
#endif

public:
 
#ifdef __USE_MPI
 /// Constructor. Creeates buffer of given size, using given communicator for packing
 CommunicationBuffer (MPI_Comm comm, int size, bool dynamic=0);
 /// Constructor. Creeates empty buffer, using given communicator for packing
 CommunicationBuffer (MPI_Comm comm, bool dynamic=0);
#endif
 /// Destructor.
 virtual ~CommunicationBuffer ();

 /**
  Resizes buffer to given size. If bufeer size is to be enlarged, then previously packed data
  are kept in new buffer. Otherwise buffer is cleared using init service.
  Current implementation only performs buffer growing, request for size decrease is ignored
  to avoid realocation if further request for groving is encountered.
  @param newSize new buffer size in bytes.
  @return nonzero if succesfull.
  */
 int resize (int newSize);
 /**
  Initializes buffer to empty state. All packed data are lost.
 */
 void init ();
 
 /// Returns the  current buffer size.
 int giveSize () {return size;}
 /// Returns current buffer position.
 int givePosition () {return curr_pos;}

#ifdef __USE_MPI
 /**
  Returns associated MPI request handle
  */
 MPI_Request giveRequest () {return this->request;}
#endif

 /**@name Methods for datatype packing/unpacking to/from buffer */
 //@{
 /** 
  Packs single integer value into buffer. 
  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
  @return nonzero if succesfull
 */
 int packInt(int value) {return packIntArray (&value, 1);}
 /** 
  Packs single double value into buffer. 
  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
  @return nonzero if succesfull
 */
 int packDouble (double value) {return packDoubleArray (&value, 1);}
 /**
  Packs array od integer value into buffer. 
  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
  @param src adress of first value in memory
  @param n number of packed integers
  @return nonzero if succesfull
  */
 int packIntArray (int* src, int n);
 /**
  Packs given IntArray  value into buffer. 
  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
  @return nonzero if succesfull
  */
 int packIntArray (const IntArray &arry);
 /**
  Packs array od double value into buffer. 
  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
  @param src adress of first value in memory
  @param n number of packed doubles
  @return nonzero if succesfull
  */
 int packDoubleArray (double* src, int n);
 /**
  Packs given FloatArray  value into buffer. 
  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
  @return nonzero if succesfull
  */
 int packFloatArray (const FloatArray &arry);
 /**
  Packs given FloatMatrix  value into buffer. 
  Buffer is enlarged if isDynamic flag is set, but it requires memory allocation and dealocation.
  @return nonzero if succesfull
  */
 int packFloatMatrix (const FloatMatrix &mtrx);

 /** 
  Unpacks single integer value from buffer. 
  @return nonzero if succesfull
 */
 int unpackInt(int& value) {return unpackIntArray (&value, 1);}
 /** 
  Unpacks single double value from buffer. 
  @return nonzero if succesfull
 */
 int unpackDouble (double& value) {return unpackDoubleArray (&value, 1);}
 /**
  unpacks array od integer value from buffer. 
  @param dest adress of first value in memory, where to store values
  @param n number of unpacked integers
  @return nonzero if succesfull
  */
 int unpackIntArray (int* dest, int n);

 /**
  Unpacks given IntArray  value from buffer. 
  @return nonzero if succesfull
  */
 int unpackIntArray (IntArray &arry);
 /**
  Packs array od double value from buffer. 
  @param dest adress of destination of unpacked  values in memory
  @param n number of unpacked doubles
  @return nonzero if succesfull
  */
 int unpackDoubleArray (double* dest, int n);
 /**
  Unpacks given FloatArray  value from buffer. 
  @return nonzero if succesfull
  */
 int unpackFloatArray (FloatArray &arry);
 /**
  Unpacks given FloatMatrix  value from buffer. 
  @return nonzero if succesfull
  */
 int unpackFloatMatrix (FloatMatrix&mtrx);
 //@}


 /**@name Methods for determining pack size of datatype to pack/unpack to/from buffer */
 //@{
 /** 
  Returns pack size required to pack integer array (c-style).
  @param array size
  @return  pack size required
  */
 int giveIntVecPackSize(int size) ;
 /** 
  Returns pack size required to pack double array (c-style).
  @param array size
  @return  pack size required
  */
 int giveDoubleVecPackSize(int size) ;
 //@}

 /**@name Services for buffer sending/receiving */
 //@{
#ifdef __USE_MPI
 /**
  Starts standart mode, nonblocking send.
  @param dest rank of destination
  @param tag message tag
  @param communicator (handle)
  @return sends MIP_Succes if ok
  */
 int iSend (int dest, int tag);
 /**
  Starts standart mode, nonblocking receive. The buffer must be large enough to receive all data.
  @param source rank of source 
  @param tag message tag
  @param count number of elements to receive (bytes). Causes receive buffer to resize to count elements.
  If zero (default value) buffer is not resized.
  @param reguest communicator request (handle)
  @return MIP_Succes if ok
  */
 int iRecv (int source, int tag, int count = 0);
 /**
  Tests if the operation identified by this->request is complete. 
  In such case, true is returned and
  if communication was initiated by nonblocking send/receive, then request handle 
  is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
  @return true if operation complete, false otherwise.
  */
 int testCompletion () ;
 /**
  Tests if the operation identified by this->request is complete. 
  In such case, true is returned and
  if communication was initiated by nonblocking send/receive, then request handle 
  is set to MPI_REQUEST_NULL. Otherwise call returns flag=false.
  @param source contain the source tag
  @param tag contain the tag of received message
  @return true if operation complete, false otherwise.
  */
 int testCompletion (int &source, int& tag) ;
 /**
  Initalizes broadcast over colaborating processes.
  The whole buffer size is broadcasted. All buffers participating in broadcast
  should have the same size.
  @param root rank of broadcast root
  @return MIP_Succes if ok
  */
 int bcast (int root);
#endif
 //@}
};

#define combuff_h
#endif
#endif
