/* $Header:$ */
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
// class CommunicationPacket
//

#ifndef compacket_h
#ifdef __PARALLEL_MODE

#include "parallel.h"
#include "combuff.h"

class IntArray;
class FloatArray;
class FloatMatrix;

#define __CommunicationPacket_DEFAULT_SIZE 10240

/**
 Class CommunicationPacket represent a data-packet, that is used to implement dynamic 
 communicator. Dynamic Communicator can pack messages into a dynamic message.
 This dynamic message is splitted into a series 
 of data packets of fixed size (this is not necessary) that are send over network.

 A special header is put at the begining of each packet buffer. This header keeps the message number 
 as well as the EOF flag indicating the last packet in message. This header is packed at the begining 
 of each packet.
*/

class CommunicationPacket : public CommunicationBuffer
{
 protected:
  int number;
  bool EOF_Flag;

public:
 
#ifdef __USE_MPI
 /// Constructor. Creeates buffer of given size, using given communicator for packing
 CommunicationPacket (MPI_Comm comm, int size);
 /// Constructor. Creeates empty buffer, using given communicator for packing
 CommunicationPacket (MPI_Comm comm);
#endif
 /// Destructor.
 ~CommunicationPacket ();

 /**
  Initializes buffer to empty state. All packed data are lost.
 */
 virtual void init ();
 
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
#endif
 //@}
 
 void setNumber(int _num) {this->number=_num;}
 void setEOFFlag() {this->EOF_Flag = true;}
 int getNumber() {return number;}
 bool isEOFFlag() {return EOF_Flag;}
 

 protected:
 /// packs packet header info at receiver beginning
 int packHeader ();
 int unpackHeader ();
};

#define compacket_h
#endif
#endif
