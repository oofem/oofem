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
// #include <stdlib.h>
#endif
#include "compiler.h"
#include "compacket.h"
#include "error.h"

#ifdef __USE_MPI
CommunicationPacket :: CommunicationPacket(MPI_Comm comm, int size, int num) : CommunicationBuffer (comm, size, false)
{
  this->EOF_Flag = false;
  this->number = num;
  // reserve space for packet header
  curr_pos += giveIntVecPackSize(2);
}


CommunicationPacket :: CommunicationPacket (MPI_Comm comm, int num) : CommunicationBuffer (comm, __CommunicationPacket_DEFAULT_SIZE, false);
{
  this->EOF_Flag = false;
  this->number = num;
  // reserve space for packet header
  curr_pos += giveIntVecPackSize(2);
}

#endif

CommunicationPacket :: ~CommunicationPacket ()
{}


void
CommunicationPacket :: init ()
{
  CommunicationBuffer::init();
  // reserve space for packet header
  curr_pos += giveIntVecPackSize(2);
  
}

int 
CommunicationPacket::iSend (int dest, int tag)
{
  this->packHeader();
  return (MPI_Isend (this->buff, this->curr_pos, MPI_PACKED, dest, tag, 
                     this->communicator, &this->request) == MPI_SUCCESS);
}


int 
CommunicationPacket::iRecv (int source, int tag, int count)
{
  if (count) {
    if (count >= this->size)  {
      // reallocate itself
      if (this->resize (count) == 0) return 0;
    }
  }
  return (MPI_Irecv (this->buff,this->size, MPI_PACKED, source, tag, 
                     this->communicator, &this->request) == MPI_SUCCESS);
  this->unpackHeader();
}

#endif

int 
CommunicationPacket::packHeader()
{
  int _arry(2);
  int _pos  =0;
  
  _arry[0]=this->number;
  _arry[1]=this->EOF_Flag;

  return (MPI_Pack (_arry, 2, MPI_INT, this->buff, size, &_pos, this->communicator) == MPI_SUCCESS);
}

int
CommunicationPacket::unpackHeader()
{
  int _arry(2);
  int _res, _pos  =0;

  _res = MPI_Unpack (this->buff, this->size, &_pos, _arry, 2, MPI_INT, this->communicator);
  this->number =   _arry[0];
  this->EOF_Flag = _arry[1];

  return (_res == MPI_SUCCESS);
}

#endif
