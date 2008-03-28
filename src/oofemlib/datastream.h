/* $Header: $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2006   Borek Patzak                                       



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

#ifndef datastream_h
#define datastream_h

#include "combuff.h"
#include "processcomm.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

/**
   The purpose of DataStream abstract class is to allow to store/restore context to different streams, 
   including file, communication buffers, etc., using the same routine. 
   This will facilitate many algorithms relying on saving/moving state of componets 
   (such as load balancing), without writting new (and wery similar) routines. 
   This  will lead to a  better consistency of code.
 */
class DataStream {
 public:
  /// Destructor
  virtual ~DataStream() {}
  /** @name Data Stream reading methods.
      These methods read "count" values from data stream into 
      array passed as the first argument. 
      All functions return nonzero if successful.
  */
  //@{
  /// Reads count integer values into array pointed by data. 
  virtual int read (int* data, const unsigned int count) = 0; 
  /// Reads count unsigned long values into array pointed by data. 
  virtual int read (unsigned long* data, const unsigned int count) = 0; 
  /// Reads count double values into array pointed by data. 
  virtual int read (double* data, const unsigned int count) = 0;
  /// Reads count char values into array pointed by data. 
  virtual int read (char* data, const unsigned int count) = 0;
  //@}

  /** @name Data Stream writting methods.
      These methods write "count" values of data into stream.
      All functions return nonzero if successful.
  */
  //@{
  /// Writes count int values from array pointed by data. 
  virtual int write (const int* data, const unsigned int count) = 0;
  /// Writes count  values from array pointed by data. 
  virtual int write (const unsigned long* data, const unsigned int count) = 0;
  /// Writes count char values from array pointed by data. 
  virtual int write (const double* data, const unsigned int count) = 0;
  /// Writes count char values from array pointed by data. 
  virtual int write (const char* data, const unsigned int count) = 0;
  //@}
};


/**
   Implementation of FileDataStream representing DataStream interface to file i/o.
   This class creates a Datastream shell around c file i/o routines. This class will 
   not provide any methods for opening/closing file. This is the responsibility of user.  
   @see DataStream class.
*/
class FileDataStream : public DataStream {
 private:
  /// FILE pointer of associated stream
  FILE* stream;
 public:
  /// Constructor, takes associated stream pointer as parameter
  FileDataStream (FILE* s) {stream = s;}
  /// Destructor (will not close stream!)
  ~FileDataStream () {}

  /** @name Data Stream reading methods.
      These methods read "count" values from data stream into 
      array passed as the first argument. 
      All functions return nonzero if successful.
  */
  //@{
  virtual int read (int* data, const unsigned int count) 
    {if (fread (data, sizeof(int), count, stream)== count) return 1; else return 0;}
  virtual int read (unsigned long* data, const unsigned int count) 
    {if (fread (data, sizeof(unsigned long), count, stream)== count) return 1; else return 0;}
  virtual int read (double* data, const unsigned int count) 
    {if (fread (data, sizeof(double), count, stream)== count) return 1; else return 0;}
  virtual int read (char* data, const unsigned int count) 
    {if (fread (data, sizeof(char), count, stream)== count) return 1; else return 0;}
  //@}
  /** @name Data Stream writting methods.
      These methods write "count" values of data into stream.
      All functions return nonzero if successful.
  */
  //@{
  virtual int write (const int* data, const unsigned int count) 
    {if (fwrite (data, sizeof(int), count, stream)== count) return 1; else return 0;}
  virtual int write (const unsigned long* data, const unsigned int count) 
    {if (fwrite (data, sizeof(unsigned long), count, stream)== count) return 1; else return 0;}
  virtual int write (const double* data, const unsigned int count) 
    {if (fwrite (data, sizeof(double), count, stream)== count) return 1; else return 0;}
  virtual int write (const char* data, const unsigned int count) 
    {if (fwrite (data, sizeof(char), count, stream)== count) return 1; else return 0;}
  //@}
};

#ifdef __PARALLEL_MODE

/**
   Implementation of ComBuffDataStream representing DataStream interface to (MPI) communication bufer i/o.
   This class creates a Datastream shell around communication buffer routines. 
   @see DataStream class.
*/

class ComBuffDataStream : public DataStream {
 private:
  /// associated communication buffer
  CommunicationBuffer* buff;

 public:
  /// Constructor, takes associated communication buffer pointer as paramemetr 
  ComBuffDataStream (CommunicationBuffer* b) {buff = b;}
  /// Destructor
  ~ComBuffDataStream () {}

  /** @name Data Stream reading methods.
      These methods read "count" values from data stream into 
      array passed as the first argument. 
      All functions return nonzero if successful.
  */
  //@{
  virtual int read (int* data, const unsigned int count) {return buff->unpackArray (data, count);}
  virtual int read (unsigned long* data, const unsigned int count) {return buff->unpackArray(data, count);}
  virtual int read (double* data, const unsigned int count) {return buff->unpackArray (data, count);}
  virtual int read (char* data, const unsigned int count) {return buff->unpackArray(data, count);}
  //@}
  /** @name Data Stream writting methods.
      These methods write "count" values of data into stream.
      All functions return nonzero if successful.
  */
  //@{
  virtual int write (const int* data, const unsigned int count) {return buff->packArray (data, count);}
  virtual int write (const unsigned long* data, const unsigned int count) {return buff->packArray (data, count);}
  virtual int write (const double* data, const unsigned int count) {return buff->packArray (data, count);}
  virtual int write (const char* data, const unsigned int count) {return buff->packArray(data, count);}
  //@}
};


/**
   Implementation of ComBuffDataStream representing DataStream interface to (MPI) process communicator.
   This class creates a Datastream shell around process communicator routines. 
   @see DataStream class.
 */
class ProcessCommDataStream : public DataStream {
 private:
  /// associated process communicator buffer
  ProcessCommunicatorBuff* pc;
 public:
  /// Constructor
  ProcessCommDataStream(ProcessCommunicatorBuff* b) {pc=b;}
  /// Destructor
  ~ProcessCommDataStream() {}
 
  /** @name Data Stream reading methods.
      These methods read "count" values from data stream into 
      array passed as the first argument. 
      All functions return nonzero if successful.
  */
  //@{
  virtual int read (int* data, const unsigned int count) {return pc->unpackArray (data, count);}
  virtual int read (unsigned long* data, const unsigned int count) {return pc->unpackArray(data, count);}
  virtual int read (double* data, const unsigned int count) {return pc->unpackArray (data, count);}
  virtual int read (char* data, const unsigned int count) {return pc->unpackArray(data, count);}
  //@}
  /** @name Data Stream writting methods.
      These methods write "count" values of data into stream.
      All functions return nonzero if successful.
  */
  //@{
  virtual int write (const int* data, const unsigned int count) {return pc->packArray (data, count);}
  virtual int write (const unsigned long* data, const unsigned int count) {return pc->packArray (data, count);}
  virtual int write (const double* data, const unsigned int count) {return pc->packArray (data, count);}
  virtual int write (const char* data, const unsigned int count) {return pc->packArray(data, count);}
  //@}
};

#endif
#endif // datastream_h
