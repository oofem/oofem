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

#ifndef datastream_h
#define datastream_h

#include "combuff.h"
#include "processcomm.h"

#include <sstream>
#include <cstdio>

namespace oofem {

/**
 * The purpose of DataStream abstract class is to allow to store/restore context to different streams,
 * including file, communication buffers, etc., using the same routine.
 * This will facilitate many algorithms relying on saving/moving state of components
 * (such as load balancing), without writing new (and very similar) routines.
 * This  will lead to a  better consistency of code.
 */
class DataStream
{
public:
    /// Destructor
    virtual ~DataStream() { }
    /**
     * @name Data Stream reading methods.
     * These methods read "count" values from data stream into
     * array passed as the first argument.
     * All functions return nonzero if successful.
     */
    //@{
    /// Reads count integer values into array pointed by data.
    virtual int read(int *data, const unsigned int count) = 0;
    /// Reads count unsigned long values into array pointed by data.
    virtual int read(unsigned long *data, const unsigned int count) = 0;
    /// Reads count long values into array pointed by data.
    virtual int read(long *data, const unsigned int count) = 0;
    /// Reads count double values into array pointed by data.
    virtual int read(double *data, const unsigned int count) = 0;
    /// Reads count char values into array pointed by data.
    virtual int read(char *data, const unsigned int count) = 0;
    //@}

    /**
     * @name Data Stream writing methods.
     * These methods write "count" values of data into stream.
     * All functions return nonzero if successful.
     */
    //@{
    /// Writes count integer values from array pointed by data.
    virtual int write(const int *data, const unsigned int count) = 0;
    /// Writes count unsigned long values from array pointed by data.
    virtual int write(const unsigned long *data, const unsigned int count) = 0;
    /// Writes count long values from array pointed by data.
    virtual int write(const long *data, const unsigned int count) = 0;
    /// Writes count double values from array pointed by data.
    virtual int write(const double *data, const unsigned int count) = 0;
    /// Writes count char values from array pointed by data.
    virtual int write(const char *data, const unsigned int count) = 0;
    //@}
};

/**
 * Implementation of StringDataStream representing DataStream interface to stringstream i/o.
 * This class creates a DataStream shell around c stream i/o routines. This class will
 * not provide any methods for opening/closing stream. This is the responsibility of user.
 * @see DataStream class.
 */
class StringDataStream : public DataStream
{
private:
    /// FILE pointer of associated stream
    std::stringstream *stream;

public:
    /// Constructor, takes associated stream pointer as parameter
    StringDataStream(std::stringstream *s) { stream = s; }
    /// Destructor (will not close stream!)
    virtual ~StringDataStream() { }

    virtual int read(int *data, const unsigned int count)
    { if ( stream->read(reinterpret_cast<char *>(data), sizeof *data) ) { return 1; } else { return 0; } }
    virtual int read(unsigned long *data, const unsigned int count)
    { if ( stream->read(reinterpret_cast<char *>(data), sizeof *data) ) { return 1; } else { return 0; } }
    virtual int read(long *data, const unsigned int count)
    { if ( stream->read(reinterpret_cast<char *>(data), sizeof *data) ) { return 1; } else { return 0; } }
    virtual int read(double *data, const unsigned int count)
    { if ( stream->read(reinterpret_cast<char *>(data), sizeof *data) ) { return 1; } else { return 0; } }
    virtual int read(char *data, const unsigned int count)
    { if ( stream->read(reinterpret_cast<char *>(data), sizeof *data) ) { return 1; } else { return 0; } }

    virtual int write(const int *data, const unsigned int count)
    { if ( stream->write(reinterpret_cast<const char *>(data), sizeof *data) ) { return 1; } else { return 0; } }
    virtual int write(const unsigned long *data, const unsigned int count)
    { if ( stream->write(reinterpret_cast<const char *>(data), sizeof *data) ) { return 1; } else { return 0; } }
    virtual int write(const long *data, const unsigned int count)
    { if ( stream->write(reinterpret_cast<const char *>(data), sizeof *data) ) { return 1; } else { return 0; } }
    virtual int write(const double *data, const unsigned int count)
    { if ( stream->write(reinterpret_cast<const char *>(data), sizeof *data) ) { return 1; } else { return 0; } }
    virtual int write(const char *data, const unsigned int count)
    { if ( stream->write(reinterpret_cast<const char *>(data), sizeof *data) ) { return 1; } else { return 0; } }
};

/**
 * Implementation of FileDataStream representing DataStream interface to file i/o.
 * This class creates a DataStream shell around c file i/o routines. This class will
 * not provide any methods for opening/closing file. This is the responsibility of user.
 * @see DataStream class.
 */
class FileDataStream : public DataStream
{
private:
    /// FILE pointer of associated stream
    FILE *stream;
public:
    /// Constructor, takes associated stream pointer as parameter
    FileDataStream(FILE *s) { stream = s; }
    /// Destructor (will not close stream!)
    virtual ~FileDataStream() { }

    virtual int read(int *data, const unsigned int count)
    { if ( fread(data, sizeof( int ), count, stream) == count ) { return 1; } else { return 0; } }
    virtual int read(unsigned long *data, const unsigned int count)
    { if ( fread(data, sizeof( unsigned long ), count, stream) == count ) { return 1; } else { return 0; } }
    virtual int read(long *data, const unsigned int count)
    { if ( fread(data, sizeof( long ), count, stream) == count ) { return 1; } else { return 0; } }
    virtual int read(double *data, const unsigned int count)
    { if ( fread(data, sizeof( double ), count, stream) == count ) { return 1; } else { return 0; } }
    virtual int read(char *data, const unsigned int count)
    { if ( fread(data, sizeof( char ), count, stream) == count ) { return 1; } else { return 0; } }

    virtual int write(const int *data, const unsigned int count)
    { if ( fwrite(data, sizeof( int ), count, stream) == count ) { return 1; } else { return 0; } }
    virtual int write(const unsigned long *data, const unsigned int count)
    { if ( fwrite(data, sizeof( unsigned long ), count, stream) == count ) { return 1; } else { return 0; } }
    virtual int write(const long *data, const unsigned int count)
    { if ( fwrite(data, sizeof( long ), count, stream) == count ) { return 1; } else { return 0; } }
    virtual int write(const double *data, const unsigned int count)
    { if ( fwrite(data, sizeof( double ), count, stream) == count ) { return 1; } else { return 0; } }
    virtual int write(const char *data, const unsigned int count)
    { if ( fwrite(data, sizeof( char ), count, stream) == count ) { return 1; } else { return 0; } }
};

#ifdef __PARALLEL_MODE

/**
 * Implementation of ComBuffDataStream representing DataStream interface to (MPI) communication buffer i/o.
 * This class creates a DataStream shell around communication buffer routines.
 * @see DataStream class.
 */

class ComBuffDataStream : public DataStream
{
private:
    /// Associated communication buffer
    CommunicationBuffer *buff;

public:
    /// Constructor, takes associated communication buffer pointer as parameter
    ComBuffDataStream(CommunicationBuffer *b) { buff = b; }
    /// Destructor
    virtual ~ComBuffDataStream() { }

    virtual int read(int *data, const unsigned int count) { return buff->unpackArray(data, count); }
    virtual int read(unsigned long *data, const unsigned int count) { return buff->unpackArray(data, count); }
    virtual int read(long *data, const unsigned int count) { return buff->unpackArray(data, count); }
    virtual int read(double *data, const unsigned int count) { return buff->unpackArray(data, count); }
    virtual int read(char *data, const unsigned int count) { return buff->unpackArray(data, count); }

    virtual int write(const int *data, const unsigned int count) { return buff->packArray(data, count); }
    virtual int write(const unsigned long *data, const unsigned int count) { return buff->packArray(data, count); }
    virtual int write(const long *data, const unsigned int count) { return buff->packArray(data, count); }
    virtual int write(const double *data, const unsigned int count) { return buff->packArray(data, count); }
    virtual int write(const char *data, const unsigned int count) { return buff->packArray(data, count); }
};


/**
 * Implementation of ComBuffDataStream representing DataStream interface to (MPI) process communicator.
 * This class creates a DataStream shell around process communicator routines.
 * @see DataStream class.
 */
class ProcessCommDataStream : public DataStream
{
private:
    /// Associated process communicator buffer
    ProcessCommunicatorBuff *pc;
public:
    /// Constructor
    ProcessCommDataStream(ProcessCommunicatorBuff *b) { pc = b; }
    /// Destructor
    virtual ~ProcessCommDataStream() { }

    virtual int read(int *data, const unsigned int count) { return pc->unpackArray(data, count); }
    virtual int read(unsigned long *data, const unsigned int count) { return pc->unpackArray(data, count); }
    virtual int read(long *data, const unsigned int count) { return pc->unpackArray(data, count); }
    virtual int read(double *data, const unsigned int count) { return pc->unpackArray(data, count); }
    virtual int read(char *data, const unsigned int count) { return pc->unpackArray(data, count); }

    virtual int write(const int *data, const unsigned int count) { return pc->packArray(data, count); }
    virtual int write(const unsigned long *data, const unsigned int count) { return pc->packArray(data, count); }
    virtual int write(const long *data, const unsigned int count) { return pc->packArray(data, count); }
    virtual int write(const double *data, const unsigned int count) { return pc->packArray(data, count); }
    virtual int write(const char *data, const unsigned int count) { return pc->packArray(data, count); }
};

#endif
} // end namespace oofem
#endif // datastream_h
