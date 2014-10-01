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
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef datastream_h
#define datastream_h

#include "oofemcfg.h"

#include <sstream>
#include <cstdio>

namespace oofem {
/**
 * The purpose of DataStream abstract class is to allow to store/restore context to different streams,
 * including file, communication buffers, etc., using the same routine.
 * This will facilitate many algorithms relying on saving/moving state of components
 * (such as load balancing), without writing new (and very similar) routines.
 * This  will lead to a  better consistency of code.
 * @todo Are the "long" and "unsigned long" functions really necessary?
 */
class OOFEM_EXPORT DataStream
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
    virtual int read(int *data, int count) = 0;
    int read(int &data) { return this->read(&data, 1); }
    /// Reads count unsigned long values into array pointed by data.
    virtual int read(unsigned long *data, int count) = 0;
    int read(unsigned long &data) { return this->read(&data, 1); }
    /// Reads count long values into array pointed by data.
    virtual int read(long *data, int count) = 0;
    int read(long &data) { return this->read(&data, 1); }
    /// Reads count double values into array pointed by data.
    virtual int read(double *data, int count) = 0;
    int read(double &data) { return this->read(&data, 1); }
    /// Reads count char values into array pointed by data.
    virtual int read(char *data, int count) = 0;
    int read(char &data) { return this->read(&data, 1); }
    /// Reads a bool value from data.
    virtual int read(bool &data) = 0;
    /// Reads a string (stored as an int for the length followed by char*).
    int read(std :: string &data);
    //@}

    /**
     * @name Data Stream writing methods.
     * These methods write "count" values of data into stream.
     * All functions return nonzero if successful.
     */
    //@{
    /// Writes count integer values from array pointed by data.
    virtual int write(const int *data, int count) = 0;
    int write(int data) { return this->write(&data, 1); }
    /// Writes count unsigned long values from array pointed by data.
    virtual int write(const unsigned long *data, int count) = 0;
    int write(unsigned long data) { return this->write(&data, 1); }
    /// Writes count long values from array pointed by data.
    virtual int write(const long *data, int count) = 0;
    int write(long data) { return this->write(&data, 1); }
    /// Writes count double values from array pointed by data.
    virtual int write(const double *data, int count) = 0;
    int write(double data) { return this->write(&data, 1); }
    /// Writes count char values from array pointed by data.
    virtual int write(const char *data, int count) = 0;
    int write(char data) { return this->write(&data, 1); }
    /// Writes a bool value.
    virtual int write(bool data) = 0;
    /// Reads a string (stored as an int for the length followed by char*).
    int write(const std :: string &data);
    /// Writes a string (wrapper needed, otherwise write(bool) is called )
    int write(const char *data) { return this->write(std :: string(data)); }
    //@}

    /**
     * @name Sizing functions.
     * These methods compute the stored size (in bytes) of an array containing "count" elements.
     */
    //@{
    virtual int givePackSizeOfInt(int count) = 0;
    virtual int givePackSizeOfDouble(int count) = 0;
    virtual int givePackSizeOfChar(int count) = 0;
    virtual int givePackSizeOfBool(int count) = 0;
    virtual int givePackSizeOfLong(int count) = 0;
    //@}
};


/**
 * Implementation of FileDataStream representing DataStream interface to file i/o.
 * This class creates a DataStream shell around c file i/o routines. This class will
 * not provide any methods for opening/closing file. This is the responsibility of user.
 * @see DataStream class.
 */
class OOFEM_EXPORT FileDataStream : public DataStream
{
private:
    /// FILE pointer of associated stream
    FILE *stream;
public:
    /// Constructor, takes associated stream pointer as parameter
    FileDataStream(FILE * s) {
        stream = s;
    }
    /// Destructor (will not close stream!)
    virtual ~FileDataStream() { }

    virtual int read(int *data, int count);
    virtual int read(unsigned long *data, int count);
    virtual int read(long *data, int count);
    virtual int read(double *data, int count);
    virtual int read(char *data, int count);
    virtual int read(bool &data);

    virtual int write(const int *data, int count);
    virtual int write(const unsigned long *data, int count);
    virtual int write(const long *data, int count);
    virtual int write(const double *data, int count);
    virtual int write(const char *data, int count);
    virtual int write(bool data);

    virtual int givePackSizeOfInt(int count);
    virtual int givePackSizeOfDouble(int count);
    virtual int givePackSizeOfChar(int count);
    virtual int givePackSizeOfBool(int count);
    virtual int givePackSizeOfLong(int count);
};

} // end namespace oofem
#endif // datastream_h
