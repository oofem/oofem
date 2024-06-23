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

#include "oofemenv.h"

#include <sstream>
#include <cstdio>
#include <exception>
#include <stdexcept>

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
    virtual int read(int *data, std::size_t count) = 0;
    int read(int &data) { return this->read(&data, 1); }
    /// Reads count unsigned long values into array pointed by data.
    virtual int read(unsigned long *data, std::size_t count) = 0;
    int read(unsigned long &data) { return this->read(&data, 1); }
#ifdef _MSC_VER
    /// Reads count unsigned std::size_t values into array pointed by data.
    virtual int read(std::size_t* data, std::size_t count) = 0;
    int read(std::size_t& data) { return this->read(&data, 1); }
#endif
    /// Reads count long values into array pointed by data.
    virtual int read(long *data, std::size_t count) = 0;
    int read(long &data) { return this->read(&data, 1); }
    /// Reads count double values into array pointed by data.
    virtual int read(double *data, std::size_t count) = 0;
    int read(double &data) { return this->read(&data, 1); }
    /// Reads count char values into array pointed by data.
    virtual int read(char *data, std::size_t count) = 0;
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
    virtual int write(const int *data, std::size_t count) = 0;
    int write(int data) { return this->write(&data, 1); }
    /// Writes count unsigned long values from array pointed by data.
    virtual int write(const unsigned long *data, std::size_t count) = 0;
    int write(unsigned long data) { return this->write(&data, 1); }
#ifdef _MSC_VER
    /// Writes count std::size_t values from array pointed by data.
    virtual int write(const std::size_t* data, std::size_t count) = 0;
    int write(std::size_t data) { return this->write(&data, 1); }
#endif
    /// Writes count long values from array pointed by data.
    virtual int write(const long *data, std::size_t count) = 0;
    int write(long data) { return this->write(&data, 1); }
    /// Writes count double values from array pointed by data.
    virtual int write(const double *data, std::size_t count) = 0;
    int write(double data) { return this->write(&data, 1); }
    /// Writes count char values from array pointed by data.
    virtual int write(const char *data, std::size_t count) = 0;
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
    virtual int givePackSizeOfInt(std::size_t count) = 0;
    virtual int givePackSizeOfDouble(std::size_t count) = 0;
    virtual int givePackSizeOfChar(std::size_t count) = 0;
    virtual int givePackSizeOfBool(std::size_t count) = 0;
    virtual int givePackSizeOfLong(std::size_t count) = 0;
    virtual int givePackSizeOfSizet(std::size_t count) = 0;
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
public:
    class CantOpen : public std::runtime_error
    {
    public:
        std::string filename;
        CantOpen(std::string file): std::runtime_error("can't open file"), filename(std::move(file)) {}
    };

private:
    /// FILE pointer of associated stream
    FILE *stream;
    /// Filename
    std :: string filename;
public:
    /// Constructor, takes associated stream pointer as parameter
    FileDataStream(std :: string filename, bool write);

    /// Destructor (will not close stream!)
    virtual ~FileDataStream();

    int read(int *data, std::size_t count) override;
    int read(unsigned long *data, std::size_t count) override;
#ifdef _MSC_VER
    int read(std::size_t *data, std::size_t count) override;
#endif
    int read(long *data, std::size_t count) override;
    int read(double *data, std::size_t count) override;
    int read(char *data, std::size_t count) override;
    int read(bool &data) override;

    int write(const int *data, std::size_t count) override;
    int write(const unsigned long *data, std::size_t count) override;
#ifdef _MSC_VER
    int write(const std::size_t* data, std::size_t count);
#endif
    int write(const long *data, std::size_t count) override;
    int write(const double *data, std::size_t count) override;
    int write(const char *data, std::size_t count) override;
    int write(bool data) override;

    int givePackSizeOfInt(std::size_t count) override;
    int givePackSizeOfDouble(std::size_t count) override;
    int givePackSizeOfChar(std::size_t count) override;
    int givePackSizeOfBool(std::size_t count) override;
    int givePackSizeOfLong(std::size_t count) override;
    int givePackSizeOfSizet(std::size_t count) override;

};

} // end namespace oofem
#endif // datastream_h
