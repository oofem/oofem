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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "datastream.h"
#include "error.h"
#include <vector>

namespace oofem
{
int DataStream :: read(std :: string &data)
{
    int n;
    std :: vector< char >str;
    if ( !this->read(& n, 1) ) {
        data = "";
        return 0;
    }
    str.resize(n);
    if ( !this->read(str.data(), n) ) {
        data = "";
        return 0;
    }
    data = std::string(str.data(), n);
    return 1;
}

int DataStream :: write(const std :: string &data)
{
    int n = ( int ) data.size();
    this->write(& n, 1);
    return this->write(data.data(), n);
}

FileDataStream :: FileDataStream(std::string filename, bool write): 
    stream(nullptr),
    filename(std::move(filename))
{
    //stream.open(filename, (write ? std::ios::out : std::ios:in) | std::ios::binary );
    this->stream = fopen(this->filename.c_str(), write ? "wb" : "rb" );
    if ( !this->stream ) {
        throw CantOpen(this->filename);
    }
}

FileDataStream :: ~FileDataStream()
{
    fclose(this->stream);
}

int FileDataStream :: read(int *data, std::size_t count)
{
    return ( fread(data, sizeof( int ), count, stream) == count );
    //this->stream.read(reinterpret_cast< char* >(data), sizeof(int)*count);
    //return this->stream.good();
}

int FileDataStream :: read(unsigned long *data, std::size_t count)
{
    return ( fread(data, sizeof( unsigned long ), count, stream) == count );
}

int FileDataStream :: read(long *data, std::size_t count)
{
    return ( fread(data, sizeof( long ), count, stream) == count );
}

#ifdef _WIN32
int FileDataStream::read(std::size_t* data, std::size_t count)
{
    return (fread(data, sizeof(std::size_t), count, stream) == count);
}
#endif

int FileDataStream :: read(double *data, std::size_t count)
{
    return ( fread(data, sizeof( double ), count, stream) == count );
}

int FileDataStream :: read(char *data, std::size_t count)
{
    return ( fread(data, sizeof( char ), count, stream) == count );
}

int FileDataStream :: read(bool &data)
{
    return ( fread(& data, sizeof( bool ), 1, stream) == 1 );
}

int FileDataStream :: write(const int *data, std::size_t count)
{
    return ( fwrite(data, sizeof( int ), count, stream) == count );
}

int FileDataStream :: write(const unsigned long *data, std::size_t count)
{
    return ( fwrite(data, sizeof( unsigned long ), count, stream) == count );
}

#ifdef _WIN32
int FileDataStream::write(const std::size_t* data, std::size_t count)
{
    return (fwrite(data, sizeof(std::size_t), count, stream) == count);
}
#endif

int FileDataStream :: write(const long *data, std::size_t count)
{
    return ( fwrite(data, sizeof( long ), count, stream) == count );
}

int FileDataStream :: write(const double *data, std::size_t count)
{
    return ( fwrite(data, sizeof( double ), count, stream) == count );
}

int FileDataStream :: write(const char *data, std::size_t count)
{
    return ( fwrite(data, sizeof( char ), count, stream) == count );
}

int FileDataStream :: write(bool data)
{
    return ( fwrite(& data, sizeof( bool ), 1, stream) == 1 );
}

int FileDataStream :: givePackSizeOfInt(std::size_t count)
{
    return (int) (sizeof(int)*count);
}

int FileDataStream :: givePackSizeOfDouble(std::size_t count)
{
    return (int) (sizeof(double)*count);
}

int FileDataStream :: givePackSizeOfChar(std::size_t count)
{
    return (int) (sizeof(char)*count);
}

int FileDataStream :: givePackSizeOfBool(std::size_t count)
{
    return (int) (sizeof(bool)*count);
}

int FileDataStream :: givePackSizeOfLong(std::size_t count)
{
    return (int) (sizeof(long)*count);
}

int FileDataStream :: givePackSizeOfSizet(std::size_t count)
{
    return (int) (sizeof(std::size_t)*count);
}

}
