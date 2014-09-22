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

#include "datastream.h"
#include "error.h"

#ifdef __PARALLEL_MODE
 #include "processcomm.h"
 #include "combuff.h"
#endif

namespace oofem
{
int DataStream :: read(std :: string &data)
{
    int n;
    char *str;
    this->read(& n, 1);
    str = new char [ n + 1 ];
    this->read(str, n);
    str [ n ] = '\0';
    data = str;
    delete [] str;
    return 1;
}

int DataStream :: write(const std :: string &data)
{
    int n = ( int ) data.size();
    this->write(& n, 1);
    return this->write(data.data(), n);
}

int FileDataStream :: read(int *data, unsigned int count)
{
    return ( fread(data, sizeof( int ), count, stream) == count );
}

int FileDataStream :: read(unsigned long *data, unsigned int count)
{
    return ( fread(data, sizeof( unsigned long ), count, stream) == count );
}

int FileDataStream :: read(long *data, unsigned int count)
{
    return ( fread(data, sizeof( long ), count, stream) == count );
}

int FileDataStream :: read(double *data, unsigned int count)
{
    return ( fread(data, sizeof( double ), count, stream) == count );
}

int FileDataStream :: read(char *data, unsigned int count)
{
    return ( fread(data, sizeof( char ), count, stream) == count );
}

int FileDataStream :: read(bool &data)
{
    return ( fread(& data, sizeof( bool ), 1, stream) == 1 );
}

int FileDataStream :: write(const int *data, unsigned int count)
{
    return ( fwrite(data, sizeof( int ), count, stream) == count );
}

int FileDataStream :: write(const unsigned long *data, unsigned int count)
{
    return ( fwrite(data, sizeof( unsigned long ), count, stream) == count );
}

int FileDataStream :: write(const long *data, unsigned int count)
{
    return ( fwrite(data, sizeof( long ), count, stream) == count );
}

int FileDataStream :: write(const double *data, unsigned int count)
{
    return ( fwrite(data, sizeof( double ), count, stream) == count );
}

int FileDataStream :: write(const char *data, unsigned int count)
{
    return ( fwrite(data, sizeof( char ), count, stream) == count );
}

int FileDataStream :: write(bool data)
{
    return ( fwrite(& data, sizeof( bool ), 1, stream) == 1 );
}

int FileDataStream :: givePackSizeOfInt(int count)
{
    return sizeof(int)*count;
}

int FileDataStream :: givePackSizeOfDouble(int count)
{
    return sizeof(double)*count;
}

int FileDataStream :: givePackSizeOfChar(int count)
{
    return sizeof(char)*count;
}

int FileDataStream :: givePackSizeOfBool(int count)
{
    return sizeof(bool)*count;
}

int FileDataStream :: givePackSizeOfLong(int count)
{
    return sizeof(int)*count;
}


#ifdef __PARALLEL_MODE

int ComBuffDataStream :: read(int *data, unsigned int count)
{
    return buff->unpackArray(data, count);
}

int ComBuffDataStream :: read(unsigned long *data, unsigned int count)
{
    return buff->unpackArray(data, count);
}

int ComBuffDataStream :: read(long *data, unsigned int count)
{
    return buff->unpackArray(data, count);
}

int ComBuffDataStream :: read(double *data, unsigned int count)
{
    return buff->unpackArray(data, count);
}

int ComBuffDataStream :: read(char *data, unsigned int count)
{
    return buff->unpackArray(data, count);
}

int ComBuffDataStream :: read(bool &data)
{
    char val;
    int ret = buff->unpackArray(& val, 1);
    data = val != 0;
    return ret;
}

int ComBuffDataStream :: write(const int *data, unsigned int count)
{
    return buff->packArray(data, count);
}

int ComBuffDataStream :: write(const unsigned long *data, unsigned int count)
{
    return buff->packArray(data, count);
}

int ComBuffDataStream :: write(const long *data, unsigned int count)
{
    return buff->packArray(data, count);
}

int ComBuffDataStream :: write(const double *data, unsigned int count)
{
    return buff->packArray(data, count);
}

int ComBuffDataStream :: write(const char *data, unsigned int count)
{
    return buff->packArray(data, count);
}

int ComBuffDataStream :: write(bool data)
{
    char val = data;
    return buff->packArray(& val, 1);
}

int ComBuffDataStream :: givePackSizeOfInt(int count)
{
    return buff->givePackSize(MPI_INT, count);
}

int ComBuffDataStream :: givePackSizeOfDouble(int count)
{
    return buff->givePackSize(MPI_DOUBLE, count);
}

int ComBuffDataStream :: givePackSizeOfChar(int count)
{
    return buff->givePackSize(MPI_CHAR, count);
}

int ComBuffDataStream :: givePackSizeOfBool(int count)
{
    return buff->givePackSize(MPI_CHAR, count);
}

int ComBuffDataStream :: givePackSizeOfLong(int count)
{
    return buff->givePackSize(MPI_LONG, count);
}



int ProcessCommDataStream :: read(int *data, unsigned int count)
{
    return pc->unpackArray(data, count);
}

int ProcessCommDataStream :: read(unsigned long *data, unsigned int count)
{
    return pc->unpackArray(data, count);
}

int ProcessCommDataStream :: read(long *data, unsigned int count)
{
    return pc->unpackArray(data, count);
}

int ProcessCommDataStream :: read(double *data, unsigned int count)
{
    return pc->unpackArray(data, count);
}

int ProcessCommDataStream :: read(char *data, unsigned int count)
{
    return pc->unpackArray(data, count);
}

int ProcessCommDataStream :: read(bool &data)
{
    char val;
    int ret = pc->unpackArray(& val, 1);
    data = val != 0;
    return ret;
}

int ProcessCommDataStream :: write(const int *data, unsigned int count)
{
    return pc->packArray(data, count);
}

int ProcessCommDataStream :: write(const unsigned long *data, unsigned int count)
{
    return pc->packArray(data, count);
}

int ProcessCommDataStream :: write(const long *data, unsigned int count)
{
    return pc->packArray(data, count);
}

int ProcessCommDataStream :: write(const double *data, unsigned int count)
{
    return pc->packArray(data, count);
}

int ProcessCommDataStream :: write(const char *data, unsigned int count)
{
    return pc->packArray(data, count);
}

int ProcessCommDataStream :: write(bool data)
{
    char val = data;
    return pc->packArray(& val, 1);
}

int ProcessCommDataStream :: givePackSizeOfInt(int count)
{
    return pc->givePackSize(MPI_INT, count);
}

int ProcessCommDataStream :: givePackSizeOfDouble(int count)
{
    return pc->givePackSize(MPI_DOUBLE, count);
}

int ProcessCommDataStream :: givePackSizeOfChar(int count)
{
    return pc->givePackSize(MPI_CHAR, count);
}

int ProcessCommDataStream :: givePackSizeOfBool(int count)
{
    return pc->givePackSize(MPI_BOOL, count);
}

int ProcessCommDataStream :: givePackSizeOfLong(int count)
{
    return pc->givePackSize(MPI_CHAR, count);
}

#endif //__PARALLEL_MODE
}
