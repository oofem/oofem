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
#endif

namespace oofem
{
int DataStream :: read(std :: string &data)
{
    int n;
    char *str;
    if ( !this->read(& n, 1) ) {
        data = "";
        return 0;
    }
    str = new char [ n + 1 ];
    if ( !this->read(str, n) ) {
        data = "";
        return 0;
    }
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

int FileDataStream :: read(int *data, int count)
{
    return ( (int)fread(data, sizeof( int ), count, stream) == count );
}

int FileDataStream :: read(unsigned long *data, int count)
{
    return ( (int)fread(data, sizeof( unsigned long ), count, stream) == count );
}

int FileDataStream :: read(long *data, int count)
{
    return ( (int)fread(data, sizeof( long ), count, stream) == count );
}

int FileDataStream :: read(double *data, int count)
{
    return ( (int)fread(data, sizeof( double ), count, stream) == count );
}

int FileDataStream :: read(char *data, int count)
{
    return ( (int)fread(data, sizeof( char ), count, stream) == count );
}

int FileDataStream :: read(bool &data)
{
    return ( (int)fread(& data, sizeof( bool ), 1, stream) == 1 );
}

int FileDataStream :: write(const int *data, int count)
{
    return ( (int)fwrite(data, sizeof( int ), count, stream) == count );
}

int FileDataStream :: write(const unsigned long *data, int count)
{
    return ( (int)fwrite(data, sizeof( unsigned long ), count, stream) == count );
}

int FileDataStream :: write(const long *data, int count)
{
    return ( (int)fwrite(data, sizeof( long ), count, stream) == count );
}

int FileDataStream :: write(const double *data, int count)
{
    return ( (int)fwrite(data, sizeof( double ), count, stream) == count );
}

int FileDataStream :: write(const char *data, int count)
{
    return ( (int)fwrite(data, sizeof( char ), count, stream) == count );
}

int FileDataStream :: write(bool data)
{
    return ( (int)fwrite(& data, sizeof( bool ), 1, stream) == 1 );
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

int ProcessCommDataStream :: read(int *data, int count)
{
    return pc->read(data, count);
}

int ProcessCommDataStream :: read(unsigned long *data, int count)
{
    return pc->read(data, count);
}

int ProcessCommDataStream :: read(long *data, int count)
{
    return pc->read(data, count);
}

int ProcessCommDataStream :: read(double *data, int count)
{
    return pc->read(data, count);
}

int ProcessCommDataStream :: read(char *data, int count)
{
    return pc->read(data, count);
}

int ProcessCommDataStream :: read(bool &data)
{
    char val;
    int ret = pc->read(& val, 1);
    data = val != 0;
    return ret;
}

int ProcessCommDataStream :: write(const int *data, int count)
{
    return pc->write(data, count);
}

int ProcessCommDataStream :: write(const unsigned long *data, int count)
{
    return pc->write(data, count);
}

int ProcessCommDataStream :: write(const long *data, int count)
{
    return pc->write(data, count);
}

int ProcessCommDataStream :: write(const double *data, int count)
{
    return pc->write(data, count);
}

int ProcessCommDataStream :: write(const char *data, int count)
{
    return pc->write(data, count);
}

int ProcessCommDataStream :: write(bool data)
{
    char val = data;
    return pc->write(& val, 1);
}

int ProcessCommDataStream :: givePackSizeOfInt(int count)
{
    OOFEM_ERROR("remove this function");
    return 0;
}

int ProcessCommDataStream :: givePackSizeOfDouble(int count)
{
    OOFEM_ERROR("remove this function");
    return 0;
}

int ProcessCommDataStream :: givePackSizeOfChar(int count)
{
    OOFEM_ERROR("remove this function");
    return 0;
}

int ProcessCommDataStream :: givePackSizeOfBool(int count)
{
    OOFEM_ERROR("remove this function");
    return 0;
}

int ProcessCommDataStream :: givePackSizeOfLong(int count)
{
    OOFEM_ERROR("remove this function");
    return 0;
}

#endif //__PARALLEL_MODE
}
