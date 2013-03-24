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

#include "datastream.h"
#include "processcomm.h"
#include "combuff.h"
#include "error.h"

namespace oofem
{

int DataStream::read ( std::string& data )
{
    int n;
    char *str;
    this->read ( &n, 1 );
    str = new char[n + 1];
    n = this->read ( str, n );
    str[n] = '\0';
    data = str;
    delete [] str;
    return n;
}

int DataStream::write ( const std::string& data )
{
    int n = ( int ) data.size();
    this->write ( &n, 1 );
    return this->write ( data.data(), n );
}

int FileDataStream::read ( int* data, unsigned int count )
{
    return ( fread ( data, sizeof ( int ), count, stream ) == count );
}

int FileDataStream::read ( long unsigned int* data, unsigned int count )
{
    return ( fread ( data, sizeof ( unsigned long ), count, stream ) == count );
}

int FileDataStream::read ( long int* data, unsigned int count )
{
    return ( fread ( data, sizeof ( long ), count, stream ) == count );
}

int FileDataStream::read ( double* data, unsigned int count )
{
    return ( fread ( data, sizeof ( double ), count, stream ) == count );
}

int FileDataStream::read ( char* data, unsigned int count )
{
    return ( fread ( data, sizeof ( char ), count, stream ) == count );
}

int FileDataStream::read ( bool* data, unsigned int count )
{
    return ( fread ( data, sizeof ( bool ), count, stream ) == count );
}

int FileDataStream::write ( const int* data, unsigned int count )
{
    return ( fwrite ( data, sizeof ( int ), count, stream ) == count );
}

int FileDataStream::write ( const long unsigned int* data, unsigned int count )
{
    return ( fwrite ( data, sizeof ( unsigned long ), count, stream ) == count );
}

int FileDataStream::write ( const long int* data, unsigned int count )
{
    return ( fwrite ( data, sizeof ( long ), count, stream ) == count );
}

int FileDataStream::write ( const double* data, unsigned int count )
{
    return ( fwrite ( data, sizeof ( double ), count, stream ) == count );
}

int FileDataStream::write ( const char* data, unsigned int count )
{
    return ( fwrite ( data, sizeof ( char ), count, stream ) == count );
}

int FileDataStream::write ( const bool* data, unsigned int count )
{
    return ( fwrite ( data, sizeof ( bool ), count, stream ) == count );
}

#ifdef __PARALLEL_MODE

int ComBuffDataStream::read ( int* data, unsigned int count )
{
    return buff->unpackArray ( data, count );
}

int ComBuffDataStream::read ( long unsigned int* data, unsigned int count )
{
    return buff->unpackArray ( data, count );
}

int ComBuffDataStream::read ( long int* data, unsigned int count )
{
    return buff->unpackArray ( data, count );
}

int ComBuffDataStream::read ( double* data, unsigned int count )
{
    return buff->unpackArray ( data, count );
}

int ComBuffDataStream::read ( char* data, unsigned int count )
{
    return buff->unpackArray ( data, count );
}

int ComBuffDataStream::read ( bool* data, unsigned int count )
{
    OOFEM_ERROR("ComBuffDataStream :: Can't read bool type");
    return 0;
}

int ComBuffDataStream::write ( const int* data, unsigned int count )
{
    return buff->packArray ( data, count );
}

int ComBuffDataStream::write ( const long unsigned int* data, unsigned int count )
{
    return buff->packArray ( data, count );
}

int ComBuffDataStream::write ( const long int* data, unsigned int count )
{
    return buff->packArray ( data, count );
}

int ComBuffDataStream::write ( const double* data, unsigned int count )
{
    return buff->packArray ( data, count );
}

int ComBuffDataStream::write ( const char* data, unsigned int count )
{
    return buff->packArray ( data, count );
}

int ComBuffDataStream::write ( const bool* data, unsigned int count )
{
    OOFEM_ERROR("ComBuffDataStream :: Can't write bool type");
    return 0;
}

int ProcessCommDataStream::read ( int* data, unsigned int count )
{
    return pc->unpackArray ( data, count );
}

int ProcessCommDataStream::read ( long unsigned int* data, unsigned int count )
{
    return pc->unpackArray ( data, count );
}

int ProcessCommDataStream::read ( long int* data, unsigned int count )
{
    return pc->unpackArray ( data, count );
}

int ProcessCommDataStream::read ( double* data, unsigned int count )
{
    return pc->unpackArray ( data, count );
}

int ProcessCommDataStream::read ( char* data, unsigned int count )
{
    return pc->unpackArray ( data, count );
}

int ProcessCommDataStream::read ( bool* data, unsigned int count )
{
    OOFEM_ERROR("ProcessCommDataStream :: Can't read bool type");
    return 0;
}

int ProcessCommDataStream::write ( const int* data, unsigned int count )
{
    return pc->packArray ( data, count );
}

int ProcessCommDataStream::write ( const long unsigned int* data, unsigned int count )
{
    return pc->packArray ( data, count );
}

int ProcessCommDataStream::write ( const long int* data, unsigned int count )
{
    return pc->packArray ( data, count );
}

int ProcessCommDataStream::write ( const double* data, unsigned int count )
{
    return pc->packArray ( data, count );
}

int ProcessCommDataStream::write ( const char* data, unsigned int count )
{
    return pc->packArray ( data, count );
}

int ProcessCommDataStream::write ( const bool* data, unsigned int count )
{
    OOFEM_ERROR("ProcessCommDataStream :: Can't write bool type");
    return 0;
}

#endif //__PARALLEL_MODE

}