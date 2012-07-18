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

#include "buffereddatareader.h"
#include "error.h"

#ifndef __MAKEDEPEND
 #include <cstdio>
 #include <cctype>
#endif

namespace oofem {
BufferedDataReader :: BufferedDataReader()
{ }

BufferedDataReader :: ~BufferedDataReader()
{ }

InputRecord *
BufferedDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
{
    char line [ OOFEM_MAX_LINE_LENGTH + 1 ];
    if ( typeId == IR_outFileRec ) {
        this->giveRawLineFromInput(line);
    } else {
        this->giveLineFromInput(line);
    }

    ir.setRecordString(line);
    return & ir;
}

void
BufferedDataReader :: finish()
{
    buffer.clear();
}

void
BufferedDataReader :: rewind()
{
    pos = buffer.begin();
}

void
BufferedDataReader :: seek(int position)
{
    int i = 0;

    pos = buffer.begin();

    while ( pos != buffer.end() ) {
        if ( (*pos)[0] != '#' ) {
            if ( ++i == position ) {
                return;
            }
        }

        ++pos;
    }

    OOFEM_ERROR("BufferedDataReader: seek: already at the end");
}


void
BufferedDataReader :: printYourself()
{
    dynaList< std :: string > :: iterator rec;

    printf( "\nBuffer with %d records\n\n", buffer.size() );

    rec = buffer.begin();
    while ( rec != buffer.end() ) {
        printf( "%s\n", ( * rec ).c_str() );
        ++rec;
    }

    printf("\n");
}


void
BufferedDataReader :: giveLine(char *line)
{
  this->giveRawLineFromInput(line);
}

void
BufferedDataReader :: writeToFile(char *fileName)
{
    dynaList< std :: string > :: iterator rec;
    FILE *dataStream;

    if ( ( dataStream = fopen(fileName, "w") ) == NULL ) {
        OOFEM_ERROR2("BufferedDataReader::writeToFile : Can't open data stream %s", fileName);
    }

    rec = buffer.begin();
    while ( rec != buffer.end() ) {
        fprintf( dataStream, "%s\n", ( * rec ).c_str() );
        ++rec;
    }

    fclose(dataStream);
}


void
BufferedDataReader :: appendInputString(std :: string &str)
{
    // Append zero char to the end of the string to prevent rubish to occur
    // at the of the line read in giveRawLineFromInput if previous string was longer !!!
    str += '\0';
    buffer.pushBack(str);
}

void
BufferedDataReader :: appendInputString(const char *line)
{
    std :: string str(line);

    // Append zero char to the end of the string to prevent rubish to occur
    // at the of the line read in giveRawLineFromInput if previous string was longer !!!
    str += '\0';
    buffer.pushBack(str);
}


void
BufferedDataReader :: giveLineFromInput(char *line)
{
    char *ptr;

    giveRawLineFromInput(line);
    // convert line to lowercase
    for ( ptr = line; ( * ptr = tolower(* ptr) ); ptr++ ) {
        ;
    }
}


void
BufferedDataReader :: giveRawLineFromInput(char *line)
{
    do {
        if ( pos == buffer.end() ) {
            OOFEM_ERROR("BufferedDataReader: giveRawLineFromInput: already at the end");
        }

        ( * pos ).copy(line, std :: string :: npos);
        ++pos;
    } while ( * line == '#' ); // skip comments

}
} // end namespace oofem
