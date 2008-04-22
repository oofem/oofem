/* $Header: /home/cvs/bp/oofem/oofemlib/src/plaintextdatareader.C,v 1.2.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <ctype.h>
#endif
#include "plaintextdatareader.h"

PlainTextDataReader :: PlainTextDataReader(char *inputfilename)
{
    if ( ( inputStream = fopen(inputfilename, "r") ) == NULL ) {
        OOFEM_ERROR2("PlainTextDataReader::PlainTextDataReader: Can't open input stream (%s)", inputfilename);
    }

    strncpy(dataSourceName, inputfilename, MAX_FILENAME_LENGTH);
    dataSourceName [ MAX_FILENAME_LENGTH - 1 ] = '\0';
}

PlainTextDataReader :: ~PlainTextDataReader()
{
    finish();
}

InputRecord *
PlainTextDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
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
PlainTextDataReader :: finish()
{
    if ( inputStream ) {
        fclose(inputStream);
    }

    inputStream = NULL;
}

void
PlainTextDataReader :: giveLineFromInput(char *line)
{
    //
    // reads one line from inputStream - for private use only.
    //
    char *ptr;

    giveRawLineFromInput(line);
    // convert line to lowercase
    for ( ptr = line; ( * ptr = tolower(* ptr) ); ptr++ ) {
        ;
    }
}

void
PlainTextDataReader :: giveRawLineFromInput(char *line)
{
    //
    // reads one line from inputStream - for private use only.
    //
    int maxchar = OOFEM_MAX_LINE_LENGTH;
    do {
        fgets(line, maxchar, inputStream);
        if ( line == NULL ) {
            OOFEM_ERROR("PlainTextDataReader::giveRawLineFromInput: Premature end  of file encountered");
        }
    } while ( * line == '#' ); // skip comments

}
