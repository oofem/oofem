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

#include "oofemtxtdatareader.h"
#include "error.h"

namespace oofem {
OOFEMTXTDataReader :: OOFEMTXTDataReader(const char *inputfilename) : DataReader(), ir()
{
    inputStream.open( inputfilename );
    if ( !inputStream.is_open() ) {
        OOFEM_ERROR2("OOFEMTXTDataReader::OOFEMTXTDataReader: Can't open input stream (%s)", inputfilename);
    }
    dataSourceName = inputfilename;
    lineNumber = 0;
}

OOFEMTXTDataReader :: ~OOFEMTXTDataReader()
{
    finish();
}

InputRecord *
OOFEMTXTDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
{
    std::string line;
    if ( typeId == IR_outFileRec ) {
        this->giveRawLineFromInput(line);
    } else {
        this->giveLineFromInput(line);
    }

    ir.setRecordString(line);
    ir.setLineNumber(this->lineNumber);
    return & ir;
}


std::string
OOFEMTXTDataReader :: giveLine()
{
    return ir.giveRecordAsString();
}

void
OOFEMTXTDataReader :: finish()
{
    if ( inputStream.is_open() )
        inputStream.close();
}

void
OOFEMTXTDataReader :: giveLineFromInput(std::string &line)
{
    // reads one line from inputStream
    // if " detected, start/stop changing to lower case characters
    bool flag = false; //0-tolower, 1-remain with capitals

    this->giveRawLineFromInput(line);

    for ( std::size_t i = 0; i < line.size(); i++ ) {
        if ( line[i] == '"' ) { //do not change to lowercase inside quotation marks
            flag = !flag; // switch flag
        }

        if ( !flag ) {
            line[i] = tolower(line[i]); // convert line to lowercase
        }
    }
}

void
OOFEMTXTDataReader :: giveRawLineFromInput(std::string &line)
{
    //
    // reads one line from inputStream - for private use only.
    //
    do {
        this->lineNumber++;
        std::getline(this->inputStream, line);
    } while ( line[0] == '#' ); // skip comments

}
} // end namespace oofem
