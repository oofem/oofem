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

#include "oofemtxtdatareader.h"
#include "error.h"

#include <string>

namespace oofem {
OOFEMTXTDataReader :: OOFEMTXTDataReader(const char *inputfilename) : DataReader(),
    ir(), inputStream(), dataSourceName(inputfilename)
{
    inputStream.open(inputfilename);
    if ( !inputStream.is_open() ) {
        OOFEM_SIMPLE_ERROR("OOFEMTXTDataReader::OOFEMTXTDataReader: Can't open input stream (%s)", inputfilename);
    }
    dataSourceName = inputfilename;
    lineNumber = 0;

    this->giveRawLineFromInput(outputFileName);
    this->giveRawLineFromInput(description);
}

OOFEMTXTDataReader :: OOFEMTXTDataReader(const OOFEMTXTDataReader &x) : DataReader(),
    ir(), inputStream(), dataSourceName(x.dataSourceName)
{
    inputStream.open( dataSourceName.c_str() );
    if ( !inputStream.is_open() ) {
        OOFEM_SIMPLE_ERROR( "OOFEMTXTDataReader::OOFEMTXTDataReader: Can't copy open input stream (%s)", dataSourceName.c_str() );
    }
    lineNumber = 0;

    this->giveRawLineFromInput(outputFileName);
    this->giveRawLineFromInput(description);
}

OOFEMTXTDataReader :: ~OOFEMTXTDataReader()
{
    finish();
}

InputRecord *
OOFEMTXTDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
{
    std :: string line;
    this->giveLineFromInput(line);

    ir.setRecordString(line);
    ir.setLineNumber(this->lineNumber);
    return & ir;
}


std :: string
OOFEMTXTDataReader :: giveLine()
{
    return ir.giveRecordAsString();
}

void
OOFEMTXTDataReader :: finish()
{
    if ( inputStream.is_open() ) {
        inputStream.close();
    }
}

void
OOFEMTXTDataReader :: giveLineFromInput(std :: string &line)
{
    // reads one line from inputStream
    // if " detected, start/stop changing to lower case characters
    bool flag = false; //0-tolower, 1-remain with capitals

    this->giveRawLineFromInput(line);

    for ( std :: size_t i = 0; i < line.size(); i++ ) {
        if ( line [ i ] == '"' ) { //do not change to lowercase inside quotation marks
            flag = !flag; // switch flag
        }

        if ( !flag ) {
            line [ i ] = tolower(line [ i ]); // convert line to lowercase
        }
    }
}

void
OOFEMTXTDataReader :: giveRawLineFromInput(std :: string &line)
{
    //
    // reads one line from inputStream - for private use only.
    //
    do {
        this->lineNumber++;
        std :: getline(this->inputStream, line);
#if __cplusplus > 199711L
        if ( line.back() == '\\' ) {
            line.pop_back();
            std :: string continuedLine;
            this->lineNumber++;
            do {
                std :: getline(this->inputStream, continuedLine);
                continuedLine.pop_back();
                line += continuedLine;
            } while ( continuedLine.back() == '\\' );
        }
#endif
    } while ( line [ 0 ] == '#' ); // skip comments
}
} // end namespace oofem
