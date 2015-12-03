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
OOFEMTXTDataReader :: OOFEMTXTDataReader(std :: string inputfilename) : DataReader(),
    dataSourceName(std :: move(inputfilename)), recordList()
{
    std :: list< std :: pair< int, std :: string > >lines;
    // Read all the lines in the main input file:
    {
        std :: ifstream inputStream(dataSourceName);
        if ( !inputStream.is_open() ) {
            OOFEM_ERROR("Can't open input stream (%s)", dataSourceName.c_str());
        }

        int lineNumber = 0;
        std :: string line;

        this->giveRawLineFromInput(inputStream, lineNumber, outputFileName);
        this->giveRawLineFromInput(inputStream, lineNumber, description);

        while (this->giveLineFromInput(inputStream, lineNumber, line)) {
            lines.emplace_back(make_pair(lineNumber, line));
        }
    }
    // Check for included files: @include "somefile"
    for ( auto it = lines.begin(); it != lines.end(); ++it ) {
        if ( it->second.compare(0, 8, "@include") == 0 ) {
            std :: string fname = it->second.substr(10, it->second.length()-11);
            OOFEM_LOG_INFO("Reading included file: %s\n", fname.c_str());

            // Remove the include line
            lines.erase(it++);
            // Add all the included lines:
            int includedLine = 0;
            std :: string line;
            std :: ifstream includedStream(fname);
            if ( !includedStream.is_open() ) {
                OOFEM_ERROR("Can't open input stream (%s)", fname.c_str());
            }
            while (this->giveLineFromInput(includedStream, includedLine, line)) {
                lines.emplace(it, make_pair(includedLine, line));
            }
        }
    }
    ///@todo This could be parallelized, but I'm not sure it is worth it 
    /// (might make debugging faulty input files harder for users as well)
    for ( auto &line: lines ) {
        //printf("line: %s\n", line.second.c_str());
        this->recordList.emplace_back(line.first, line.second);
    }
    this->it = this->recordList.begin();
}

OOFEMTXTDataReader :: OOFEMTXTDataReader(const OOFEMTXTDataReader &x) : OOFEMTXTDataReader(x.dataSourceName) {}

OOFEMTXTDataReader :: ~OOFEMTXTDataReader()
{
}

InputRecord *
OOFEMTXTDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
{
    if ( this->it == this->recordList.end() ) {
        OOFEM_ERROR("Out of input records, file contents must be missing");
    }
    return &(*this->it++);
}

bool
OOFEMTXTDataReader :: peakNext(const std :: string &keyword)
{
    std :: string nextKey;
    this->it->giveRecordKeywordField(nextKey);
    return keyword.compare( nextKey ) == 0;
}

void
OOFEMTXTDataReader :: finish()
{
    if ( this->it != this->recordList.end() ) {
        OOFEM_WARNING("There are unread lines in the input file\n"
            "The most common cause are missing entries in the domain record, e.g. 'nset'");
    }
    this->recordList.clear();
}

bool
OOFEMTXTDataReader :: giveLineFromInput(std :: ifstream &stream, int &lineNum, std :: string &line)
{
    // reads one line from inputStream
    // if " detected, start/stop changing to lower case characters
    bool flag = false; //0-tolower, 1-remain with capitals

    bool read = this->giveRawLineFromInput(stream, lineNum, line);
    if ( !read ) {
        return false;
    }

    for ( auto &c: line ) {
        if ( c == '"' ) { //do not change to lowercase inside quotation marks
            flag = !flag; // switch flag
        }

        if ( !flag ) {
            c = (char)tolower(c); // convert line to lowercase
        }
    }
    return true;
}

bool
OOFEMTXTDataReader :: giveRawLineFromInput(std :: ifstream &stream, int &lineNum, std :: string &line)
{
    //
    // reads one line from inputStream - for private use only.
    //
    do {
        lineNum++;
        std :: getline(stream, line);
        if ( !stream ) {
            return false;
        } if ( line.length() > 0 ) {
            if ( line.back() == '\\' ) {
                std :: string continuedLine;
                do {
                    lineNum++;
                    std :: getline(stream, continuedLine);
                    if ( !stream ) {
                        return false;
                    }
                    line.pop_back();
                    line += continuedLine;
                } while ( continuedLine.back() == '\\' );
            }
        }
    } while ( line.length() == 0 || line [ 0 ] == '#' ); // skip comments
    return true;
}
} // end namespace oofem
