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

#ifndef buffereddatareader_h
#define buffereddatareader_h

#include <string>

#include "datareader.h"
#include "oofemtxtinputrecord.h"
#include "dynalist.h"

namespace oofem {
/**
 * Class representing the implementation of plain text data reader.
 * It reads a sequential sequence of input records from data file
 * and creates the corresponding input records.
 * There is no check for record type requested, it is assumed that records are
 * written in correct order, which determined by the coded sequence of
 * component initialization and described in input manual.
 */
class BufferedDataReader : public DataReader
{
protected:
    OOFEMTXTInputRecord ir;
    dynaList< std :: string > buffer;
    dynaList< std :: string > :: iterator pos;

public:
    /** Constructor. */
    BufferedDataReader();
    virtual ~BufferedDataReader();

    virtual InputRecord *giveInputRecord(InputRecordType, int recordId);

    /// @return The name (shortened) of data source.
    virtual const char *giveDataSourceName() const { return "BufferedDataReader source"; }

    /**
     * Appends string to end of buffer.
     * @param str String to add to append
     */
    void appendInputString(std :: string &str);
    /**
     * Appends string to end of buffer.
     * @param line Null terminated string to append.
     */
    void appendInputString(const char *line);
    /// Puts buffer iterator back at start.
    void rewind();
    /**
     * Puts buffer iterator at given position.
     * @param position Position to seek to.
     */
    void seek(int position);

    virtual void finish();

    /// Prints buffer to stdout.
    void printYourself();

    virtual std::string giveLine();

    void setOutputFileName(const std::string &fname) { this->outputFileName = fname; }
    void setDescription(const std::string &desc) { this->description = desc; }

    /**
     * Dumps buffer to file.
     * @param fileName Name of file to dump data to.
     */
    void writeToFile(char *fileName);

protected:
    void giveRawLineFromInput(std::string &line);
    void giveLineFromInput(std::string &line);
};
} // end namespace oofem
#endif // buffereddatareader_h
