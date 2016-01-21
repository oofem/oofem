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

#ifndef dynamicdatareader_h
#define dynamicdatareader_h

#include "datareader.h"
#include <list>
#include <memory>

namespace oofem {
class InputRecord;

/**
 * Class representing the implementation of a dynamic data reader for in-code use.
 * Despite its name, it can also write OOFEM text-input files.
 *
 * @see DynamicInputRecord It is the intended complement for in-code generation of FE-problem intialization.
 * @author Mikael Öhman
 * @todo InputRecordType is ignored. It shouldn't be too difficult to respect it, but it's not necessary.
 */
class OOFEM_EXPORT DynamicDataReader : public DataReader
{
protected:
    /// Keeps track of the current position in the list
    std :: list< std :: unique_ptr< InputRecord > > :: iterator it;
    /// All record types will be appended to this list, no split in terms of InputRecordType is implemented yet.
    std :: list< std :: unique_ptr< InputRecord > >recordList;

public:
    /// Constructor.
    DynamicDataReader();
    virtual ~DynamicDataReader();

    DynamicDataReader(const DynamicDataReader& src) = delete;
    DynamicDataReader &operator = (const DynamicDataReader &src) = delete;

    /**
     * Main purpose of this class it the possibility to add new input records in code.
     * The input records can be any implementation, but the intended use would be together with the DynamicInputRecord class.
     * @param type Currently ignored, but left here for consistency with giveInputRecord
     * @param record New record to be added at the end. New input records have to be added in the same order as the text input files do it.
     */
    void insertInputRecord(InputRecordType type, InputRecord *record);

    /**
     * Sets the output file name. Used for writing input files.
     */
    void setOutputFileName(const std :: string &outputFileName) { this->outputFileName = outputFileName; }
    /**
     * Sets the description line. Used for writing input files.
     */
    void setDescription(const std :: string &description) { this->description = description; }

    virtual InputRecord *giveInputRecord(InputRecordType, int recordId);
    virtual bool peakNext(const std :: string &keyword);
    virtual void finish();
    virtual const char *giveDataSourceName() const { return ""; }

    /**
     * Writes all input records to file.
     * @param fileName Name of file to dump data to.
     */
    void writeToFile(const char *fileName);
};
} // end namespace oofem
#endif // dynamicdatareader_h
