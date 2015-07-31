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

#include "dynamicdatareader.h"
#include "inputrecord.h"
#include "error.h"

#include <fstream>

namespace oofem {
DynamicDataReader :: DynamicDataReader() : DataReader()
{
    this->it = recordList.end();
}

DynamicDataReader :: ~DynamicDataReader()
{
}

void
DynamicDataReader :: insertInputRecord(InputRecordType type, InputRecord *record)
{
    // Should care about the record type, but its a hassle.
    this->recordList.emplace_back(record);
}

InputRecord *
DynamicDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
{
    // Ignores recordId in favor of having a dynamic list (just incremental access). typeId could be supported, but its a hassle.
    // The txt data reader makes the same assumptions.
    if ( this->it == this->recordList.end() ) {
        this->it = this->recordList.begin();
    } else {
        ++this->it;
    }
    return this->it->get();
}

bool
DynamicDataReader :: peakNext(const std :: string &keyword)
{
    std :: string nextKey;
    (*this->it)->giveRecordKeywordField(nextKey);
    return keyword.compare( nextKey ) == 0;
}

void
DynamicDataReader :: finish()
{
    this->recordList.clear();
}

void
DynamicDataReader :: writeToFile(const char *fileName)
{
    std :: ofstream fout(fileName);

    fout << this->outputFileName << '\n';
    fout << this->description << '\n';
    for ( auto &rec: this->recordList ) {
        fout << rec->giveRecordAsString() << "\n";
    }
    fout.close();
}
} // end namespace oofem
