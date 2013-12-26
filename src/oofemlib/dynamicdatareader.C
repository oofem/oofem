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
#include "oofemtxtinputrecord.h"
#include "dynamicinputrecord.h"
#include "error.h"

#include <fstream>

namespace oofem {
DynamicDataReader :: DynamicDataReader() : DataReader()
{
    this->it = recordList.end();
}

DynamicDataReader :: ~DynamicDataReader()
{
    this->finish();
}

void
DynamicDataReader :: insertInputRecord(InputRecordType type, InputRecord *record)
{
    // Should care about the record type, but its a hassle.
    this->recordList.push_back(record);
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
    return * ( this->it );
}

void
DynamicDataReader :: finish()
{
    // Not sure if i need to do this;
    for ( std :: list< InputRecord * > :: iterator tempit = this->recordList.begin(); tempit != this->recordList.end(); ++tempit ) {
        delete * tempit;
    }
    this->recordList.clear();
}

void
DynamicDataReader :: writeToFile(const char *fileName)
{
    std :: ofstream fout(fileName);

    fout << this->outputFileName << '\n';
    fout << this->description << '\n';
    for ( std :: list< InputRecord * > :: iterator it = this->recordList.begin(); it != this->recordList.end(); ++it ) {
        DynamicInputRecord *dyn;
        OOFEMTXTInputRecord *txt;
        if ( ( dyn = dynamic_cast< DynamicInputRecord * >( * it ) ) ) {
            fout << dyn->giveRecordAsString() << "\n";
        } else if ( ( txt = dynamic_cast< OOFEMTXTInputRecord * >( * it ) ) ) {
            fout << txt->giveRecordAsString() << '\n';
        } else {
            OOFEM_ERROR("DynamicDataReader :: writeToFile - A non-text or dynamic input record found, can't be printed to file\n");
        }
    }
    fout.close();
}
} // end namespace oofem
