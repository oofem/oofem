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

#include "dynamicdatareader.h"
#include "inputrecord.h"

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
    return *(this->it);
}

void
DynamicDataReader :: finish()
{
    // Not sure if i need to do this;
    for (std::list<InputRecord*>::iterator tempit = this->recordList.begin(); tempit != this->recordList.end(); ++tempit) {
        delete *tempit;
    }
    this->recordList.clear();
}

} // end namespace oofem
