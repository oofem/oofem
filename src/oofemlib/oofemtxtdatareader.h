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

#ifndef oofemtxtdatareader_h
#define oofemtxtdatareader_h

#include "datareader.h"
#include "oofemtxtinputrecord.h"
#include "oofem_limits.h"

namespace oofem {
/**
 * Class representing the implementation of plain text date reader.
 * It reads a sequence of input records from data file
 * and creates the corresponding input records.
 * There is no check for record type requested, it is assumed that records are
 * written in correct order, which determined by the coded sequence of
 * component initialization and described in input manual.
 */
class OOFEMTXTDataReader : public DataReader
{
protected:
    OOFEMTXTInputRecord ir;
    FILE *inputStream;
    std::string dataSourceName;

public:
    /// Constructor.
    OOFEMTXTDataReader(const char *inputfilename);
    virtual ~OOFEMTXTDataReader();

    virtual InputRecord *giveInputRecord(InputRecordType, int recordId);
    virtual void giveLine(char *line);
    virtual void finish();
    virtual const char *giveDataSourceName() const { return dataSourceName.c_str(); }

protected:
    void giveLineFromInput(char *line);
    void giveRawLineFromInput(char *line);
};
} // end namespace oofem
#endif // oofemtxtdatareader_h
