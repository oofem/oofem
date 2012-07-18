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

#ifndef datareader_h
#define datareader_h


#include "inputrecord.h"

namespace oofem {
/**
 * Class representing the abstraction for input data source.
 * Its role is to provide input records for particular components.
 * The input records are identified by record type and component number.
 * The order of input records is in fact determined by the coded sequence of
 * component initialization. The input record identification facilitates the
 * implementation of database readers with direct or random access.
 */
class DataReader
{
public:
    /// Determines the type of input record.
    enum InputRecordType {
        IR_outFileRec, IR_jobRec, IR_domainRec, IR_outManRec, IR_domainCompRec, IR_geometryRec, IR_gbpmRec,
        IR_emodelRec, IR_mstepRec, IR_expModuleRec, IR_dofmanRec, IR_elemRec,
        IR_crosssectRec, IR_matRec, IR_nlocBarRec, IR_bcRec, IR_icRec, IR_ltfRec,
        IR_nRandomFieldGenRec, IR_xfemManRec, IR_enrichFuncRec, IR_geoRec, IR_enrichItemRec
    };

    DataReader() { }
    virtual ~DataReader() { }

    /**
     * Returns input record corresponding to given InputRecordType value and its record_id.
     * The returned InputRecord reference is valid only until the next call.
     * @param irType Determines type of record to be returned.
     * @param recordId Determines the record  number corresponding to component number.
     */
    virtual InputRecord *giveInputRecord(InputRecordType irType, int recordId) = 0;

    /**
     * Reads the whole line
     */
    virtual void giveLine(char *line) { };
    
    /// Return a line number, which is helpful for tracking errors.
    virtual int giveLineNumber() { return lineNumber; };
    
    /**
     * Allows to detach all data connections.
     */
    virtual void finish() = 0;

    /// Prints the name (shortened) of data source.
    virtual const char *giveDataSourceName() const = 0;
    /// Prints the error message.
    void report_error(const char *_class, const char *proc, const char *kwd, IRResultType result, const char *file, int line);
    
protected:
    /// Keep track of read line from stream. Used for error repors.
    int lineNumber;
};
} // end namespace oofem
#endif // datareader_h
