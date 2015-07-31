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

#ifndef datareader_h
#define datareader_h

#include "oofemcfg.h"
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
class OOFEM_EXPORT DataReader
{
protected:
    /// Output file name (first line in OOFEM input files).
    std :: string outputFileName;
    /// Description line (second line in OOFEM input files).
    std :: string description;

public:
    /// Determines the type of input record.
    enum InputRecordType {
        IR_domainRec, IR_outManRec, IR_domainCompRec, IR_geometryRec, IR_gbpmRec,
        IR_emodelRec, IR_mstepRec, IR_expModuleRec, IR_dofmanRec, IR_elemRec,
        IR_crosssectRec, IR_matRec, IR_nlocBarRec, IR_bcRec, IR_icRec, IR_funcRec, IR_setRec,
        IR_xfemManRec, IR_enrichFuncRec, IR_geoRec, IR_enrichItemRec,
        IR_enrichFrontRec, IR_propagationLawRec, IR_fracManRec, IR_failCritRec, 
        IR_contactManRec, IR_contactDefRec
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
     * Peak in advance into the record list.
     * @return True if next keyword is a set.
     */
    virtual bool peakNext(const std :: string &keyword) { return false; }

    /**
     * Allows to detach all data connections.
     */
    virtual void finish() = 0;

    /// Gives the output file name
    std :: string giveOutputFileName() { return this->outputFileName; }
    /// Gives the problem description
    std :: string giveDescription() { return this->description; }

    /// Prints the name (shortened) of data source.
    virtual const char *giveDataSourceName() const = 0;
    /// Prints the error message.
    void report_error(const char *_class, const char *proc, const char *kwd, IRResultType result, const char *file, int line);
};
} // end namespace oofem
#endif // datareader_h
