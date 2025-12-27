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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "oofemenv.h"
#include "inputrecord.h"
#include "error.h"

#include<iostream>
#include<fstream>

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
        IR_enrichFrontRec, IR_propagationLawRec, IR_crackNucleationRec, IR_fracManRec, IR_failCritRec,
        IR_contactSurfaceRec, IR_fieldRec, 
        // MPM specific
        IR_mpmVarRec, IR_mpmTermRec, IR_mpmIntegralRec,
        IR_unspecified // internal use only, signifies error in setting record type
    };
    /* XML tags corresponding to record types; those with "" are just enumeration group where arbitrary tags may be used */
    static constexpr const char* InputRecordTags[]={
        "Domain","OutputManager","DomainComp","Geometry","GBPM",
        "Analysis","Metastep",/*ExportModule*/"",/*Node*/"",/*Element*/"",
        /*CrossSection*/"",/*Material*/"",/*"NonlocalBarrier"*/"",/*BoundaryCondition*/"","InitialCondition",/*TimeFunction*/"","Set",
        "XFemManager","EnrichmentFunction","EnrichmentGeometry",/*EnrichmentItem*/"",
        /*EnrichmentFront*/"","PropagationLaw","CrackNucleation","FractureManager","FailCriterion",
        /*ContactSurface*/"","Field",
        "MPMVariable",/*"MPMTerm"*/"","MPMIntegral",
        "UNSPECIFIED"
    };

    DataReader() { }
    virtual ~DataReader() { }

    /**
     * Returns input record corresponding to given InputRecordType value and its record_id.
     * The returned InputRecord reference is valid only until the next call.
     * @param irType Determines type of record to be returned.
     * @param recordId Determines the record  number corresponding to component number.
     */
    virtual InputRecord &giveInputRecord(InputRecordType irType, int recordId) = 0;
    /**
     * Returns top input record, for readers which support it; others return empty pointer
     */
    virtual InputRecord* giveTopInputRecord(){ return nullptr; }

    /**
     * Peak in advance into the record list.
     * @return True if next keyword is a set.
     */
    virtual bool peekNext(const std :: string &keyword) { return false; }

    /**
     * Allows to detach all data connections.
     */
    virtual void finish() = 0;

    /// Gives the reference file name (e.g. file name)
    virtual std :: string giveReferenceName() const = 0;
    /// Gives the output file name
    std :: string giveOutputFileName() { return this->outputFileName; }
    /// Gives the problem description
    std :: string giveDescription() { return this->description; }

    virtual bool hasFlattenedStructure() { return false; }

    virtual void enterGroup(const std::string& name) {};
    virtual void leaveGroup(const std::string& name) {};
    virtual void enterRecord(InputRecord* rec) {};
    virtual void leaveRecord(InputRecord* rec) {};

    /// RAII guard for DataReader::enterRecord and DataReader::leaveRecord.
    class RecordGuard{
        DataReader& reader;
        InputRecord* rec;
    public:
        RecordGuard(DataReader& reader_, InputRecord* rec_): reader(reader_), rec(rec_) { reader.enterRecord(rec); }
        ~RecordGuard() { reader.leaveRecord(rec); }
    };

    /// Internal range-like class, return type for giveGroupRecords methods
    class GroupRecords {
        DataReader& dr;
        std::string group;
        InputRecordType irType;
        int size_;
    public:
        class Iterator {
            DataReader& dr;
            std::string group;
            InputRecordType irType;
            int size;
            int index;
            InputRecord* irPtr=nullptr;
            bool entered=false;
        public:
            Iterator( DataReader &dr_, const std::string &group_, InputRecordType irType_, int size_, int index_ );
            Iterator &operator++();
            InputRecord& operator*() { return *irPtr; }
            bool operator!=(const Iterator& other){ return this->index!=other.index; }
            int index1() const { return index+1; }
        };
        GroupRecords( DataReader &dr_, const std::string &group_, InputRecordType irType_, int size );
        Iterator begin(){ return Iterator(dr,group,irType,size_,0); }
        Iterator end(){ return Iterator(dr,group,irType,size_,std::max(size_,0)); }
        int size(){ return size_; }
    };

    static const int NoSuchGroup=-1;
    /**
     * Return number of items in a named group
     * @param name name of the group; if empty, current group is used
     * @return Number of items in the named group; returns 0 if the group is empty, and -1 (DatReader::NoSuchGroup) if the group does not exist at all.
     */
    virtual int giveGroupCount(const std::string& name){ return NoSuchGroup; }
    /// Predicate whether a named group exists (it can still be empty).
    bool hasGroup(const std::string& name){ return giveGroupCount(name)>=0; }
    /**
     * Give range (provides begin(), end(), size()) to iterate over records.
     * Some readers (text) specify the number of records in the InputRecord with the given input field type.
     * Other readers (XML) find the number of records as number of elements in the subgroup enclosed with *name* tag.
     * @param ir Input record which may hold the number of subsequent entries to be read from the stream
     * @param ift Field type in *ir* record specifying the number of records.
     * @param name Subgroup name for tag-based readers.
     * @param irType Expected input record type for all records to be read.
     * @param optional If not optional and the number of records is not given, fail with error. Otherwise assume 0-sized subgroup.
     * @return Object providing begin(), end() iterators and size().
     */
    GroupRecords giveGroupRecords(const std::shared_ptr<InputRecord> &ir, InputFieldType ift, const std::string &name, InputRecordType irType, bool optional );
    /**
     * Give range to iterate over records within a named group
     * @param name Subgroup name; if not given, give records within the current group
     * @param irType Input record type for records to be read
     * @param numRequired if non-negative, this number is checked against number of records present (for readers which can determine that), and mismatch error is thrown when they are different
     * @return Object providing being(), end() iterators and size().
     */
    GroupRecords giveGroupRecords(const std::string& name, InputRecordType irType, int numRequired=-1);
    /// Return pointer to subrecord of given type (must be exactly one); if not present, returns nullptr.
    InputRecord *giveChildRecord( const std::shared_ptr<InputRecord> &ir, InputFieldType ift, const std::string &name, InputRecordType irType, bool optional );


};
} // end namespace oofem
#endif // datareader_h
