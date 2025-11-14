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

#ifndef xmldatareader_h
#define xmldatareader_h

#include "datareader.h"
#include "xmlinputrecord.h"
#include <pugixml.hpp>
#include <set>




namespace oofem {
/**
 * Class representing the implementation of XML reader.
 * It reads a sequence of input records from data file and creates the corresponding input records.
 */
class OOFEM_EXPORT XMLDataReader : public DataReader
{
protected:
    friend XMLInputRecord;
    std::string xmlFile;
    std::vector<size_t> newlines;
    pugi::xml_document doc;
    struct StackItem{
        pugi::xml_node parent;
        pugi::xml_node curr;
        std::unique_ptr<XMLInputRecord> lastRecord;
        std::set<pugi::xml_node> seen;
    };
    std::vector<StackItem> stack;
    std::string giveStackPath(); // string representation
    std::tuple<size_t,size_t> offset2lc(size_t offset);
    std::string loc();
    std::string loc(const pugi::xml_node&);
    std::unique_ptr<InputRecord> topRecord;
public:
    XMLDataReader(const std::string& xmlFile);
    // XMLDataReader(const XMLDataReader&, pugi::xml_node subRoot);
    virtual ~XMLDataReader(){};
    bool hasFlattenedStructure() override { return true; }

    //! guess whether given file is XML
    static bool canRead(const std::string& xmlFile);
    InputRecord &giveInputRecord(InputRecordType, int recordId) override;
    InputRecord* giveTopInputRecord() override { return topRecord.get(); }
    bool peekNext(const std :: string &keyword) override { return false; } /* no peeking, it is used for hacks only */
    void finish() override;
    std::string giveReferenceName() const override { return xmlFile; }
    void enterGroup(const std::string& name) override;
    void leaveGroup(const std::string& name) override;
    void enterRecord(InputRecord*) override;
    void leaveRecord(InputRecord*) override;

    int giveGroupCount(const std::string& name) override;
    int giveCurrentGroupCount();
};
} // end namespace oofem
#endif // xmldatareader_h
