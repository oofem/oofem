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

#ifndef xmldatareader_h
#define xmldatareader_h

#include "datareader.h"
#include "xmlinputrecord.h"
#include <pugixml.hpp>
#include <set>
#include <map>




namespace oofem {
/**
 * Class representing the implementation of XML reader.
 * It reads a sequence of input records from data file and creates the corresponding input records.
 */
class OOFEM_EXPORT XMLDataReader : public DataReader
{
protected:
    friend XMLInputRecord;
    std::string topXmlFile;
    /* map parent xml node (which is either empty for top-level file for xi:include node for included files) to filenames */
    std::map<pugi::xml_node,std::string> xmlFiles;
    /* map document node to ordered list of newline offsets (used to compute line:column from offset in messages) */
    std::map<pugi::xml_node,std::vector<size_t>> newlines;
    /* map parent xml node (which is either empty for top-level file for xi:include node for included files) to (sub)document */
    std::map<pugi::xml_node,pugi::xml_document> docs;
    struct StackItem{
        pugi::xml_node parent;
        pugi::xml_node curr;
        std::shared_ptr<XMLInputRecord> lastRecord;
        std::set<pugi::xml_node> seen;
        int lastRecId=0;
    };
    std::vector<StackItem> stack;
    std::string giveStackPath(); // string representation
    pugi::xml_document& loadXml(pugi::xml_node parent, const std::string& xml);
    pugi::xml_node resolveXiInclude(const pugi::xml_node& n);
    std::tuple<size_t,size_t> offset2lc(const std::vector<size_t>& nl, size_t offset) const;
    std::string loc() const ;
    std::string loc(const pugi::xml_node&) const;
    std::shared_ptr<InputRecord> topRecord;
    pugi::xml_node giveNamedChild(const pugi::xml_node& parent, const std::string& name);
    const std::string XiIncludeTag="xi:include";
    int setRecId(int lastRecId);
public:
    XMLDataReader(const std::string& xmlFile);
    virtual ~XMLDataReader(){};
    bool hasFlattenedStructure() override { return true; }

    //! guess whether given file is XML
    static bool canRead(const std::string& xmlFile);
    InputRecord& giveInputRecord(InputRecordType, int recordId) override;
    InputRecord* giveTopInputRecord() override { return topRecord.get(); }
    bool peekNext(const std :: string &keyword) override { return false; } /* no peeking, it is used for hacks only */
    void finish() override;
    std::string giveReferenceName() const override { return topXmlFile; }
    void enterGroup(const std::string& name) override;
    void leaveGroup(const std::string& name) override;
    void enterRecord(InputRecord*) override;
    void leaveRecord(InputRecord*) override;

    int giveGroupCount(const std::string& name) override;
    int giveCurrentGroupCount();
};
} // end namespace oofem
#endif // xmldatareader_h
