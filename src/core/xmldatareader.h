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
#include "xmlutil.h"
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
    static constexpr int FormatLowest=1;
    static constexpr int FormatHighest=2;
    int formatVersion=FormatLowest;
    std::map<pugi::xml_node,std::shared_ptr<xmlutil::XmlDoc>> docs;
    struct StackItem{
        pugi::xml_node parent;
        pugi::xml_node curr;
        std::shared_ptr<XMLInputRecord> lastRecord;
        std::set<pugi::xml_node> seen;
        int lastRecId=0;
    };
    std::vector<StackItem> stack;
    std::string giveStackPath(); // string representation
    xmlutil::XmlDoc& loadXml(pugi::xml_node parent, const std::string& xml);
    pugi::xml_node resolveXiInclude(pugi::xml_node& n);

    std::string loc() const ;
    std::string loc(const pugi::xml_node&) const;
    std::shared_ptr<InputRecord> topRecord;
    pugi::xml_node giveNamedChild(const pugi::xml_node& parent, const std::string& name);
    const std::string XiIncludeTag="xi:include";
    int setRecId(int lastRecId);
public:
    XMLDataReader(const std::string& xmlFile);
    virtual ~XMLDataReader(){};
    bool hasFeature(FormatFeature f) override;

    //! guess whether given file is XML
    static bool canRead(const std::string& xmlFile);
    std::shared_ptr<InputRecord> giveNextInputRecord(InputRecordType) override;
    std::shared_ptr<InputRecord> giveTopInputRecord() override;
    bool peekNext(const std :: string &keyword) override { return false; } /* no peeking, it is used for hacks only */
    void finish() override;
    std::string giveReferenceName() const override { return topXmlFile; }
    void enterGroup(const std::string& name) override;
    void leaveGroup(const std::string& name) override;
    void enterRecord(const std::shared_ptr<InputRecord>) override;
    void leaveRecord(const std::shared_ptr<InputRecord>) override;

    int giveGroupCount(const std::string& name) override;
    int giveCurrentGroupCount();
};
} // end namespace oofem
#endif // xmldatareader_h
