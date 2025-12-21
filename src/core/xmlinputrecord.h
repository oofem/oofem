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

#ifndef xmlinputrecord_h
#define xmlinputrecord_h

#include "inputrecord.h"

#include <pugixml.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <set>

#define _XML_NI std::cerr<<__PRETTY_FUNCTION__<<": not yet implemented."<<std::endl; abort();

namespace oofem {
class XMLDataReader;
/**
 * Class representing the Input Record for OOFEM XML input file format.
 * It is represented as XML node with attributes, nested nodes, but without text.
 */
class OOFEM_EXPORT XMLInputRecord : public InputRecord
{
    pugi::xml_node node;
    friend XMLDataReader;
    std::set<std::string> attrQueried;
    std::set<std::string> attrRead;
    int recId=-1;
    XMLDataReader* _reader() const { return (XMLDataReader*)(this->giveReader()); }
    static std::string xmlizeAttrName(const std::string& s);
public:
    std::string _attr_traced_read(const char* name){ return std::get<0>(_attr_traced_read_with_node(name)); }
    std::tuple<std::string,pugi::xml_node> _attr_traced_read_with_node(const char* name);

    XMLInputRecord(XMLDataReader* reader_, const pugi::xml_node& node_);
    std::shared_ptr<InputRecord> clone() const override { return std::make_shared<XMLInputRecord>(*this); }

    void finish(bool wrn = true) override;

    std::string loc() const { return loc(node); }
    std::string loc(const pugi::xml_node& node) const ;
    int setRecId(int lastRecId);


    void giveRecordKeywordField(std :: string &answer, int &value) override;
    void giveRecordKeywordField(std :: string &answer) override;
    void giveField(int &answer, InputFieldType id) override;
    void giveField(double &answer, InputFieldType id) override;
    void giveField(bool &answer, InputFieldType id) override;
    void giveField(std :: string &answer, InputFieldType id) override;
    void giveField(FloatArray &answer, InputFieldType id) override;
    void giveField(IntArray &answer, InputFieldType id) override;
    void giveField(FloatMatrix &answer, InputFieldType id) override;
    void giveField(std :: vector< std :: string > &answer, InputFieldType id) override { _XML_NI; }
    void giveField(Dictionary &answer, InputFieldType id) override;
    void giveField(std :: list< Range > &answer, InputFieldType id) override;
    void giveField(ScalarFunction &answer, InputFieldType id) override;

    int giveGroupCount(InputFieldType id, const std::string& name, bool optional) override;
    bool hasChild(InputFieldType id, const std::string& name, bool optional) override;

    bool hasField(InputFieldType id) override;
    void printYourself() override { _XML_NI; }

    std::string giveRecordAsString() const override;
    std::string giveRecordInTXTFormat() const override;
    std::string giveLocation() const override { return loc(); }


    static bool node_seen_get(const pugi::xml_node& n);
    static void node_seen_set(pugi::xml_node& n, bool seen);
    static constexpr char SeenMark[]="__seen";
};
} // end namespace oofem
#endif // xmlinputrecord_h
