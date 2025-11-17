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

#include "xmldatareader.h"
#include "error.h"

#include <string>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <limits>
#include <pugixml.hpp>


namespace oofem {
    pugi::xml_document& XMLDataReader::loadXml(const std::string& xml){
        // trasverse just to get line breaks
        std::ifstream i(xml.c_str(),std::ifstream::in|std::ios::binary);
        std::list<size_t> br;
        size_t pos=0;
        while(i.good()){ char c=i.get(); pos++; if(c=='\n') br.push_back(pos); }
        std::vector<size_t> nl;
        nl.reserve(br.size());
        if(br.size()==std::numeric_limits<size_t>().max()) abort(); // workaround GCC warning: iteration 2305843009213693952 invokes undefined behavior [-Waggressive-loop-optimizations]
        nl.assign(br.begin(),br.end());
        docs[xml]=pugi::xml_document();
        auto& doc=docs[xml];
        pugi::xml_parse_result result=doc.load_file(xml.c_str());
        if(!result){
            auto [line,col]=offset2lc(nl,result.offset);
            OOFEM_ERROR("Error parsing %s:%d:%d: %s",xml.c_str(),line,col,result.description());
        }
        newlines[doc]=std::move(nl);
        xmlFiles[doc]=xml;
        return doc;
    }


    XMLDataReader::XMLDataReader(const std::string& xmlFile_): topXmlFile(xmlFile_) {
        pugi::xml_document& doc=loadXml(topXmlFile);
        stack.push_back(StackItem{doc,doc});
        enterGroup("oofem");
        pugi::xml_node root=stack.back().parent;
        pugi::xml_node output_node=root.child("Output"), description_node=root.child("Description");
        if(!output_node) OOFEM_ERROR("Error reading %s: <Output> node missing in %s.",topXmlFile.c_str(),root.name());
        if(!description_node) OOFEM_ERROR("Error reading %s: <Description> node missing.",topXmlFile.c_str());
        XMLInputRecord::node_seen_set(output_node,true);
        XMLInputRecord::node_seen_set(description_node,true);
        outputFileName=output_node.text().as_string();
        description=description_node.text().as_string();
        topRecord=std::unique_ptr<InputRecord>(new XMLInputRecord(this,root,-1));
    }

    bool XMLDataReader :: canRead(const std::string& xml){
        return std::filesystem::path(xml).extension() == ".xml";
    }

    std::tuple<size_t,size_t> XMLDataReader::offset2lc(const std::vector<size_t>& nl, size_t offset){
        if(nl.empty()) return std::make_tuple(0,0);
        size_t ix=std::distance(nl.begin(),std::lower_bound(nl.begin(),nl.end(),offset));
        // std::cerr<<"offset="<<offset<<",ix="<<ix<<std::endl;
        return std::make_tuple(ix,offset-(ix==0?0:nl[ix-1]));
    }
    std::string XMLDataReader::loc(){ return loc(stack.back().curr?stack.back().curr:stack.back().parent); }
    std::string XMLDataReader::loc(const pugi::xml_node& n){
        std::ostringstream oss;
        auto [line,col]=offset2lc(newlines[n.root()],n.offset_debug());
        oss<<xmlFiles[n.root()]<<":"<<line<<":"<<col;
        return oss.str();
    }

    int XMLDataReader::giveGroupCount(const std::string& name) {
        if(name.empty()) return giveCurrentGroupCount();
        enterGroup(name);
        int ret=this->giveCurrentGroupCount();
        leaveGroup(name);
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<giveStackPath()<<" // "<<name<<": size is "<<ret);
        return ret;
    }

    int XMLDataReader::giveCurrentGroupCount(){
        int ret=0;
        for([[maybe_unused]] pugi::xml_node n: stack.back().parent.children()) ret++;
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<giveStackPath()<<": size is "<<ret);
        return ret;
    }

    std::string XMLDataReader::giveStackPath(){
        std::ostringstream oss;
        for(const StackItem& i: stack) oss<<"/"<<i.parent.name();
        return oss.str();
    }
    void XMLDataReader::enterGroup(const std::string& name) {
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<xmlFile<<"::"<<giveStackPath()<<": enter '"<<name<<"'");
        pugi::xml_node grp=stack.back().parent.child(name.c_str());
        XMLInputRecord::node_seen_set(grp,true);
        if(!grp){ std::cerr<<"No "<<giveStackPath()<<" // "<<name<<std::endl; abort(); OOFEM_ERROR("Error reading %s: %s has no child %s",loc(grp).c_str(),giveStackPath().c_str(),name.c_str()); }
        pugi::xml_node ch1=grp.first_child();
        if(!ch1) OOFEM_ERROR("Error reading %s: %s: %s has no child nodes",loc().c_str(),giveStackPath().c_str(),name.c_str());
        stack.push_back(StackItem{grp,ch1});
    }
    void XMLDataReader::leaveGroup(const std::string& name) {
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<xmlFile<<"::"<<giveStackPath()<<": leave ");
        if(stack.back().parent.name()!=name) OOFEM_ERROR("Error reading %s: %s: bottom-most group is %s, request to leave %s.",loc().c_str(),giveStackPath().c_str(),stack.back().parent.name(),name.c_str());
        stack.pop_back();
    }
    void XMLDataReader::enterRecord(InputRecord* rec) {
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<xmlFile<<"::"<<giveStackPath());
        XMLInputRecord* r=dynamic_cast<XMLInputRecord*>(rec);
        if(!r) OOFEM_ERROR("Error reading %s: input record is not a XMLInputRecord?",loc().c_str());
        _XML_DEBUG("   entering '"<<r->node.name()<<"' @ "<<(void*)(&r->node));
        stack.push_back(StackItem{r->node,r->node.first_child()});
    }
    void XMLDataReader::leaveRecord(InputRecord* rec) {
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<xmlFile<<"::"<<giveStackPath());
        XMLInputRecord* r=dynamic_cast<XMLInputRecord*>(rec);
        if(!r) OOFEM_ERROR("Error reading %s: input record is not a XMLInputRecord?",loc().c_str());
        _XML_DEBUG("   leaving '"<<r->node.name()<<"' @ "<<(void*)(&r->node));
        if(r->node!=stack.back().parent) OOFEM_ERROR("Error reading %s: %s: bottom-most node is %s @ %x, request to leave %s @ %x",loc().c_str(),giveStackPath().c_str(),stack.back().parent.name(),(void*)&(stack.back().parent),r->node.name(),(void*)&(r->node));
        // leave "<<std::endl;
        stack.pop_back();
    }

    void XMLDataReader::finish(){
        leaveGroup("oofem");
        // doc.print(std::cerr,"  ");
        if(stack.size()>1) OOFEM_WARNING("XML stack not popped (%d entries)",stack.size());
        pugi::xml_node n;
        for(const auto& [xml,doc]: docs){
            while((n=doc.find_node([](const pugi::xml_node& n)->bool{ return !XMLInputRecord::node_seen_get(n); }))){
                // handle Output and Description nodes (PCDATA), provided the parent tag was __seen (PCDATA node cannot be assigned attribute)
                if(n.type()==pugi::xml_node_type::node_pcdata && XMLInputRecord::node_seen_get(n.parent())){ n.parent().remove_child(n); continue; }
                std::ostringstream oss;
                oss<<"Unprocessed XML fragment "<<loc(n)<<"\n";
                n.print(oss,"  ");
                OOFEM_WARNING(oss.str().c_str());
                n.parent().remove_child(n);
            }
}
    }

    InputRecord &
    XMLDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
    {
        std::string tag=DataReader::InputRecordTags[typeId];
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<xmlFile<<"::"<<giveStackPath()<<": tag '"<<tag<<"'");
        StackItem& tip(stack.back());
        if(tag.empty()){
            // tip.curr=tip.curr.next_sibling();
            while(tip.curr && tip.seen.count(tip.curr)>0){
                _XML_DEBUG("  --- "<<tip.curr.name());
                tip.curr=tip.curr.next_sibling();
            }
            if(!tip.curr) OOFEM_ERROR("Error reading %s: %s: no more entries to read (arbitrary tag).",loc().c_str(),giveStackPath().c_str());
        } else {
            tip.curr=tip.parent.first_child();
            while(tip.curr && (tip.curr.name()!=tag || tip.seen.count(tip.curr)>0)){
                _XML_DEBUG("  --- "<<tip.curr.name());
                tip.curr=tip.curr.next_sibling();
            }
            if(!tip.curr) OOFEM_ERROR("Error reading %s: %s: no tag %s left.",loc().c_str(),giveStackPath().c_str(),tag.c_str());
        }
        _XML_DEBUG("  ==> "<<tip.curr.name());
        tip.seen.insert(tip.curr);
        tip.lastRecord=std::make_unique<XMLInputRecord>(this,tip.curr,/*tag.empty()?tip.seen.size():-1*/tip.seen.size());
        _XML_DEBUG("   tip.curr="<<tip.curr.name()<<": "<<XMLDataReader::node_seen_get(tip.curr));
        return *tip.lastRecord;
    }
} // end namespace oofem
