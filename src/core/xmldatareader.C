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

#include "xmldatareader.h"
#include "error.h"

#include <string>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <limits>
#include <cstring>
#include <filesystem>
#include <pugixml.hpp>
#include <bits/stdc++.h>
#include <string.h>

#ifdef _MSC_VER
#define strcasecmp _stricmp
#endif


// #define _XML_DEBUG(m) std::cerr<<std::string(__FUNCTION__).substr(0,20)<<": "<<m<<std::endl;
#define _XML_DEBUG(m)


namespace oofem {
    pugi::xml_document& XMLDataReader::loadXml(pugi::xml_node parent, const std::string& xml_){
        std::filesystem::path xml=(xml_);
        if(xml.is_relative() && parent){
            std::filesystem::path parentXml(xmlFiles[parent.root()]);
            xml=parentXml.parent_path()/xml;
            _XML_DEBUG("Relative path "<<xml_<<" resolved to "<<xml<<" (parent is "<<parentXml<<")");
        }
        // trasverse just to get line breaks
        std::ifstream i(xml,std::ifstream::in|std::ios::binary);
        if(!i.is_open()) OOFEM_ERROR("Error opening %s: %s",xml.c_str(),std::strerror(errno));
        std::list<size_t> br;
        size_t pos=0;
        while(i.good()){ char c=i.get(); pos++; if(c=='\n') br.push_back(pos); }
        std::vector<size_t> nl;
        nl.reserve(br.size());
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Waggressive-loop-optimizations"
            nl.assign(br.begin(),br.end());
        #pragma GCC diagnostic pop
        _XML_DEBUG(xml<<": "<< nl.size()<<" newlines found, about to parse the XML...");
        docs[parent]=pugi::xml_document();
        if(parent) XMLInputRecord::node_seen_set(parent,true);
        auto& doc=docs[parent];
        pugi::xml_parse_result result=doc.load_file(xml.c_str());
        if(!result){
            auto [line,col]=offset2lc(nl,result.offset);
            OOFEM_ERROR("Error parsing %s:%ld:%ld: %s",xml.c_str(),line,col,result.description());
        }
        newlines[doc]=std::move(nl);
        xmlFiles[doc]=xml;
        size_t nChildren=std::distance(doc.begin(),doc.end());
        if(nChildren!=1) OOFEM_ERROR("Error reading %s: must have exactly 1 top-level tag (not %ld)",xml.c_str(),nChildren);
        _XML_DEBUG(xml<<": "<<" parsing done, root is '"<<doc.first_child().name()<<"'");
        return doc;
    }

    pugi::xml_node XMLDataReader::resolveXiInclude(const pugi::xml_node& n){
        if(n.name()!=XiIncludeTag) OOFEM_ERROR("Error resolving %s: must be xi:include (not %s)",loc().c_str(),n.name());
        if(auto xi=docs.find(n); xi!=docs.end()) return xi->second.first_child();
        pugi::xml_attribute href=n.attribute("href");
        if(!href) OOFEM_ERROR("%s: xi:include must have 'href' attribute.",loc().c_str());
        _XML_DEBUG(loc()<<":"<<giveStackPath()<<": reading xi:include, href='"<<href.value()<<"'");
        return loadXml(n,href.value()).first_child();
    }

    pugi::xml_node XMLDataReader::giveNamedChild(const pugi::xml_node& parent, const std::string& name){
        _XML_DEBUG("finding child named '"<<name<<"'");
        for(pugi::xml_node c=parent.first_child(); c; c=c.next_sibling()){
            if(c.name()==name){
                _XML_DEBUG("returning direct child "<<loc(c));
                return c;
            }
            if(strcasecmp(c.name(),name.c_str())==0){
                std::cerr<<loc(parent)<<": case-insensitive XML tag match ('"<<c.name()<<"', requested '"<<name<<"')"<<std::endl;
                break;
            }
            if(c.name()==XiIncludeTag){
                pugi::xml_node xi=resolveXiInclude(c);
                if(xi.name()==name){
                    _XML_DEBUG("returning xi:include "<<loc(c));
                    return xi;
                }
            }
        }
        return pugi::xml_node();
    }

    XMLDataReader::XMLDataReader(const std::string& xmlFile_): topXmlFile(xmlFile_) {
        pugi::xml_document& doc=loadXml(pugi::xml_node(),topXmlFile);
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
        topRecord=std::shared_ptr<InputRecord>(new XMLInputRecord(this,root));
    }

    bool XMLDataReader :: canRead(const std::string& xml){
        return std::filesystem::path(xml).extension() == ".xml";
    }

    std::tuple<size_t,size_t> XMLDataReader::offset2lc(const std::vector<size_t>& nl, size_t offset) const {
        if(nl.empty()) return std::make_tuple(0,0);
        size_t ix=std::distance(nl.begin(),std::lower_bound(nl.begin(),nl.end(),offset));
        // _XML_DEBUG("offset="<<offset<<",ix="<<ix); // comment out so that it does not damage most other _XML_DEBUG messages
        return std::make_tuple(ix+1,offset-(ix==0?0:nl[ix-1]));
    }
    std::string XMLDataReader::loc() const { return loc(stack.back().curr?stack.back().curr:stack.back().parent); }
    std::string XMLDataReader::loc(const pugi::xml_node& n) const {
        std::ostringstream oss;
        auto [line,col]=offset2lc(newlines.at(n.root()),n.offset_debug());
        oss<<xmlFiles.at(n.root())<<":"<<line<<":"<<col;
        return oss.str();
    }

    int XMLDataReader::giveGroupCount(const std::string& name) {
        if(name.empty()) return giveCurrentGroupCount();
        enterGroup(name);
        int ret=this->giveCurrentGroupCount();
        leaveGroup(name);
        _XML_DEBUG(giveStackPath()<<" // "<<name<<": size is "<<ret);
        return ret;
    }

    int XMLDataReader::giveCurrentGroupCount(){
        int ret=0;
        for([[maybe_unused]] pugi::xml_node n: stack.back().parent.children()) ret++;
        _XML_DEBUG(giveStackPath()<<": size is "<<ret);
        return ret;
    }

    std::string XMLDataReader::giveStackPath(){
        std::ostringstream oss;
        for(const StackItem& i: stack) oss<<"/"<<i.parent.name();
        return oss.str();
    }
    void XMLDataReader::enterGroup(const std::string& name) {
        _XML_DEBUG("::"<<giveStackPath()<<": enter '"<<name<<"'");
        pugi::xml_node grp=giveNamedChild(stack.back().parent,name);
        XMLInputRecord::node_seen_set(grp,true);
        if(!grp){ std::cerr<<"No "<<giveStackPath()<<" // "<<name<<std::endl; OOFEM_ERROR("Error reading %s: %s has no child %s",loc(grp).c_str(),giveStackPath().c_str(),name.c_str()); }
        // stack.back().seen.insert(grp);
        pugi::xml_node ch1=grp.first_child();
        if(!ch1) OOFEM_ERROR("Error reading %s: %s: %s has no child nodes",loc().c_str(),giveStackPath().c_str(),name.c_str());
        stack.push_back(StackItem{grp,ch1});
    }
    void XMLDataReader::leaveGroup(const std::string& name) {
        _XML_DEBUG(loc()<<"::"<<giveStackPath()<<": leave "<<name);
        // if an exception is being propagated, this error would hide it
        if(!std::uncaught_exceptions()){
            if(stack.back().parent.name()!=name) OOFEM_ERROR("Error reading %s: %s: bottom-most group is %s, request to leave %s.",loc().c_str(),giveStackPath().c_str(),stack.back().parent.name(),name.c_str());
        }
        stack.pop_back();
    }
    void XMLDataReader::enterRecord(InputRecord* rec) {
        _XML_DEBUG(loc()<<"::"<<giveStackPath());
        XMLInputRecord* r=dynamic_cast<XMLInputRecord*>(rec);
        if(!r) OOFEM_ERROR("Error reading %s: input record is not a XMLInputRecord?",loc().c_str());
        _XML_DEBUG("   entering '"<<r->node.name()<<"' @ "<<(void*)(&r->node));
        stack.push_back(StackItem{r->node,r->node.first_child()});
    }
    void XMLDataReader::leaveRecord(InputRecord* rec) {
        _XML_DEBUG(loc()<<"::"<<giveStackPath());
        XMLInputRecord* r=dynamic_cast<XMLInputRecord*>(rec);
        if(!r) OOFEM_ERROR("Error reading %s: input record is not a XMLInputRecord?",loc().c_str());
        _XML_DEBUG("   leaving '"<<r->node.name()<<"' @ "<<(void*)(&r->node));
        // if an exception is being propagated, and enterRecord/leaveRecord were not used via RAII-guard (RecordGuard),
        // our would abort immediately (with misleading error message) and the exception would never make it to the handler
        if(!std::uncaught_exceptions()){
            if(r->node!=stack.back().parent) OOFEM_ERROR("Error reading %s: %s: bottom-most node is %s @ %p, request to leave %s @ %p",loc().c_str(),giveStackPath().c_str(),stack.back().parent.name(),(void*)&(stack.back().parent),r->node.name(),(void*)&(r->node));
        }
        stack.pop_back();
    }

    void XMLDataReader::finish(){
        leaveGroup("oofem");
        // doc.print(std::cerr,"  ");
        if(stack.size()>1) OOFEM_WARNING("XML stack not popped (%ld entries)",stack.size());
        pugi::xml_node n;
        for(const auto& [xml,doc]: docs){
            while((n=doc.find_node([](const pugi::xml_node& n)->bool{ return !XMLInputRecord::node_seen_get(n); }))){
                // handle Output and Description nodes (PCDATA), provided the parent tag was __seen (PCDATA node cannot be assigned attribute)
                if(n.type()==pugi::xml_node_type::node_pcdata && XMLInputRecord::node_seen_get(n.parent())){ n.parent().remove_child(n); continue; }
                std::ostringstream oss;
                oss<<"Unprocessed XML fragment "<<loc(n)<<"\n";
                n.print(oss,"  ");
                OOFEM_WARNING("%s",oss.str().c_str());
                n.parent().remove_child(n);
            }
        }
    }

    InputRecord &
    XMLDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
    {
        std::string tag=DataReader::InputRecordTags[typeId];
        _XML_DEBUG(loc()<<"::"<<giveStackPath()<<": tag '"<<tag<<"'");
        StackItem& tip(stack.back());
        if(tag.empty()){
            tip.curr=tip.parent.first_child();
            while(tip.curr && tip.seen.count(tip.curr)>0){
                _XML_DEBUG("  --- "<<tip.curr.name());
                tip.curr=tip.curr.next_sibling();
            }
            if(!tip.curr) OOFEM_ERROR("Error reading %s: %s: no more entries to read (arbitrary tag).",loc().c_str(),giveStackPath().c_str());
        } else {
            tip.curr=tip.parent.first_child();
            while(true){
                if(!tip.curr) break;
                if(tip.seen.count(tip.curr)==0){
                    if(tip.curr.name()==XiIncludeTag && resolveXiInclude(tip.curr).name()==tag) break;
                    else if(tip.curr.name()==tag) break;
                    else if(strcasecmp(tip.curr.name(),tag.c_str())==0){
                        std::cerr<<loc(tip.curr)<<": case-insensitive XML tag match ('"<<tip.curr.name()<<"', requested '"<<tag<<"')"<<std::endl;
                        break;
                    }
                }
                tip.curr=tip.curr.next_sibling();
            }
            if(!tip.curr) OOFEM_ERROR("Error reading %s: %s: no tag %s left.",loc().c_str(),giveStackPath().c_str(),tag.c_str());
        }
        _XML_DEBUG("  ==> "<<tip.curr.name());
        tip.seen.insert(tip.curr);
        pugi::xml_node n=((tip.curr.name()==XiIncludeTag) ? resolveXiInclude(tip.curr) : tip.curr);
        // empty tag means all nodes are traversed; in that case it is plus minus safe to assume that previous record was already processed
        if(tip.lastRecord && tag.empty()) tip.lastRecord->finish(); //for checking everything has been read before moving onto the next one
        tip.lastRecord=std::make_shared<XMLInputRecord>(this,n);
        tip.lastRecId=tip.lastRecord->setRecId(tip.lastRecId);
        _XML_DEBUG("   tip.curr="<<tip.curr.name()<<": "<<XMLInputRecord::node_seen_get(tip.curr));
        return *tip.lastRecord;
    }
} // end namespace oofem
