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

#include "xmlinputrecord.h"
#include "xmldatareader.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dictionary.h"
#include "range.h"
#include "scalarfunction.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <ostream>
#include <sstream>
#include <iostream>
#include <regex>
#include <charconv>
#include <cassert>
#include <string.h>

#ifdef _MSC_VER
    #define strcasecmp _stricmp
#endif


// #define _XML_DEBUG(m) std::cerr<<__FUNCTION__<<": "<<m<<std::endl;
#define _XML_DEBUG(m)


namespace oofem {

    bool XMLInputRecord::node_seen_get(const pugi::xml_node& n){ return !!n.attribute(SeenMark); }

    void XMLInputRecord::node_seen_set(pugi::xml_node& n, bool seen){
        if(seen && !n.attribute(SeenMark)) n.append_attribute(SeenMark);
        else if(!seen && n.attribute(SeenMark)) n.remove_attribute(SeenMark);
    }

    std::string XMLInputRecord::loc(const pugi::xml_node& node) const {
        return _reader()->loc(node);
    }

    template<typename T>
    T string_to(const std::string& s, std::function<std::string()> where){
        T val;
        const char* last=s.data()+s.size();
        auto [p,e]=std::from_chars(s.data(),last,val);
        if(p!=last) OOFEM_ERROR("%s: error parsing '%s' as typeid '%s' (leftover chars)",where().c_str(),s.c_str(),typeid(T).name());
        if(e==std::errc()) return val;
        OOFEM_ERROR("%s: error parsing '%s' as typeid '%s' (std::from_chars error).",where().c_str(),s.c_str(),typeid(T).name());
    }

    template<>
    Range string_to(const std::string& s, std::function<std::string()> where){
        if(std::regex_match(s,std::regex("[0-9]+"))){
            return Range(std::atoi(s.c_str()));
        }
        std::smatch match;
        if(std::regex_match(s,match,std::regex("\\s*([0-9]+)\\s*(-|[.]{2,3})\\s*([0-9]+)"))){
            assert(match.size()==4);
            return Range(std::atoi(match[1].str().c_str()),std::atoi(match[3].str().c_str()));
        }
        OOFEM_ERROR("%s: error parsing '%s' as range (single integer or range between two integers separated with -, .., ...).",where().c_str(),s.c_str());
    };
    template<>
    bool string_to(const std::string& s, std::function<std::string()> where){
        if(s=="0" || s=="n" || s=="N" || s=="no"  || s=="No"  || s=="NO" ){ return false; }
        if(s=="1" || s=="y" || s=="Y" || s=="yes" || s=="Yes" || s=="YES"){ return true;  }
        OOFEM_ERROR("%s: error parsing '%s' as bool (alllowed values: 0, n, N, no, No, NO; 1, y, Y, yes, Yes, YES).",where().c_str(),s.c_str())
    }


    // trim string from both sides
    std::string _lrtrim(std::string &s) { s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int c) {return !std::isspace(c);})); return s; }

    struct Tokens{
        std::string attr;
        std::string str;
        pugi::xml_node node;
        std::vector<std::string> toks;
        std::function<std::string()> loc;
        Tokens(const std::string& attr_, XMLInputRecord* rec, const char* sep_regex="\\s+"): attr(attr_) {
            std::tie(str,node)=rec->_attr_traced_read_with_node(attr.c_str());
            //std::cerr<<"attr_="<<attr_<<", attr="<<attr<<", offset="<<node.offset_debug()<<std::endl;
            loc=[rec,this](){ return rec->loc(this->node); };
            _makeToks(sep_regex);
        };
        Tokens(const std::string& attr_, const std::string& str_, std::function<std::string()> loc_, const char* sep_regex="\\s+"): attr(attr_), str(str_), loc(loc_){
            str=_lrtrim(str);
            _makeToks(sep_regex);
        }
        void _makeToks(const char* sep_regex){
            if(str.empty()) return; // this would create spurious empty token, stay on zero size instead
            const std::regex sep_re(sep_regex);
            std::copy(std::sregex_token_iterator(str.begin(), str.end(), sep_re, -1),
                      std::sregex_token_iterator(),
                      std::back_inserter(toks)
            );
        };
        size_t size() { return toks.size(); }
        void assertSize(size_t req){ if(size()!=req) OOFEM_ERROR("%s: attribute %s: length mismatch (%d items, %d required)",loc().c_str(),attr.c_str(),(int)size(),(int)req); }
        template<typename T> T as(size_t ix){
            if(ix>size()) OOFEM_ERROR("%s: attribute %s: invalid index %d in sequence of length %d: %s",loc().c_str(),attr.c_str(),(int)ix,(int)size(),str.c_str());
            const std::string& s(toks[ix]);
            return string_to<T>(s,[this,ix](){ return loc()+": attribute "+attr+" ["+std::to_string(ix)+"]"; });
        }
    };

    XMLInputRecord :: XMLInputRecord(XMLDataReader* reader_, const pugi::xml_node& node_): InputRecord((DataReader*)reader_), node(node_) {
        node_seen_set(node,true);
        _XML_DEBUG(loc()<<": node.name()="<<node.name());
    }

    int XMLInputRecord::setRecId(int lastRecId){
        this->recId=lastRecId+1;
        giveOptionalField(this->recId,"id");
        if(lastRecId>0 && this->recId<=lastRecId) OOFEM_WARNING("%s: descencing ids (previous %d, now %d)",loc().c_str(),lastRecId,recId);
        return this->recId;
    };

    int XMLInputRecord::giveGroupCount(InputFieldType id, const std::string& name, bool optional){
        _XML_DEBUG("id="<<id<<", name="<<name<<", optional="<<optional);
        _XML_DEBUG(loc(node)<<", node.name()="<<node.name());
        pugi::xml_node ch=_reader()->giveNamedChild(node,name);
        if(!ch){
            if(optional) return 0; // return DataReader::NoSuchGroup;
            OOFEM_ERROR("%s: %s has no child node %s.",loc().c_str(),node.name(),name.c_str());
        }
        int ret=0;
        node_seen_set(ch,true);
        for([[maybe_unused]] pugi::xml_node n: ch.children()) ret++;
        _XML_DEBUG(node.name()<<" has "<<ret<<" children.");
        return ret;
    }
    bool XMLInputRecord::hasChild(InputFieldType id, const std::string& name, bool optional){
        bool has=!!_reader()->giveNamedChild(node,name);
        if(!has && !optional) OOFEM_ERROR("%s: %s has no no child node '%s' (required).",loc().c_str(),node.name(),name.c_str());
        return has;
    }

    void
    XMLInputRecord :: giveRecordKeywordField(std :: string &answer){
        _XML_DEBUG("node.name()="<<node.name());
        if(node.attribute("type")) answer=_attr_traced_read("type");
        else answer=node.name();
    }
    void XMLInputRecord::giveRecordKeywordField(std::string& answer, int& value){
        _XML_DEBUG("node.name()="<<node.name());
        answer=node.name();
        value=recId;
    }

    bool XMLInputRecord::hasField(InputFieldType id0){
        std::string id=xmlizeAttrName(id0);
        pugi::xml_attribute att=node.attribute(id.c_str());
        if(!att){            // retry case-insensitive
            for(att=node.first_attribute(); att; att=att.next_attribute()){
                if(strcasecmp(att.name(),id.c_str())==0){
                    std::cerr<<"XML: case-insensitive hasField ('"<<att.name()<<"', requested '"<<id<<"')"<<std::endl;
                    break;
                }
            }
        }
        if(att){ attrQueried.insert(att.name()); }
        return !!att;
    }
    std::string XMLInputRecord::xmlizeAttrName(const std::string& s){
        std::string n2(s);
        for(size_t i=0; i<n2.size(); i++) if(n2[i]=='(' || n2[i]==')' || n2[i]=='/') n2[i]='_';
        return n2;
    }

    std::tuple<std::string,pugi::xml_node> XMLInputRecord::_attr_traced_read_with_node(const char* name){
        std::string n2=xmlizeAttrName(std::string(name));
        pugi::xml_attribute att=node.attribute(n2.c_str());
        if(!att){
            // retry case-insensitive
            for (att=node.first_attribute(); att; att=att.next_attribute()){
                #ifdef _MSC_VER
                    #define strcasecmp _stricmp
                #endif
                if(strcasecmp(att.name(),name)==0){
                    std::cerr<<"XML: case-insensitive match ('"<<att.name()<<"', requested '"<<name<<"')"<<std::endl;
                    break;
                }
            }
        }
        if(!att) OOFEM_ERROR("%s: no such attribute: %s",loc().c_str(),n2.c_str());
        std::string ret=att.as_string();
        attrRead.insert(att.name());
        return std::make_tuple(ret,node);
    }

    void XMLInputRecord::finish(bool wrn) {
        _XML_DEBUG(loc());
        pugi::xml_document tmp;
        pugi::xml_node xNotseen=tmp.append_child(node.name());
        pugi::xml_node xNotempty=tmp.append_child(node.name());
        int nNotseen=0, nNotempty=0;
        for([[maybe_unused]] const pugi::xml_attribute& a: node.attributes()) {
            if(std::string(a.name())==SeenMark) continue;
            bool queried(attrQueried.count(a.name())), read(attrRead.count(a.name())), hasValue=(a.value()[0]!=0);
            if(read) continue;
            if(!queried){ xNotseen.append_attribute(a.name())=a.value(); nNotseen++; continue; }
            if(/* !read && queried && */ hasValue){ xNotempty.append_attribute(a.name())=a.value(); nNotempty++; }
        }
        if(nNotseen==0 && nNotempty==0) return;
        if(!wrn) return;
        std::ostringstream oss; oss<<loc();
        if(nNotseen){ oss<<"\n   attribute"<<(nNotseen>1?"s":"")<<" ignored:\n      "; xNotseen.print(oss); }
        if(nNotempty){ oss<<"\n   attribute"<<(nNotempty>1?"s":"")<<" with ignored non-empty value:\n      "; xNotempty.print(oss); }
        OOFEM_WARNING("%s",oss.str().c_str());
    }
    std::string XMLInputRecord::giveRecordAsString() const {
        pugi::xml_document tmp;
        pugi::xml_node n=tmp.append_child(node.name());
        for([[maybe_unused]] const pugi::xml_attribute& a: node.attributes()) {
            if(std::string(a.name())==XMLInputRecord::SeenMark) continue;
            n.append_attribute(a.name())=a.value();
        }
        std::ostringstream oss; n.print(oss,"",pugi::format_raw);
        return oss.str();
    }
    std::string XMLInputRecord::giveRecordInTXTFormat() const {
        OOFEM_ERROR("%s: not (yet?) implemented.",__PRETTY_FUNCTION__);
    }

    void XMLInputRecord::giveField(std::string& answer, InputFieldType id){
        pugi::xml_node node;
        std::tie(answer,node)=_attr_traced_read_with_node(id);
        _XML_DEBUG(loc(node)<<": "<<node.name()<<"::"<<id<<"="<<answer);
    }
    void XMLInputRecord::giveField(FloatArray &answer, InputFieldType id){
        Tokens tt(id,this);
        // std::cerr<<id<<": FloatArray from '"<<tt.str<<"' ("<<tt.size()<<" items)"<<std::endl;
        answer.resize(tt.size());
        for(size_t i=0; i<tt.size(); i++) answer[i]=tt.as<double>(i);
        _XML_DEBUG(tt.loc()<<": parsed attribute "<<id<<" as "<<answer);
    }
    void XMLInputRecord::giveField(FloatMatrix &answer, InputFieldType id){
        Tokens trows(id,this,"\\s*;\\s*");
        int rows=trows.size();
        for(int row=0; row<rows; row++){
            Tokens tcols(id,trows.toks[row],[this,trows](){ return this->loc(trows.node); },"\\s+");
            if(row==0){ answer.resize(rows,tcols.size()); }
            else if((int)answer.cols()!=(int)tcols.size()) OOFEM_ERROR("%s: row %d has inconsistent number of columns (%d != %d)",loc().c_str(),rows,(int)answer.cols(),(int)tcols.size());
            for(int col=0; col<answer.cols(); col++){ answer(row,col)=tcols.as<double>(col); }
        }
        _XML_DEBUG(tt.loc()<<": parsed attribute "<<id<<" as "<<answer);
    }
    void XMLInputRecord::giveField(IntArray &answer, InputFieldType id){
        Tokens tt(id,this);
        //std::cerr<<id<<": IntArray from '"<<tt.str<<" ("<<tt.size()<<" items"<<std::endl;
        answer.resize(tt.size());
        for(size_t i=0; i<tt.size(); i++) answer[i]=tt.as<int>(i);
        _XML_DEBUG(tt.loc()<<": parsed attribute "<<id<<" as "<<answer);
    }
    void XMLInputRecord::giveField(double& answer, InputFieldType id){
        std::string s; pugi::xml_node n;
        std::tie(s,n)=_attr_traced_read_with_node(id);
        answer=string_to<double>(s,[this,n,id](){ return loc(n)+": attribute '"+id+"'"; });
    }
    void XMLInputRecord::giveField(int& answer, InputFieldType id){
        std::string s; pugi::xml_node n;
        std::tie(s,n)=_attr_traced_read_with_node(id);
        answer=string_to<int>(s,[this,n,id](){ return loc(n)+": attribute '"+id+"'"; });
    }
    void XMLInputRecord::giveField(bool& answer, InputFieldType id){
        if(!hasField(id)){ answer=false; return; }
        std::string s; pugi::xml_node n;
        std::tie(s,n)=_attr_traced_read_with_node(id);
        answer=string_to<bool>(s,[this,n,id](){ return loc(n)+": attribute '"+id+"'"; });
    }
    void XMLInputRecord::giveField(std::list<Range>& answer, InputFieldType id){
        Tokens tt(id,this,"\\s*,\\s*");
        for(size_t i=0; i<tt.size(); i++) answer.push_back(tt.as<Range>(i));
    }
    void XMLInputRecord::giveField(ScalarFunction& answer, InputFieldType id){
        std::string s; pugi::xml_node n;
        std::tie(s,n)=_attr_traced_read_with_node(id);
        auto where=[this,n,id](){ return loc(n)+": attribute '"+id+"'"; };
        if(s[0]=='@') answer.setReference(string_to<int>(s.substr(1,s.size()-1),where));
        else if(s[0]=='$'){ std::string s2=s.substr(1,s.size()-2); answer.setSimpleExpression(s2); }
        else answer.setValue(string_to<double>(s,where));
    }
    void XMLInputRecord::giveField(Dictionary &answer, InputFieldType id){
        std::string s; pugi::xml_node n;
        std::tie(s,n)=_attr_traced_read_with_node(id);
        Tokens items(id,this,"\\s*;\\s*");
        for(const std::string& tok: items.toks){
            Tokens kv(id,tok,[this,items](){ return this->loc(items.node); },/*sep_regex*/"\\s+");
            if(kv.size()!=2) OOFEM_ERROR("%s: dictionary items must have 2 whitespace-separated fields (%d tokens found)",loc(n).c_str(),(int)kv.size());
            int key;
            if(std::regex_match(kv.toks[0],std::regex("[0-9]+"))){ key=std::atoi(kv.toks[0].c_str()); }
            else if(kv.toks[0].size()==1){ key=(char)kv.toks[0][0]; }
            else OOFEM_ERROR("%s: dictionary key must be integer or single letter (not '%s')",loc().c_str(),kv.toks[0].c_str());
            answer.add(key,kv.as<double>(1));
        }

    }

#if 0
    void
    XMLInputRecord :: giveField(std :: vector< std :: string > &answer, InputFieldType id)
    {
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            int size;
            setReadFlag(indx);
            auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
            if ( ptr == nullptr || *ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }
            answer.reserve(size);
            setReadFlag(indx);
            for ( int i = 1; i <= size; i++ ) {
                answer.push_back( tokenizer.giveToken(indx + i) );
                setReadFlag(indx + i);
            }

        } else {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }

#endif
} // end namespace oofem
