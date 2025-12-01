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
#include <cassert>
#include <cctype>
#include <ostream>
#include <sstream>
#include <iostream>
#include <regex>
#include <charconv>

// #define _XML_DEBUG(m) std::cerr<<__FUNCTION__<<": "<<m<<std::endl;
#define _XML_DEBUG(m)


namespace oofem {


    bool XMLInputRecord::node_seen_get(const pugi::xml_node& n){ return !!n.attribute(SeenMark); }

    void XMLInputRecord::node_seen_set(pugi::xml_node& n, bool seen){
        if(seen && !n.attribute(SeenMark)) n.append_attribute(SeenMark);
        else if(!seen && n.attribute(SeenMark)) n.remove_attribute(SeenMark);
    }

    std::string XMLInputRecord::loc(const pugi::xml_node& node){
        return _reader()->loc(node);
    }

    template<typename T>
    T string_to(const std::string& s, std::function<std::string()> where){
        T val;
        const char* last=s.data()+s.size();
        auto [p,e]=std::from_chars(s.data(),last,val);
        if(p!=last) OOFEM_ERROR("%s: error parsing %s as %s (leftover chars)",where().c_str(),s.c_str(),typeid(T).name());
        if(e==std::errc()) return val;
        OOFEM_ERROR("%s: error parsing %s (from_chars<%s> error)",where().c_str(),s.c_str(),typeid(T).name());
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
        OOFEM_ERROR("%s: error parsing %s as range (single integer or range between two integers separated with -, .., ...)",where().c_str(),s.c_str());
    };



    struct Tokens{
        std::string attr;
        std::string str;
        std::vector<std::string> toks;
        std::function<std::string()> loc;
        Tokens(const std::string& attr_, XMLInputRecord* rec, const char* sep_regex="\\s+"): attr(attr_) {
            pugi::xml_node node;
            std::tie(str,node)=rec->_attr_traced_read_with_node(attr.c_str());
            //std::cerr<<"attr_="<<attr_<<", attr="<<attr<<", offset="<<node.offset_debug()<<std::endl;
            loc=[rec,node](){ return rec->loc(node); };
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

    XMLInputRecord :: XMLInputRecord(XMLDataReader* reader_, const pugi::xml_node& node_, int ordinal_): InputRecord((DataReader*)reader_), node(node_), ordinal(ordinal_) {
        node_seen_set(node,true);
        _XML_DEBUG(loc()<<": node.name()="<<node.name());
    }

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
        if(!has && !optional) OOFEM_ERROR("...");
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
        pugi::xml_attribute id=node.attribute("id");
        if(!id){
            if(ordinal<0) OOFEM_ERROR("%s: node %s: id='...' not specified (and not tracked automatically).",loc().c_str(),node.name());
            value=ordinal;
        } else {
            std::string val; pugi::xml_node n;
            std::tie(val,n)=_attr_traced_read_with_node("id");
            value=string_to<int>(val,[this,n](){ return loc(n)+" attribute 'id'"; });
            if(ordinal>0 && ordinal!=value) OOFEM_ERROR("%s: node %s: id='%d' but node position is %d",loc(n).c_str(),node.name(),value,ordinal);
        }
    }

    bool XMLInputRecord::hasField(InputFieldType id){
        pugi::xml_attribute attr=node.attribute(id);
        if(attr){ attrSeen.insert(id); }
        return !!attr;
    }
    std::tuple<std::string,pugi::xml_node> XMLInputRecord::_attr_traced_read_with_node(const char* name){
        std::string n2(name);
        for(size_t i=0; i<n2.size(); i++) if(n2[i]=='(' || n2[i]==')') n2[i]='_';
        pugi::xml_attribute att=node.attribute(n2.c_str());
        if(!att) OOFEM_ERROR("%s: no such attribute: %s",loc().c_str(),n2.c_str());
        std::string ret=att.as_string();
        attrSeen.insert(n2);
        return std::make_tuple(ret,node);
    }

    void XMLInputRecord::finish(bool wrn) {
        int aLeft=0;
        for([[maybe_unused]] const pugi::xml_attribute& a: node.attributes()) {
            if(std::string(a.name())==SeenMark) continue;
            if(attrSeen.count(a.name())==0) aLeft++;
        }
        for(const std::string& n: attrSeen) node.remove_attribute(n.c_str());
        if(!wrn) return;
        if(aLeft==0) return;
        node.remove_children();
        std::ostringstream oss;
        oss<<" "<<loc()<<": unprocessed XML attributes:\n";
        node_seen_set(node,false);
        node.print(oss,"  ",pugi::format_default,pugi::encoding_auto,1);
        node_seen_set(node,true);
        OOFEM_WARNING("%s",oss.str().c_str());
    }

    void XMLInputRecord::giveField(std::string& answer, InputFieldType id){
        pugi::xml_node node;
        std::tie(answer,node)=_attr_traced_read_with_node(id);
        _XML_DEBUG(loc(node)<<": "<<node.name()<<"::"<<id<<"="<<answer);
    }
    void XMLInputRecord::giveField(FloatArray &answer, InputFieldType id){
        Tokens tt(id,this);
        answer.resize(tt.size());
        for(size_t i=0; i<tt.size(); i++) answer[i]=tt.as<double>(i);
        _XML_DEBUG(tt.loc()<<": parsed attribute "<<id<<" as "<<answer);
    }
    void XMLInputRecord::giveField(IntArray &answer, InputFieldType id){
        Tokens tt(id,this);
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
        if(s=="no" || s=="0" || s=="n" || s=="N" || s=="No" || s=="NO") OOFEM_WARNING("%s: %s='%s' is interpreted as TRUE (omit the attribute for false)",loc(n).c_str(),id,s.c_str());
        answer=true;
    }
    void XMLInputRecord::giveField(std::list<Range>& answer, InputFieldType id){
        Tokens tt(id,this,"\\s*,\\s*");
        for(size_t i=0; i<tt.size(); i++) answer.push_back(tt.as<Range>(i));
    }



#if 0
    void
    XMLInputRecord :: giveField(FloatMatrix &answer, InputFieldType id)
    {
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            int nrows, ncols;
            setReadFlag(indx);

            auto ptr = scanInteger(tokenizer.giveToken(++indx), nrows);
            if ( ptr == nullptr || *ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);
            ptr = scanInteger(tokenizer.giveToken(++indx), ncols);
            if ( ptr == nullptr || *ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);

            if ( readMatrix(tokenizer.giveToken(++indx), nrows, ncols, answer) == 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);
        } else {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }

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

    void
    XMLInputRecord :: giveField(Dictionary &answer, InputFieldType id)
    {
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            setReadFlag(indx);
            int size;
            auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
            if ( ptr == nullptr || *ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);

            answer.clear();
            for ( int i = 1; i <= size; i++ ) {
                int key = 0;
                const char * token = tokenizer.giveToken(++indx);
                auto ptr1 = scanInteger( token, key );
                if ( ptr1 == nullptr || *ptr1 != 0 ) {
                    key = token [ 0 ];
                    // throw BadFormatInputException(*this, id, lineNumber);
                }
                double value;
                setReadFlag(indx);
                auto ptr = scanDouble(tokenizer.giveToken(++indx), value);
                if ( ptr == nullptr || *ptr != 0 ) {
                    throw BadFormatInputException(*this, id, lineNumber);
                }

                setReadFlag(indx);
                answer.add(key, value);
            }
        } else {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }

    void
    XMLInputRecord :: giveField(ScalarFunction &answer, InputFieldType id)
    {
        const char *rec;
        int indx = this->giveKeywordIndx(id);

        if ( indx ) {
            setReadFlag(indx);
            rec = tokenizer.giveToken(++indx);
            if ( * rec == '@' ) {
                // reference to function
                int refVal;
                auto ptr = scanInteger(rec + 1, refVal);
                if ( ptr == nullptr || *ptr != 0 ) {
                    throw BadFormatInputException(*this, id, lineNumber);
                }
                setReadFlag(indx);
                answer.setReference(refVal);
            } else if ( * rec == '$' ) {
                // simple expression
                std :: string expr;

                expr = std :: string( tokenizer.giveToken(indx) );
                setReadFlag(indx);
                std :: string _v = expr.substr(1, expr.size() - 2);

                answer.setSimpleExpression(_v); // get rid of enclosing '$'
            } else {
                double val;
                auto ptr = scanDouble(tokenizer.giveToken(indx), val);
                if ( ptr == nullptr || *ptr != 0 ) {
                    throw BadFormatInputException(*this, id, lineNumber);
                }

                setReadFlag(indx);
                answer.setValue(val);
            }
        } else {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }
#endif
} // end namespace oofem
