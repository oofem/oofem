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
#include <cctype>
#include <ostream>
#include <sstream>
#include <iostream>
#include <regex>
#include <charconv>

namespace oofem {


    bool XMLInputRecord::node_seen_get(const pugi::xml_node& n){ return !!n.attribute(SeenMark); }

    void XMLInputRecord::node_seen_set(pugi::xml_node& n, bool seen){
        if(seen && !n.attribute(SeenMark)) n.append_attribute(SeenMark);
        else if(!seen && n.attribute(SeenMark)) n.remove_attribute(SeenMark);
    }

    std::string XMLInputRecord::loc(const pugi::xml_node& node){
        return ((XMLDataReader*)this->giveReader())->loc(node);
    }

    template<typename T>
    T string_to(const std::string& s, std::function<std::string()> where){
        T val;
        const char* last=s.data()+s.size();
        auto [p,e]=std::from_chars(s.data(),last,val);
        if(p!=last) OOFEM_ERROR("%s: error parsing %s as %s (leftover chars)",where().c_str(),s.c_str(),typeid(T).name());
        if(e==std::errc()) return val;
        OOFEM_ERROR("%s: error parsing %s (from_chars error)",where().c_str(),s.c_str(),typeid(T).name());
    }


    struct Tokens{
        std::string attr;
        std::string str;
        std::vector<std::string> toks;
        std::function<std::string()> loc;
        Tokens(const std::string& attr_, XMLInputRecord* rec): attr(attr_) {
            pugi::xml_node node;
            std::tie(str,node)=rec->_attr_traced_read_with_node(attr.c_str());
            //std::cerr<<"attr_="<<attr_<<", attr="<<attr<<", offset="<<node.offset_debug()<<std::endl;
            loc=[rec,node](){ return rec->loc(node); };
            const std::regex ws_re("\\s+");
            std::copy(std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1),
                      std::sregex_token_iterator(),
                      std::back_inserter(toks)
            );
        };
        size_t size() { return toks.size(); }
        void assertSize(size_t req){ if(size()!=req) OOFEM_ERROR("%s: attribute %s: length mismatch (%d items, %d required)",loc().c_str(),attr.c_str(),size(),req); }
        template<typename T> T as(size_t ix){
            if(ix>size()) OOFEM_ERROR("%s: attribute %s: invalid index %d in sequence of length %d: %s",loc().c_str(),attr.c_str(),ix,size(),str.c_str());
            const std::string& s(toks[ix]);
            return string_to<T>(s,[this,ix](){ return loc()+": attribute "+attr+" ["+std::to_string(ix)+"]"; });
        }
    };

    XMLInputRecord :: XMLInputRecord(XMLDataReader* reader_, const pugi::xml_node& node_, int ordinal_): InputRecord((DataReader*)reader_), node(node_), ordinal(ordinal_) {
        node_seen_set(node,true);
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<loc()<<": node.name()="<<node.name());
    }

    int XMLInputRecord::giveGroupCount(InputFieldType id, const std::string& name, bool optional){
        pugi::xml_node ch=node.child(name.c_str());
        if(!ch){
            _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<node.name()<<" has no children (optional="<<optional<<".");
            if(optional) return 0; // return DataReader::NoSuchGroup;
            OOFEM_ERROR("%s: %s has no child node %s.",loc().c_str(),node.name(),name.c_str());
        }
        int ret=0;
        node_seen_set(ch,true);
        for([[maybe_unused]] pugi::xml_node n: ch.children()) ret++;
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<node.name()<<" has "<<ret<<" children.");
        return ret;
    }
    bool XMLInputRecord::hasChild(InputFieldType id, const std::string& name, bool optional){
        bool has=!!node.child(name.c_str());
        if(!has && !optional) OOFEM_ERROR("...");
        return has;
    }

    void
    XMLInputRecord :: giveRecordKeywordField(std :: string &answer){
        _XML_DEBUG(__PRETTY_FUNCTION__<<": node.name()="<<node.name());
        if(node.attribute("type")) answer=_attr_traced_read("type");
        else answer=node.name();
    }
    void XMLInputRecord::giveRecordKeywordField(std::string& answer, int& value){
        _XML_DEBUG(__PRETTY_FUNCTION__<<": node.name()="<<node.name());
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
        // node.remove_attribute(att); std::cerr<<" -- <"<<node.name()<<" "<<name<<">: destroyed"<<std::endl;
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
        std::ostringstream oss;
        oss<<"Unprocessed XML attributes:\n";
        node_seen_set(node,false);
        node.print(oss,"  ");
        node_seen_set(node,true);
        OOFEM_WARNING(oss.str().c_str());
    }

    void XMLInputRecord::giveField(std::string& answer, InputFieldType id){
        pugi::xml_node node;
        std::tie(answer,node)=_attr_traced_read_with_node(id);
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<loc(*node)<<": "<node.name()<<"::"<<id<<"="<<answer);
    }


    void XMLInputRecord::giveField(FloatArray &answer, InputFieldType id){
        Tokens tt(id,this);
        size_t len=tt.as<int>(0);
        tt.assertSize(len+1);
        answer.resize(len);
        for(size_t i=0; i<len; i++) answer[i]=tt.as<double>(i+1);
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<tt.loc()<<": parsed attribute "<<id<<" as "<<answer);
    }
    void XMLInputRecord::giveField(IntArray &answer, InputFieldType id){
        Tokens tt(id,this);
        size_t len=tt.as<int>(0);
        tt.assertSize(len+1);
        answer.resize(len);
        for(size_t i=0; i<len; i++) answer[i]=tt.as<int>(i+1);
        _XML_DEBUG(__PRETTY_FUNCTION__<<": "<<tt.loc()<<": parsed attribute "<<id<<" as "<<answer);
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



#if 0
    void
    XMLInputRecord :: giveRecordKeywordField(std :: string &answer)
    {
        std::cerr<<__PRETTY_FUNCTION__<<": node.name()="<<node.name()<<std::endl;
        answer=node.name();

        if ( tokenizer.giveNumberOfTokens() > 0 ) {
            answer = std :: string( tokenizer.giveToken(1) );
            setReadFlag(1);
        } else {
            throw BadFormatInputException(*this, "RecordID", lineNumber);
        }
    }
#endif
#if 0
    void
    XMLInputRecord :: giveField(int &answer, InputFieldType id)
    {
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            auto ptr = scanInteger(tokenizer.giveToken(indx + 1), answer);
            if ( ptr == nullptr || *ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);
            setReadFlag(indx + 1);
        } else {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }

    void
    XMLInputRecord :: giveField(double &answer, InputFieldType id)
    {
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            auto ptr = scanDouble(tokenizer.giveToken(indx + 1), answer);
            if ( ptr == nullptr || *ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);
            setReadFlag(indx + 1);
        } else {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }

    void
    XMLInputRecord :: giveField(bool &answer, InputFieldType id)
    {
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            int val;
            auto ptr = scanInteger(tokenizer.giveToken(indx + 1), val);
            if ( ptr == nullptr || *ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);
            setReadFlag(indx + 1);
            answer = val != 0;
        } else {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }

    void
    XMLInputRecord :: giveField(std :: string &answer, InputFieldType id)
    {
        int indx = 0;
        if ( id ) {
            if ( ( indx = this->giveKeywordIndx(id) ) == 0 ) {
                throw MissingKeywordInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);
            indx++;
        } else {
            indx = 1;
        }

        const char *_token = tokenizer.giveToken(indx);
        if ( _token ) {
            answer = std :: string( tokenizer.giveToken(indx) );
            setReadFlag(indx);
        } else {
            answer = "";
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }

    void
    XMLInputRecord :: giveField(IntArray &answer, InputFieldType id)
    {
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            int size;
            setReadFlag(indx);
            auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
            if ( ptr == nullptr || *ptr != 0) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            answer.resize(size);
            setReadFlag(indx);

            for ( int i = 1; i <= size; i++ ) {
                int value;
                ptr = scanInteger(tokenizer.giveToken(indx + i), value);
                if ( ptr == nullptr || *ptr != 0 ) {
                    throw BadFormatInputException(*this, id, lineNumber);
                }

                answer.at(i) = value;
                setReadFlag(indx + i);
            }

        } else {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }

    void
    XMLInputRecord :: giveField(FloatArray &answer, InputFieldType id)
    {
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            int size;
            setReadFlag(indx);
            auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
            if ( ptr == nullptr || *ptr != 0 ) {
                throw BadFormatInputException(*this, id, lineNumber);
            }

            answer.resize(size);
            setReadFlag(indx);

            for ( int i = 1; i <= size; i++ ) {
                double value;
                auto ptr = scanDouble(tokenizer.giveToken(indx + i), value);
                if ( ptr == nullptr || *ptr != 0 ) {
                    throw BadFormatInputException(*this, id, lineNumber);
                }

                answer.at(i) = value;
                setReadFlag(indx + i);
            }

        } else {
            throw MissingKeywordInputException(*this, id, lineNumber);
        }
    }

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
    XMLInputRecord :: giveField(std :: list< Range > &list, InputFieldType id)
    {
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            int li, hi;
            setReadFlag(indx);
            const char *rec = tokenizer.giveToken(++indx);
            if ( * rec != '{' ) {
                OOFEM_WARNING("missing left '{'");
                list.clear();
                throw BadFormatInputException(*this, id, lineNumber);
            }

            setReadFlag(indx);
            rec++;
            // read ranges
            while ( readRange(& rec, li, hi) ) {
                Range range(li, hi);
                list.push_back(range);
            }

            // skip whitespaces after last range
            while ( isspace(* rec) ) {
                rec++;
            }

            // test for enclosing bracket
            if ( * rec != '}' ) {
                OOFEM_WARNING("missing end '}'");
                list.clear();
                throw BadFormatInputException(*this, id, lineNumber);
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

    bool
    XMLInputRecord :: hasField(InputFieldType id)
    {
        //returns nonzero if id is present in source
        int indx = this->giveKeywordIndx(id);
        if ( indx ) {
            setReadFlag(indx);
        }

        return ( indx > 0 ) ? true : false;
    }

    void
    XMLInputRecord :: printYourself()
    {
        printf( "%s", this->record.c_str() );
    }

    const char *
    XMLInputRecord :: scanInteger(const char *source, int &value)
    {
        //
        // reads integer value from source, returns pointer to char after this number
        //
        char *endptr;

        if ( source == nullptr ) {
            value = 0;
            return nullptr;
        }

        value = strtol(source, & endptr, 10);
        return endptr;
    }

    const char *
    XMLInputRecord :: scanDouble(const char *source, double &value)
    {
        //
        // reads double value from source, returns pointer to char after this number
        //
        char *endptr;

        if ( source == nullptr ) {
            value = 0;
            return nullptr;
        }

        value = strtod(source, & endptr);
        return endptr;
    }

    int
    XMLInputRecord :: giveKeywordIndx(const char *kwd)
    {
        int ntokens = tokenizer.giveNumberOfTokens();
        for ( int i = 1; i <= ntokens; i++ ) {
            if ( strcmp( kwd, tokenizer.giveToken(i) ) == 0 ) {
                return i;
            }
        }

        return 0;
    }

    void
    XMLInputRecord :: finish(bool wrn)
    {
        if ( !wrn ) {
            return;
        }

        std :: ostringstream buff;
        bool pf = true, wf = false;
        int ntokens = tokenizer.giveNumberOfTokens();
        for ( int i = 0; i < ntokens; i++ ) {
            //fprintf (stderr, "[%s] ", tokenizer.giveToken(i+1));
            if ( !readFlag [ i ] ) {
                if ( pf ) {
                    buff << "Unread token(s) detected in the following record\n\"";
                    for ( int j = 0; j < 40; j++ ) {
                        if ( this->record [ j ] == '\n' || this->record [ j ] == '\0' ) {
                            break;
                        } else {
                            buff << this->record [ j ];
                        }
                    }
                    if ( this->record.size() > 41 ) {
                        buff << "...";
                    }
                    buff << "\":\n";

                    pf = false;
                    wf = true;
                }

                buff << "[" << tokenizer.giveToken(i + 1) << "]";
            }
        }

        if ( wf ) {
            OOFEM_WARNING( buff.str().c_str() );
        }
    }

    int
    XMLInputRecord :: readRange(const char **helpSource, int &li, int &hi)
    {
        char *endptr;
        // skip whitespaces
        while ( isspace(* * helpSource) ) {
            ( * helpSource )++;
        }

        // test if character is digit
        if ( isdigit(* * helpSource) ) {
            // digit character - read one value range
            li = hi = strtol(* helpSource, & endptr, 10);
            * helpSource = endptr;
            return 1;
        } else if ( * * helpSource == '(' ) {
            // range left parenthesis found
            ( * helpSource )++;
            // read lower index
            li = strtol(* helpSource, & endptr, 10);
            * helpSource = endptr;
            // test whitespaces
            if ( * * helpSource != ' ' && * * helpSource != '\t' ) {
                OOFEM_WARNING("unexpected token while reading range value");
                return 0;
            }

            // read end index
            hi = strtol(* helpSource, & endptr, 10);
            * helpSource = endptr;
            // skip whitespaces
            while ( isspace(* * helpSource) ) {
                ( * helpSource )++;
            }

            // test for enclosing bracket
            if ( * * helpSource == ')' ) {
                ( * helpSource )++;
                return 1;
            } else {
                OOFEM_WARNING("end ')' missing while parsing range value");
                return 0;
            }
        }

        return 0;
    }

    int
    XMLInputRecord :: readMatrix(const char *helpSource, int r, int c, FloatMatrix &ans)
    {
        const char *endptr = helpSource;

        if ( helpSource == NULL ) {
            ans.clear();
            return 0;
        }

        ans.resize(r, c);
        // skip whitespaces
        while ( isspace(* endptr) ) {
            ( endptr )++;
        }

        if ( * endptr == '{' ) {
            // range left parenthesis found
            ( endptr )++;
            // read row by row separated by semicolon
            for ( int i = 1; i <= r; i++ ) {
                for ( int j = 1; j <= c; j++ ) {
                    endptr = scanDouble( endptr, ans.at(i, j) );
                }

                if ( i < r ) {
                    // skip whitespaces
                    while ( isspace(* endptr) ) {
                        ( endptr )++;
                    }

                    // test for row terminating semicolon
                    if ( * endptr == ';' ) {
                        ( endptr )++;
                    } else {
                        OOFEM_WARNING("missing row terminating semicolon");
                        return 0;
                    }
                }
            }

            // skip whitespaces
            while ( isspace(* endptr) ) {
                ( endptr )++;
            }

            // test for enclosing bracket
            if ( * endptr == '}' ) {
                return 1;
            } else {
                OOFEM_WARNING("end '}' missing while parsing matrix value");
                return 0;
            }
        } else {
            return 0;
        }

    }
#endif
} // end namespace oofem
