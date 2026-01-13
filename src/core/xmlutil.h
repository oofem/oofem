#pragma once
#include<charconv>
#include<cstring>
#include<pugixml.hpp>
#include<functional>
#include"error.h"
#include"range.h"
#include<regex>
#include<cassert>

namespace oofem{
   namespace xmlutil{
      template<typename T>
      T string_to(const std::string& s, std::function<std::string()> where){
         T val;
         const char* last=s.data()+s.size();
         auto [p,e]=std::from_chars(s.data(),last,val);
         if(p!=last) OOFEM_ERROR("%s: error parsing '%s' as typeid '%s' (leftover chars)",where().c_str(),s.c_str(),typeid(T).name());
         if(e==std::errc()) return val;
         OOFEM_ERROR("%s: error parsing '%s' as typeid '%s' (std::from_chars error).",where().c_str(),s.c_str(),typeid(T).name());
      }
      template<> Range string_to(const std::string& s, std::function<std::string()> where);
      template<> bool string_to(const std::string& s, std::function<std::string()> where);
      template<> std::string string_to(const std::string& s, std::function<std::string()> where);

      /* helper class (derived from pugi::xml_document) which additionally holds mapping from offset to line/column, and source filename.
       * plus defines walk_get<T> utility method for quick retrieval of values from XML
      */
      struct XmlDoc: public pugi::xml_document {
         std::string filename;
         std::vector<size_t> newlines;
         std::tuple<size_t,size_t> offset2lc(size_t offset) const;
         std::string offset2loc(size_t offset) const;
         XmlDoc(){}
         XmlDoc(const std::string& xml);

         template<typename T>
         T walk_get(const std::vector<std::string>& path, const char* attr=""){
            pugi::xml_node curr=*this;
            for(const auto& p: path){
               pugi::xml_node prev=curr;
               curr=prev.child(p);
               if(!curr) OOFEM_ERROR("%s: no child element <%s>",offset2loc(prev.offset_debug()).c_str(),p.c_str());
            }
            if(attr[0]=='\0') return string_to<T>(curr.text().get(),[this,&curr](){ return (this->offset2loc(curr.offset_debug())+": "+curr.name()); });
            pugi::xml_attribute a=curr.attribute(attr);
            if(!a) OOFEM_ERROR("%s: no attribute '%s'",offset2loc(curr.offset_debug()).c_str(),a.name());
            return string_to<T>(a.value(),[this,&curr,&attr](){ return (this->offset2loc(curr.offset_debug())+": attribute '"+attr+"'"); });
         }
      };
   }
}
