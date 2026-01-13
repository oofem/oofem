#include"xmlutil.h"
#include<fstream>
#include<list>
namespace oofem::xmlutil {

   template<>
   Range string_to(const std::string& s, std::function<std::string()> where) {
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
   template<>
   std::string string_to(const std::string& s, std::function<std::string()> where) { return s; }

   std::tuple<size_t,size_t> XmlDoc::offset2lc(size_t offset) const {
      if(newlines.empty()) return std::make_tuple(0,0);
      size_t ix=std::distance(newlines.begin(),std::lower_bound(newlines.begin(),newlines.end(),offset));
      return std::make_tuple(ix+1,offset-(ix==0?0:newlines[ix-1]));
   }
   std::string XmlDoc::offset2loc(size_t offset) const {
      auto [line,col]=offset2lc(offset);
      return filename+":"+std::to_string(line)+":"+std::to_string(col);
   }
   XmlDoc::XmlDoc(const std::string& xml){
      filename=xml;
      std::ifstream i(xml,std::ifstream::in|std::ios::binary);
      if(!i.is_open()) OOFEM_ERROR("Error opening %s: %s",xml.c_str(),std::strerror(errno));
      std::list<size_t> br;
      size_t pos=0;
      while(i.good()){ char c=i.get(); pos++; if(c=='\n') br.push_back(pos); }
      newlines.reserve(br.size());
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Waggressive-loop-optimizations"
         newlines.assign(br.begin(),br.end());
      #pragma GCC diagnostic pop
      // _XML_DEBUG(xml<<": "<< newlines.size()<<" newlines found, about to parse the XML...");
      pugi::xml_parse_result result=load_file(xml.c_str());
      if(!result) OOFEM_ERROR("Error parsing %s: %s",offset2loc(result.offset).c_str(),result.description());
}

};
