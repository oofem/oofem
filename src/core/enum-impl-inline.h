template<typename T> struct EnumTraits;

#if !defined(ENUM_TYPE) || !defined(ENUM_DEF)
   #error both ENUM_TYPE and ENUM_DEF macros must be defined before including this header
#endif

// leaked from enumitem.hpp perhaps
#ifdef ENUM_ITEM
   #undef ENUM_ITEM
   #undef ENUM_ITEM_WITH_VALUE
   #define _MUST_RESTORE_ENUM_ITEM
#endif

#define _str(s) #s
#define _xstr(s) _str(s)

// define the enum itself
#define ENUM_ITEM(i) i,
#define ENUM_ITEM_WITH_VALUE(i,v) i=v,
enum
   #ifdef ENUM_CLASS
      class
   #endif
   ENUM_TYPE{ ENUM_DEF };
#undef ENUM_ITEM_WITH_VALUE
#undef ENUM_ITEM

#ifndef ENUM_PREFIX
   #define ENUM_PREFIX ""
#endif

// conversions
template<> struct EnumTraits<ENUM_TYPE>{
   static constexpr const char* enum_name = _xstr(ENUM_TYPE);

   #define ENUM_ITEM_WITH_VALUE(i,v) {ENUM_TYPE::i,#i},
   #define ENUM_ITEM(i) {ENUM_TYPE::i,#i},
   static constexpr std::initializer_list<std::pair<ENUM_TYPE,const char*>> value_to_name={ ENUM_DEF };
   #undef ENUM_ITEM_WITH_VALUE
   #undef ENUM_ITEM

   static constexpr const char opt_prefix[]=ENUM_PREFIX;
   static constexpr size_t opt_prefix_len=sizeof(opt_prefix)-1;
   static const bool _starts_with_prefix(const char* s){ return opt_prefix_len>0 && std::strlen(s)>opt_prefix_len && std::strncmp(s,opt_prefix,opt_prefix_len)==0; }

   static std::optional<const char*> name(ENUM_TYPE v){ for(const auto& vn: EnumTraits<ENUM_TYPE>::value_to_name) if(vn.first==v) return vn.second; return {}; };
   static std::optional<ENUM_TYPE> value(int i){ for(const auto& vn: EnumTraits<ENUM_TYPE>::value_to_name) if((int)vn.first==i) return vn.first; return {}; };
   static std::optional<ENUM_TYPE> value(const char* n){
      for(const auto& vn: EnumTraits<ENUM_TYPE>::value_to_name){
         if(std::strcmp(vn.second,n)==0) return vn.first;
         if(_starts_with_prefix(vn.second) && std::strcmp(vn.second+opt_prefix_len,n)==0) return vn.first;
      }
      return {};
   };
   static std::vector<std::string> all_names(){
      std::vector<std::string> ret;
      ret.reserve((opt_prefix_len==0?1:2)*value_to_name.size());
      for(const auto& vn: EnumTraits<ENUM_TYPE>::value_to_name){
         ret.push_back(vn.second);
         if(_starts_with_prefix(vn.second)) ret.push_back(vn.second+opt_prefix_len);
      }
      return ret;
   }
};

// compatibility with the current functions
#define __TO_STRING(T) inline const char* __ ## T ## ToString(ENUM_TYPE value){ return EnumTraits<ENUM_TYPE>::name(value).value(); }
#define _TO_STRING(T) __TO_STRING(T)
   _TO_STRING(ENUM_TYPE)
#undef _TO_STRING
#undef __TO_STRING

// don't leak macros out
#undef ENUM_TYPE
#undef ENUM_DEF
#undef ENUM_NAME
#undef ENUM_PREFIX
#ifdef ENUM_CLASS
   #undef ENUM_CLASS
#endif
#undef _str
#undef _xstr

// restore
#ifdef _MUST_RESTORE_ENUM_ITEM
   #define ENUM_ITEM(element) element,
   #define ENUM_ITEM_WITH_VALUE(element, val) element = val,
   #undef _MUST_RESTORE_ENUM_ITEM
#endif
