/*
This header is to be included in-place to define an enum.

The following macros must be defined prior to inclusion (will be #undef'd)

- ENUM_TYPE: name of the enum
- ENUM_DEF: sequence of ENUM_ITEM(name) or ENUM_ITEM_WITH_VALUE(name,value)

The following may be defined prior to includion (will be #undef'd)

- ENUM_PREFIX: if defined, name to value conversion will allow dropping prefix from the name (e.g. saying "Total" instead of "VM_Total", the prefix being "VM_")
- ENUM_CLASS: if defined, define the enum as "enum class"; otherwise, a plain enum

*/

#ifndef enum_h
   #error enum.h MUST be included prior to including this header
#endif

#if !defined(ENUM_TYPE) || !defined(ENUM_DEF)
   #error both ENUM_TYPE and ENUM_DEF macros must be defined before including this header
#endif


// leaked from enumitem.hpp perhaps, so undefine to avoid warnings and then re-define again
#ifdef ENUM_ITEM
   #undef ENUM_ITEM
   #undef ENUM_ITEM_WITH_VALUE
   #define _MUST_RESTORE_ENUM_ITEM
#endif

// helper macros
#define _str(s) #s
#define _xstr(s) _str(s)

// define the enum itself
#define ENUM_ITEM(i) i,
#define ENUM_ITEM_WITH_VALUE(i,v) i=v,
enum
   #ifdef ENUM_CLASS
      class
   #endif
   ENUM_TYPE { ENUM_DEF };
#undef ENUM_ITEM_WITH_VALUE
#undef ENUM_ITEM

#ifndef ENUM_PREFIX
   #define ENUM_PREFIX ""
#endif

/* minimal class holding metadata about the enum: values, names, opt_prefix */
template<> struct EnumData<ENUM_TYPE>{
   typedef ENUM_TYPE Enum;
   static constexpr const char* enum_name = _xstr(ENUM_TYPE);
   static constexpr const char opt_prefix[]=ENUM_PREFIX;
   static constexpr size_t opt_prefix_len=sizeof(opt_prefix)-1;

   struct EnumItem { Enum value; const char* name; };
   #define ENUM_ITEM_WITH_VALUE(i,v) {ENUM_TYPE::i,#i},
   #define ENUM_ITEM(i) {ENUM_TYPE::i,#i},
   static constexpr std::initializer_list<EnumItem> value_to_name={ ENUM_DEF };
   #undef ENUM_ITEM_WITH_VALUE
   #undef ENUM_ITEM
};

// conversion function for compatibility
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
