/*

Copied from https://github.com/therocode/meta_enum, MIT-license

Adapted for oofem:

- Central registry for enumerations at oofem::MetaEnumMeta (so that one can look up metadata via oofem::MetaEnumMeta<MyEnum>::..., as done in InputRecord handling)
- This however restricts the enums to be declared in oofem:: namespace only (not inside classes).
- To make an enum inside a class, declare it outside and then typedef it inside the class as well.

Better enumerations (with reflection) are coming with c++26.
*/

#pragma once
#include <array>
#include <string_view>
#include <optional>

namespace oofem {
    template<typename T> struct MetaEnumMeta{ };
}

template <typename EnumType>
struct MetaEnumMember
{
    EnumType value = {};
    std::string_view name;
    std::string_view string;
    size_t index = {};
};

template <typename EnumType, typename UnderlyingTypeIn, size_t size>
struct MetaEnum
{
    using UnderlyingType = UnderlyingTypeIn;
    std::string_view string;
    std::array<MetaEnumMember<EnumType>, size> members = {};
};

namespace meta_enum_internal
{


constexpr bool isNested (size_t brackets, bool quote)
{
    return brackets != 0 || quote;
}

constexpr size_t nextEnumCommaOrEnd(size_t start, std::string_view enumString)
{
    size_t brackets = 0; //()[]{}
    bool quote = false; //""
    char lastChar = '\0';
    char nextChar = '\0';

    auto feedCounters = [&brackets, &quote, &lastChar, &nextChar] (char c)
    {
        if(quote)
        {
            if(lastChar != '\\' && c == '"') //ignore " if they are backslashed
                quote = false;
            return;
        }

        switch(c)
        {
            case '"':
                if(lastChar != '\\') //ignore " if they are backslashed
                    quote = true;
                break;
            case '(':
            case '<':
                if(lastChar == '<' || nextChar == '<')
                    break;
                [[fallthrough]];
            case '{':
                ++brackets;
                break;
            case ')':
            case '>':
                if(lastChar == '>' || nextChar == '>')
                    break;
                [[fallthrough]];
            case '}':
                --brackets;
                break;
            default:
                break;
        }
    };

    size_t current = start;
    for(; current < enumString.size() && (isNested(brackets, quote) || (enumString[current] != ',')); ++current)
    {
        feedCounters(enumString[current]);
        lastChar = enumString[current];
        nextChar = current + 2 < enumString.size() ? enumString[current + 2] : '\0';
    }

    return current;
}

constexpr bool isAllowedIdentifierChar(char c)
{
    return (c >= 'a' && c <= 'z') ||
           (c >= 'A' && c <= 'Z') ||
           (c >= '0' && c <= '9') ||
           c == '_';
}

constexpr std::string_view parseEnumMemberName(std::string_view memberString)
{
    size_t nameStart = 0;
    while(!isAllowedIdentifierChar(memberString[nameStart]))
    {
        ++nameStart;
    }

    size_t nameSize = 0;

    while((nameStart + nameSize) < memberString.size() && isAllowedIdentifierChar(memberString[nameStart + nameSize]))
    {
        ++nameSize;
    }

    return std::string_view(memberString.data() + nameStart, nameSize);
}

template <typename EnumType, typename UnderlyingType, size_t size>
constexpr MetaEnum<EnumType, UnderlyingType, size> parseMetaEnum(std::string_view in, const std::array<EnumType, size>& values)
{
    MetaEnum<EnumType, UnderlyingType, size> result;
    result.string = in;

    std::array<std::string_view, size> memberStrings;
    size_t amountFilled = 0;

    size_t currentStringStart = 0;

    while(amountFilled < size)
    {
        size_t currentStringEnd = nextEnumCommaOrEnd(currentStringStart + 1, in);
        size_t currentStringSize = currentStringEnd - currentStringStart;

        if(currentStringStart != 0)
        {
            ++currentStringStart;
            --currentStringSize;
        }

        memberStrings[amountFilled] = std::string_view(in.data() + currentStringStart, currentStringSize);
        ++amountFilled;
        currentStringStart = currentStringEnd;
    }

    for(size_t i = 0; i < memberStrings.size(); ++i)
    {
        result.members[i].name = parseEnumMemberName(memberStrings[i]);
        result.members[i].string = memberStrings[i];
        result.members[i].value = values[i];
        result.members[i].index = i;
    }

    return result;
}

template <typename EnumUnderlyingType>
struct IntWrapper
{
    constexpr IntWrapper(): value(0), empty(true)
    {
    }
    constexpr IntWrapper(EnumUnderlyingType in): value(in), empty(false)
    {
    }
    constexpr IntWrapper operator=(EnumUnderlyingType in)
    {
        value = in;
        empty = false;
        return *this;
    }
    EnumUnderlyingType value;
    bool empty;
};

template <typename EnumType, typename EnumUnderlyingType, size_t size>
constexpr std::array<EnumType, size> resolveEnumValuesArray(const std::initializer_list<IntWrapper<EnumUnderlyingType>>& in)
{
    std::array<EnumType, size> result{};

    EnumUnderlyingType nextValue = 0;
    for(size_t i = 0; i < size; ++i)
    {
        auto wrapper = *(in.begin() + i);
        EnumUnderlyingType newValue = wrapper.empty ? nextValue : wrapper.value;
        nextValue = newValue + 1;
        result[i] = static_cast<EnumType>(newValue);
    }

    return result;
}
}

#define meta_enum(Type, UnderlyingType, ...)\
    enum Type : UnderlyingType { __VA_ARGS__};\
    constexpr static auto Type##_internal_size = []  () constexpr\
    {\
        using IntWrapperType = meta_enum_internal::IntWrapper<UnderlyingType>;\
        IntWrapperType __VA_ARGS__;\
        return std::initializer_list<IntWrapperType>{__VA_ARGS__}.size();\
    };\
    constexpr static auto Type##_meta = meta_enum_internal::parseMetaEnum<Type, UnderlyingType, Type##_internal_size()>(#__VA_ARGS__, []() {\
        using IntWrapperType = meta_enum_internal::IntWrapper<UnderlyingType>;\
        IntWrapperType __VA_ARGS__;\
        return meta_enum_internal::resolveEnumValuesArray<Type, UnderlyingType, Type##_internal_size()>({__VA_ARGS__});\
    }());\
    constexpr static auto Type##_value_to_string = [](Type e)\
    {\
        for(const auto& member : Type##_meta.members)\
        {\
            if(member.value == e)\
                return member.name;\
        }\
        return std::string_view("__INVALID_ENUM_VAL__");\
    \
    };\
    constexpr static auto Type##_meta_from_name = [](std::string_view s) -> std::optional<MetaEnumMember<Type>>\
    {\
        for(const auto& member : Type##_meta.members)\
        {\
            if(member.name == s)\
                return member;\
        }\
        return std::nullopt;\
    \
    };\
    constexpr static auto Type##_meta_from_value = [] (Type v) -> std::optional<MetaEnumMember<Type>>\
    {\
        for(const auto& member : Type##_meta.members)\
        {\
            if(member.value == v)\
                return member;\
        }\
        return std::nullopt;\
    \
    };\
    constexpr static auto Type##_meta_from_index = [] (size_t i)\
    {\
        std::optional<MetaEnumMember<Type>> result;\
        if(i < Type##_meta.members.size())\
            result = Type##_meta.members[i];\
        return result;\
    \
    };\
    template<> struct MetaEnumMeta<Type>{\
        constexpr static auto members = Type##_meta.members; \
        constexpr static auto value_to_string = Type##_value_to_string; \
        constexpr static auto meta_from_name = Type##_meta_from_name;\
        constexpr static auto meta_from_value = Type##_meta_from_value;\
        constexpr static auto meta_from_index = Type##_meta_from_index;\
        constexpr static const char* enum_name = #Type; \
    };\


#define meta_enum_class(Type, UnderlyingType, ...)\
    enum class Type : UnderlyingType { __VA_ARGS__};\
    constexpr static auto Type##_internal_size = [] () constexpr\
    {\
        using IntWrapperType = meta_enum_internal::IntWrapper<UnderlyingType>;\
        IntWrapperType __VA_ARGS__;\
        return std::initializer_list<IntWrapperType>{__VA_ARGS__}.size();\
    };\
    constexpr static auto Type##_meta = meta_enum_internal::parseMetaEnum<Type, UnderlyingType, Type##_internal_size()>(#__VA_ARGS__, []() {\
        using IntWrapperType = meta_enum_internal::IntWrapper<UnderlyingType>;\
        IntWrapperType __VA_ARGS__;\
        return meta_enum_internal::resolveEnumValuesArray<Type, UnderlyingType, Type##_internal_size()>({__VA_ARGS__});\
    }());\
    constexpr static auto Type##_value_to_string = [](Type e)\
    {\
        for(const auto& member : Type##_meta.members)\
        {\
            if(member.value == e)\
                return member.name;\
        }\
        return std::string_view("__INVALID_ENUM_VAL__");\
    \
    };\
    constexpr static auto Type##_meta_from_name = [](std::string_view s) -> std::optional<MetaEnumMember<Type>>\
    {\
        for(const auto& member : Type##_meta.members)\
        {\
            if(member.name == s)\
                return member;\
        }\
        return std::nullopt;\
    \
    };\
    constexpr static auto Type##_meta_from_value = [] (Type v) -> std::optional<MetaEnumMember<Type>>\
    {\
        for(const auto& member : Type##_meta.members)\
        {\
            if(member.value == v)\
                return member;\
        }\
        return std::nullopt;\
    \
    };\
    constexpr static auto Type##_meta_from_index = [] (size_t i)\
    {\
        std::optional<MetaEnumMember<Type>> result;\
        if(i < Type##_meta.members.size())\
            result = Type##_meta.members[i];\
        return result;\
    \
    };\
    template<> struct MetaEnumMeta<Type>{\
        constexpr static auto members = Type##_meta.members; \
        constexpr static auto value_to_string = Type##_value_to_string; \
        constexpr static auto meta_from_name = Type##_meta_from_name;\
        constexpr static auto meta_from_value = Type##_meta_from_value;\
        constexpr static auto meta_from_index = Type##_meta_from_index;\
        constexpr static const char* enum_name = #Type; \
    };\

