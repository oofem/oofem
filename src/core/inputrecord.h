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

#ifndef inputrecord_h
#define inputrecord_h

#include <vector>
#include <list>
#include <string>
#include <memory>
#include <exception>
#include <regex>

#include "logger.h" // for missing __func__ in MSC
#include "oofemenv.h"
#include "enum.h"

namespace oofem {
class IntArray;
class FloatArray;
class FloatMatrix;
class Dictionary;
class Range;
class ScalarFunction;
class DataReader;

/// Identifier of fields in input records.
typedef const char *InputFieldType;

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the compulsory
 * field identified by __kwd and stores the  result into __value parameter.
 * Includes also the error reporting.
 */
#define IR_GIVE_FIELD(__ir, __value, __id) (__ir).giveField(__value, __id);

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the optional
 * field identified by __kwd and stores the  result into __value parameter.
 * Includes also the error reporting.
 */
#define IR_GIVE_OPTIONAL_FIELD(__ir, __value, __id) (__ir).giveOptionalField(__value, __id);

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the compulsory record keyword (__kwd)
 * and its number (__value param). Includes also the error reporting.
 */
#define IR_GIVE_RECORD_KEYWORD_FIELD(__ir, __name, __value) \
    (__ir).giveRecordKeywordField(__name, __value);


// #define _INPUTRECORD_OPTIONAL_OLD

/**
 * Class representing the general Input Record. The input record consists of several fields.
 * Provides several requesting functions for reading field values. The derived classes of
 * Input record can represent database records or text file records, allowing the transparent
 * input operations.
 * The input record after init phase should "contain" all relevant data, so the input record should
 * resolve all dependencies. This allows to create a copy of input record instance for later use
 * without the need to re-open input files (used for metasteps).
 */
class OOFEM_EXPORT InputRecord: public std::enable_shared_from_this<InputRecord>
{
    DataReader* reader = nullptr;
public:

    InputRecord() {}
    InputRecord(DataReader* reader_);
    /// Destructor
    virtual ~InputRecord() = default;

    /** Creates a newly allocated copy of the receiver */
    virtual std::shared_ptr<InputRecord> clone() const = 0;
    std::shared_ptr<InputRecord> ptr();

    /// Returns string representation of record in OOFEMs text format.
    virtual std :: string giveRecordAsString() const = 0;
    virtual std :: string giveRecordInTXTFormat() const = 0;
    virtual std :: string giveLocation() const = 0;

    /**@name Compulsory field extraction methods
     * Reads the field value identified by keyword
     * @param answer contains result
     * @param id field keyword
     */
    //@{
    /// Reads the record id field  (type of record) and its corresponding number.
    virtual void giveRecordKeywordField(std :: string &answer, int &value) = 0;
    /// Reads the record id field  (type of record).
    virtual void giveRecordKeywordField(std :: string &answer) = 0;
    /// Reads the integer field value.
    virtual void giveField(int &answer, InputFieldType id) = 0;
    /// Reads the double field value.
    virtual void giveField(double &answer, InputFieldType id) = 0;
    /// Reads the bool field value.
    virtual void giveField(bool &answer, InputFieldType id) = 0;
    /// Reads the string field value.
    virtual void giveField(std :: string &answer, InputFieldType id) = 0;
    /// Reads the FloatArray field value.
    virtual void giveField(FloatArray &answer, InputFieldType id) = 0;
    /// Reads the IntArray field value.
    virtual void giveField(IntArray &answer, InputFieldType id) = 0;
    /// Reads the FloatMatrix field value.
    virtual void giveField(FloatMatrix &answer, InputFieldType id) = 0;
    /// Reads the vector of strings.
    virtual void giveField(std :: vector< std :: string > &answer, InputFieldType id) = 0;
    /// Reads the Dictionary field value.
    virtual void giveField(Dictionary &answer, InputFieldType id) = 0;
    /// Reads the std::list<Range> field value.
    virtual void giveField(std :: list< Range > &answer, InputFieldType id) = 0;
    /// Reads the ScalarFunction field value.
    virtual void giveField(ScalarFunction &function, InputFieldType id) = 0;

    static std::string error_msg_with_hints(const std::string& val, const std::map<int,std::vector<std::string>>& v2nn);
    static int giveLevenshteinDist(const std::string& word1, const std::string& word2);

    /// Reads enumeration (must be defined via enum-impl.hpp) directly
    template<typename AnEnum>
    void giveField(AnEnum& answer, InputFieldType id){
        typedef EnumTraits<AnEnum> Traits;
        std::string s;
        giveField(s,id);
        #ifdef _USE_TRACE_FIELDS
            if(InputRecord::TraceFields::active){
                traceEnum(Traits::enum_name,Traits::all_values_to_names());
                traceField(id,(std::string("enum:")+Traits::enum_name).c_str());
            }
        #endif
        if(std::regex_match(s,std::regex("\\s*[0-9]+\\s*"))){
            int val=std::atoi(s.c_str());
            auto v=Traits::value(val);
            if(!v) OOFEM_ERROR("%s: %s (enum %s): invalid index '%d'%s",giveLocation().c_str(),id,Traits::enum_name,val,error_msg_with_hints("",Traits::all_values_to_names()).c_str());
            answer=v.value();
        } else {
            auto v=Traits::value(s.c_str());
            if(!v){ OOFEM_ERROR("%s: %s (enum %s): unrecognized name '%s'%s",giveLocation().c_str(),id,Traits::enum_name,s.c_str(),error_msg_with_hints(s,Traits::all_values_to_names()).c_str()); }
            answer=v.value();
        }
    }
    #ifdef _USE_TRACE_FIELDS
        // field access tracing variables, set at startup from main()
        struct TraceFields {
            static bool active;
            static std::ofstream out;
            static void write(const std::string& s);
        };
        static void traceEnum(const std::string& name, const std::map<int,std::vector<std::string>>& val2names);
        virtual void traceField(InputFieldType id, const char* type) {};
    #else
        void traceField(InputFieldType id, const char* type) const { };
    #endif
    /**@name Optional field extraction methods
     * Reads the field value identified by keyword
     * @param answer contains result
     * @param id field keyword
     */
    template<typename T>
    void giveOptionalField(T& answer, InputFieldType id){ if(hasField(id)) giveField(answer,id); }
    //@}

    /**@name Child reader methods
     */
    //@{
    // return count of nested group; if not optional, the empty group is required
    virtual int giveGroupCount(InputFieldType id, const std::string& name, bool optional) = 0;
    // return whether a single child of given type exists
    virtual bool hasChild(InputFieldType id, const std::string& name, bool optional) = 0;
    //@}
#ifdef _INPUTRECORD_OPTIONAL_OLD
    //@{
    /// Reads the integer field value.
    void giveOptionalField(int &answer, InputFieldType id);
    /// Reads the double field value.
    void giveOptionalField(double &answer, InputFieldType id);
    /// Reads the bool field value.
    void giveOptionalField(bool &answer, InputFieldType id);
    /// Reads the string field value.
    void giveOptionalField(std :: string &answer, InputFieldType id);
    /// Reads the FloatArray field value.
    void giveOptionalField(FloatArray &answer, InputFieldType id);
    /// Reads the IntArray field value.
    void giveOptionalField(IntArray &answer, InputFieldType id);
    /// Reads the FloatMatrix field value.
    void giveOptionalField(FloatMatrix &answer, InputFieldType id);
    /// Reads the vector of strings.
    void giveOptionalField(std :: vector< std :: string > &answer, InputFieldType id);
    /// Reads the Dictionary field value.
    void giveOptionalField(Dictionary &answer, InputFieldType id);
    /// Reads the std::list<Range> field value.
    void giveOptionalField(std :: list< Range > &answer, InputFieldType id);
    /// Reads the ScalarFunction field value.
    void giveOptionalField(ScalarFunction &function, InputFieldType id);
    //@}
#endif
    /// Returns true if record contains field identified by idString keyword.
    virtual bool hasField(InputFieldType id) = 0;

    /// Print input record.
    virtual void printYourself() = 0;

    /// Terminates the current record session and if the flag is true, warning is printed for unscanned tokens.
    virtual void finish(bool wrn = true) = 0;

    DataReader* giveReader() const;
};


class InputException : public std::exception
{
public:
    std::string record;
    std::string keyword;
    int number;
    InputException(const InputRecord &ir, std::string keyword, int number);
};


class MissingKeywordInputException : public InputException
{
protected:
    std::string msg;

public:
    MissingKeywordInputException(const InputRecord &ir, std::string keyword, int number);
    const char* what() const noexcept override;
};


class BadFormatInputException : public InputException
{
protected:
    std::string msg;

public:
    BadFormatInputException(const InputRecord &ir, std::string keyword, int number);
    const char* what() const noexcept override;
};


class ValueInputException : public InputException
{
protected:
    std::string msg;

public:
    ValueInputException(const InputRecord &ir, std::string keyword, const std::string &reason);
    const char* what() const noexcept override;
};

class ComponentInputException : public std::exception {
public:
    enum ComponentType 
    {
        ctElement,
        ctDofManager
    };

protected:
    std::string msg;
public:
    ComponentInputException(const std::string keyword, ComponentType ct, int number, const std::string &reason);
    ComponentInputException(ComponentType ct, int number, const std::string &reason);

    const char* what() const noexcept override;
};

} // end namespace oofem
#endif // inputrecord_h
