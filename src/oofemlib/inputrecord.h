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

#ifndef inputrecord_h
#define inputrecord_h

#include <vector>
#include <list>
#include <string>

#include "logger.h" // for missing __func__ in MSC
#include "oofemcfg.h"
#include "irresulttype.h"

namespace oofem {
class IntArray;
class FloatArray;
class FloatMatrix;
class Dictionary;
class Range;
class ScalarFunction;

/// Identifier of fields in input records.
typedef const char *InputFieldType;

/**
 * Macro simplifying the error reporting.
 */
#define IR_IOERR(__keyword, __ir, __result) \
    __ir->report_error(this->giveClassName(), __func__, __keyword, __result, __FILE__, __LINE__);

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the compulsory
 * field identified by __kwd and stores the  result into __value parameter.
 * Includes also the error reporting.
 */
#define IR_GIVE_FIELD(__ir, __value, __id) result = __ir->giveField(__value, __id); \
    if ( result != IRRT_OK ) { IR_IOERR(__id, __ir, result); }

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the optional
 * field identified by __kwd and stores the  result into __value parameter.
 * Includes also the error reporting.
 */
#define IR_GIVE_OPTIONAL_FIELD(__ir, __value, __id) result = __ir->giveOptionalField(__value, __id); \
    if ( result != IRRT_OK ) { IR_IOERR(__id, __ir, result); }

/**
 * Macro facilitating the use of input record reading methods.
 * uses the given input record (__ir parameter) and reads the compulsory record keyword (__kwd)
 * and its number (__value param). Includes also the error reporting.
 */
#define IR_GIVE_RECORD_KEYWORD_FIELD(__ir, __name, __value) \
    result = __ir->giveRecordKeywordField(__name, __value); \
    if ( result != IRRT_OK ) { IR_IOERR("RecordIDField", __ir, result); }



/**
 * Class representing the general Input Record. The input record consists of several fields.
 * Provides several requesting functions for reading field values. The derived classes of
 * Input record can represent database records or text file records, allowing the transparent
 * input operations.
 * The input record after init phase should "contain" all relevant data, so the input record should
 * resolve all dependencies. This allows to create a copy of input record instance for later use
 * without the need to re-open input files (used for metasteps).
 */
class OOFEM_EXPORT InputRecord
{
public:
    /// Constructor. Creates an empty input record.
    InputRecord();
    /// Copy constructor
    InputRecord(const InputRecord &);
    /// Destructor
    virtual ~InputRecord() { }
    /// Assignment operator.
    InputRecord &operator = ( const InputRecord & );

    /** Creates a newly allocated copy of the receiver */
    virtual InputRecord *GiveCopy() = 0;

    /// Returns string representation of record in OOFEMs text format.
    virtual std :: string giveRecordAsString() const = 0;

    /**@name Compulsory field extraction methods
     * Reads the field value identified by keyword
     * @param answer contains result
     * @param idString field keyword
     * @return IRResultType
     */
    //@{
    /// Reads the record id field  (type of record) and its corresponding number.
    virtual IRResultType giveRecordKeywordField(std :: string &answer, int &value) = 0;
    /// Reads the record id field  (type of record).
    virtual IRResultType giveRecordKeywordField(std :: string &answer) = 0;
    /// Reads the integer field value.
    virtual IRResultType giveField(int &answer, InputFieldType id) = 0;
    /// Reads the double field value.
    virtual IRResultType giveField(double &answer, InputFieldType id) = 0;
    /// Reads the bool field value.
    virtual IRResultType giveField(bool &answer, InputFieldType id) = 0;
    /// Reads the string field value.
    virtual IRResultType giveField(std :: string &answer, InputFieldType id) = 0;
    /// Reads the FloatArray field value.
    virtual IRResultType giveField(FloatArray &answer, InputFieldType id) = 0;
    /// Reads the IntArray field value.
    virtual IRResultType giveField(IntArray &answer, InputFieldType id) = 0;
    /// Reads the FloatMatrix field value.
    virtual IRResultType giveField(FloatMatrix &answer, InputFieldType id) = 0;
    /// Reads the vector of strings.
    virtual IRResultType giveField(std :: vector< std :: string > &answer, InputFieldType id) = 0;
    /// Reads the Dictionary field value.
    virtual IRResultType giveField(Dictionary &answer, InputFieldType id) = 0;
    /// Reads the std::list<Range> field value.
    virtual IRResultType giveField(std :: list< Range > &answer, InputFieldType id) = 0;
    /// Reads the ScalarFunction field value.
    virtual IRResultType giveField(ScalarFunction &function, InputFieldType id) = 0;
    //@}

    /**@name Optional field extraction methods
     * Reads the field value identified by keyword
     * @param answer contains result
     * @param idString field keyword
     * @return IRResultType
     */
    //@{
    /// Reads the integer field value.
    IRResultType giveOptionalField(int &answer, InputFieldType id);
    /// Reads the double field value.
    IRResultType giveOptionalField(double &answer, InputFieldType id);
    /// Reads the bool field value.
    IRResultType giveOptionalField(bool &answer, InputFieldType id);
    /// Reads the string field value.
    IRResultType giveOptionalField(std :: string &answer, InputFieldType id);
    /// Reads the FloatArray field value.
    IRResultType giveOptionalField(FloatArray &answer, InputFieldType id);
    /// Reads the IntArray field value.
    IRResultType giveOptionalField(IntArray &answer, InputFieldType id);
    /// Reads the FloatMatrix field value.
    IRResultType giveOptionalField(FloatMatrix &answer, InputFieldType id);
    /// Reads the vector of strings.
    IRResultType giveOptionalField(std :: vector< std :: string > &answer, InputFieldType id);
    /// Reads the Dictionary field value.
    IRResultType giveOptionalField(Dictionary &answer, InputFieldType id);
    /// Reads the std::list<Range> field value.
    IRResultType giveOptionalField(std :: list< Range > &answer, InputFieldType id);
    /// Reads the ScalarFunction field value.
    IRResultType giveOptionalField(ScalarFunction &function, InputFieldType id);
    //@}

    /// Returns true if record contains field identified by idString keyword.
    virtual bool hasField(InputFieldType id) = 0;

    /// Returns error string corresponding to given value of IRResultType type.
    const char *strerror(IRResultType);
    /// Print input record.
    virtual void printYourself() = 0;

    /// Prints the error message.
    virtual void report_error(const char *_class, const char *proc, InputFieldType id,
                              IRResultType result, const char *file, int line) = 0;

    /// Terminates the current record session and if the flag is true, warning is printed for unscanned tokens.
    virtual void finish(bool wrn = true) = 0;
};
} // end namespace oofem
#endif // inputrecord_h
