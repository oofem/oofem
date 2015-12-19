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

#include "dynamicinputrecord.h"
#include "femcmpnn.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dictionary.h"
#include "range.h"
#include "node.h"
#include "element.h"
#include "scalarfunction.h"

#include <sstream>

namespace oofem {
DynamicInputRecord *CreateNodeIR(int i, InputFieldType nodeType, FloatArray coord)
{
    DynamicInputRecord *result = new DynamicInputRecord(nodeType, i);
    result->setField(std :: move(coord), _IFT_Node_coords);
    return result;
}

DynamicInputRecord *CreateElementIR(int i, InputFieldType elementType, IntArray nodes, int cs)
{
    DynamicInputRecord *result = new DynamicInputRecord(elementType, i);
    result->setField(std :: move(nodes), _IFT_Element_nodes);
    if ( cs != 0 ) {
        result->setField(cs, _IFT_Element_crosssect);
    }
    return result;
}


DynamicInputRecord :: DynamicInputRecord(std :: string keyword, int value) : InputRecord(),
    recordKeyword(keyword), ///@todo std::move here
    recordNumber(value),
    emptyRecord(),
    intRecord(),
    doubleRecord(),
    boolRecord(),
    stringRecord(),
    floatArrayRecord(),
    intArrayRecord(),
    matrixRecord(),
    stringListRecord(),
    dictionaryRecord(),
    rangeRecord()
{ }

DynamicInputRecord :: DynamicInputRecord(FEMComponent &femc) : InputRecord(),
    recordKeyword(),
    recordNumber(0),
    emptyRecord(),
    intRecord(),
    doubleRecord(),
    boolRecord(),
    stringRecord(),
    floatArrayRecord(),
    intArrayRecord(),
    matrixRecord(),
    stringListRecord(),
    dictionaryRecord(),
    rangeRecord()
{
    femc.giveInputRecord(* this);
}

DynamicInputRecord :: DynamicInputRecord(const DynamicInputRecord &src) : InputRecord(src),
    recordKeyword(src.recordKeyword),
    recordNumber(src.recordNumber),
    emptyRecord(src.emptyRecord),
    intRecord(src.intRecord),
    doubleRecord(src.doubleRecord),
    boolRecord(src.boolRecord),
    stringRecord(src.stringRecord),
    floatArrayRecord(src.floatArrayRecord),
    intArrayRecord(src.intArrayRecord),
    matrixRecord(src.matrixRecord),
    stringListRecord(src.stringListRecord),
    dictionaryRecord(src.dictionaryRecord),
    rangeRecord(src.rangeRecord)
{ }

DynamicInputRecord :: ~DynamicInputRecord()
{ }

DynamicInputRecord &DynamicInputRecord :: operator = ( const DynamicInputRecord & src )
{
    this->recordKeyword = src.recordKeyword;
    this->recordNumber = src.recordNumber;
    this->emptyRecord = src.emptyRecord;
    this->intRecord = src.intRecord;
    this->doubleRecord = src.doubleRecord;
    this->boolRecord = src.boolRecord;
    this->stringRecord = src.stringRecord;
    this->floatArrayRecord = src.floatArrayRecord;
    this->intArrayRecord = src.intArrayRecord;
    this->matrixRecord = src.matrixRecord;
    this->stringListRecord = src.stringListRecord;
    this->dictionaryRecord = src.dictionaryRecord;
    this->rangeRecord = src.rangeRecord;

    return * this;
}

void DynamicInputRecord :: finish(bool wrn)
{
    ///@todo Implement warning about unread entries
}

IRResultType DynamicInputRecord :: giveRecordKeywordField(std :: string &answer, int &value)
{
    answer = this->recordKeyword;
    value = this->recordNumber;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveRecordKeywordField(std :: string &answer)
{
    answer = this->recordKeyword;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(int &answer, InputFieldType id)
{
    std :: map< std :: string, int > :: iterator it = this->intRecord.find(id);
    if ( it == this->intRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(double &answer, InputFieldType id)
{
    std :: map< std :: string, double > :: iterator it = this->doubleRecord.find(id);
    if ( it == this->doubleRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(bool &answer, InputFieldType id)
{
    std :: map< std :: string, bool > :: iterator it = this->boolRecord.find(id);
    if ( it == this->boolRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(std :: string &answer, InputFieldType id)
{
    std :: map< std :: string, std :: string > :: iterator it = this->stringRecord.find(id);
    if ( it == this->stringRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(FloatArray &answer, InputFieldType id)
{
    std :: map< std :: string, FloatArray > :: iterator it = this->floatArrayRecord.find(id);
    if ( it == this->floatArrayRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(IntArray &answer, InputFieldType id)
{
    std :: map< std :: string, IntArray > :: iterator it = this->intArrayRecord.find(id);
    if ( it == this->intArrayRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(FloatMatrix &answer, InputFieldType id)
{
    std :: map< std :: string, FloatMatrix > :: iterator it = this->matrixRecord.find(id);
    if ( it == this->matrixRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(std :: vector< std :: string > &answer, InputFieldType id)
{
    std :: map< std :: string, std :: vector< std :: string > > :: iterator it = this->stringListRecord.find(id);
    if ( it == this->stringListRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(Dictionary &answer, InputFieldType id)
{
    std :: map< std :: string, Dictionary > :: iterator it = this->dictionaryRecord.find(id);
    if ( it == this->dictionaryRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(std :: list< Range > &answer, InputFieldType id)
{
    std :: map< std :: string, std :: list< Range > > :: iterator it = this->rangeRecord.find(id);
    if ( it == this->rangeRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(ScalarFunction &answer, InputFieldType id)
{
    std :: map< std :: string, ScalarFunction > :: iterator it = this->scalarFunctionRecord.find(id);
    if ( it == this->scalarFunctionRecord.end() ) {
        return IRRT_NOTFOUND;
    }
    answer = it->second;
    return IRRT_OK;
}

bool DynamicInputRecord :: hasField(InputFieldType id)
{
    return this->emptyRecord.find(id) != this->emptyRecord.end() ||
           this->intRecord.find(id) != this->intRecord.end() ||
           this->doubleRecord.find(id) != this->doubleRecord.end() ||
           this->boolRecord.find(id) != this->boolRecord.end() ||
           this->floatArrayRecord.find(id) != this->floatArrayRecord.end() ||
           this->intArrayRecord.find(id) != this->intArrayRecord.end() ||
           this->matrixRecord.find(id) != this->matrixRecord.end() ||
           this->stringListRecord.find(id) != this->stringListRecord.end() ||
           this->dictionaryRecord.find(id) != this->dictionaryRecord.end() ||
           this->rangeRecord.find(id) != this->rangeRecord.end();
}

void
DynamicInputRecord :: printYourself()
{
    //printf( "%s", this->record );
}

// Setters
void DynamicInputRecord :: setRecordKeywordField(std :: string keyword, int value)
{
    this->recordKeyword = std :: move(keyword);
    this->recordNumber = value;
}

void DynamicInputRecord :: setRecordKeywordNumber(int value)
{
    this->recordNumber = value;
}

void DynamicInputRecord :: setField(int item, InputFieldType id)
{
    this->intRecord [ id ] = item;
}

void DynamicInputRecord :: setField(double item, InputFieldType id)
{
    this->doubleRecord [ id ] = item;
}

void DynamicInputRecord :: setField(bool item, InputFieldType id)
{
    this->boolRecord [ id ] = item;
}

void DynamicInputRecord :: setField(std :: string item, InputFieldType id)
{
    this->stringRecord [ id ] = item;
}

void DynamicInputRecord :: setField(FloatArray item, InputFieldType id)
{
    this->floatArrayRecord [id] = std :: move(item);
}

void DynamicInputRecord :: setField(IntArray item, InputFieldType id)
{
    this->intArrayRecord [id] = std :: move(item);
}

void DynamicInputRecord :: setField(FloatMatrix item, InputFieldType id)
{
    this->matrixRecord [id] = std :: move(item);
}

void DynamicInputRecord :: setField(std :: vector< std :: string > item, InputFieldType id)
{
    this->stringListRecord [id] = std :: move(item);
}

void DynamicInputRecord :: setField(const Dictionary &item, InputFieldType id)
{
    this->dictionaryRecord [ id ] = item;
}

void DynamicInputRecord :: setField(const std :: list< Range > &item, InputFieldType id)
{
    this->rangeRecord [ id ] = item;
}

void DynamicInputRecord :: setField(const ScalarFunction &item, InputFieldType id)
{
    this->scalarFunctionRecord [ id ] = item;
}

void DynamicInputRecord :: setField(InputFieldType id)
{
    this->emptyRecord.insert(id);
}

void DynamicInputRecord :: unsetField(InputFieldType id)
{
    this->emptyRecord.erase(id);
    this->intRecord.erase(id);
    this->doubleRecord.erase(id);
    this->boolRecord.erase(id);
    this->stringRecord.erase(id);
    this->floatArrayRecord.erase(id);
    this->intArrayRecord.erase(id);
    this->matrixRecord.erase(id);
    this->stringListRecord.erase(id);
    this->dictionaryRecord.erase(id);
    this->rangeRecord.erase(id);
    this->scalarFunctionRecord.erase(id);
}

void
DynamicInputRecord :: report_error(const char *_class, const char *proc, InputFieldType id,
                                   IRResultType result, const char *file, int line)
{
    oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, NULL, file, line,
                              "Input error: \"%s\", field keyword \"%s\"\n%s::%s",
                              strerror(result), id, _class, proc);
    OOFEM_EXIT(1); ///@todo We should never directly exit when dealing with user input.
}

// Helpful macro since we have so many separate records
#define forRecord(name) \
    for ( const auto &x: name ) { \
        rec << " " << x.first << " " << x.second; \
    }

std :: string DynamicInputRecord :: giveRecordAsString() const
{
    std :: ostringstream rec;
    rec << this->recordKeyword;
    // Some records aren't numbered, here we assume that if no number is set (default 0) then it's not printed.
    // Though, technically, some things *could* be numbered arbitrarily (even negative), this is never actually done in practice.
    if ( this->recordNumber > 0 ) {
        rec << " " << this->recordNumber;
    }

    // Empty fields
    for ( const auto &x: emptyRecord ) {
        rec << " " << x;
    }

    // Standard fields;
    forRecord(intRecord);
    forRecord(doubleRecord);
    forRecord(boolRecord);
    forRecord(stringRecord);
    forRecord(floatArrayRecord);
    forRecord(intArrayRecord);
    forRecord(matrixRecord);
    forRecord(dictionaryRecord);
    forRecord(scalarFunctionRecord);

    // Have to write special code for std::vector and std::list
    for ( const auto &x: stringListRecord ) {
        rec << " " << x.first;
        const std :: vector< std :: string > &list = x.second;
        rec << " " << list.size();
        for ( const auto &y: list ) {
            rec << " " << y;
        }
    }

    for ( const auto &x: rangeRecord ) {
        rec << " " << x.first;
        const std :: list< Range > &list = x.second;
        rec << " " << list.size();
        for ( const auto &y: list ) {
            rec << " " << y;
        }
    }

    return rec.str();
}
};
