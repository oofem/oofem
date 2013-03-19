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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "dynamicinputrecord.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "dictionr.h"
#include "range.h"

namespace oofem {

DynamicInputRecord :: DynamicInputRecord() : InputRecord(),
    recordKeyword(),
    recordNumber(0),
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
}

DynamicInputRecord :: DynamicInputRecord(const DynamicInputRecord &src) :
    recordKeyword(src.recordKeyword),
    recordNumber(src.recordNumber),
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
{
}

DynamicInputRecord :: ~DynamicInputRecord()
{
}

DynamicInputRecord & DynamicInputRecord :: operator=(const DynamicInputRecord &src)
{
    this->recordKeyword = src.recordKeyword;
    this->recordNumber = src.recordNumber;
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

IRResultType DynamicInputRecord :: giveRecordKeywordField(std::string &answer, int &value)
{
    answer = this->recordKeyword;
    value = this->recordNumber;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveRecordKeywordField(std::string &answer)
{
    answer = this->recordKeyword;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(int &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, int>::iterator it = this->intRecord.find(id);
    if (it == this->intRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(double &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, double>::iterator it = this->doubleRecord.find(id);
    if (it == this->doubleRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(bool &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, bool>::iterator it = this->boolRecord.find(id);
    if (it == this->boolRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(std::string &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, std::string>::iterator it = this->stringRecord.find(id);
    if (it == this->stringRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(FloatArray &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, FloatArray>::iterator it = this->floatArrayRecord.find(id);
    if (it == this->floatArrayRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(IntArray &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, IntArray>::iterator it = this->intArrayRecord.find(id);
    if (it == this->intArrayRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(FloatMatrix &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, FloatMatrix>::iterator it = this->matrixRecord.find(id);
    if (it == this->matrixRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(std::vector< std::string > &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, std::vector< std::string> >::iterator it = this->stringListRecord.find(id);
    if (it == this->stringListRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(Dictionary &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, Dictionary>::iterator it = this->dictionaryRecord.find(id);
    if (it == this->dictionaryRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(std::list< Range > &answer, InputFieldType fieldID, const char *id)
{
    std::map< std::string, std::list< Range > >::iterator it = this->rangeRecord.find(id);
    if (it == this->rangeRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

bool DynamicInputRecord :: hasField(InputFieldType fieldID, const char *id)
{
    return this->intRecord.find(id) != this->intRecord.end() ||
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
void DynamicInputRecord :: setRecordKeywordField(const std::string &keyword, int value)
{
    this->recordKeyword = keyword;
    this->recordNumber = value;
}

void DynamicInputRecord :: setField(int item, const char *id)
{
    this->intRecord[id] = item;
}

void DynamicInputRecord :: setField(double item, const char *id)
{
    this->doubleRecord[id] = item;
}

void DynamicInputRecord :: setField(bool item, const char *id)
{
    this->boolRecord[id] = item;
}

void DynamicInputRecord :: setField(const std::string &item, const char *id)
{
    this->stringRecord[id] = item;
}

void DynamicInputRecord :: setField(const FloatArray &item, const char *id)
{
    this->floatArrayRecord[id] = item;
}

void DynamicInputRecord :: setField(const IntArray &item, const char *id)
{
    this->intArrayRecord[id] = item;
}

void DynamicInputRecord :: setField(const FloatMatrix &item, const char *id)
{
    this->matrixRecord[id] = item;
}

void DynamicInputRecord :: setField(const std::vector< std::string > &item, const char *id)
{
    this->stringListRecord[id] = item;
}

void DynamicInputRecord :: setField(const Dictionary &item, const char *id)
{
    this->dictionaryRecord[id] = item;
}

void DynamicInputRecord :: setField(const std::list< Range > &item, const char *id)
{
    this->rangeRecord[id] = item;
}

};
