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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
#include "dynalist.h"
#include "range.h"

namespace oofem {

DynamicInputRecord :: DynamicInputRecord() : InputRecord()
{
    this->recordKeyword = "";
    this->recordNumber = 0;
}

DynamicInputRecord :: DynamicInputRecord(const DynamicInputRecord &src)
{
    this->recordKeyword = src.recordKeyword;
    this->recordNumber = src.recordNumber;
    this->intRecord = src.intRecord;
    this->doubleRecord = src.doubleRecord;
    this->stringRecord = src.stringRecord;
    this->floatArrayRecord = src.floatArrayRecord;
    this->intArrayRecord = src.intArrayRecord;
    this->matrixRecord = src.matrixRecord;
    this->stringListRecord = src.stringListRecord;
    this->dictionaryRecord = src.dictionaryRecord;
    this->rangeRecord = src.rangeRecord;
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

IRResultType DynamicInputRecord :: giveField(int &answer, InputFieldType fieldID, const char *idString)
{
    std::map<InputFieldType, int>::iterator it = this->intRecord.find(fieldID);
    if (it == this->intRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(double &answer, InputFieldType fieldID, const char *idString)
{
    std::map<InputFieldType, double>::iterator it = this->doubleRecord.find(fieldID);
    if (it == this->doubleRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(std::string &answer, InputFieldType fieldID, const char *idString)
{
    std::map<InputFieldType, std::string>::iterator it = this->stringRecord.find(fieldID);
    if (it == this->stringRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(FloatArray &answer, InputFieldType fieldID, const char *idString)
{
    std::map<InputFieldType, FloatArray>::iterator it = this->floatArrayRecord.find(fieldID);
    if (it == this->floatArrayRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(IntArray &answer, InputFieldType fieldID, const char *idString)
{
    std::map<InputFieldType, IntArray>::iterator it = this->intArrayRecord.find(fieldID);
    if (it == this->intArrayRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(FloatMatrix &answer, InputFieldType fieldID, const char *idString)
{
    std::map<InputFieldType, FloatMatrix>::iterator it = this->matrixRecord.find(fieldID);
    if (it == this->matrixRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(std::vector< std::string > &answer, InputFieldType fieldID, const char *idString)
{
    std::map<InputFieldType, std::vector< std::string> >::iterator it = this->stringListRecord.find(fieldID);
    if (it == this->stringListRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(Dictionary &answer, InputFieldType fieldID, const char *idString)
{
    std::map<InputFieldType, Dictionary>::iterator it = this->dictionaryRecord.find(fieldID);
    if (it == this->dictionaryRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

IRResultType DynamicInputRecord :: giveField(dynaList< Range > &answer, InputFieldType fieldID, const char *idString)
{
    std::map<InputFieldType, dynaList< Range > >::iterator it = this->rangeRecord.find(fieldID);
    if (it == this->rangeRecord.end())
        return IRRT_NOTFOUND;
    answer = it->second;
    return IRRT_OK;
}

bool DynamicInputRecord :: hasField(InputFieldType fieldID, const char *idString)
{
    return this->intRecord.find(fieldID) != this->intRecord.end() ||
           this->doubleRecord.find(fieldID) != this->doubleRecord.end() ||
           this->floatArrayRecord.find(fieldID) != this->floatArrayRecord.end() ||
           this->intArrayRecord.find(fieldID) != this->intArrayRecord.end() ||
           this->matrixRecord.find(fieldID) != this->matrixRecord.end() ||
           this->stringListRecord.find(fieldID) != this->stringListRecord.end() ||
           this->dictionaryRecord.find(fieldID) != this->dictionaryRecord.end() ||
           this->rangeRecord.find(fieldID) != this->rangeRecord.end();
}

// Setters
void DynamicInputRecord :: setRecordKeywordField(const std::string &keyword, int value)
{
    this->recordKeyword = keyword;
    this->recordNumber = value;
}

void DynamicInputRecord :: setField(int item, InputFieldType fieldID)
{
    this->intRecord[fieldID] = item;
}

void DynamicInputRecord :: setField(double item, InputFieldType fieldID)
{
    this->doubleRecord[fieldID] = item;
}

void DynamicInputRecord :: setField(const std::string &item, InputFieldType fieldID)
{
    this->stringRecord[fieldID] = item;
}

void DynamicInputRecord :: setField(const FloatArray &item, InputFieldType fieldID)
{
    this->floatArrayRecord[fieldID] = item;
}

void DynamicInputRecord :: setField(const IntArray &item, InputFieldType fieldID)
{
    this->intArrayRecord[fieldID] = item;
}

void DynamicInputRecord :: setField(const FloatMatrix &item, InputFieldType fieldID)
{
    this->matrixRecord[fieldID] = item;
}

void DynamicInputRecord :: setField(const std::vector< std::string > &item, InputFieldType fieldID)
{
    this->stringListRecord[fieldID] = item;
}

void DynamicInputRecord :: setField(const Dictionary &item, InputFieldType fieldID)
{
    this->dictionaryRecord[fieldID] = item;
}

void DynamicInputRecord :: setField(const dynaList< Range > &item, InputFieldType fieldID)
{
    this->rangeRecord[fieldID] = item;
}

};
