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

#ifndef dynamicinputrecord_h
#define dynamicinputrecord_h

#include "inputrecord.h"
#include <map>

namespace oofem {

/**
 * Class representing the a dynamic Input Record.
 * The input record is represented as a list of fields.
 * This is intended for internal usage, where new input records and such are created dynamically.
 * @author Mikael Öhman
 */
class DynamicInputRecord : public InputRecord
{
protected:
    std::string recordKeyword;
    int recordNumber;

    // Record representation.
    //std::map<InputFieldType, void*> record;
    std::map<InputFieldType, int>  intRecord;
    std::map<InputFieldType, double> doubleRecord;
    std::map<InputFieldType, std::string> stringRecord;
    std::map<InputFieldType, FloatArray> floatArrayRecord;
    std::map<InputFieldType, IntArray> intArrayRecord;
    std::map<InputFieldType, FloatMatrix> matrixRecord;
    std::map<InputFieldType, std::vector< std::string> > stringListRecord;
    std::map<InputFieldType, Dictionary> dictionaryRecord;
    std::map<InputFieldType, dynaList< Range > > rangeRecord;

public:
    /// Constructor. Creates an empty input record.
    DynamicInputRecord();
    /// Copy constructor.
    DynamicInputRecord(const DynamicInputRecord &);
    /// Destructor.
    virtual ~DynamicInputRecord();
    /// Assignment operator.
    DynamicInputRecord & operator=(const DynamicInputRecord &);

    virtual InputRecord *GiveCopy() { return new DynamicInputRecord(* this); }
    virtual void finish(bool wrn = true);

    virtual IRResultType giveRecordKeywordField(std::string &answer, int &value);
    virtual IRResultType giveRecordKeywordField(std::string &answer);
    virtual IRResultType giveField(int &answer, InputFieldType fieldID, const char *idString);
    virtual IRResultType giveField(double &answer, InputFieldType fieldID, const char *idString);
    virtual IRResultType giveField(std::string &answer, InputFieldType fieldID, const char *idString);
    virtual IRResultType giveField(FloatArray &answer, InputFieldType fieldID, const char *idString);
    virtual IRResultType giveField(IntArray &answer, InputFieldType fieldID, const char *idString);
    virtual IRResultType giveField(FloatMatrix &answer, InputFieldType fieldID, const char *idString);
    virtual IRResultType giveField(std::vector< std::string > &answer, InputFieldType fieldID, const char *idString);
    virtual IRResultType giveField(Dictionary &answer, InputFieldType fieldID, const char *idString);
    virtual IRResultType giveField(dynaList< Range > &answer, InputFieldType fieldID, const char *idString);
    virtual IRResultType giveField(double &answer, int tokenNumber);
    
    virtual bool hasField(InputFieldType fieldID, const char *idString);
    virtual void printYourself();
    // Setters, unique for the dynamic input record
    virtual void setRecordKeywordField(const std::string &keyword, int number);
    virtual void setField(int item, InputFieldType fieldID);
    virtual void setField(double item, InputFieldType fieldID);
    virtual void setField(const std::string &item, InputFieldType fieldID);
    virtual void setField(const FloatArray &item, InputFieldType fieldID);
    virtual void setField(const IntArray &item, InputFieldType fieldID);
    virtual void setField(const FloatMatrix &item, InputFieldType fieldID);
    virtual void setField(const std::vector< std::string > &item, InputFieldType fieldID);
    virtual void setField(const Dictionary &item, InputFieldType fieldID);
    virtual void setField(const dynaList< Range > &item, InputFieldType fieldID);

};
} // end namespace oofem
#endif // dynamicinputrecord_h
