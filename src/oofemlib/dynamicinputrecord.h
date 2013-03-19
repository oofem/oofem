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
    //std::map<std::string, void*> record;
    std::map<std::string, int>  intRecord;
    std::map<std::string, double> doubleRecord;
    std::map<std::string, bool>  boolRecord;
    std::map<std::string, std::string> stringRecord;
    std::map<std::string, FloatArray> floatArrayRecord;
    std::map<std::string, IntArray> intArrayRecord;
    std::map<std::string, FloatMatrix> matrixRecord;
    std::map<std::string, std::vector< std::string> > stringListRecord;
    std::map<std::string, Dictionary> dictionaryRecord;
    std::map<std::string, std::list< Range > > rangeRecord;

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
    virtual IRResultType giveField(int &answer, InputFieldType fieldID, const char *id);
    virtual IRResultType giveField(double &answer, InputFieldType fieldID, const char *id);
    virtual IRResultType giveField(bool &answer, InputFieldType fieldID, const char *id);
    virtual IRResultType giveField(std::string &answer, InputFieldType fieldID, const char *id);
    virtual IRResultType giveField(FloatArray &answer, InputFieldType fieldID, const char *id);
    virtual IRResultType giveField(IntArray &answer, InputFieldType fieldID, const char *id);
    virtual IRResultType giveField(FloatMatrix &answer, InputFieldType fieldID, const char *id);
    virtual IRResultType giveField(std::vector< std::string > &answer, InputFieldType fieldID, const char *id);
    virtual IRResultType giveField(Dictionary &answer, InputFieldType fieldID, const char *id);
    virtual IRResultType giveField(std::list< Range > &answer, InputFieldType fieldID, const char *id);
    
    virtual bool hasField(InputFieldType fieldID, const char *id);
    virtual void printYourself();
    // Setters, unique for the dynamic input record
    virtual void setRecordKeywordField(const std::string &keyword, int number);
    virtual void setField(int item, const char *id);
    virtual void setField(double item, const char *id);
    virtual void setField(bool item, const char *id);
    virtual void setField(const std::string &item, const char *id);
    virtual void setField(const FloatArray &item, const char *id);
    virtual void setField(const IntArray &item, const char *id);
    virtual void setField(const FloatMatrix &item, const char *id);
    virtual void setField(const std::vector< std::string > &item, const char *id);
    virtual void setField(const Dictionary &item, const char *id);
    virtual void setField(const std::list< Range > &item, const char *id);

};
} // end namespace oofem
#endif // dynamicinputrecord_h
