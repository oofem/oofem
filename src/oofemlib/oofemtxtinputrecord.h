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

#ifndef oofemtxtinputrecord_h
#define oofemtxtinputrecord_h

#include "inputrecord.h"
#include "tokenizer.h"

#include <string>
#include <vector>

namespace oofem {
/**
 * Class representing the Input Record for OOFEM txt input file format.
 * The input record is represented as string consisting of several fields.
 */
class OOFEM_EXPORT OOFEMTXTInputRecord : public InputRecord
{
protected:
    /**
     * Tokenizer is used to parse records.
     * This also enables to perform additional check for input records, since
     * unread fields can be detected
     */
    Tokenizer tokenizer;
    std :: vector< bool >readFlag;

    /// Record representation.
    std :: string record;

    int lineNumber;

public:
    /// Constructor. Creates an empty input record.
    OOFEMTXTInputRecord();
    /// Constructor. Creates the input record corresponding to given string.
    OOFEMTXTInputRecord(int linenumber, std :: string source);
    /// Copy constructor.
    OOFEMTXTInputRecord(const OOFEMTXTInputRecord &);
    /// Assignment operator.
    OOFEMTXTInputRecord &operator = ( const OOFEMTXTInputRecord & );

    std::unique_ptr<InputRecord> clone() const override { return std::make_unique<OOFEMTXTInputRecord>(*this); }

public:
    /// Sets the record string.
    void setRecordString(std :: string newStr);
    /// Returns record string.
    std :: string giveRecordAsString() const override { return this->record; }

    void finish(bool wrn = true) override;

public:
    void giveRecordKeywordField(std :: string &answer, int &value) override;
    void giveRecordKeywordField(std :: string &answer) override;
    void giveField(int &answer, InputFieldType id) override;
    void giveField(double &answer, InputFieldType id) override;
    void giveField(bool &answer, InputFieldType id) override;
    void giveField(std :: string &answer, InputFieldType id) override;
    void giveField(FloatArray &answer, InputFieldType id) override;
    void giveField(IntArray &answer, InputFieldType id) override;
    void giveField(FloatMatrix &answer, InputFieldType id) override;
    void giveField(std :: vector< std :: string > &answer, InputFieldType id) override;
    void giveField(Dictionary &answer, InputFieldType id) override;
    void giveField(std :: list< Range > &answer, InputFieldType id) override;
    void giveField(ScalarFunction &answer, InputFieldType id) override;

    bool hasField(InputFieldType id) override;
    void printYourself() override;

    void setLineNumber(int num) { this->lineNumber = num; }

protected:
    int giveKeywordIndx(const char *kwd);
    const char *scanInteger(const char *source, int &value);
    const char *scanDouble(const char *source, double &value);
    void setReadFlag(int itok) { readFlag [ itok - 1 ] = true; }

    /**
     * Reads single range record from input record represented by *helpSource  string.
     * @param helpSource Pointer to current string position, on return helpSource points
     * to next character after reading range record.
     * @param li Starting range index.
     * @param hi End range index.
     * @return Nonzero on success.
     */
    int readRange(const char **helpSource, int &li, int &hi);
    /**
     * Reads single matrix record from input record represented by *helpSource  string.
     * @param helpSource Pointer to current string position, on return helpSource points
     * to next character after reading range record.
     * @param r Matrix rows.
     * @param c Matrix columns.
     * @param ans Float matrix.
     * @return Nonzero on success.
     */
    int readMatrix(const char *helpSource, int r, int c, FloatMatrix &ans);
};
} // end namespace oofem
#endif // oofemtxtinputrecord_h
