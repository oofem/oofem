/* $Header: /home/cvs/bp/oofem/oofemlib/src/inputrecord.h,v 1.4 2003/05/19 13:03:57 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#ifndef oofemtxtinputrecord_h
#define oofemtxtinputrecord_h

#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#endif
#include "inputrecord.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "dictionr.h"
#include "dynalist.h"
#include "range.h"

#define IR_MAX_ERROR_LENGTH 100

/** Tokenizer is used to parse records.
 *  This also enables to perform additional check for input records, since
 *  unread fields can be detected
 */
#include "tokenizer.h"




/**
 * Class representing the Input Record for OOFEM txt input file format. The input record is represented as string consisting of several fields.
 */
class OOFEMTXTInputRecord : public InputRecord
{
protected:

    Tokenizer tokenizer;
    bool readFlag [ MAX_TOKENS ];

    // record representation
    char record [ OOFEM_MAX_LINE_LENGTH ];

public:
    /// Constructor. Creates an empty input record.
    OOFEMTXTInputRecord();
    /// Constructor. Creates the input record corresponding to given string
    OOFEMTXTInputRecord(char *source);
    /// Copy constructor
    OOFEMTXTInputRecord(const OOFEMTXTInputRecord &);
    /// Destructor
    ~OOFEMTXTInputRecord() { }
    /// Assingnment operator
    OOFEMTXTInputRecord &operator=(const OOFEMTXTInputRecord &);

    /** Creates a newly allocated copy of the receiver */
    virtual InputRecord *GiveCopy() { return new OOFEMTXTInputRecord(* this); }

public:
    /// Sets the record string
    void setRecordString(char *);
    /// Returns record string
    char *giveRecordAsString() { return this->record; }

    /** terminates the current record session and if flag is true warnin is printed for unscanned tokens */
    void finish(bool wrn = true);



public:
    /**@name Compulsory field extraction methods
     * Reads the field value identified by keyword
     * @param answer contains result
     * @param idString field keyword
     * @return IRResultType
     */
    //@{
    /// Reads the record id field  (type of record) and its corresponding number
    virtual IRResultType giveRecordKeywordField(char *answer, int &value, int maxchar);
    /// Reads the record id field  (type of record)
    virtual IRResultType giveRecordKeywordField(char *answer, int maxchar);
    /// Reads the integer field value
    virtual IRResultType giveField(int &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the double field value
    virtual IRResultType giveField(double &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the char* field value
    virtual IRResultType giveField(char *answer, int maxchar, const InputFieldType fieldI, const char *idString);
    /// Reads the FloatArray field value
    virtual IRResultType giveField(FloatArray &answer, const InputFieldType fieldI, const char *idString);
    /// Reads the IntArray field value
    virtual IRResultType giveField(IntArray &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the FloatMatrix field value
    virtual IRResultType giveField(FloatMatrix &answer, const InputFieldType fieldI, const char *idString);
    /// Reads the Dictionary field value
    virtual IRResultType giveField(Dictionary &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the dynaList<Range> field value
    virtual IRResultType giveField(dynaList< Range > &answer, const InputFieldType fieldID, const char *idString);
    //@}

    /**@name Optional field extraction methods
     * Reads the field value identified by keyword
     * @param answer contains result
     * @param idString field keyword
     * @return IRResultType
     */
    //@{
    /// Reads the integer field value
    virtual IRResultType giveOptionalField(int &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the double field value
    virtual IRResultType giveOptionalField(double &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the char* field value
    virtual IRResultType giveOptionalField(char *answer, int maxchar, const InputFieldType fieldID, const char *idString);
    /// Reads the FloatArray field value
    virtual IRResultType giveOptionalField(FloatArray &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the IntArray field value
    virtual IRResultType giveOptionalField(IntArray &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the FloatMatrix field value
    virtual IRResultType giveOptionalField(FloatMatrix &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the Dictionary field value
    virtual IRResultType giveOptionalField(Dictionary &answer, const InputFieldType fieldID, const char *idString);
    /// Reads the dynaList<Range> field value
    virtual IRResultType giveOptionalField(dynaList< Range > &answer, const InputFieldType fieldID, const char *idString);
    //@}

    /// Returns true if record contains field identified by idString keyword
    virtual bool         hasField(const InputFieldType fieldID, const char *idString);

protected:

    int giveKeywordIndx(const char *kwd);
    int scanInteger(const char *source, int &value);
    int scanDouble(const char *source, double &value);
    void setReadFlag(int itok) { readFlag [ itok - 1 ] = true; }


    const char *__getPosAfter(const char *, const char *);
    const char *__scanInteger(const char *source, int *value);
    const char *__scanDouble(const char *source, double *value);
    const char *__skipNextWord(const char *src);

    char *__readSimpleString(const char *source, char *simpleString, int maxchar, const char **remain);
    const char *__readKeyAndVal(const char *source, char *key, int *val, int maxchar, const char **remain);
    const char *__readKeyAndVal(const char *source, char *key, double *val, int maxchar, const char **remain);

    /**
     * Reads single range record from input record represented by *helpSource  string.
     * @param helpSource pointer to current string possition, on return helpSource points
     * to next charcter after reading range record.
     * @param li starting range index
     * @param hi end range index
     * @return on success nonzero valur returned
     */
    int    readRange(const char **helpSource, int &li, int &hi);
    /**
     * Reads single matrix record from input record represented by *helpSource  string.
     * @param helpSource pointer to current string possition, on return helpSource points
     * to next charcter after reading range record.
     * @param r,c matrix dimensions
     * @param ans float matrix
     * @return on success nonzero valur returned
     */
    int    readMatrix(const char *helpSource, int r, int c, FloatMatrix &ans);
};

#endif // oofemtxtinputrecord_h
