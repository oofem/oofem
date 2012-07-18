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

#include "oofemtxtinputrecord.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "dictionr.h"
#include "dynalist.h"
#include "range.h"

#ifndef __MAKEDEPEND
 #include <cstdlib>
 #include <cstdio>
 #include <cstring>
 #include <cctype>
 #include <ostream>
 #include <sstream>
#endif

namespace oofem {
OOFEMTXTInputRecord :: OOFEMTXTInputRecord() : InputRecord(), tokenizer(), record()
{
}

OOFEMTXTInputRecord :: OOFEMTXTInputRecord(const OOFEMTXTInputRecord &src) : InputRecord(src), tokenizer()
{
    this->record = src.record;
    tokenizer.tokenizeLine(this->record.c_str());
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = src.readFlag [ i ];
    }
}

OOFEMTXTInputRecord :: OOFEMTXTInputRecord(const char *source) : InputRecord(), tokenizer()
{
    this->record = source;
    tokenizer.tokenizeLine(this->record.c_str());
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = false;
    }
}

OOFEMTXTInputRecord &
OOFEMTXTInputRecord :: operator = ( const OOFEMTXTInputRecord & src )
{
    this->record = src.record;
    tokenizer.tokenizeLine(this->record.c_str());
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = src.readFlag [ i ];
    }

    return * this;
}

void
OOFEMTXTInputRecord :: setRecordString(const std::string &newRec)
{
    this->record = newRec;
    tokenizer.tokenizeLine(this->record.c_str());
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = false;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveRecordKeywordField(std::string &answer, int &value)
{
    if ( tokenizer.giveNumberOfTokens() > 0 ) {
        answer = std::string(tokenizer.giveToken(1));
        setReadFlag(1);
        if ( scanInteger(tokenizer.giveToken(2), value) == 0 ) {
            return IRRT_BAD_FORMAT;
        }
        setReadFlag(2);

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveRecordKeywordField(std::string &answer)
{
    if ( tokenizer.giveNumberOfTokens() > 0 ) {
        answer = std::string(tokenizer.giveToken(1));
        setReadFlag(1);

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(int &answer, InputFieldType fieldID, const char *idString)
{
    int indx = this->giveKeywordIndx(idString);
    if ( indx ) {
        if ( scanInteger(tokenizer.giveToken(indx + 1), answer) == 0 ) {
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);
        setReadFlag(indx + 1);
        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(double &answer, InputFieldType fieldID, const char *idString)
{
    int indx = this->giveKeywordIndx(idString);
    if ( indx ) {
        if ( scanDouble(tokenizer.giveToken(indx + 1), answer) == 0 ) {
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);
        setReadFlag(indx + 1);
        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(std::string &answer, InputFieldType fieldI, const char *idString)
{
    int indx = 0;
    if ( idString ) {
        if ( ( indx = this->giveKeywordIndx(idString) ) == 0 ) {
            return IRRT_NOTFOUND;
        }

        setReadFlag(indx);
        indx++;
    } else {
        indx = 1;
    }

    const char *_token = tokenizer.giveToken(indx);
    if ( _token ) {
        answer = std::string(tokenizer.giveToken(indx));
        setReadFlag(indx);
        return IRRT_OK;
    } else {
        answer = "";
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(IntArray &answer, InputFieldType fieldID, const char *idString)
{
    int value, size;
    int indx = this->giveKeywordIndx(idString);
    if ( indx ) {
        setReadFlag(indx);
        if ( scanInteger(tokenizer.giveToken(++indx), size) == 0 ) {
            return IRRT_BAD_FORMAT;
        }

        answer.resize(size);
        setReadFlag(indx);

        for ( int i = 1; i <= size; i++ ) {
            if ( scanInteger(tokenizer.giveToken(indx + i), value) == 0 ) {
                return IRRT_BAD_FORMAT;
            }

            answer.at(i) = value;
            setReadFlag(indx + i);
        }

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}


IRResultType
OOFEMTXTInputRecord :: giveField(FloatArray &answer, InputFieldType fieldID, const char *idString)
{
    double value;
    int size;
    int indx = this->giveKeywordIndx(idString);
    if ( indx ) {
        setReadFlag(indx);
        if ( scanInteger(tokenizer.giveToken(++indx), size) == 0 ) {
            return IRRT_BAD_FORMAT;
        }

        answer.resize(size);
        setReadFlag(indx);

        for ( int i = 1; i <= size; i++ ) {
            if ( scanDouble(tokenizer.giveToken(indx + i), value) == 0 ) {
                return IRRT_BAD_FORMAT;
            }

            answer.at(i) = value;
            setReadFlag(indx + i);
        }

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}


IRResultType
OOFEMTXTInputRecord :: giveField(FloatMatrix &answer, InputFieldType fieldID, const char *idString)
{
    int nrows, ncols;
    int indx = this->giveKeywordIndx(idString);
    if ( indx ) {
        setReadFlag(indx);

        if ( scanInteger(tokenizer.giveToken(++indx), nrows) == 0 ) {
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);
        if ( scanInteger(tokenizer.giveToken(++indx), ncols) == 0 ) {
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);

        if ( readMatrix(tokenizer.giveToken(++indx), nrows, ncols, answer) == 0 ) {
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);
        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}


IRResultType
OOFEMTXTInputRecord :: giveField(std::vector< std::string > &answer, InputFieldType fieldID, const char *idString)
{
    int size;
    int indx = this->giveKeywordIndx(idString);
    if ( indx ) {
        setReadFlag(indx);
        if ( scanInteger(tokenizer.giveToken(++indx), size) == 0 ) {
            return IRRT_BAD_FORMAT;
        }
        answer.reserve(size);
        setReadFlag(indx);
        for ( int i = 1; i <= size; i++ ) {
            answer.push_back(tokenizer.giveToken(indx + i));
            setReadFlag(indx + i);
        }

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}


IRResultType
OOFEMTXTInputRecord :: giveField(Dictionary &answer, InputFieldType fieldID, const char *idString)
{
    double value;
    int size;
    char key;
    int indx = this->giveKeywordIndx(idString);
    if ( indx ) {
        setReadFlag(indx);
        if ( scanInteger(tokenizer.giveToken(++indx), size) == 0 ) {
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);

        answer.clear();
        for ( int i = 1; i <= size; i++ ) {
            key = ( tokenizer.giveToken(++indx) ) [ 0 ];
            setReadFlag(indx);
            if ( scanDouble(tokenizer.giveToken(++indx), value) == 0 ) {
                return IRRT_BAD_FORMAT;
            }

            setReadFlag(indx);
            answer.add(key, value);
        }

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(dynaList< Range > &list, InputFieldType fieldID, const char *idString)
{
    int li, hi;
    const char *rec;
    int indx = this->giveKeywordIndx(idString);
    if ( indx ) {
        setReadFlag(indx);
        rec = tokenizer.giveToken(++indx);
        if ( * rec != '{' ) {
            OOFEM_WARNING("OOFEMTXTInputRecord::readRangeList: parse error - missing left '{'");
            list.clear();
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);
        rec++;
        // read ranges
        while ( readRange(& rec, li, hi) ) {
            Range range(li, hi);
            list.pushBack(range);
        }

        // skip whitespaces after last range
        while ( isspace(* rec) ) {
            rec++;
        }

        // test for enclosing bracket
        if ( * rec != '}' ) {
            OOFEM_WARNING("OOFEMTXTInputRecord::readRangeList: parse error - missing end '}'");
            list.clear();
            return IRRT_BAD_FORMAT;
        }

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(double &answer, int tokenNumber){
    if ( scanDouble(tokenizer.giveToken(tokenNumber), answer) == 0 ) {
        OOFEM_ERROR4("Double not found on line %d in token number %d, string %s", this->lineNumber, tokenNumber, this->record.c_str() );
        return IRRT_BAD_FORMAT;
    }

    setReadFlag(tokenNumber);
    return IRRT_OK;
}


bool
OOFEMTXTInputRecord :: hasField(InputFieldType fieldID, const char *idString)
{
    //returns nonzero if idString is present in source
    int indx = this->giveKeywordIndx(idString);
    if ( indx ) {
        setReadFlag(indx);
    }

    return ( indx > 0 ) ? true : false;
}

void
OOFEMTXTInputRecord :: printYourself()
{
  printf( "%s", this->record.c_str() );
}


int
OOFEMTXTInputRecord :: scanInteger(const char *source, int &value)
{
    //
    // reads integer value from source, returns nonzero if o.k.
    //
    int i;
    if ( source == NULL ) {
        value = 0;
        return 0;
    }

    i = sscanf(source, "%d", & value);
    if ( ( i == EOF ) || ( i == 0 ) ) {
        value = 0;
        return 0;
    }

    return 1;
}


int
OOFEMTXTInputRecord :: scanDouble(const char *source, double &value)
{
    //
    // reads integer value from source, returns pointer to char after this number
    //
    int i;

    if ( source == NULL ) {
        value = 0.0;
        return 0;
    }

    i = sscanf(source, "%lf", & value);
    if ( ( i == EOF ) || ( i == 0 ) ) {
        value = 0.0;
        return 0;
    }

    return 1;
}

int
OOFEMTXTInputRecord :: giveKeywordIndx(const char *kwd)
{
    int ntokens = tokenizer.giveNumberOfTokens();
    for ( int i = 1; i <= ntokens; i++ ) {
        if ( strcmp( kwd, tokenizer.giveToken(i) ) == 0 ) {
            return i;
        }
    }

    return 0;
}

void
OOFEMTXTInputRecord :: finish(bool wrn)
{
    if ( !wrn ) {
        return;
    }

    std::ostringstream buff;
    bool pf = true, wf = false;
    int ntokens = tokenizer.giveNumberOfTokens();
    for ( int i = 0; i < ntokens; i++ ) {
        //fprintf (stderr, "[%s] ", tokenizer.giveToken(i+1));
        if ( !readFlag [ i ] ) {
            if ( pf ) {
                buff << "Unread token(s) detected in the following record\n\"";
                for ( int j = 0; j < 40; j++ ) {
                    if ( this->record [ j ] == '\n' || this->record [ j ] == '\0' ) {
                        break;
                    } else {
                        buff << this->record [ j ];
                    }
                }
                if ( this->record.size() > 41 ) {
                    buff << "...";
                }
                buff << "\":\n";

                pf = false;
                wf = true;
            }

            buff << "[" << tokenizer.giveToken(i + 1) << "]";
        }
    }

    if ( wf ) {
        OOFEM_WARNING(buff.str().c_str());
    }
}

char *
OOFEMTXTInputRecord :: __readSimpleString(const char *source, char *simpleString, int maxchar, const char **remain)
// reads Simple string from source according to following rules:
// at beginning skips whitespace (blank, tab)
// read string terminated by whitespace or end-of-line
// remain is unread remain of source string.
// maximum of maxchar (including terminating '\0') is copied into simpleString.
{
    const char *curr = source;
    char *ss = simpleString;
    int count = 0;

    if ( source == NULL ) {
        * remain = NULL;
        return NULL;
    }

    while ( isspace(* curr) || !* curr ) {
        curr++;
    }

    if ( !curr ) {
        OOFEM_ERROR("OOFEMTXTInputRecord::readSimpleString : unexpected end-of-line");
    }

    while ( ( !( isspace(* curr) || !* curr ) ) && ( ++count < maxchar ) ) {
        * ss++ = * curr++;
    }

    * ss = '\0';
    * remain = curr;
    return simpleString;
}

const char *
OOFEMTXTInputRecord :: __readKeyAndVal(const char *source, char *key, int *val, int maxchar, const char **remain)
//
//
//
{
    key = __readSimpleString(source, key, maxchar, remain);
    * remain = __scanInteger(* remain, val);
    return * remain;
}

const char *
OOFEMTXTInputRecord :: __readKeyAndVal(const char *source, char *key, double *val, int maxchar, const char **remain)
//
//
//
{
    key = __readSimpleString(source, key, maxchar, remain);
    * remain = __scanDouble(* remain, val);
    return * remain;
}

const char *
OOFEMTXTInputRecord :: __getPosAfter(const char *source, const char *idString)
//
// returns position of substring idString in source
// return value pointer at the end of occurrence idString in source
// (idString must be separated from rest by blank or by tabulator
// if string not found, returns NULL
//
{
    const char *str1;
    const char *helpSource = source;
    int len = strlen(idString);
    int whitespaceBefore, whitespaceAfter;

    do {
        if ( ( str1 = strstr(helpSource, idString) ) == NULL ) {
            return NULL;
        }

        helpSource = str1 + 1;
        whitespaceAfter = isspace( * ( helpSource + len - 1 ) );
        if ( str1 == source ) {
            whitespaceBefore = 1;
        } else {
            whitespaceBefore = isspace( * ( str1 - 1 ) );
        }
    } while ( !( whitespaceBefore && whitespaceAfter ) );

    return str1 + len;
}

const char *
OOFEMTXTInputRecord :: __skipNextWord(const char *src)
//
// skips next word in src ; returns pointer after it
//
{
    while ( isspace(* src) || !* src ) {
        src++;
    }

    // skips whitespaces if any
    while ( !( isspace(* src) || !* src ) ) {
        src++;
    }

    // skips one word
    return src;
}

const char *
OOFEMTXTInputRecord :: __scanInteger(const char *source, int *value)
{
    //
    // reads integer value from source, returns pointer to char after this number
    //
    char *endptr;

    if ( source == NULL ) {
        * value = 0;
        return NULL;
    }

    * value = ( int ) strtol(source, & endptr, 10);
    //  i = sscanf (source,"%d",value);
    // if (i == EOF ){ *value =0 ; return NULL ;}
    return endptr;
}

const char *
OOFEMTXTInputRecord :: __scanDouble(const char *source, double *value)
{
    //
    // reads integer value from source, returns pointer to char after this number
    //
    char *endptr;

    if ( source == NULL ) {
        * value = 0;
        return NULL;
    }

    * value = strtod(source, & endptr);
    //i = sscanf (source,"%lf",value);
    //if (i == EOF ){ *value =0 ; return NULL ;}
    //return __skipNextWord(source);
    return endptr;
}

int
OOFEMTXTInputRecord :: readRange(const char **helpSource, int &li, int &hi)
{
    char *endptr;
    // skip whitespaces
    while ( isspace(* * helpSource) ) {
        ( * helpSource )++;
    }

    // test if character is digit
    if ( isdigit(* * helpSource) ) {
        // digit character - read one value range
        li = hi = strtol(* helpSource, & endptr, 10);
        * helpSource = endptr;
        return 1;
    } else if ( * * helpSource == '(' ) {
        // range left parenthesis found
        ( * helpSource )++;
        // read lower index
        li = strtol(* helpSource, & endptr, 10);
        * helpSource = endptr;
        // test whitespaces
        if ( * * helpSource != ' ' && * * helpSource != '\t' ) {
            OOFEM_WARNING("OOFEMTXTInputRecord::readRange: unexpected token while reading range value");
            return 0;
        }

        // read end index
        hi = strtol(* helpSource, & endptr, 10);
        * helpSource = endptr;
        // skip whitespaces
        while ( isspace(* * helpSource) ) {
            ( * helpSource )++;
        }

        // test for enclosing bracket
        if ( * * helpSource == ')' ) {
            ( * helpSource )++;
            return 1;
        } else {
            OOFEM_WARNING("OOFEMTXTInputRecord::readRange: end ')' missing while parsing range value");
            return 0;
        }
    }

    return 0;
}

int
OOFEMTXTInputRecord :: readMatrix(const char *helpSource, int r, int c, FloatMatrix &ans)
{
    const char *endptr = helpSource;
    int i, j;

    if ( helpSource == NULL ) {
        ans.resize(0, 0);
        return 0;
    }

    ans.resize(r, c);
    // skip whitespaces
    while ( isspace(* endptr) ) {
        ( endptr )++;
    }

    if ( * endptr == '{' ) {
        // range left parenthesis found
        ( endptr )++;
        // read row by row separated by semicolon
        for ( i = 1; i <= r; i++ ) {
            for ( j = 1; j <= c; j++ ) {
                endptr = __scanDouble( endptr, & ans.at(i, j) );
            }

            if ( i < r ) {
                // skip whitespaces
                while ( isspace(* endptr) ) {
                    ( endptr )++;
                }

                // test for row terminating semicolon
                if ( * endptr == ';' ) {
                    ( endptr )++;
                } else {
                    OOFEM_WARNING("OOFEMTXTInputRecord::readMatrix: missing row terminating semicolon");
                    return 0;
                }
            }
        }

        // skip whitespaces
        while ( isspace(* endptr) ) {
            ( endptr )++;
        }

        // test for enclosing bracket
        if ( * endptr == '}' ) {
            ( endptr )++;
            helpSource = endptr;
            return 1;
        } else {
            OOFEM_WARNING("OOFEMTXTInputRecord::readMatrix: end '}' missing while parsing matrix value");
            return 0;
        }
    } else {
        return 0;
    }

    return 0;
}
} // end namespace oofem
