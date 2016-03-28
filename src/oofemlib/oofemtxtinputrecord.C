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

#include "oofemtxtinputrecord.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dictionary.h"
#include "range.h"
#include "scalarfunction.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <ostream>
#include <sstream>

namespace oofem {
OOFEMTXTInputRecord :: OOFEMTXTInputRecord() : InputRecord(), tokenizer(), record()
{ }

OOFEMTXTInputRecord :: OOFEMTXTInputRecord(const OOFEMTXTInputRecord &src) : InputRecord(src), tokenizer(),
    record(src.record), lineNumber(src.lineNumber)
{
    tokenizer.tokenizeLine( this->record );
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = src.readFlag [ i ];
    }
}

OOFEMTXTInputRecord :: OOFEMTXTInputRecord(int linenumber, std :: string source) : InputRecord(), tokenizer(),
    record(std :: move(source)), lineNumber(linenumber)
{
    tokenizer.tokenizeLine( this->record );
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
    tokenizer.tokenizeLine( this->record );
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = src.readFlag [ i ];
    }

    return * this;
}

void
OOFEMTXTInputRecord :: setRecordString(std :: string newRec)
{
    this->record = std :: move(newRec);
    tokenizer.tokenizeLine( this->record );
    int ntok = tokenizer.giveNumberOfTokens();
    readFlag.resize(ntok);
    for ( int i = 0; i < ntok; i++ ) {
        readFlag [ i ] = false;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveRecordKeywordField(std :: string &answer, int &value)
{
    if ( tokenizer.giveNumberOfTokens() > 0 ) {
        answer = std :: string( tokenizer.giveToken(1) );
        setReadFlag(1);
        auto ptr = scanInteger(tokenizer.giveToken(2), value);
        if ( ptr == NULL || *ptr != 0 ) {
            return IRRT_BAD_FORMAT;
        }
        setReadFlag(2);

        return IRRT_OK;
    } else {
        return IRRT_BAD_FORMAT;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveRecordKeywordField(std :: string &answer)
{
    if ( tokenizer.giveNumberOfTokens() > 0 ) {
        answer = std :: string( tokenizer.giveToken(1) );
        setReadFlag(1);

        return IRRT_OK;
    } else {
        return IRRT_BAD_FORMAT;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(int &answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        auto ptr = scanInteger(tokenizer.giveToken(indx + 1), answer);
        if ( ptr == NULL || *ptr != 0 ) {
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
OOFEMTXTInputRecord :: giveField(double &answer, InputFieldType id)
{
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        auto ptr = scanDouble(tokenizer.giveToken(indx + 1), answer);
        if ( ptr == NULL || *ptr != 0 ) {
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
OOFEMTXTInputRecord :: giveField(bool &answer, InputFieldType id)
{
    int val;
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        auto ptr = scanInteger(tokenizer.giveToken(indx + 1), val);
        if ( ptr == NULL || *ptr != 0 ) {
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);
        setReadFlag(indx + 1);
        answer = val != 0;
        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(std :: string &answer, InputFieldType id)
{
    int indx = 0;
    if ( id ) {
        if ( ( indx = this->giveKeywordIndx(id) ) == 0 ) {
            return IRRT_NOTFOUND;
        }

        setReadFlag(indx);
        indx++;
    } else {
        indx = 1;
    }

    const char *_token = tokenizer.giveToken(indx);
    if ( _token ) {
        answer = std :: string( tokenizer.giveToken(indx) );
        setReadFlag(indx);
        return IRRT_OK;
    } else {
        answer = "";
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(IntArray &answer, InputFieldType id)
{
    int value, size;
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        setReadFlag(indx);
        auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
        if ( ptr == NULL || *ptr != 0) {
            return IRRT_BAD_FORMAT;
        }

        answer.resize(size);
        setReadFlag(indx);

        for ( int i = 1; i <= size; i++ ) {
            ptr = scanInteger(tokenizer.giveToken(indx + i), value);
            if ( ptr == NULL || *ptr != 0 ) {
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
OOFEMTXTInputRecord :: giveField(FloatArray &answer, InputFieldType id)
{
    double value;
    int size;
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        setReadFlag(indx);
        auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
        if ( ptr == NULL || *ptr != 0 ) {
            return IRRT_BAD_FORMAT;
        }

        answer.resize(size);
        setReadFlag(indx);

        for ( int i = 1; i <= size; i++ ) {
            auto ptr = scanDouble(tokenizer.giveToken(indx + i), value);
            if ( ptr == NULL || *ptr != 0 ) {
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
OOFEMTXTInputRecord :: giveField(FloatMatrix &answer, InputFieldType id)
{
    int nrows, ncols;
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        setReadFlag(indx);

        auto ptr = scanInteger(tokenizer.giveToken(++indx), nrows);
        if ( ptr == NULL || *ptr != 0 ) {
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);
        ptr = scanInteger(tokenizer.giveToken(++indx), ncols);
        if ( ptr == NULL || *ptr != 0 ) {
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
OOFEMTXTInputRecord :: giveField(std :: vector< std :: string > &answer, InputFieldType id)
{
    int size;
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        setReadFlag(indx);
        auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
        if ( ptr == NULL || *ptr != 0 ) {
            return IRRT_BAD_FORMAT;
        }
        answer.reserve(size);
        setReadFlag(indx);
        for ( int i = 1; i <= size; i++ ) {
            answer.push_back( tokenizer.giveToken(indx + i) );
            setReadFlag(indx + i);
        }

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(Dictionary &answer, InputFieldType id)
{
    double value;
    int size;
    char key;
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        setReadFlag(indx);
        auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
        if ( ptr == NULL || *ptr != 0 ) {
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);

        answer.clear();
        for ( int i = 1; i <= size; i++ ) {
            key = ( tokenizer.giveToken(++indx) ) [ 0 ];
            setReadFlag(indx);
            auto ptr = scanDouble(tokenizer.giveToken(++indx), value);
            if ( ptr == NULL || *ptr != 0 ) {
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
OOFEMTXTInputRecord :: giveField(std :: list< Range > &list, InputFieldType id)
{
    int li, hi;
    const char *rec;
    int indx = this->giveKeywordIndx(id);
    if ( indx ) {
        setReadFlag(indx);
        rec = tokenizer.giveToken(++indx);
        if ( * rec != '{' ) {
            OOFEM_WARNING("missing left '{'");
            list.clear();
            return IRRT_BAD_FORMAT;
        }

        setReadFlag(indx);
        rec++;
        // read ranges
        while ( readRange(& rec, li, hi) ) {
            Range range(li, hi);
            list.push_back(range);
        }

        // skip whitespaces after last range
        while ( isspace(* rec) ) {
            rec++;
        }

        // test for enclosing bracket
        if ( * rec != '}' ) {
            OOFEM_WARNING("missing end '}'");
            list.clear();
            return IRRT_BAD_FORMAT;
        }

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

IRResultType
OOFEMTXTInputRecord :: giveField(ScalarFunction &answer, InputFieldType id)
{
    const char *rec;
    int indx = this->giveKeywordIndx(id);

    if ( indx ) {
        setReadFlag(indx);
        rec = tokenizer.giveToken(++indx);
        if ( * rec == '@' ) {
            // reference to function
            int refVal;
            auto ptr = scanInteger(rec + 1, refVal);
            if ( ptr == NULL || *ptr != 0 ) {
                return IRRT_BAD_FORMAT;
            }
            setReadFlag(indx);
            answer.setReference(refVal);
        } else if ( * rec == '$' ) {
            // simple expression
            std :: string expr;

            expr = std :: string( tokenizer.giveToken(indx) );
            setReadFlag(indx);
            std :: string _v = expr.substr(1, expr.size() - 2);

            answer.setSimpleExpression(_v); // get rid of enclosing '$'
        } else {
            double val;
            auto ptr = scanDouble(tokenizer.giveToken(indx), val);
            if ( ptr == NULL || *ptr != 0 ) {
                return IRRT_BAD_FORMAT;
            }

            setReadFlag(indx);
            answer.setValue(val);
        }

        return IRRT_OK;
    } else {
        return IRRT_NOTFOUND;
    }
}

bool
OOFEMTXTInputRecord :: hasField(InputFieldType id)
{
    //returns nonzero if id is present in source
    int indx = this->giveKeywordIndx(id);
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

const char *
OOFEMTXTInputRecord :: scanInteger(const char *source, int &value)
{
    //
    // reads integer value from source, returns pointer to char after this number
    //
    char *endptr;

    if ( source == NULL ) {
        value = 0;
        return NULL;
    }

    value = strtol(source, & endptr, 10);
    return endptr;
}

const char *
OOFEMTXTInputRecord :: scanDouble(const char *source, double &value)
{
    //
    // reads double value from source, returns pointer to char after this number
    //
    char *endptr;

    if ( source == NULL ) {
        value = 0;
        return NULL;
    }

    value = strtod(source, & endptr);
    return endptr;
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

    std :: ostringstream buff;
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
        OOFEM_WARNING( buff.str().c_str() );
    }
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
            OOFEM_WARNING("unexpected token while reading range value");
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
            OOFEM_WARNING("end ')' missing while parsing range value");
            return 0;
        }
    }

    return 0;
}

int
OOFEMTXTInputRecord :: readMatrix(const char *helpSource, int r, int c, FloatMatrix &ans)
{
    const char *endptr = helpSource;

    if ( helpSource == NULL ) {
        ans.clear();
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
        for ( int i = 1; i <= r; i++ ) {
            for ( int j = 1; j <= c; j++ ) {
                endptr = scanDouble( endptr, ans.at(i, j) );
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
                    OOFEM_WARNING("missing row terminating semicolon");
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
            return 1;
        } else {
            OOFEM_WARNING("end '}' missing while parsing matrix value");
            return 0;
        }
    } else {
        return 0;
    }

}


void
OOFEMTXTInputRecord :: report_error(const char *_class, const char *proc, InputFieldType id,
                                    IRResultType result, const char *file, int line)
{
    oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, NULL, file, line,
                              "Input error on line %d: \"%s\", field keyword \"%s\"\nIn function %s::%s\nRecord:\"%s\"",
                              lineNumber, strerror(result), id, _class, proc, this->giveRecordAsString().c_str());
    OOFEM_EXIT(1); ///@todo We should never directly exit when dealing with user input.
}
} // end namespace oofem
