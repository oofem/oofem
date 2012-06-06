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

#include "strreader.h"

#include "range.h"
#include "oofem_limits.h"

#ifndef __MAKEDEPEND
 #include <cctype>
#endif

namespace oofem {
const char *StringReader :: getPosAfter(const char *source, const char *idString)
//
// returns position of substring idString in source
// return value pointer at the end of occurrence idString in source
// (idString must be separated from rest by blank or by tabulator
//
{
    const char *str1, *helpSource = source;
    int len = strlen(idString);
    int whitespaceBefore, whitespaceAfter;

    do {
        if ( ( str1 = strstr(helpSource, idString) ) == NULL ) {
            return NULL;
        }

        helpSource = str1 + 1;
        whitespaceAfter = ( isspace( * ( helpSource + len - 1 ) ) );
        if ( str1 == source ) {
            whitespaceBefore = 1;
        } else {
            whitespaceBefore = ( isspace( * ( str1 - 1 ) ) );
        }
    } while ( !( whitespaceBefore && whitespaceAfter ) );

    return str1 + len;
}

const char *StringReader :: skipNextWord(const char *src)
//
// skips next word in src ; returns pointer after it
//
{
    while ( isspace(* src) || !* src ) {
        src++;
    }

    // skips white spaces if any
    while ( !( isspace(* src) || !* src ) ) {
        src++;
    }

    // skips one word
    return src;
}


const char *StringReader :: scanInteger(const char *source, int *value)
{
    //
    // reads integer value from source, returns pointer to char after this number
    //
    int i;

#ifdef __OOFEM_DO_NOT_PARSE_NULL
    if ( source == NULL ) {
        * value = 0;
        return NULL;
    }

#endif

    i = sscanf(source, "%d", value);
    if ( i == EOF ) {
        * value = 0;
        return NULL;
    }

    return skipNextWord(source);
}


const char *StringReader :: scanDouble(const char *source, double *value)
{
    //
    // reads integer value from source, returns pointer to char after this number
    //
    int i;

#ifdef __OOFEM_DO_NOT_PARSE_NULL
    if ( source == NULL ) {
        * value = 0;
        return NULL;
    }

#endif

    i = sscanf(source, "%lf", value);
    if ( i == EOF ) {
        * value = 0;
        return NULL;
    }

    return skipNextWord(source);
}


int StringReader :: readInteger(const char *source, const char *idString)
{
    const char *str = getPosAfter(source, idString);

#ifdef __OOFEM_DO_NOT_PARSE_NULL
    if ( str == NULL ) {
        return 0;
    }

#endif
    return atoi(str);
}

double StringReader :: readDouble(const char *source, const char *idString)
{
    const char *str = getPosAfter(source, idString);

#ifdef __OOFEM_DO_NOT_PARSE_NULL
    if ( str == NULL ) {
        return 0;
    }

#endif

    return strtod(str, NULL);
}


const char *StringReader :: readString(const char *source, const char *idString, char *string, int maxchar)
{
    const char *s;

    return readSimpleString(getPosAfter(source, idString), string, maxchar, & s);
}


const char *
StringReader :: readQuotedString(const char *source, const char *idString, char *string, int maxchar)
{
    const char *curr = getPosAfter(source, idString);
    const char *result = string;
    int len = 0;

    if ( !curr ) {
        string [ 0 ] = '\0';
        return result;
    }

    // skip whitespaces
    while ( isspace(* curr) || !* curr ) {
        curr++;
    }

    if ( !curr ) {
        fprintf(stderr, "End-of-line encountered\a\n");
        exit(1);
    }

    if ( * curr++ == '"' ) {
        while ( !( * curr == '"' ) ) {
            if ( ( * curr == '\n' ) || !* curr ) {
                fprintf(stderr, "readQuotedString: final \" expected\a\n");
                exit(1);
            }

            if ( ++len == ( maxchar ) ) {
                * string = 0;
                break;
            }

            * string++ = * curr++;
        }

        * string = 0;
    } else {
        fprintf(stderr, "readQuotedString: expected initial  \"\a\n");
        exit(1);
    }

    return result;
}


IntArray *StringReader :: ReadIntArray(const char *source, const char *idString)
{
    const char *str1 = source;
    int value, size;
    IntArray *arry;


    str1 = getPosAfter(source, idString);
    str1 = scanInteger(str1, & size);
    //  if(size ==0) return NULL;   // if commented returns new array of size = 0 for size = 0
    if ( str1 == NULL ) {
        size = 0;
    }

    arry = new IntArray(size);
    for ( int i = 1; i <= size; i++ ) {
        str1 = scanInteger(str1, & value);
        arry->at(i) = value;
    }

    return arry;
}


FloatArray *StringReader :: ReadFloatArray(const char *source, const char *idString)
{
    const char *str1 = source;
    double value;
    int size;
    FloatArray *arry;

    str1 = getPosAfter(source, idString);
    str1 = scanInteger(str1, & size);
    //   if(size ==0) return NULL;  // if commented returns new array of size 0 for size = 0
    if ( str1 == NULL ) {
        size = 0;
    }

    arry = new FloatArray(size);
    for ( int i = 1; i <= size; i++ ) {
        str1 = scanDouble(str1, & value);
        arry->at(i) = value;
    }

    return arry;
}



Dictionary *StringReader :: ReadDictionary(const char *source, const char *idString)
{
    const char *str1 = source;
    double value;
    int size;
    char key [ 100 + 1 ]; // 'key' is eventually of size 1, but some words that are
                                      // read in the meantime in the data file can be larger !
    Dictionary *dict;

    dict = new Dictionary();
    str1 = getPosAfter(source, idString); // move after idString
    str1 = scanInteger(str1, & size);    // read number of conditions
    if ( size == 0 ) {
        return dict;                     // no - conditions - return void dictionary
    }

    if ( str1 == NULL ) {
        return dict;                     // end of line - return void dictionary
    }

    for ( int i = 1; i <= size; i++ ) {
        readSimpleString(str1, key, 100 + 1, & str1);
        str1 = scanDouble(str1, & value);
        dict->add(key [ 0 ], value);
    }

    return dict;
}



char *StringReader :: readSimpleString(const char *source, char *simpleString, int maxchar, const char **remain)
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
        fprintf(stderr, "End-of-line encountered\a\n");
        exit(1);
    }

    while ( ( !( isspace(* curr) || !* curr ) ) && ( ++count < maxchar ) ) {
        * ss++ = * curr++;
    }

    * ss = '\0';
    * remain = curr;
    return simpleString;
}


const char *StringReader :: readKeyAndVal(const char *source, char *key, int *val, int maxchar, const char **remain)
{
    key = readSimpleString(source, key, maxchar, remain);
    * remain = scanInteger(* remain, val);
    return * remain;
}


const char *StringReader :: readKeyAndVal(const char *source, char *key, double *val, int maxchar, const char **remain)
{
    key = readSimpleString(source, key, maxchar, remain);
    * remain = scanDouble(* remain, val);
    return * remain;
}


bool StringReader :: hasString(const char *source, const char *idString)
{
    //returns nonzero if idString is present in source
    const char *str = strstr(source, idString);
    return str != NULL;
}


void StringReader :: readRangeList(dynaList< Range > &list, const char *source, const char *idString)
{
    int li, hi;
    const char *str1, *helpSource = source;
    // Range* range;

    // find first valid occurrence of idString
    int len = strlen(idString);
    do {
        if ( ( str1 = strstr(helpSource, idString) ) == NULL ) {
            return;
        }

        helpSource = str1 + 1;
    } while ( !( isspace( * ( helpSource + len - 1 ) ) ) );


    helpSource = str1 + len;
    // find first non whitespace character
    // skip whitespaces
    while ( isspace(* helpSource) ) {
        helpSource++;
    }

    // test if list left bracketed found
    if ( * helpSource != '{' ) {
        OOFEM_WARNING("StringReader::readRangeList: parse error - missing left '{'");
        list.clear();
        return;
    }

    helpSource++;
    // read ranges
    while ( readRange(& helpSource, li, hi) ) {
        Range range(li, hi);
        list.pushBack(range);
    }

    // skip whitespaces after last range
    while ( isspace(* helpSource) ) {
        helpSource++;
    }

    // test for enclosing bracket
    if ( * helpSource != '}' ) {
        OOFEM_WARNING("StringReader::readRangeList: parse error - missing end '}'");
        list.clear();
        return;
    }
}


int StringReader :: readRange(const char **helpSource, int &li, int &hi)
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
        if ( isspace(* * helpSource) ) {
            OOFEM_WARNING("StringReader::readRange: unexpected token while reading range value");
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
            OOFEM_WARNING("StringReader::readRange: end ')' missing while parsing range value");
            return 0;
        }
    }

    return 0;
}
} // end namespace oofem
