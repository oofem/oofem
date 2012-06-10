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

#ifndef strreader_h
#define strreader_h

#include "intarray.h"
#include "flotarry.h"
#include "dictionr.h"
#include "dynalist.h"

namespace oofem {
class Range;

/**
 * This class implements simple string reader.
 * Reads number (integer or double ) or array (integer or double) depicted by string identifier.
 * The domain class converts all record characters to lower case.
 * But the StringReader is case sensitive.
 */
class StringReader
{
public:
    StringReader() { }
    ~StringReader() { }

    int readInteger(const char *source, const char *idString);
    double readDouble(const char *source, const char *idString);
    const char *readString(const char *source, const char *idString, char *string, int maxchar);
    /**
     * Reads quoted string identified by keyword.
     * The read string must be bounded by " characters.
     * @param source Source record.
     * @param idString Keyword identifying the value.
     * @param string Buffer to store result. Should be at least of maxchar size.
     * @param maxchar Maximum character to be read.
     * @return Pointer to string parameter.
     */
    const char *readQuotedString(const char *source, const char *idString, char *string, int maxchar);
    IntArray *ReadIntArray(const char *source, const char *idString);
    FloatArray *ReadFloatArray(const char *source, const char *idString);
    Dictionary *ReadDictionary(const char *source, const char *idString);
    char *readSimpleString(const char *source, char *simpleString, int maxchar, const char **remain);
    const char *readKeyAndVal(const char *source, char *key, int *val, int maxchar, const char **remain);
    const char *readKeyAndVal(const char *source, char *key, double *val, int maxchar, const char **remain);
    /**
     * Checks for keywords in strings.
     * @param idString Keyword identifying the value.
     * @param source Source record.
     * @return Nonzero if idString is present in source.
     */
    bool hasString(const char *source, const char *idString);

    /**
     * Reads range list. The range syntax is "range_list_name { number, (start end) }",
     * where range_list_name is name of list in source record. It is necessary to enclose list definition into {}.
     * The single number is shortcut for one value range. The range is specified using two numbers in parenthesis.
     * @param list List to be read. Actually items are added into this list.
     * @param source Source string with corresponding record.
     * @param idString String containing range_list_name.
     */
    void readRangeList(dynaList< Range > &list, const char *source, const char *idString);
    /**
     * Reads single range record from input record represented by helpSource string.
     * @param helpSource Pointer to current string position, on return helpSource points
     * to next character after reading range record.
     * @param li Starting range index.
     * @param hi End range index.
     * @return Nonzero on success.
     */
    int readRange(const char **helpSource, int &li, int &hi);

private:
    //
    // I don't expect, that You will use following functions
    //

    const char *getPosAfter(const char *, const char *);
    const char *scanInteger(const char *source, int *value);
    const char *scanDouble(const char *source, double *value);
    const char *skipNextWord(const char *src);
};
} // end namespace oofem
#endif // strreader_h
