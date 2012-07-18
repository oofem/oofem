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
#ifndef tokenizer_h
#define tokenizer_h

#include "oofem_limits.h"

namespace oofem {
/**
 * Tokenizer class. This class splits given record (represented as string) to
 * particular tokens, which are usually separated by white spaces (but other separators
 * are also supported). Tokenizer recognizes quoted strings and structured tokens that are
 * bounded by '{' '}' pairs, can be nested and represent single token.
 *
 * The implementation uses continuous token buffer, from which space for particular tokens is allocated
 * and initial position of each token is kept in array to provide quick access.
 */
class Tokenizer
{
private:
    /// Common token buffer, tokens are terminated by '\0'
    char tokens [ OOFEM_MAX_TOKENS_LENGTH ];
    /// Array of pointers to token buffer. The i-th pointer points to the position of i-th token in token buffer.
    char *tokenPosition [ OOFEM_MAX_TOKENS ];
    /// Token separator; zero represents whitespace chars
    char separator;
    /// Number of tokens.
    int nTokens;

public:
    /// Constructor. Creates tokenizer with given character as separator.
    Tokenizer(char separator = 0);
    /// Tokenizes given record (string).
    int tokenizeLine(const char *line);
    /// returns the number of tokens.
    int giveNumberOfTokens();
    /// Returns pointer to i-th token.
    const char *giveToken(int);

protected:
    /**
     * Reads next token.
     * @param bpos Position (index) in token buffer.
     * @param line Pointer pointer to record from which token is parsed.
     * @param token Pointer to next free token buffer position.
     * @param sep Separator.
     */
    void readToken(int &bpos, const char * &line, char * &token, char sep);
    /**
     * Reads next structured token (bounded by '{' '}' pairs, possibly nested).
     * @param bpos Position (index) in token buffer.
     * @param line Pointer pointer to record from which token is parsed.
     * @param token Pointer to next free token buffer position.
     */
    void readStructToken(int &bpos, const char * &line, char * &token);
    /**
     * Reads next string token (quoted).
     * @param bpos Position (index) in token buffer.
     * @param line Pointer pointer to record from which token is parsed.
     * @param token Pointer to next free token buffer position.
     */
    void readStringToken(int &bpos, const char * &line, char * &token);
};
} // end namespace oofem
#endif // tokenizer_h
