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

#include "tokenizer.h"
#include "error.h"
#ifndef __MAKEDEPEND
 #include <ctype.h>
#endif

namespace oofem {
Tokenizer :: Tokenizer(char separator)
{
    this->separator = separator;
    nTokens = -1;
}

void
Tokenizer :: readStringToken(int &bpos, const char * &line, char * &token)
{
    ++line; // skip leading '"'
    this->readToken(bpos, line, token, '"'); // read everything up to terminating '"'
    if ( * line == '"' ) {
        line++;            // check if terminating '"' was found
    } else {
        OOFEM_WARNING("Tokenizer::readStringToken : Missing closing separator (\")");
    }
}

void
Tokenizer :: readStructToken(int &bpos, const char * &line, char * &token)
{
    int refLevel = 0;
    char c;

    do {
        c = * line++;
        if ( bpos < OOFEM_MAX_TOKENS_LENGTH - 1 ) {
            * token = c;
            token++;
            bpos++;
        }

        if ( c == '{' ) {
            refLevel++;
        } else if ( c == '}' ) {
            refLevel--;
        }
    } while ( refLevel && c != '\0' && c != '\n' );

    if ( refLevel ) {
        OOFEM_WARNING("Tokenizer::readStructToken : Missing closing separator (})");
    }
}


void
Tokenizer :: readToken(int &bpos, const char * &line, char * &token, char sep)
{
    char c;
    if ( sep == 0 ) {
        for ( ; ( c = * line ) != '\0' && !isspace(c); line++ ) {
            if ( bpos < OOFEM_MAX_TOKENS_LENGTH - 1 ) {
                * token = c;
                token++;
                bpos++;
            }
        }
    } else {
        for ( ; ( c = * line ) != '\0' && c != '\n' && c != sep; line++ ) {
            if ( bpos < OOFEM_MAX_TOKENS_LENGTH - 1 ) {
                * token = c;
                token++;
                bpos++;
            }
        }
    }
}


int Tokenizer :: tokenizeLine(const char *currentLine)
{
    /*
     *       Parses input line, tokens are separated by separator.
     *       if separator == 0, tokens are separated by whitespaces.
     *       return code:
     * 0 - ok
     * 1 - too many tokens
     */
    bool overflow = false;
    int c = 0;
    const char *ptr = currentLine;
    char *token;

    int bpos = 0;
    nTokens = 0;

    while ( nTokens < OOFEM_MAX_TOKENS && * ptr != '\n' && * ptr != '\0' ) {
        tokenPosition [ nTokens ] = token = tokens + bpos;
        if ( separator == 0 ) {
            /* skip whitespaces */
            while ( ( c = * ptr ) != '\0' && c != '\n' && isspace(c) ) {
                ptr++;
            }

            if ( c == '"' ) {
                readStringToken(bpos, ptr, token);
            } else if ( c == '{' ) {
                readStructToken(bpos, ptr, token);
            } else {
                this->readToken(bpos, ptr, token, 0);
            }
        } else {
            this->readToken(bpos, ptr, token, separator);
            if ( c == separator ) {
                ptr++;                /* skip separator */
            }
        }

        if ( token - tokenPosition [ nTokens ] ) {  // if some chars were valid
            * token = '\0';
            if ( bpos < OOFEM_MAX_TOKENS_LENGTH - 1 ) {
                token++;
                bpos++;
                nTokens++;
            } else {
                overflow = true;
            }
        }
    }

    if ( overflow ) {
        OOFEM_WARNING("Tokenizer :: tokenizeLine: overflow detected, increase token buffer");
    }

    if ( nTokens >= OOFEM_MAX_TOKENS ) {
        OOFEM_WARNING("Tokenizer :: tokenizeLine: overflow detected, increase number of tokens");
    }

    if ( * ptr == '\n' || * ptr == '\0' ) {
        return 0;
    } else {
        return 1;
    }
}

int Tokenizer :: giveNumberOfTokens()
{
    // if EOF currentTokens == -1
    return nTokens;
}

const char *Tokenizer :: giveToken(int i)
{
    // tokens are numbered from 1

    if ( i <= nTokens ) {
        return tokenPosition [ i - 1 ];
    } else {
        return NULL;
    }
}
} // end namespace oofem
