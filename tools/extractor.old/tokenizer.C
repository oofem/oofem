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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

Tokenizer :: Tokenizer(FILE *inFile, char separator)
{
    inputStream = inFile;
    this->separator = separator;
    currentTokens = -1;
    isEOFflag = 0;
}

char *Tokenizer :: giveLineFromInput()
/*
 * reads one line from inputStream - for private use only.
 * sets EOF flag if EOF reached.
 * returns NULL if EOF reached and zero length string parsed.
 */
{
    int i, c;

    for ( i = 0; i < MAX_LINE_LENGTH - 1 && ( c = getc(inputStream) ) != '\n' && c != EOF; ++i ) {
        currentLine [ i ] = c;
    }

    currentLine [ i ] = '\0';

    if ( c == EOF ) {
        isEOFflag = 1;
        if ( i == 0 ) {
            currentTokens = -1;
            return NULL;
        }
    }

    this->tokenizeLine();
    return currentLine;
}


int Tokenizer :: tokenizeLine()
{
    /*
     * Parses input line, tokens are separated by seperator.
     * if separator == 0, tokens are separated by whitespaces.
     * total number of tokens parsed is return value.
     */
    int i, c, currToken = 0;
    char *ptr = currentLine;
    char *token;

    while ( currToken < MAX_TOKENS && * ptr != '\n' && * ptr != '\0' ) {
        token = tokens [ currToken ];
        if ( separator == 0 ) {
            /* skip whitespaces */
            while ( ( c = * ptr ) != '\0' && c != '\n' && ( c == ' ' || c == '\t' ) ) {
                ptr++;
            }

            for ( i = 0; i < MAX_TOKEN_LENGTH - 1 && ( c = * ptr ) != '\0' && c != '\n' && c != ' ' && c != '\t'; i++, ptr++ ) {
                * token = c;
                token++;
            }

            * token = '\0';
            token++;
            /* skip whitespaces */
            while ( ( c = * ptr ) != '\0' && c != '\n' && ( c == ' ' || c == '\t' ) ) {
                ptr++;
            }
        } else {
            for ( i = 0; i < MAX_TOKEN_LENGTH - 1 && ( c = * ptr ) != separator; i++, ptr++ ) {
                * token = c;
                token++;
            }

            * token = '\0';
            token++;
            ptr++; /* skip separator */
        }

        currToken++;
    }

    currentTokens = currToken;
    return currToken;
}

int Tokenizer :: giveNumberOfTokens()
{
    // if EOF currentTokens == -1
    return currentTokens;
}

char *Tokenizer :: giveToken(int i)
{
    // tokens are numbered from 1

    if ( i <= currentTokens ) {
        return tokens [ i - 1 ];
    } else {
        return NULL;
    }
}

int Tokenizer :: hasToken(const char *str)
{
    int i;
    for ( i = 1; i <= currentTokens; i++ ) {
        if ( !strncmp(tokens [ i - 1 ], str, strlen(str) - 1) ) {
            return i;
        }
    }

    return 0;
}



