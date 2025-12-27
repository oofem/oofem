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

#ifndef tokenizer_h
#define tokenizer_h

#include <cstdio>
#include <string.h>

#define MAX_LINE_LENGTH 1024
#define MAX_TOKEN_LENGTH 80
#define MAX_TOKENS 20

typedef char tokenType [ MAX_TOKENS ] [ MAX_TOKEN_LENGTH ];

char *giveLineFromInput(char *line);
int tokenizeLine(char separator, char *line, tokenType *tokens);


class Tokenizer
{
private:
    tokenType tokens;
    FILE *inputStream;
    char separator;

    int currentTokens, isEOFflag;
    char currentLine [ MAX_LINE_LENGTH ];
    int tokenizeLine();

public:
    Tokenizer(FILE * inFile, char separator = 0);
    char *giveLineFromInput();
    int   giveNumberOfTokens();
    int   isEOF() { return isEOFflag; }
    char *giveToken(int);
    char *giveCurrentLine() { return currentLine; }
    // return token number if token present; zero otherwise
    int   hasToken(const char *str);
};


#endif // tokenizer_h
