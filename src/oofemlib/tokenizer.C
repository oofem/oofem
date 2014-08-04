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

#include "tokenizer.h"
#include "error.h"

#include <cctype>
#include <list>
#include <iterator>

namespace oofem {
Tokenizer :: Tokenizer() :
    tokens()
{ }


std :: string
Tokenizer :: readStringToken(std :: size_t &pos, const std :: string &line)
{
    pos++;
    std :: string x = this->readToken(pos, line, '"'); // read everything up to terminating '"' (or to the end of the string)
    if ( line [ pos ] == '"' ) {
        pos++;            // check if terminating '"' was found
    } else {
        OOFEM_WARNING("Missing closing separator (\") inserted at end of line");
    }
    return x;
}


std :: string
Tokenizer :: readStructToken(std :: size_t &pos, const std :: string &line)
{
    std :: string x = this->readToken(pos, line, '}'); // read everything up to terminating '}' (or to the end of the string)
    if ( line [ pos ] == '}' ) {
        pos++;            // check if terminating '}' was found
    } else {
        OOFEM_WARNING("Missing closing separator (}) inserted at end of line");
    }
    return x + '}'; // structs are left with surrounding brackets, unlike strings ""
}

std :: string
Tokenizer :: readSimpleExpressionToken(std :: size_t &pos, const std :: string &line)
{
    pos++;
    std :: string x = this->readToken(pos, line, '$'); // read everything up to terminating '$' (or to the end of the string)
    if ( line [ pos ] == '$' ) {
        pos++;            // check if terminating '"' was found
    } else {
        OOFEM_WARNING("Missing closing separator (\"$\") inserted at end of line");
    }
    return '$' + x + '$'; // simple expressions are left with surrounding '$";
}


std :: string
Tokenizer :: readSimpleToken(std :: size_t &pos, const std :: string &line)
{
    std :: size_t startpos = pos;
    while ( pos < line.size() && !isspace(line [ pos ]) ) {
        pos++;
    }
    return line.substr(startpos, pos - startpos);
}


std :: string
Tokenizer :: readToken(std :: size_t &pos, const std :: string &line, char sep)
{
    std :: size_t startpos = pos;
    while ( pos < line.size() && line [ pos ] != sep ) {
        pos++;
    }
    return line.substr(startpos, pos - startpos);
}


void Tokenizer :: tokenizeLine(const std :: string &currentLine)
{
    std :: list< std :: string >sList;
    std :: size_t bpos = 0;
    char c = 0;
    int nTokens = 0;

    while ( bpos < currentLine.size() ) {
        c = currentLine [ bpos ];

        if ( isspace(c) ) {
            bpos++;
            continue;
        } else if ( c == '"' ) {
            sList.push_back( this->readStringToken(bpos, currentLine) );
        } else if ( c == '{' ) {
            sList.push_back( this->readStructToken(bpos, currentLine) );
        } else if ( c == '$' ) {
            sList.push_back( this->readSimpleExpressionToken(bpos, currentLine) );
        } else {
            sList.push_back( this->readSimpleToken(bpos, currentLine) );
        }
    }

    // Clear the old stuff and copy list to vector
    this->tokens.clear();
    this->tokens.reserve(nTokens);
    std :: copy( sList.begin(), sList.end(), std :: back_inserter(tokens) );
}

int Tokenizer :: giveNumberOfTokens()
{
    // if EOF currentTokens == -1
    return ( int ) tokens.size();
}

const char *Tokenizer :: giveToken(int i)
{
    // tokens are numbered from 1

    if ( i <= ( int ) tokens.size() ) {
        return tokens [ i - 1 ].c_str();
    } else {
        return NULL;
    }
}
} // end namespace oofem
