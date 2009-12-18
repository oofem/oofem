/* $Header: /home/cvs/bp/oofem/oofemlib/src/parser.C,v 1.6.4.1 2004/04/05 15:19:43 bp Exp $ */
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


#include "parser.h"
#include "error.h"
#ifndef __MAKEDEPEND
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#endif
#include "compiler.h"

namespace oofem {

double Parser :: expr(int get)
{
    // get indicates whether there is need to to call get_token() to get next token.

    // epression:
    //     expression + term
    //     expression - term
    //
    double left = term(get);

    for ( ; ; ) { // forever
        switch ( curr_tok ) {
        case PLUS:
            left += term(TRUE);
            break;
        case MINUS:
            left -= term(TRUE);
            break;
        default:
            return left;
        }
    }
}

double Parser :: term(int get) // multiply and divide
{
    double d, left = prim(get);

    for ( ; ; ) { // forever
        switch ( curr_tok ) {
        case BOOL_EQ:
            left = ( left == prim(TRUE) );
            break;
        case BOOL_LE:
            left = ( left <= prim(TRUE) );
            break;
        case BOOL_LT:
            left = ( left < prim(TRUE) );
            break;
        case BOOL_GE:
            left = ( left >= prim(TRUE) );
            break;
        case BOOL_GT:
            left = ( left > prim(TRUE) );
            break;
        case MUL:
            left *= prim(TRUE);
            break;
        case DIV:
            if ( ( d = prim(TRUE) ) ) {
                left /= d;
                break;
            }

            return error("divide by 0");

        default:
            return left;
        }
    }
}



double Parser :: prim(int get) // handle primaries
{
    if ( get ) {
        get_token();
    }

    switch ( curr_tok ) {
    case NUMBER:
    {
        double v = number_value;
        get_token();
        return v;
    }
    case NAME:
    {
        //   double &v = table[string_value];
        //   if (get_token() == ASSIGN) v = expr (TRUE);
        //   return v;
        if ( get_token() == ASSIGN ) {
            name *n = insert(string_value);
            n->value = expr(TRUE);
            return n->value;
        }

        return look(string_value)->value;
    }
    case MINUS:  // unary minus
        return -prim(TRUE);

    case LP:
    {
        double e = expr(TRUE);
        if ( curr_tok != RP ) {
            return error(") expected");
        }

        get_token(); // eat ')'
        return e;
    }
    case SQRT_FUNC:
    {
        double e = agr(TRUE);
        return sqrt(e);
    }
    case SIN_FUNC:
    {
        double e = agr(TRUE);
        return sin(e);
    }
    case COS_FUNC:
    {
        double e = agr(TRUE);
        return cos(e);
    }
    case TAN_FUNC:
    {
        double e = agr(TRUE);
        return tan(e);
    }
    case ATAN_FUNC:
    {
        double e = agr(TRUE);
        return atan(e);
    }
    case ASIN_FUNC:
    {
        double e = agr(TRUE);
        return asin(e);
    }
    case ACOS_FUNC:
    {
        double e = agr(TRUE);
        return acos(e);
    }
    case EXP_FUNC:
    {
        double e = agr(TRUE);
        return exp(e);
    }
    default:
        return error("primary expected");
    }
}

double Parser :: agr(int get)
{
    if ( get ) {
        get_token();
    }

    switch ( curr_tok ) {
    case LP:
    {
        double e = expr(TRUE);
        if ( curr_tok != RP ) {
            return error(") expected");
        }

        get_token(); // eat ')'
        return e;
    }
    default:
        return error("function argument expected");
    }
}


Parser :: Token_value Parser :: get_token()
{
    char ch = 0;
    int len;

    do { // skip whitespaces except '\n'
        //  if (! input->get(ch)) return curr_tok = END;
        if ( !( ch = * ( parsedLine++ ) ) ) {
            return curr_tok = END;
        }
    } while ( ch != '\n' && isspace(ch) );

    switch ( ch ) {
    case 0:
    case '\n':
        return curr_tok = END;

    case ';':
        return curr_tok = PRINT;

    case '*':
    case '/':
    case '+':
    case '-':
    case '(':
    case ')':
        return curr_tok = Token_value(ch);

    case '=':
        if ( ( ch = * ( parsedLine++ ) ) == '=' ) {
            return curr_tok = BOOL_EQ;
        } else {
            parsedLine--;
            return curr_tok = ASSIGN;
        }

    case '<':
        if ( ( ch = * ( parsedLine++ ) ) == '=' ) {
            return curr_tok = BOOL_LE;
        } else {
            parsedLine--;
            return curr_tok = BOOL_LT;
        }

    case '>':
        if ( ( ch = * ( parsedLine++ ) ) == '=' ) {
            return curr_tok = BOOL_GE;
        } else {
            parsedLine--;
            return curr_tok = BOOL_GT;
        }

    case '0': case '1': case '2': case '3': case '4': case '5':
    case '6': case '7': case '8': case '9': case '.':
        //input->putback(ch);
        parsedLine--;
        //*input >> number_value;
        char *endParse;
        number_value = strtod(parsedLine, & endParse);
        parsedLine = endParse;

        return curr_tok = NUMBER;

    default:
        if ( isalpha(ch) ) {
            //   string_value = ch;
            //   while (input->get(ch) && isalnum (ch)) string_value += ch;
            //   input->putback (ch);
            char *p = string_value;
            * p++ = ch;
            len = 1;
            //   while (input->get(ch) && isalnum (ch)) *p++ = ch;
            while ( ( ch = * ( parsedLine++ ) ) && isalnum(ch) ) {
                * p++ = ch;
                if ( len++ >= Parser_CMD_LENGTH ) {
                    error("command too long");
                }
            }

            * p = 0;
            //   input->putback(ch);
            parsedLine--;

            if ( !strncmp(string_value, "sqrt", 4) ) {
                return curr_tok = SQRT_FUNC;
            } else if ( !strncmp(string_value, "sin", 3) )  {
                return curr_tok = SIN_FUNC;
            } else if ( !strncmp(string_value, "cos", 3) )  {
                return curr_tok = COS_FUNC;
            } else if ( !strncmp(string_value, "tan", 3) )  {
                return curr_tok = TAN_FUNC;
            } else if ( !strncmp(string_value, "atan", 4) )  {
                return curr_tok = ATAN_FUNC;
            } else if ( !strncmp(string_value, "asin", 4) )  {
                return curr_tok = ASIN_FUNC;
            } else if ( !strncmp(string_value, "acos", 4) )  {
                return curr_tok = ACOS_FUNC;
            } else if ( !strncmp(string_value, "exp", 4) )  {
                return curr_tok = EXP_FUNC;
            } else {
                return curr_tok = NAME;
            }
        }

        error("bad token");
        return curr_tok = PRINT;
    }
}


Parser :: name *Parser :: look(const char *p, int ins) {
    int ii = 0;                               // hash
    const char *pp = p;
    while ( * pp ) {
        ii = ii << 1 ^ * pp++;
    }

    if ( ii < 0 ) {
        ii = -ii;
    }

    ii %= Parser_TBLSZ;

    for ( name *n = table [ ii ]; n; n = n->next ) { // search
        if ( strcmp(p, n->string) == 0 ) {
            return n;
        }
    }

    if ( ins == 0 ) {
        error("name not found");
    }

    name *nn = new name;
    nn->string = new char [ strlen(p) + 1 ];
    strcpy(nn->string, p);
    nn->value = 0.;
    nn->next = table [ ii ];
    table [ ii ] = nn;
    return nn;
}


double Parser :: error(const char *s)
{
    no_of_errors++;
    OOFEM_WARNING2("Parser: erorr: %s", s);
    return 1;
}

double Parser :: eval(const char *string, int &err) {
    parsedLine = string;
    double result;
    no_of_errors = 0;
    do {
        result = expr(TRUE);
    } while ( curr_tok != END );

    err = no_of_errors;
    return result;
}

void Parser :: reset() {
    // empty Parser table
    name *entry, *next;
    for ( int i = 0; i < Parser_TBLSZ; i++ ) {
        if ( ( entry = table [ i ] ) ) {
            do {
                next = entry->next;
                delete entry;
            } while ( ( entry = next ) );

            table [ i ] = 0;
        }
    }
}

} // end namespace oofem
