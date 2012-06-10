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

#include "parser.h"
#include "error.h"
#include "mathfem.h"

#ifndef __MAKEDEPEND
 #include <cctype>
 #include <cstdlib>
 #include <cstring>
#endif

namespace oofem {
double Parser :: expr(bool get)
{
    // get indicates whether there is need to to call get_token() to get next token.

    // expression:
    //     expression + term
    //     expression - term
    //
    double left = term(get);

    for ( ; ; ) { // forever
        switch ( curr_tok ) {
        case PLUS:
            left += term(true);
            break;
        case MINUS:
            left -= term(true);
            break;
        default:
            return left;
        }
    }
}

double Parser :: term(bool get) // multiply and divide
{
    double d, left = prim(get);

    for ( ; ; ) { // forever
        switch ( curr_tok ) {
        case BOOL_EQ:
            left = ( left == prim(true) );
            break;
        case BOOL_LE:
            left = ( left <= prim(true) );
            break;
        case BOOL_LT:
            left = ( left < prim(true) );
            break;
        case BOOL_GE:
            left = ( left >= prim(true) );
            break;
        case BOOL_GT:
            left = ( left > prim(true) );
            break;
        case MUL:
            left *= prim(true);
            break;
        case DIV:
            if ( ( d = prim(true) ) ) {
                left /= d;
                break;
            }
            error("divide by 0");
            return 1;
        case POW:
            left = __OOFEM_POW(left,prim(true));
            break;
        default:
            return left;
        }
    }
}

double Parser :: prim(bool get) // handle primaries
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
        //   if (get_token() == ASSIGN) v = expr (true);
        //   return v;
        if ( get_token() == ASSIGN ) {
            name *n = insert(string_value);
            n->value = expr(true);
            return n->value;
        }

        return look(string_value)->value;
    }
    case MINUS:  // unary minus
        return -prim(true);

    case LP:
    {
        double e = expr(true);
        if ( curr_tok != RP ) {
            error(") expected");
            return 1;
        }

        get_token(); // eat ')'
        return e;
    }
    case SQRT_FUNC:
    {
        double e = agr(true);
        return sqrt(e);
    }
    case SIN_FUNC:
    {
        double e = agr(true);
        return sin(e);
    }
    case COS_FUNC:
    {
        double e = agr(true);
        return cos(e);
    }
    case TAN_FUNC:
    {
        double e = agr(true);
        return tan(e);
    }
    case ATAN_FUNC:
    {
        double e = agr(true);
        return atan(e);
    }
    case ASIN_FUNC:
    {
        double e = agr(true);
        return asin(e);
    }
    case ACOS_FUNC:
    {
        double e = agr(true);
        return acos(e);
    }
    case EXP_FUNC:
    {
        double e = agr(true);
        return exp(e);
    }
    case HEAVISIDE_FUNC://Heaviside function
    {
        double time = look("t")->value;
        double e = agr(true);

        return time<e?0:1;
    }

    default:
        error("primary expected");
        return 1;
    }
}

double Parser :: agr(bool get)
{
    if ( get ) {
        get_token();
    }

    switch ( curr_tok ) {
    case LP:
    {
        double e = expr(true);
        if ( curr_tok != RP ) {
            error(") expected");
            return 1;
        }

        get_token(); // eat ')'
        return e;
    }
    default:
        error("function argument expected");
        return 1;
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
    case '^':
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
            } else if ( !strncmp(string_value, "sin", 3) ) {
                return curr_tok = SIN_FUNC;
            } else if ( !strncmp(string_value, "cos", 3) ) {
                return curr_tok = COS_FUNC;
            } else if ( !strncmp(string_value, "tan", 3) ) {
                return curr_tok = TAN_FUNC;
            } else if ( !strncmp(string_value, "atan", 4) ) {
                return curr_tok = ATAN_FUNC;
            } else if ( !strncmp(string_value, "asin", 4) ) {
                return curr_tok = ASIN_FUNC;
            } else if ( !strncmp(string_value, "acos", 4) ) {
                return curr_tok = ACOS_FUNC;
            } else if ( !strncmp(string_value, "exp", 4) ) {
                return curr_tok = EXP_FUNC;
            } else if ( !strncmp(string_value, "h", 1) ) {
                return curr_tok = HEAVISIDE_FUNC;
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


void Parser :: error(const char *s)
{
    no_of_errors++;
    OOFEM_WARNING2("Parser :: error: %s", s);
}

double Parser :: eval(const char *string, int &err) {
    parsedLine = string;
    double result;
    no_of_errors = 0;
    do {
        result = expr(true);
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
