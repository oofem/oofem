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

#ifndef parser_h
#define parser_h

namespace oofem {
#define Parser_CMD_LENGTH 1024
#define Parser_TBLSZ 23

/**
 * Class for evaluating mathematical expressions in strings.
 * Strings should be in MATLAB syntax. The parser understands variable names with values set by "x=expression;"
 * The following binary operators are recognized: < > == + - * /
 * and the following functions: sqrt, sin, cos, tan, atan, asin, acos, exp,
 * as well as parentheses ( )
 *
 * Example string:
 * x=3;y=7;sqrt(x*(x/y+3))
 */
class Parser
{
public:
    Parser() {
        curr_tok = PRINT;
        no_of_errors = 0;
        for ( int i = 0; i < Parser_TBLSZ; i++ ) { table [ i ] = 0; } }
    ~Parser() { reset(); }

    double eval(const char *string, int &err);
    void   reset();

private:
    enum Token_value {
        NAME, NUMBER, END,
        SQRT_FUNC, SIN_FUNC, COS_FUNC, TAN_FUNC, ATAN_FUNC, ASIN_FUNC, ACOS_FUNC, EXP_FUNC, HEAVISIDE_FUNC,
        PLUS='+', MINUS='-', MUL='*', DIV='/', POW='^', BOOL_EQ, BOOL_LE, BOOL_LT, BOOL_GE, BOOL_GT,
        PRINT=';', ASSIGN='=', LP='(', RP=')'
    };

    int no_of_errors;
    Token_value curr_tok;
    struct name {
        char *string;
        name *next;
        double value;
    };
    name *table [ Parser_TBLSZ ];
    double number_value;
    char string_value [ Parser_CMD_LENGTH ];
    const char *parsedLine;

    name *look(const char *p, int ins = 0);
    inline name *insert(const char *s) { return look(s, 1); }
    void error(const char *s);
    double expr(bool get);
    double term(bool get);
    double prim(bool get);
    double agr(bool get);
    Token_value get_token();
};
} // end namespace oofem
#endif // parser_h
