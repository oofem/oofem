/* $Header: /home/cvs/bp/oofem/oofemlib/src/parser.h,v 1.5 2003/04/06 14:08:25 bp Exp $ */
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

#ifndef parser_h
#define parser_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <math.h>
#endif

#define Parser_CMD_LENGTH 1024
#define Parser_TBLSZ 23

class Parser
{
public:
    Parser() { curr_tok = PRINT;
               no_of_errors = 0;
               for ( int i = 0; i < Parser_TBLSZ; i++ ) { table [ i ] = 0; } }
    ~Parser() { reset(); }

    double eval(const char *string, int &err);
    void   reset();

private:
    enum Token_value {
        NAME, NUMBER, END,
        SQRT_FUNC, SIN_FUNC, COS_FUNC, TAN_FUNC, ATAN_FUNC, ASIN_FUNC, ACOS_FUNC, EXP_FUNC,
        PLUS='+', MINUS='-', MUL='*', DIV='/', BOOL_EQ, BOOL_LE, BOOL_LT, BOOL_GE, BOOL_GT,
        PRINT=';', ASSIGN='=', LP='(', RP=')'
    };

    int no_of_errors;
    Token_value curr_tok;
    struct name { char *string;
                  name *next;
                  double value; };
    name *table [ Parser_TBLSZ ];
    double number_value;
    char string_value [ Parser_CMD_LENGTH ];
    const char *parsedLine;

    name *look(const char *p, int ins = 0);
    inline name *insert(const char *s) { return look(s, 1); }
    double error(const char *s);
    double expr(int get);
    double term(int get);
    double prim(int get);
    double agr(int get);
    Token_value get_token();
};


#endif // parser_h
