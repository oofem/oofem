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

/*
 * This code is part of OOFEM source tree.
 * Copyright (c) Borek Patzak, 1999.
 */



#include "seek.h"
#if defined ( __NetBSD__ ) || defined ( __FreeBSD__ )
 #include <libgen.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
// #define SILICON_GRAPHICS

/*
 * CHECKER_MODE ... reads check record description in input file and
 * compares exact results in check record with computed results
 * extracted from output file. prints check records which fails as weel as
 * final result.
 *
 * SYSYTEM_STUDENT_EXTRACTOR_MODE ...reads check record description in input file
 * and prints corresponding extracted values using post_expr format
 * to stdout. Used by system student (Homeworks on Internet project).
 *
 * EXTRACTOR_MODE ... works like grMaker. The grMaker will not be further maintained.
 * In extractor mode the listed records are extracted for each time step,
 * and printed into rows (one row per solution step) and in coolumns the required record values
 * for corresponding step are printed in order in which they heve been specified.
 * The new EXTRACTOR_MODE specific keyword #TIME is added and it will cause the current time step time
 * beiing printed.
 */

//#define COLOR_MODE
#define DEFAULT_CHECK_TOLERANCE 0.0001
#define MAX_CHECK_ITEMS 80
#define MAX_NAME_LENGTH 20
#define MAX_COMMENT_LENGTH 80
#define MAX_EXPR_LENGTH 80
#define MAX_FRMT_LENGTH 120
#define MAX_KEYWORD_LENGTH 20

#ifdef COLOR_MODE
 #define  NORMAL "\033[0;39m"
 #define  GREEN  "\033[1;32m"
 #define  RED    "\033[1;31m"
#else
 #define  NORMAL ""
 #define  GREEN  ""
 #define  RED    ""
#endif

//extern stateType currState;

//enum typeOfRequest {Eigval = 1, LoadLevel = 0, NodeSideVal = 2, ElementVal  = 3, BeamElementVal = 4, ReactionForce = 5};
enum typeOfRequest {
#ifdef EXTRACTOR_MODE
    TimeVal = 0, // introduced to extract current time step value
#endif
    LoadLevel = 1,
    NumberOfIter = 2,
    QuasiReaction = 3,
    NodeSideVal = 4,
    ElementVal  = 5,
    BeamElementVal = 6,
    ReactionForce = 7,
    Eigval = 8,
    EigNodeSideVal = 9,
    EigElementVal=10,
    EigBeamElementVal=11,
    SolutionStepUserTime=12
};
//
// Note if user wants to extract eigen value node or element values, then separate keyworks like
// Eig_NodeSideVal must be introduced and should be appended after Eigval typeOfRequest value.
//
struct checkRequestType {
    double timeStep;
    int eigValNum;
    typeOfRequest type;
    int elNumber, dofNumber, compNumber, stressFlag;
    int gpNum;
    int recNumber;
    double tolerance;
    char keyword [ MAX_KEYWORD_LENGTH + 1 ];
#ifdef CHECKER_MODE
    double checkValue;
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
    char value_name [ MAX_NAME_LENGTH + 1 ];
    char comment   [ MAX_COMMENT_LENGTH + 1 ];
    /* post_expr is in form of format of printf like function
     * like "-1.*%lf", where there are following restrictions:
     * - using this format extracted value is printed out - %lf should be used;
     * - this expression will be evaluated by parser program - the expressin must be inside ` `;
     * - if extracted value must be used more than once times  in expression, an assingment to
     * temporary variable can be used - for example `temp=%lf;temp*temp`. The parser evaluates
     * all expressions (must be seperated by ;), the value of multiple expression returned
     * is the result from last parsed expression. See parser documentation for details;
     *
     * if post_exp is missing (expr keyword not found) then extracted value is directly printed out.
     */
    char post_expr [ MAX_EXPR_LENGTH + 1 ];
    int flag;   // indicates, if student system will not check value or print value in teacher mode
    double extractedValue;
#endif
#ifdef EXTRACTOR_MODE
    double extractedValue;
#endif
    // int checked;
};

checkRequestType checkRequestArry [ MAX_CHECK_ITEMS ];
int checkRequestOrderArry [ MAX_CHECK_ITEMS ];
int nChecks;
double tolerance;

void printHelp();
//int popNextCheckRecord () ;
int checkValue(double firstVal, double secondVal, double tolerance);
void printCheckError(const char *str);
void printCheckErrori1(const char *str, int val);
void printCheckErrors1(const char *str, const char *str2);
void printCheckErrors2(const char *str, const char *str2, const char *str3);
void printCheckRuleError(int val);
char *getPosAfter(char *source, const char *idString);
int readInteger(char *source, const char *idString, int compulsory = 1);
double readDouble(char *source, const char *idString, int compulsory = 1);
char readChar(char *source, const char *idString);
void readString(char *string, char *source, const char *idString, int maxLen, int compulsory = 1);
int hasString(char *source, const char *idString);
void printStudentSystemExtractedVariable(double checkVal, char *name, char *expr, double tolerance, int flag, char *comment);
void printStudentSystemExtractorError(int indx);
void printExtractorError(int indx);
int compareRecords(int i, int j);
void sortCheckRecords();
void quickSortArry(int *map, int l, int r);
int quickSortPartition(int *map, int l, int r);

char inputStreamName [ 80 ], outFileName [ 80 ];


int main(int argc, char *argv[]) {
    int result, i, indx, flag = 0;
    double checkVal, dummy;
    FILE *inputStream, *outputFileStream;
    Tokenizer *it, *outFileTokenizer;
#ifdef EXTRACTOR_MODE
    stateType currState;
#endif

    // parse command line arguments
    if ( argc == 1 ) {
        printHelp();
        exit(1);
    }

    for ( i = 1; i < argc; i++ ) {
        if ( !strncmp(argv [ i ], "-f", 2) && ( i + 1 < argc ) ) {
            strcpy(inputStreamName, argv [ i + 1 ]);
            i++;
            flag = 1;
        } else if ( !strncmp(argv [ i ], "-help", 5) ) {
            printHelp();
        } else {
            printHelp();
        }
    }

    // check if input file parameter has been read
    if ( flag == 0 ) {
        printf("\n\n");
        printCheckError("no input file given, exiting");
        exit(1);
    }

    // open file
    if ( ( inputStream = fopen(inputStreamName, "r") ) != NULL ) {
        it = new Tokenizer(inputStream);
    } else {
        printCheckError("can't open input file, exiting");
        exit(1);
    }

    // read Output file, where results are stored and check section.
    int ReadOutFileFlag = 0, ReadCheckSection = 0;
    it->giveLineFromInput();
    while ( !( ( it->isEOF() ) || ( ReadOutFileFlag != 0 && ReadCheckSection != 0 ) ) ) {
        if ( !ReadOutFileFlag && it->giveToken(1) [ 0 ] != '#' ) {
            // first uncommented line contains outputFileName
            strcpy( outFileName, it->giveToken(1) );
            ReadOutFileFlag = 1;
        }

        if ( it->giveNumberOfTokens() >= 1 && !strncmp(it->giveToken(1), "#%BEGIN_CHECK%", 14) ) {
            // check section found
            // first  check if default tolerance is specified
            if ( getPosAfter(it->giveCurrentLine(), "tolerance") != NULL ) {
                tolerance = readDouble(it->giveCurrentLine(), "tolerance");
            } else {
                tolerance = DEFAULT_CHECK_TOLERANCE; // if no use default value
            }

            nChecks = 0;

            do {
                it->giveLineFromInput();
                if ( it->giveToken(1) [ 1 ] == '#' ) {
                    continue;
                }

                if ( !strncmp(it->giveToken(1), "#%END_CHECK%", 12) ) {
                    ReadCheckSection = 1;
                    break;
                } else if ( !strncmp(it->giveToken(1), "#NODE", 5) ||
                           !strncmp(it->giveToken(1), "#SIDE", 5) ||
                           !strncmp(it->giveToken(1), "#EIG_NODE", 9) ||
                           !strncmp(it->giveToken(1), "#EIG_SIDE", 9) ) {
                    /*     if (it->giveNumberOfTokens() != 6) {
                     * printCheckErrors1 ("Missing arguments on line: %s",it->giveCurrentLine());
                     * exit(1);
                     * } */
                    if ( !strncmp(it->giveToken(1), "#EIG_NODE", 9) || !strncmp(it->giveToken(1), "#EIG_SIDE", 9) ) {
                        checkRequestArry [ nChecks ].type       = EigNodeSideVal;
                        checkRequestArry [ nChecks ].eigValNum   = readInteger(it->giveCurrentLine(), "EigNum"); // eig Val Num
                    } else {
                        checkRequestArry [ nChecks ].type       = NodeSideVal;
                        checkRequestArry [ nChecks ].eigValNum   = 0;
                    }

#ifndef EXTRACTOR_MODE
                    // in EXTRACTOR_MODE tstep value has no meaning
                    checkRequestArry [ nChecks ].timeStep   = readDouble(it->giveCurrentLine(), "tStep"); // time step
#endif
                    checkRequestArry [ nChecks ].elNumber   = readInteger(it->giveCurrentLine(), "number"); // nodeSide number
                    checkRequestArry [ nChecks ].dofNumber  = readInteger(it->giveCurrentLine(), "dof"); // dof number
                    checkRequestArry [ nChecks ].compNumber = readChar(it->giveCurrentLine(), "unknown"); // unknown type
                    if ( hasString(it->giveCurrentLine(), "tolerance") ) {
                        checkRequestArry [ nChecks ].tolerance = readDouble(it->giveCurrentLine(), "tolerance"); // check tolerance
                    } else {
                        checkRequestArry [ nChecks ].tolerance = tolerance;
                    }

#ifdef CHECKER_MODE
                    checkRequestArry [ nChecks ].checkValue = readDouble(it->giveCurrentLine(), "value"); // check Value
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
                    // check value name
                    readString(checkRequestArry [ nChecks ].value_name, it->giveCurrentLine(), "name", MAX_NAME_LENGTH);
                    readString(checkRequestArry [ nChecks ].comment, it->giveCurrentLine(), "comment", MAX_COMMENT_LENGTH, 0);
                    readString(checkRequestArry [ nChecks ].post_expr, it->giveCurrentLine(), "expr", MAX_EXPR_LENGTH, 0);
                    checkRequestArry [ nChecks ].flag = readInteger(it->giveCurrentLine(), "noCheckFlag", 0);
#endif
                    //checkRequestArry[nChecks].checked    = 0;
                    nChecks++;
                } else if ( !strncmp(it->giveToken(1), "#ELEMENT", 8) ||
                           !strncmp(it->giveToken(1), "#EIG_ELEMENT", 12) ) {
                    /* if (it->giveNumberOfTokens() != 7) {
                     * printCheckErrors1 ("Missing arguments on line: %s",it->giveCurrentLine());
                     * exit(1);
                     * }*/
                    if ( !strncmp(it->giveToken(1), "#EIG_ELEMENT", 12) ) {
                        checkRequestArry [ nChecks ].type       = EigElementVal;
                        checkRequestArry [ nChecks ].eigValNum   = readInteger(it->giveCurrentLine(), "EigNum"); // eig Val Num
                    } else {
                        checkRequestArry [ nChecks ].type       = ElementVal;
                        checkRequestArry [ nChecks ].eigValNum   = 0;
                    }

#ifndef EXTRACTOR_MODE
                    // in EXTRACTOR_MODE tstep value has no meaning
                    checkRequestArry [ nChecks ].timeStep   = readDouble(it->giveCurrentLine(), "tStep");
#endif
                    checkRequestArry [ nChecks ].elNumber   = readInteger(it->giveCurrentLine(), "number");
                    checkRequestArry [ nChecks ].gpNum      = readInteger(it->giveCurrentLine(), "gp"); // gp number
                    checkRequestArry [ nChecks ].recNumber  = readInteger(it->giveCurrentLine(), "record"); // element record num.
                    checkRequestArry [ nChecks ].compNumber = readInteger(it->giveCurrentLine(), "component"); // component number
                    readString(checkRequestArry [ nChecks ].keyword, it->giveCurrentLine(), "keyword", MAX_KEYWORD_LENGTH, 0); // read optional keyword
                    if ( hasString(it->giveCurrentLine(), "tolerance") ) {
                        checkRequestArry [ nChecks ].tolerance = readDouble(it->giveCurrentLine(), "tolerance"); // check tolerance
                    } else {
                        checkRequestArry [ nChecks ].tolerance = tolerance;
                    }

#ifdef CHECKER_MODE
                    checkRequestArry [ nChecks ].checkValue = readDouble(it->giveCurrentLine(), "value"); // check Value
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
                    // check value name
                    readString(checkRequestArry [ nChecks ].value_name, it->giveCurrentLine(), "name", MAX_NAME_LENGTH);
                    readString(checkRequestArry [ nChecks ].comment, it->giveCurrentLine(), "comment", MAX_COMMENT_LENGTH, 0);
                    readString(checkRequestArry [ nChecks ].post_expr, it->giveCurrentLine(), "expr", MAX_EXPR_LENGTH, 0);
                    checkRequestArry [ nChecks ].flag = readInteger(it->giveCurrentLine(), "noCheckFlag", 0);
#endif
                    //checkRequestArry[nChecks].checked    = 0;
                    nChecks++;
                } else if ( !strncmp(it->giveToken(1), "#BEAM_ELEMENT", 13) ||
                           !strncmp(it->giveToken(1), "#EIG_BEAM_ELEMENT", 17) ) {
                    /* if (it->giveNumberOfTokens() != 7) {
                     * printCheckErrors1 ("Missing arguments on line: %s",it->giveCurrentLine());
                     * exit(1);
                     * }*/
                    if ( !strncmp(it->giveToken(1), "#EIG_BEAM_ELEMENT", 17) ) {
                        checkRequestArry [ nChecks ].type       = EigNodeSideVal;
                        checkRequestArry [ nChecks ].eigValNum   = readInteger(it->giveCurrentLine(), "EigNum"); // eig Val Num
                    } else {
                        checkRequestArry [ nChecks ].type       = BeamElementVal;
                        checkRequestArry [ nChecks ].eigValNum   = 0;
                    }

#ifndef EXTRACTOR_MODE
                    // in EXTRACTOR_MODE tstep value has no meaning
                    checkRequestArry [ nChecks ].timeStep   = readDouble(it->giveCurrentLine(), "tStep");
#endif
                    checkRequestArry [ nChecks ].elNumber   = readInteger(it->giveCurrentLine(), "number");
                    //checkRequestArry[nChecks].gpNum      = readInteger (it->giveCurrentLine(), "gp"); // gp number
                    checkRequestArry [ nChecks ].gpNum      = 1;
                    checkRequestArry [ nChecks ].recNumber  = readInteger(it->giveCurrentLine(), "record"); // element record num.
                    checkRequestArry [ nChecks ].compNumber = readInteger(it->giveCurrentLine(), "component"); // component number
                    if ( hasString(it->giveCurrentLine(), "tolerance") ) {
                        checkRequestArry [ nChecks ].tolerance = readDouble(it->giveCurrentLine(), "tolerance"); // check tolerance
                    } else {
                        checkRequestArry [ nChecks ].tolerance = tolerance;
                    }

#ifdef CHECKER_MODE
                    checkRequestArry [ nChecks ].checkValue = readDouble(it->giveCurrentLine(), "value"); // check Value
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
                    // check value name
                    readString(checkRequestArry [ nChecks ].value_name, it->giveCurrentLine(), "name", MAX_NAME_LENGTH);
                    readString(checkRequestArry [ nChecks ].comment, it->giveCurrentLine(), "comment", MAX_COMMENT_LENGTH, 0);
                    readString(checkRequestArry [ nChecks ].post_expr, it->giveCurrentLine(), "expr", MAX_EXPR_LENGTH, 0);
                    checkRequestArry [ nChecks ].flag = readInteger(it->giveCurrentLine(), "noCheckFlag", 0);
#endif
                    //checkRequestArry[nChecks].checked    = 0;
                    nChecks++;
                } else if ( !strncmp(it->giveToken(1), "#REACTION", 9) ) {
                    /*if (it->giveNumberOfTokens() != 5) {
                     * printCheckErrors1 ("Missing arguments on line: %s",it->giveCurrentLine());
                     * exit(1);
                     * }*/
                    checkRequestArry [ nChecks ].type       = ReactionForce;
#ifndef EXTRACTOR_MODE
                    // in EXTRACTOR_MODE tstep value has no meaning
                    checkRequestArry [ nChecks ].timeStep   = readDouble(it->giveCurrentLine(), "tStep");
#endif
                    checkRequestArry [ nChecks ].elNumber   = readInteger(it->giveCurrentLine(), "comNumber"); // node or side number
                    checkRequestArry [ nChecks ].dofNumber  = readInteger(it->giveCurrentLine(), "dof"); // dof number
                    if ( hasString(it->giveCurrentLine(), "tolerance") ) {
                        checkRequestArry [ nChecks ].tolerance = readDouble(it->giveCurrentLine(), "tolerance"); // check tolerance
                    } else {
                        checkRequestArry [ nChecks ].tolerance = tolerance;
                    }

#ifdef CHECKER_MODE
                    checkRequestArry [ nChecks ].checkValue = readDouble(it->giveCurrentLine(), "value"); // check Value
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
                    // check value name
                    readString(checkRequestArry [ nChecks ].value_name, it->giveCurrentLine(), "name", MAX_NAME_LENGTH);
                    readString(checkRequestArry [ nChecks ].comment, it->giveCurrentLine(), "comment", MAX_COMMENT_LENGTH, 0);
                    readString(checkRequestArry [ nChecks ].post_expr, it->giveCurrentLine(), "expr", MAX_EXPR_LENGTH, 0);
                    checkRequestArry [ nChecks ].flag = readInteger(it->giveCurrentLine(), "noCheckFlag", 0);
#endif
                    //checkRequestArry[nChecks].checked    = 0;
                    nChecks++;
                } else if ( !strncmp(it->giveToken(1), "#QUASI-REACTION", 9) ) {
                    /*if (it->giveNumberOfTokens() != 5) {
                     * printCheckErrors1 ("Missing arguments on line: %s",it->giveCurrentLine());
                     * exit(1);
                     * }*/
                    checkRequestArry [ nChecks ].type       = QuasiReaction;
#ifndef EXTRACTOR_MODE
                    // in EXTRACTOR_MODE tstep value has no meaning
                    checkRequestArry [ nChecks ].timeStep   = readDouble(it->giveCurrentLine(), "tStep");
#endif
                    checkRequestArry [ nChecks ].elNumber   = readInteger(it->giveCurrentLine(), "comNumber"); // node or side number
                    checkRequestArry [ nChecks ].dofNumber  = readInteger(it->giveCurrentLine(), "dof"); // dof number
                    if ( hasString(it->giveCurrentLine(), "tolerance") ) {
                        checkRequestArry [ nChecks ].tolerance = readDouble(it->giveCurrentLine(), "tolerance"); // check tolerance
                    } else {
                        checkRequestArry [ nChecks ].tolerance = tolerance;
                    }

#ifdef CHECKER_MODE
                    checkRequestArry [ nChecks ].checkValue = readDouble(it->giveCurrentLine(), "value"); // check Value
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
                    // check value name
                    readString(checkRequestArry [ nChecks ].value_name, it->giveCurrentLine(), "name", MAX_NAME_LENGTH);
                    readString(checkRequestArry [ nChecks ].comment, it->giveCurrentLine(), "comment", MAX_COMMENT_LENGTH, 0);
                    readString(checkRequestArry [ nChecks ].post_expr, it->giveCurrentLine(), "expr", MAX_EXPR_LENGTH, 0);
                    checkRequestArry [ nChecks ].flag = readInteger(it->giveCurrentLine(), "noCheckFlag", 0);
#endif
                    //checkRequestArry[nChecks].checked    = 0;
                    nChecks++;
                } else if ( !strncmp(it->giveToken(1), "#EIGVAL", 7) ) {
                    /*if (it->giveNumberOfTokens() != 3) {
                     * printCheckErrors1 ("Missing arguments on line: %s",it->giveCurrentLine());
                     * exit(1);
                     * }*/
                    checkRequestArry [ nChecks ].type       = Eigval;
#ifndef EXTRACTOR_MODE
                    // in EXTRACTOR_MODE tstep value has no meaning
                    checkRequestArry [ nChecks ].timeStep   = readDouble(it->giveCurrentLine(), "tStep");
#endif
                    checkRequestArry [ nChecks ].eigValNum   = readInteger(it->giveCurrentLine(), "EigNum"); // eig Val Num
                    if ( hasString(it->giveCurrentLine(), "tolerance") ) {
                        checkRequestArry [ nChecks ].tolerance = readDouble(it->giveCurrentLine(), "tolerance"); // check tolerance
                    } else {
                        checkRequestArry [ nChecks ].tolerance = tolerance;
                    }

#ifdef CHECKER_MODE
                    checkRequestArry [ nChecks ].checkValue = readDouble(it->giveCurrentLine(), "value"); // check Value
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
                    // check value name
                    readString(checkRequestArry [ nChecks ].value_name, it->giveCurrentLine(), "name", MAX_NAME_LENGTH);
                    readString(checkRequestArry [ nChecks ].comment, it->giveCurrentLine(), "comment", MAX_COMMENT_LENGTH, 0);
                    readString(checkRequestArry [ nChecks ].post_expr, it->giveCurrentLine(), "expr", MAX_EXPR_LENGTH, 0);
                    checkRequestArry [ nChecks ].flag = readInteger(it->giveCurrentLine(), "noCheckFlag", 0);
#endif
                    //checkRequestArry[nChecks].checked    = 0;
                    nChecks++;
                } else if ( !strncmp(it->giveToken(1), "#LOADLEVEL", 10) ) {
                    /*if (it->giveNumberOfTokens() != 3) {
                     * printCheckErrors1 ("Missing arguments on line: %s",it->giveCurrentLine());
                     * exit(1);
                     * }*/
                    checkRequestArry [ nChecks ].type       = LoadLevel;
#ifndef EXTRACTOR_MODE
                    // in EXTRACTOR_MODE tstep value has no meaning
                    checkRequestArry [ nChecks ].timeStep   = readDouble(it->giveCurrentLine(), "tStep");
#endif
                    if ( hasString(it->giveCurrentLine(), "tolerance") ) {
                        checkRequestArry [ nChecks ].tolerance = readDouble(it->giveCurrentLine(), "tolerance"); // check tolerance
                    } else {
                        checkRequestArry [ nChecks ].tolerance = tolerance;
                    }

#ifdef CHECKER_MODE
                    checkRequestArry [ nChecks ].checkValue = readDouble(it->giveCurrentLine(), "value"); // check Value
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
                    // check value name
                    readString(checkRequestArry [ nChecks ].value_name, it->giveCurrentLine(), "name", MAX_NAME_LENGTH);
                    readString(checkRequestArry [ nChecks ].comment, it->giveCurrentLine(), "comment", MAX_COMMENT_LENGTH, 0);
                    readString(checkRequestArry [ nChecks ].post_expr, it->giveCurrentLine(), "expr", MAX_EXPR_LENGTH, 0);
                    checkRequestArry [ nChecks ].flag = readInteger(it->giveCurrentLine(), "noCheckFlag", 0);
#endif
                    //checkRequestArry[nChecks].checked    = 0;
                    nChecks++;
                } else if ( !strncmp(it->giveToken(1), "#NITE", 5) ) {
                    /*if (it->giveNumberOfTokens() != 3) {
                     * printCheckErrors1 ("Missing arguments on line: %s",it->giveCurrentLine());
                     * exit(1);
                     * }*/
                    checkRequestArry [ nChecks ].type       = NumberOfIter;
#ifndef EXTRACTOR_MODE
                    // in EXTRACTOR_MODE tstep value has no meaning
                    checkRequestArry [ nChecks ].timeStep   = readDouble(it->giveCurrentLine(), "tStep");
#endif
                    if ( hasString(it->giveCurrentLine(), "tolerance") ) {
                        checkRequestArry [ nChecks ].tolerance = readDouble(it->giveCurrentLine(), "tolerance"); // check tolerance
                    } else {
                        checkRequestArry [ nChecks ].tolerance = tolerance;
                    }

#ifdef CHECKER_MODE
                    checkRequestArry [ nChecks ].checkValue = readDouble(it->giveCurrentLine(), "value"); // check Value
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
                    // check value name
                    readString(checkRequestArry [ nChecks ].value_name, it->giveCurrentLine(), "name", MAX_NAME_LENGTH);
                    readString(checkRequestArry [ nChecks ].comment, it->giveCurrentLine(), "comment", MAX_COMMENT_LENGTH, 0);
                    readString(checkRequestArry [ nChecks ].post_expr, it->giveCurrentLine(), "expr", MAX_EXPR_LENGTH, 0);
                    checkRequestArry [ nChecks ].flag = readInteger(it->giveCurrentLine(), "noCheckFlag", 0);
#endif
                    //checkRequestArry[nChecks].checked    = 0;
                    nChecks++;

#ifdef EXTRACTOR_MODE
                } else if ( !strncmp(it->giveToken(1), "#TIME", 5) ) {
                    checkRequestArry [ nChecks ].type       = TimeVal;
                    nChecks++;
                } else if ( !strncmp(it->giveToken(1), "#UTIME", 6) ) {
                    checkRequestArry [ nChecks ].type       = SolutionStepUserTime;
                    nChecks++;
#endif
                } else {
                    printCheckErrors1( "Check section syntax error on line: %s", it->giveCurrentLine() );
                    exit(1);
                }

                if ( nChecks >= MAX_CHECK_ITEMS ) {
                    printCheckError("Too many check records");
                    exit(1);
                }

                // move to next line
            }  while ( !it->isEOF() );

            it->giveLineFromInput();
        }

        it->giveLineFromInput();
    }

    if ( ReadOutFileFlag && ReadCheckSection ) {
        // input file succesfully parsed, out file section and check sections found
    } else if ( ReadCheckSection == 0 ) {
        printCheckError("Check section not found, exiting");
        exit(1);
    } else if ( ReadOutFileFlag == 0 ) {
        printCheckError("Output file record not found, exiting");
        exit(1);
    }

    // open ouput file
    if ( ( outputFileStream = fopen(outFileName, "r") ) != NULL ) {
        outFileTokenizer = new Tokenizer(outputFileStream);
    } else {
        printCheckError("Can't open input file, exiting");
        exit(1);
    }


    //
    // begin check
    //

    int errorFlag = 0; // no error yet
#ifdef EXTRACTOR_MODE
    // ensure timeSteps are equal to be sorted not regarding these values
    // and timeStep or eigValNum equal to -1 means any step when seeking
    // for these values.
    for ( i = 0; i < nChecks; i++ ) {
        checkRequestArry [ i ].timeStep = -1.0;
        checkRequestArry [ i ].eigValNum = -1;
    }

#endif
    // initaialize check order array
    for ( i = 0; i < nChecks; i++ ) {
        checkRequestOrderArry [ i ] = i;
    }

    // sort check records
    sortCheckRecords();

#ifdef EXTRACTOR_MODE
    // do infinite loop over records (records are extracted for each time step)
    do {
        // seek solotion step only once
        result =  seekSolutionStep(outFileTokenizer, -1);
#endif
    // pop records, seek to appropriate position and check value
    for ( i = 0; i < nChecks; i++ ) {
        indx = checkRequestOrderArry [ i ];
        // pop next check record, seek to appropriate position and check value
        switch ( checkRequestArry [ indx ].type ) {
#ifdef EXTRACTOR_MODE
        case TimeVal:
            //result = seekSolutionStep (outFileTokenizer, checkRequestArry[indx].timeStep);
            getCurrState(& currState);
            checkRequestArry [ indx ].extractedValue = currState.currStep;
            break;
#endif



        case Eigval:
#ifndef EXTRACTOR_MODE
            result =  seekSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep);
#endif
            result =  result && seekEigvalSolutionStep(outFileTokenizer, checkRequestArry [ indx ].eigValNum, checkVal);
#ifdef CHECKER_MODE
            result = result && checkValue(checkVal, checkRequestArry [ indx ].checkValue, checkRequestArry [ indx ].tolerance);
            if ( result <= 0 ) {
                printCheckRuleError(indx);
                errorFlag = 1;
            }

#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
            if ( result <= 0 ) {
                printStudentSystemExtractorError(indx);
                errorFlag = 1;
            }

            printStudentSystemExtractedVariable(checkVal, checkRequestArry [ indx ].value_name, checkRequestArry [ indx ].post_expr,
                                                checkRequestArry [ indx ].tolerance,
                                                checkRequestArry [ indx ].flag, checkRequestArry [ indx ].comment);
#endif
#ifdef EXTRACTOR_MODE
            checkRequestArry [ indx ].extractedValue = checkVal;
#endif
            break;


        case LoadLevel:
#ifndef EXTRACTOR_MODE
            result =  seekSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep);
#endif
            result = result && seekNLSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep, checkVal, dummy);
#ifdef CHECKER_MODE
            result = result && checkValue(checkVal, checkRequestArry [ indx ].checkValue, checkRequestArry [ indx ].tolerance);
            if ( result <= 0 ) {
                printCheckRuleError(indx);
                errorFlag = 1;
            }

#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
            if ( result <= 0 ) {
                printStudentSystemExtractorError(indx);
                errorFlag = 1;
            }

            printStudentSystemExtractedVariable(checkVal, checkRequestArry [ indx ].value_name, checkRequestArry [ indx ].post_expr,
                                                checkRequestArry [ indx ].tolerance,
                                                checkRequestArry [ indx ].flag, checkRequestArry [ indx ].comment);
#endif
#ifdef EXTRACTOR_MODE
            checkRequestArry [ indx ].extractedValue = checkVal;
#endif

            break;

        case NumberOfIter:
#ifndef EXTRACTOR_MODE
            result =  seekSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep);
#endif
            result = result && seekNLSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep, dummy, checkVal);
#ifdef CHECKER_MODE
            result = result && checkValue(checkVal, checkRequestArry [ indx ].checkValue, checkRequestArry [ indx ].tolerance);
            if ( result <= 0 ) {
                printCheckRuleError(indx);
                errorFlag = 1;
            }

#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
            if ( result <= 0 ) {
                printStudentSystemExtractorError(indx);
                errorFlag = 1;
            }

            printStudentSystemExtractedVariable(checkVal, checkRequestArry [ indx ].value_name, checkRequestArry [ indx ].post_expr,
                                                checkRequestArry [ indx ].tolerance,
                                                checkRequestArry [ indx ].flag, checkRequestArry [ indx ].comment);
#endif
#ifdef EXTRACTOR_MODE
            checkRequestArry [ indx ].extractedValue = checkVal;
#endif

            break;


        case NodeSideVal:
        case EigNodeSideVal:
#ifndef EXTRACTOR_MODE
            result =  seekSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep);
#endif
            if ( checkRequestArry [ indx ].type == EigNodeSideVal ) {
                result =  result && seekEigvalSolutionStep(outFileTokenizer, checkRequestArry [ indx ].eigValNum, checkVal);
            }

            result = result && seekNodeRecord(outFileTokenizer, checkRequestArry [ indx ].elNumber);
            result = result && seekDofRecord(outFileTokenizer, checkRequestArry [ indx ].dofNumber,
                                             checkRequestArry [ indx ].compNumber, checkVal);
#ifdef CHECKER_MODE
            result = result && checkValue(checkVal, checkRequestArry [ indx ].checkValue, checkRequestArry [ indx ].tolerance);
            if ( result <= 0 ) {
                printCheckRuleError(indx);
                errorFlag = 1;
            }

#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
            if ( result <= 0 ) {
                printStudentSystemExtractorError(indx);
                errorFlag = 1;
            }

            printStudentSystemExtractedVariable(checkVal, checkRequestArry [ indx ].value_name, checkRequestArry [ indx ].post_expr,
                                                checkRequestArry [ indx ].tolerance,
                                                checkRequestArry [ indx ].flag, checkRequestArry [ indx ].comment);
#endif
#ifdef EXTRACTOR_MODE
            checkRequestArry [ indx ].extractedValue = checkVal;
#endif
            break;


        case ElementVal:
        case EigElementVal:
#ifndef EXTRACTOR_MODE
            result =  seekSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep);
#endif
            if ( checkRequestArry [ indx ].type == EigElementVal ) {
                result =  result && seekEigvalSolutionStep(outFileTokenizer, checkRequestArry [ indx ].eigValNum, checkVal);
            }

            result = result && seekElementRecord(outFileTokenizer, checkRequestArry [ indx ].elNumber);
            result = result && seekGPRecord(outFileTokenizer, checkRequestArry [ indx ].gpNum);
            if ( checkRequestArry [ indx ].recNumber == 0 ) { // strain
                result = result && seekStressStrainGPRecord(outFileTokenizer, 0, 1,
                                                            checkRequestArry [ indx ].compNumber,
                                                            checkVal);
            } else if ( checkRequestArry [ indx ].recNumber == 1 ) { // stress
                result = result && seekStressStrainGPRecord(outFileTokenizer, 1, 0,
                                                            checkRequestArry [ indx ].compNumber,
                                                            checkVal);
            } else if ( checkRequestArry [ indx ].recNumber == 2 ) { // status
                result = result && seekAndParseMaterialStatusRecord(outFileTokenizer, checkRequestArry [ indx ].keyword,
                                                                    checkRequestArry [ indx ].compNumber,
                                                                    checkVal);
            }

#ifdef CHECKER_MODE
            result = result && checkValue(checkVal, checkRequestArry [ indx ].checkValue, checkRequestArry [ indx ].tolerance);
            if ( result <= 0 ) {
                printCheckRuleError(indx);
                errorFlag = 1;
            }

#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
            if ( result <= 0 ) {
                printStudentSystemExtractorError(indx);
                errorFlag = 1;
            }

            printStudentSystemExtractedVariable(checkVal, checkRequestArry [ indx ].value_name, checkRequestArry [ indx ].post_expr,
                                                checkRequestArry [ indx ].tolerance,
                                                checkRequestArry [ indx ].flag, checkRequestArry [ indx ].comment);
#endif
#ifdef EXTRACTOR_MODE
            checkRequestArry [ indx ].extractedValue = checkVal;
#endif
            break;



        case BeamElementVal:
        case EigBeamElementVal:
#ifndef EXTRACTOR_MODE
            result =  seekSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep);
#endif
            if ( checkRequestArry [ indx ].type == EigBeamElementVal ) {
                result =  result && seekEigvalSolutionStep(outFileTokenizer, checkRequestArry [ indx ].eigValNum, checkVal);
            }

            result = result && seekElementRecord(outFileTokenizer, checkRequestArry [ indx ].elNumber);
            if ( checkRequestArry [ indx ].recNumber == 0 ) { // strain
                result = result && seekBeamRecord(outFileTokenizer, 0, 1,
                                                  checkRequestArry [ indx ].compNumber,
                                                  checkVal);
            } else {
                result = result && seekBeamRecord(outFileTokenizer, 1, 0,
                                                  checkRequestArry [ indx ].compNumber,
                                                  checkVal);
            }

#ifdef CHECKER_MODE
            result = result && checkValue(checkVal, checkRequestArry [ indx ].checkValue, checkRequestArry [ indx ].tolerance);
            if ( result <= 0 ) {
                printCheckRuleError(indx);
                errorFlag = 1;
            }

#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
            if ( result <= 0 ) {
                printStudentSystemExtractorError(indx);
                errorFlag = 1;
            }

            printStudentSystemExtractedVariable(checkVal, checkRequestArry [ indx ].value_name, checkRequestArry [ indx ].post_expr,
                                                checkRequestArry [ indx ].tolerance,
                                                checkRequestArry [ indx ].flag, checkRequestArry [ indx ].comment);
#endif
#ifdef EXTRACTOR_MODE
            checkRequestArry [ indx ].extractedValue = checkVal;
#endif
            break;
        case ReactionForce:
#ifndef EXTRACTOR_MODE
            result =   seekSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep);
#endif
            result = result && seekFirstReactionsRecord(outFileTokenizer);
            result = result && seekReactionRecord(outFileTokenizer, checkRequestArry [ indx ].elNumber,
                                                  checkRequestArry [ indx ].dofNumber, checkVal);
#ifdef CHECKER_MODE
            result = result && checkValue(checkVal, checkRequestArry [ indx ].checkValue, checkRequestArry [ indx ].tolerance);
            if ( result <= 0 ) {
                printCheckRuleError(indx);
                errorFlag = 1;
            }

#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
            if ( result <= 0 ) {
                printStudentSystemExtractorError(indx);
                errorFlag = 1;
            }

            printStudentSystemExtractedVariable(checkVal, checkRequestArry [ indx ].value_name, checkRequestArry [ indx ].post_expr,
                                                checkRequestArry [ indx ].tolerance,
                                                checkRequestArry [ indx ].flag, checkRequestArry [ indx ].comment);
#endif
#ifdef EXTRACTOR_MODE
            checkRequestArry [ indx ].extractedValue = checkVal;
#endif
            break;
        case QuasiReaction:
#ifndef EXTRACTOR_MODE
            result =   seekSolutionStep(outFileTokenizer, checkRequestArry [ indx ].timeStep);
#endif
            result = result && seekFirstQuasiReactionRecord(outFileTokenizer);
            result = result && seekQuasiReactionRecord(outFileTokenizer, checkRequestArry [ indx ].elNumber,
                                                       checkRequestArry [ indx ].dofNumber, checkVal);
            if ( result <= 0 ) {
                checkVal = 0.0;
                result = 1;
            }

#ifdef CHECKER_MODE
            result = result && checkValue(checkVal, checkRequestArry [ indx ].checkValue, checkRequestArry [ indx ].tolerance);
            if ( result <= 0 ) {
                printCheckRuleError(indx);
                errorFlag = 1;
            }

#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
            if ( result <= 0 ) {
                printStudentSystemExtractorError(indx);
                errorFlag = 1;
            }

            printStudentSystemExtractedVariable(checkVal, checkRequestArry [ indx ].value_name, checkRequestArry [ indx ].post_expr,
                                                checkRequestArry [ indx ].tolerance,
                                                checkRequestArry [ indx ].flag, checkRequestArry [ indx ].comment);
#endif
#ifdef EXTRACTOR_MODE
            checkRequestArry [ indx ].extractedValue = checkVal;
#endif
            break;

#ifdef EXTRACTOR_MODE
        case SolutionStepUserTime:
            //result = seekSolutionStep (outFileTokenizer, checkRequestArry[indx].timeStep);
            result = seekStepUserTime(outFileTokenizer, checkVal);
            checkRequestArry [ indx ].extractedValue = checkVal;
            break;
#endif


        default:
            printCheckError("Unknown check type encountered - internal error");
        }

        //checkRequestArry[indx].checked = 1;
    } // end loop over records

#ifdef EXTRACTOR_MODE

    if ( result <= 0 ) {
        if ( outFileTokenizer->isEOF() ) {
            errorFlag = 1;
        } else {
            printExtractorError(indx);
            errorFlag = 1;
        }
    }


    // print extracted values into columns, preserving input order
    if ( !errorFlag ) {
        for ( i = 0; i < nChecks; i++ ) {
            //indx = checkRequestOrderArry[i];
            printf("%10le ", checkRequestArry [ i ].extractedValue);
        }

        printf("\n");
    }

    // end of infinite loop until errorFlag
}
while ( !errorFlag ) {
    ;
}

#endif


    if ( errorFlag ) {
#ifndef EXTRACTOR_MODE
        //printCheckErrors1 ("Error(s) encountered while parsing %s file",inputStreamName);
        printf("---Checker: %50s    [%sFAILED%s]\n", basename(inputStreamName), RED, NORMAL);
#endif
    } else {
#ifdef CHECKER_MODE
        printf("---Checker: %50s    [  %sOK%s  ]\n", basename(inputStreamName), GREEN, NORMAL);
#endif
#ifdef SYSYTEM_STUDENT_EXTRACTOR_MODE
        printf("---StudentSystemExtractor: %20s file succesfully processed\n", inputStreamName);
#endif
    }

    delete it;
    delete outFileTokenizer;
    fclose(inputStream);
    fclose(outputFileStream);
}



void
printHelp() {
#ifdef EXTRACTOR_MODE
    printf("Usage: extractor -f inputStream\n\n");
    printf("where inputStream is input file name path, containing description of\n");
    printf("exracted values. The format of input file name is following:\n");
    printf("\toofemOutputFileName\n");
    printf("\t#%%BEGIN_CHECK%%\n");
    printf("\t[#NODE number # dof # unknown #]\n");
    printf("\t[#SIDE number # dof # unknown #]\n");
    printf("\t[#EIG_NODE number # dof # unknown #]\n");
    printf("\t[#EIG_SIDE number # dof # unknown #]\n");
    printf("\t[#ELEMENT number # gp # record # keyword # component #]\n");
    printf("\t[#EIG_ELEMENT number # gp # record # component #]\n");
    printf("\t[#BEAM_ELEMENT number # record # component #]\n");
    printf("\t[#EIG_BEAM_ELEMENT number # record # component #]\n");
    printf("\t[#QUASI-REACTION comNumber # dof # ]\n");
    printf("\t[#REACTION comNumber # dof # ]\n");
    printf("\t[#EIGVAL]\n");
    printf("\t[#LOADLEVEL]\n");
    printf("\t[#NITE]\n");
    printf("\t[#TIME]\n");
    printf("\t[#UTIME]\n");
    printf("\t#%%END_CHECK%%\n\n");
    printf("The %%BEGIN_CHECK%% and %%END_CHECK%% records are compulsory,\n");
    printf("the records in [ ] are optional, and can be used once or more times.\n");
    printf("values expected are integer numbers, except unknown value, where \n");
    printf("single char value is expected (i.e. 'd' for displacement, 'v' for velocity,\n");
    printf("and 'a' for acceleration) and keyword, where quoted string is expected.\n");
    printf("The keyword is used only when element status value is requested,\n");
    printf("i.e., when 'record' is equal to 2\n");
    printf("The output is printed to stdout, one row per each solution step,\n");
    printf("In particular columns, the extracted values are printed, preserving\n");
    printf("their order in input file\n\n");
    printf("(c) Borek Patzak 1999\n");
#endif
}

/*
 * int
 * popNextCheckRecord ()
 * {
 *
 * // returns index of next check message, which follows current reached state in input file.
 * // if all records have been checked, negative value is returned.
 * // if an error has been encountered, then negative number is returned.
 * //
 * stateType currState;
 * int i, nearestIndx = -1;
 * // find first unchecked record
 * for (i=0; i < nChecks; i++) { // loop over all records
 * if (!checkRequestArry[i].checked) {
 * nearestIndx = i;
 * break;
 * }
 * }
 *
 * if (nearestIndx == -1) return -1; // all records have been already checked
 *
 * getCurrState (&currState);
 *
 * for (i=1; i < nChecks; i++) { // loop over all records
 * if (checkRequestArry[i].checked) continue; // already checked, move to next record
 * if ((checkRequestArry[i].timeStep < checkRequestArry[nearestIndx].timeStep) &&
 * (checkRequestArry[i].timeStep >= currState.currStep)) {
 * nearestIndx = i;
 * continue;
 * } else if ((checkRequestArry[i].timeStep == checkRequestArry[nearestIndx].timeStep) &&
 *    (checkRequestArry[i].timeStep >= currState.currStep)) {
 * // same time step -> check further details
 * if (checkRequestArry[i].eigValNum < checkRequestArry[nearestIndx].eigValNum) {
 *  nearestIndx = i;
 *  continue;
 * } else if (checkRequestArry[i].eigValNum > checkRequestArry[nearestIndx].eigValNum) {
 *  continue;
 * } else if (checkRequestArry[i].type < checkRequestArry[nearestIndx].type) {
 *  nearestIndx = i;
 *  continue;
 * } else if ((checkRequestArry[i].type == Eigval)  || (checkRequestArry[i].type == LoadLevel)) {
 *   // multiple check for Eigval or LoadLevel data
 *   continue;
 * } else if (checkRequestArry[i].type > checkRequestArry[nearestIndx].type) {
 *  continue;
 * } else if ((checkRequestArry[i].type == NodeSideVal) || (checkRequestArry[i].type == EigNodeSideVal)) {
 *   // node side records
 *   if (checkRequestArry[i].elNumber < checkRequestArry[nearestIndx].elNumber) {
 *   nearestIndx = i;
 *   continue;
 *   } else if (checkRequestArry[i].elNumber == checkRequestArry[nearestIndx].elNumber) {
 *   // same nodeSide Record - check dofs
 *   if (checkRequestArry[i].dofNumber < checkRequestArry[nearestIndx].dofNumber) {
 *   nearestIndx = i;
 *   continue;
 *   } else continue; // sane dof record, checks with possible different unknowns
 *   }
 * } else if ((checkRequestArry[i].type == ElementVal) || (checkRequestArry[i].type == BeamElementVal) ||
 *       (checkRequestArry[i].type == EigElementVal) || (checkRequestArry[i].type == EigBeamElementVal)) {
 *   if (checkRequestArry[i].elNumber < checkRequestArry[nearestIndx].elNumber) {
 *   nearestIndx = i;
 *   continue;
 *   } else if (checkRequestArry[i].elNumber == checkRequestArry[nearestIndx].elNumber) {
 *   // same elements, check further
 *   if (checkRequestArry[i].gpNum < checkRequestArry[nearestIndx].gpNum) {
 *   nearestIndx = i;
 *   continue;
 *   } else if (checkRequestArry[i].gpNum == checkRequestArry[nearestIndx].gpNum) {
 *   // same gp -> check for component to check
 *   if (checkRequestArry[i].recNumber < checkRequestArry[nearestIndx].recNumber) {
 *     nearestIndx = i;
 *     continue;
 *   } else continue; // same element record numbers,
 *   // only different componenets of same record checked.
 *   }
 *   }
 * } else if (checkRequestArry[i].type == ReactionForce) {
 *   // node side records
 *   if (checkRequestArry[i].elNumber < checkRequestArry[nearestIndx].elNumber) {
 *   nearestIndx = i;
 *   continue;
 *   } else if (checkRequestArry[i].elNumber == checkRequestArry[nearestIndx].elNumber) {
 *   // same nodeSide Record - check dofs
 *   if (checkRequestArry[i].dofNumber < checkRequestArry[nearestIndx].dofNumber) {
 *   nearestIndx = i;
 *   continue;
 *   } else continue; // sane dof record, checks with possible different unknowns
 *  }
 * }
 * }
 * }
 * return nearestIndx;
 * }
 */


int
compareRecords(int i, int j)
{
    /*
     * returns result < 0 if i-th record is before j-th record,
     * zero if records are equal,
     * otherwise returns > 0.
     */
    if ( i == j ) {
        return 0;
    }

    if ( ( checkRequestArry [ i ].timeStep < checkRequestArry [ j ].timeStep ) ) {
        return -1;
    } else if ( checkRequestArry [ i ].timeStep > checkRequestArry [ j ].timeStep ) {
        return 1;
        // same time step -> check further details
        // compare eigen valur number (like time step)
    } else if ( checkRequestArry [ i ].eigValNum < checkRequestArry [ j ].eigValNum ) {
        return -1;
    } else if ( checkRequestArry [ i ].eigValNum > checkRequestArry [ j ].eigValNum ) {
        return 1;
        // same eigval num -> check further details
    } else  if ( checkRequestArry [ i ].type < checkRequestArry [ j ].type ) {
        return -1;
    } else if ( checkRequestArry [ i ].type > checkRequestArry [ j ].type ) {
        return 1;
    } else {
        // same types -> check further details
        if ( checkRequestArry [ i ].type == LoadLevel ) {
            // multiple check for LoadLevel data
            return 0;
        } else  if ( checkRequestArry [ i ].type == Eigval ) {
            // multiple check for Eigval data -> compare eigValNum
            if ( checkRequestArry [ i ].eigValNum < checkRequestArry [ j ].eigValNum ) {
                return -1;
            } else if ( checkRequestArry [ i ].eigValNum > checkRequestArry [ j ].eigValNum ) {
                return 1;
            } else {
                return 0;
            }
        } else if ( ( checkRequestArry [ i ].type == NodeSideVal ) || ( checkRequestArry [ i ].type == EigNodeSideVal ) ) {
            // node side records
            if ( checkRequestArry [ i ].elNumber < checkRequestArry [ j ].elNumber ) {
                return -1;
            } else if ( checkRequestArry [ i ].elNumber > checkRequestArry [ j ].elNumber ) {
                return 1;
            } else {
                // same nodeSide Record - check dofs
                if ( checkRequestArry [ i ].dofNumber < checkRequestArry [ j ].dofNumber ) {
                    return -1;
                } else if ( checkRequestArry [ i ].dofNumber > checkRequestArry [ j ].dofNumber ) {
                    return 1;
                } else {
                    return 0; // sane dof record, checks with possible different unknowns
                }
            }
        } else if ( ( checkRequestArry [ i ].type == ElementVal ) || ( checkRequestArry [ i ].type == BeamElementVal ) ||
                   ( checkRequestArry [ i ].type == EigElementVal ) || ( checkRequestArry [ i ].type == EigBeamElementVal ) ) {
            // element records
            if ( checkRequestArry [ i ].elNumber < checkRequestArry [ j ].elNumber ) {
                return -1;
            } else if ( checkRequestArry [ i ].elNumber > checkRequestArry [ j ].elNumber ) {
                return 1;
            } else {
                // same elements, check further
                if ( checkRequestArry [ i ].gpNum < checkRequestArry [ j ].gpNum ) {
                    return -1;
                } else if ( checkRequestArry [ i ].gpNum > checkRequestArry [ j ].gpNum ) {
                    return 1;
                } else {
                    // same gp -> check for component to check
                    if ( checkRequestArry [ i ].recNumber < checkRequestArry [ j ].recNumber ) {
                        return -1;
                    } else if ( checkRequestArry [ i ].recNumber > checkRequestArry [ j ].recNumber ) {
                        return 1;
                    } else {
                        // same element record numbers,
                        // only different componenets of same record checked.
                        return 0;
                    }
                }
            }
        } else if ( ( checkRequestArry [ i ].type == ReactionForce ) || ( checkRequestArry [ i ].type == QuasiReaction ) ) {
            // reaction node side records
            if ( checkRequestArry [ i ].elNumber < checkRequestArry [ j ].elNumber ) {
                return -1;
            } else if ( checkRequestArry [ i ].elNumber > checkRequestArry [ j ].elNumber ) {
                return 1;
            } else {
                // same nodeSide Record - check dofs
                if ( checkRequestArry [ i ].dofNumber < checkRequestArry [ j ].dofNumber ) {
                    return -1;
                } else if ( checkRequestArry [ i ].dofNumber > checkRequestArry [ j ].dofNumber ) {
                    return 1;
                } else {
                    return 0; // sane dof record, checks with possible different unknowns
                }
            }
        }
    }

    return -1;
}




int checkValue(double firstVal, double secondVal, double tolerance)
{
    // checks two values for equality, with given tolerance.
    // returns nonzero if two given values are equal, zero
    // if not.

    return ( fabs(firstVal - secondVal) <= tolerance );
}

void printCheckError(const char *str) {
    printf("***Checker Error: %s\n", str);
}

void printCheckErrori1(const char *str, int val) {
    printf("***Checker Error:");
    printf(str, val);
    printf("\n");
}

void printCheckErrors1(const char *str, const char *str2) {
    printf("***Checker Error:");
    printf(str, str2);
    printf("\n");
}

void printCheckErrors2(const char *str, const char *str2, const char *str3) {
    printf("***Checker Error:");
    printf(str, str2, str3);
    printf("\n");
}

void printCheckRuleError(int val) {
    printf("*Error, when checking rule %d", val + 1);
    printf("\n");
}

char *getPosAfter(char *source, const char *idString)
//
// returns possition of substring idString in source
// return value pointer at the end of occurence idString in source
// (idString must be separated from rest by blank or by tabulator
//
{
    char *str1, *helpSource = source;
    int len = strlen(idString);
    do {
        if ( ( str1 = strstr(helpSource, idString) ) == NULL ) {
            return NULL;
        }

        helpSource = str1 + 1;
    } while ( !( ( * ( helpSource + len - 1 ) == ' ' ) || ( * ( helpSource + len - 1 ) == '\t' ) ) );

    return str1 + len;
}

int readInteger(char *source, const char *idString, int compulsory)
{
    char *str = getPosAfter(source, idString);


    if ( str == NULL ) {
        if ( compulsory ) {
            printCheckErrors2("Missing keyword \"%s\", while parsing source line\n\t%s", idString, source);
            printCheckErrors1("Errors encountered while parsing %s file", inputStreamName);
            exit(1);
        } else {
            return 0;
        }
    }

    return atoi(str);
}

double readDouble(char *source, const char *idString, int compulsory)
{
    char *str = getPosAfter(source, idString);

    if ( str == NULL ) {
        if ( compulsory ) {
            printCheckErrors2("Missing keyword \"%s\", while parsing source line\n\t%s", idString, source);
            printCheckErrors1("Errors encountered while parsing %s file", inputStreamName);
            exit(1);
        } else {
            return 0.;
        }
    }

    return strtod(str, NULL);
}

int hasString(char *source, const char *idString)
{
    char *str = getPosAfter(source, idString);

    if ( str == NULL ) {
        return 0;
    }

    return 1;
}

void readString(char *string, char *source, const char *idString, int maxLen, int compulsory)
{
    char *str = getPosAfter(source, idString);
    if ( str == NULL ) {
        if ( compulsory ) {
            printCheckErrors2("Missing keyword \"%s\", while parsing source line \n\t%s", idString, source);
            printCheckErrors1("Errors encountered while parsing %s file", inputStreamName);
            exit(1);
        } else {
            * string = 0;
            return;
        }
    }

    char *curr = str;
    int curr_len = 0;

    // skip whitespaces
    while ( ( * curr == ' ' ) || ( * curr == '\n' ) || ( * curr == '\t' ) || ! * curr ) {
        curr++;
    }

    if ( * curr == '\"' ) {
        curr++; // skip quotes
        while ( ( * curr != '\"' ) && ( * curr != '\n' ) && ( * curr != 0 ) && ( curr_len < maxLen ) ) {
            string [ curr_len++ ] = * curr;
            curr++;
        }

        if ( * curr != '\"' ) {
            printCheckErrors1("Missing end \" or max name length exceeded, while parsing source line \n\t%s", source);
            printCheckErrors1("Errors encountered while parsing %s file", inputStreamName);
            exit(1);
        }
    } else {
        while ( isalnum(* curr) && ( curr_len < maxLen ) ) {
            string [ curr_len++ ] = * curr;
            curr++;
        }
    }
}

char readChar(char *source, const char *idString)
{
    char *str = getPosAfter(source, idString);

    if ( str == NULL ) {
        printCheckErrors2("Missing keyword \"%s\", while parsing source line \n\t%s", idString, source);
        printCheckErrors1("Errors encountered while parsing %s file", inputStreamName);
        exit(1);
    }

    char *curr = str;
    // skip whitespaces
    while ( ( * curr == ' ' ) || ( * curr == '\n' ) || ( * curr == '\t' ) || ! * curr ) {
        curr++;
    }

    return * curr;
}

void printStudentSystemExtractedVariable(double checkVal, char *name, char *expr_frmt, double tolerance, int flag, char *comment)
{
    if ( * expr_frmt == 0 ) {
        printf("%s:%le:%le:%d:\"%s\":\n", name, tolerance, checkVal, flag, comment);
    } else {
        //char frmt[MAX_FRMT_LENGTH];
        //sprintf (frmt, "\%s:\%lf:%s:\%d:\"\%s\":\n", expr_frmt);
        printf("%s:%le:", name, tolerance);
        printf(expr_frmt, checkVal);
        printf(":%d:\"%s\":\n", flag, comment);
    }
}

void printStudentSystemExtractorError(int indx)
{
    printf("***StudentSystemExtractor Error in processing rule %d", indx + 1);
}

void
printExtractorError(int indx)
{
    //
    if ( indx != 0 ) {
        printf("***Extractor Error in processing rule %d", indx + 1);
    }
}



void
sortCheckRecords()
{
    fprintf(stderr, "Sorting ...");
    quickSortArry(checkRequestOrderArry, 0, nChecks - 1);
    fprintf(stderr, "done\n");
    for ( int i = 0; i < nChecks; i++ ) {
        fprintf(stderr, "%d ", checkRequestOrderArry [ i ]);
    }

    fprintf(stderr, "\n");
}


void
quickSortArry(int *map, int l, int r)
{
    if ( r <= l ) {
        return;
    }

    int i = quickSortPartition(map, l, r);
    quickSortArry(map, l, i - 1);
    quickSortArry(map, i + 1, r);
}

int
quickSortPartition(int *map, int l, int r)
{
    int i = l - 1, j = r;
    int v = map [ r ];
    int swap;

    for ( ; ; ) {
        while ( compareRecords(map [ ++i ], v) < 0 ) {
            ;
        }

        while ( compareRecords(v, map [ --j ]) < 0 ) {
            if ( j == l ) {
                break;
            }
        }

        if ( i >= j ) {
            break;
        }

        swap = map [ i ];
        map [ i ] = map [ j ];
        map [ j ] = swap;
    }

    swap = map [ i ];
    map [ i ] = map [ r ];
    map [ r ] = swap;
    return i;
}
