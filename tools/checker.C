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



#include "seek.h"
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
// #define SILICON_GRAPHICS

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
            } else if ( checkRequestArry [ i ].eigValNum > checkRequestArry [ j ].eigValNum )       {
                return 1;
            } else                                                                                             {
                return 0;
            }
        } else if ( ( checkRequestArry [ i ].type == NodeSideVal ) || ( checkRequestArry [ i ].type == EigNodeSideVal ) ) {
            // node side records
            if ( checkRequestArry [ i ].elNumber < checkRequestArry [ j ].elNumber ) {
                return -1;
            } else if ( checkRequestArry [ i ].elNumber > checkRequestArry [ j ].elNumber )       {
                return 1;
            } else                                                                                           {
                // same nodeSide Record - check dofs
                if ( checkRequestArry [ i ].dofNumber < checkRequestArry [ j ].dofNumber ) {
                    return -1;
                } else if ( checkRequestArry [ i ].dofNumber > checkRequestArry [ j ].dofNumber )      {
                    return 1;
                } else                                                                                            {
                    return 0; // sane dof record, checks with possible different unknowns
                }
            }
        } else if ( ( checkRequestArry [ i ].type == ElementVal ) || ( checkRequestArry [ i ].type == BeamElementVal ) ||
                   ( checkRequestArry [ i ].type == EigElementVal ) || ( checkRequestArry [ i ].type == EigBeamElementVal ) ) {
            // element records
            if ( checkRequestArry [ i ].elNumber < checkRequestArry [ j ].elNumber ) {
                return -1;
            } else if ( checkRequestArry [ i ].elNumber > checkRequestArry [ j ].elNumber )       {
                return 1;
            } else                                                                                           {
                // same elements, check further
                if ( checkRequestArry [ i ].gpNum < checkRequestArry [ j ].gpNum ) {
                    return -1;
                } else if ( checkRequestArry [ i ].gpNum > checkRequestArry [ j ].gpNum )      {
                    return 1;
                } else                                                                                    {
                    // same gp -> check for component to check
                    if ( checkRequestArry [ i ].recNumber < checkRequestArry [ j ].recNumber ) {
                        return -1;
                    } else if ( checkRequestArry [ i ].recNumber > checkRequestArry [ j ].recNumber )     {
                        return 1;
                    } else                                                                                           {
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
            } else if ( checkRequestArry [ i ].elNumber > checkRequestArry [ j ].elNumber )       {
                return 1;
            } else                                                                                           {
                // same nodeSide Record - check dofs
                if ( checkRequestArry [ i ].dofNumber < checkRequestArry [ j ].dofNumber ) {
                    return -1;
                } else if ( checkRequestArry [ i ].dofNumber > checkRequestArry [ j ].dofNumber )      {
                    return 1;
                } else                                                                                            {
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

void printCheckError(char *str) {
    printf("***Checker Error: %s\n", str);
}

void printCheckErrori1(char *str, int val) {
    printf("***Checker Error:");
    printf(str, val);
    printf("\n");
}

void printCheckErrors1(char *str, char *str2) {
    printf("***Checker Error:");
    printf(str, str2);
    printf("\n");
}

void printCheckErrors2(char *str, char *str2, char *str3) {
    printf("***Checker Error:");
    printf(str, str2, str3);
    printf("\n");
}

void printCheckRuleError(int val) {
    printf("*Error, when checking rule %d", val + 1);
    printf("\n");
}

char *getPosAfter(char *source, char *idString)
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

int readInteger(char *source, char *idString, int compulsory)
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

double readDouble(char *source, char *idString, int compulsory)
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

int hasString(char *source, char *idString)
{
    char *str = getPosAfter(source, idString);

    if ( str == NULL ) {
        return 0;
    }

    return 1;
}

void readString(char *string, char *source, char *idString, int maxLen, int compulsory)
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

char readChar(char *source, char *idString)
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
