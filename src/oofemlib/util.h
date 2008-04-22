/* $Header: /home/cvs/bp/oofem/oofemlib/src/util.h,v 1.7 2003/05/19 13:03:58 bp Exp $ */
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


#ifndef util_h
#define util_h

/**
 * Reads the input line from file and converts all characters to lower case
 * @param line contains read line
 * @param len determines max number of read characters including terminationg '\0'
 */
char *giveLineFromInput(FILE *inputStream, char *line, int len);
/// Reads the input line from file
char *giveRawLineFromInput(FILE *inputStream, char *line, int len);

// Returns the name of the file containing the data of the problem.
char *giveInputDataFileName(char *dataInputFileName, int maxlen);

//Instanciates the new problem
EngngModel *InstanciateProblem(DataReader *dr, problemMode mode, int contextFlag, EngngModel *master = NULL);

#endif // util_h
