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

#ifndef util_h
#define util_h

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif
#include "problemmode.h"

namespace oofem {
class DataReader;
class EngngModel;

/**
 * Reads the input line from file and converts all characters to lower case.
 * @param inputStream Stream to read from.
 * @param line Contains read line.
 * @param len Determines max number of read characters including terminating '\0'.
 * @return Newly allocated pointer containing string.
 */
char *giveLineFromInput(FILE *inputStream, char *line, int len);
/**
 * Reads the input line from file.
 * @param inputStream Stream to read from.
 * @param line Contains read line.
 * @param len Determines max number of read characters including terminating '\0'.
 * @return Newly allocated pointer containing string.
 */
char *giveRawLineFromInput(FILE *inputStream, char *line, int len);

/**
 * Returns the name of the file containing the data of the problem.
 * @param dataInputFileName Char buffer at least the size maxlen.
 * @param maxlen Maximum length of file name.
 * @return Same as dataInputFileName.
 */
char *giveInputDataFileName(char *dataInputFileName, int maxlen);

/**
 * Instanciates the new problem.
 * @param dr DataReader containing the problem data.
 * @param mode Mode determining macro or micro problem.
 * @param master Master problem in case of multiscale computations.
 * @param parallelFlag Determines if the problem should be run in parallel or not.
 * @todo Document the contextFlag input.
 */
EngngModel *InstanciateProblem(DataReader *dr, problemMode mode, int contextFlag, EngngModel *master = 0, bool parallelFlag = false);

/**
 * Static storage for temporary strings to solve compiler warnings about conversion from string constant to char*.
 * @param src Constant data to be copied over.
 * @return Pointer to static array of data, overwritten with src.
 */
char *oofem_tmpstr(const char *src);
} // end namespace oofem
#endif // util_h
