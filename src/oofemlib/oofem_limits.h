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

//
// FILE: oofem_limits.h
//

#ifndef oofem_limits_h
#define oofem_limits_h

namespace oofem {
/// Maximum input line length read.
#define OOFEM_MAX_LINE_LENGTH 32768
// size of token buffer, which holds tokens
#define OOFEM_MAX_TOKENS_LENGTH 32768
/// Maximum keyword name string length.
#define MAX_NAME_LENGTH 40
/// Maximum file name path string length.
#define MAX_FILENAME_LENGTH 120
/// Maximum class name string length
#define MAX_CLASSNAME_LENGTH 50
/// max number of tokens
#define OOFEM_MAX_TOKENS 8000
/// max legth of error message
#define MAX_ERROR_MSG_LENGTH 2048
} // end namespace oofem
#endif // oofem_limits_h

