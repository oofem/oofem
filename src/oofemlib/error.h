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

//   ************************************
//   *** Error macros                 ***
//   ************************************

#ifndef error_h
#define error_h

#include "logger.h"

#include <string>
#include <cstdio>
#include <cstdlib>

namespace oofem {
/** Cause oofem program termination by calling exit. */
#define OOFEM_EXIT(code) \
    oofem_logger.printStatistics(); \
    fprintf(stderr, "oofem exit code %d\n", code); \
    exit(code);

/**
 * Macros for printing errors.
 * This macro can be used only within classes that implement errorInfo function.
 */
//@{
#define OOFEM_FATAL(...) { oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, errorInfo(__func__).c_str(), __FILE__, __LINE__, __VA_ARGS__); OOFEM_EXIT(1); }
#define OOFEM_ERROR(...) { oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, errorInfo(__func__).c_str(), __FILE__, __LINE__, __VA_ARGS__); OOFEM_EXIT(1); }
#define OOFEM_WARNING(...) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, errorInfo(__func__).c_str(), __FILE__, __LINE__, __VA_ARGS__)
#define OOFEM_SERROR(...) { oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __func__, __FILE__, __LINE__, __VA_ARGS__); OOFEM_EXIT(1); }
#define OOFEM_SWARNING(...) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __func__, __FILE__, __LINE__, __VA_ARGS__)
//@}

OOFEM_EXPORT std::string errorInfo(const char *func);

} // end namespace oofem
#endif // error_h
