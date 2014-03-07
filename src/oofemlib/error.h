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

namespace oofem {
/** Cause oofem program termination by calling exit. */
OOFEM_EXPORT void oofem_exit(int code);

/**
 * Macros for printing errors.
 * This macro can be used only within classes that implement errorInfo function.
 */
//@{
#define OOFEM_FATAL(...) { oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, this->errorInfo(__func__).c_str(), __FILE__, __LINE__, __VA_ARGS__); oofem_exit(1); }
#define OOFEM_ERROR(...) { oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, this->errorInfo(__func__).c_str(), __FILE__, __LINE__, __VA_ARGS__); oofem_exit(1); }
#define OOFEM_WARNING(...) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, this->errorInfo(__func__).c_str(), __FILE__, __LINE__, __VA_ARGS__)
//@}

/**
 * Log-family macros that use OOFEM loggers (and exit for errors and fatals).
 */
//@{
#define OOFEM_SIMPLE_FATAL(...) { OOFEM_LOG_FATAL(__VA_ARGS__); oofem_exit(1); }
#define OOFEM_SIMPLE_ERROR(...) { OOFEM_LOG_ERROR(__VA_ARGS__); oofem_exit(1); }
#define OOFEM_SIMPLE_WARNING(...) OOFEM_LOG_WARNING(__VA_ARGS__)
//@}

} // end namespace oofem
#endif // error_h
