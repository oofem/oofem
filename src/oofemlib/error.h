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
 * Macros calling the error/warning service, which is defined for all classes
 * derived form FEMComponent and some others classes.
 * It inserts automatically file and line informations into call.
 * This macro can be used only within classes that implement error service.
 * The corresponding error service is assumed to have fist argument the filename, second argument line number,
 * followed by format string and optional arguments (as in printf family of functions). It typically formats the
 * final message and uses oofem default loggers to report message.
 */
//@{
#define OOFEM_CLASS_ERROR(...) error(__FILE__, __LINE__, __VA_ARGS__);
#define _error(...) OOFEM_CLASS_ERROR(__VA_ARGS__)
#define OOFEM_CLASS_WARNING(...) warning(__FILE__, __LINE__, __VA_ARGS__);
#define _warning(...) OOFEM_CLASS_WARNING(__VA_ARGS__)




#define _error1(_1) OOFEM_CLASS_ERROR(_1)
#define _error2(_1, _2) OOFEM_CLASS_ERROR(_1, _2)
#define _error3(_1, _2, _3) OOFEM_CLASS_ERROR(_1, _2, _3)
#define _error4(_1, _2, _3, _4) OOFEM_CLASS_ERROR(_1, _2, _3, _4)
#define _error5(_1, _2, _3, _4, _5) OOFEM_CLASS_ERROR(_1, _2, _3, _4, _5)
#define _warning1(_1) OOFEM_CLASS_WARNING(_1)
#define _warning2(_1, _2) OOFEM_CLASS_WARNING(_1, _2)
#define _warning3(_1, _2, _3) OOFEM_CLASS_WARNING(_1, _2, _3)
#define _warning4(_1, _2, _3, _4) OOFEM_CLASS_WARNING(_1, _2, _3, _4)
//@}

/**
 * Log-family macros that use OOFEM loggers (and exit for errors and fatals).
 */
//@{
#define OOFEM_FATAL(...) { OOFEM_LOG_FATAL(__VA_ARGS__); oofem_exit(1); }
#define OOFEM_ERROR(...) { OOFEM_LOG_ERROR(__VA_ARGS__); oofem_exit(1); }
#define OOFEM_WARNING(...) { OOFEM_LOG_WARNING(__VA_ARGS__); }




#define OOFEM_FATAL1(_1) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, _1); oofem_exit(1); }
#define OOFEM_FATAL2(_1, _2) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, _1, _2); oofem_exit(1); }
#define OOFEM_FATAL3(_1, _2, _3) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, _1, _2, _3); oofem_exit(1); }

#define OOFEM_ERROR1(_1) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1); oofem_exit(1); }
#define OOFEM_ERROR2(_1, _2) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1, _2); oofem_exit(1); }
#define OOFEM_ERROR3(_1, _2, _3) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1, _2, _3); oofem_exit(1); }
#define OOFEM_ERROR4(_1, _2, _3, _4) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1, _2, _3, _4); oofem_exit(1); }
#define OOFEM_ERROR5(_1, _2, _3, _4, _5) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1, _2, _3, _4, _5); oofem_exit(1); }
#define OOFEM_ERROR6(_1, _2, _3, _4, _5, _6) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1, _2, _3, _4, _5, _6); oofem_exit(1); }

#define OOFEM_WARNING1(_1) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1); }
#define OOFEM_WARNING2(_1, _2) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1, _2); }
#define OOFEM_WARNING3(_1, _2, _3) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1, _2, _3); }
#define OOFEM_WARNING4(_1, _2, _3, _4) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1, _2, _3, _4); }
#define OOFEM_WARNING5(_1, _2, _3, _4, _5) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1, _2, _3, _4, _5); }
#define OOFEM_WARNING6(_1, _2, _3, _4, _5, _6) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1, _2, _3, _4, _5, _6); }
//@}

/**
 * Log-family macros that allow to pass file and line info. They use OOFEM loggers (and exit for fatals and errors).
 */
//@{
#define __OOFEM_FATAL(_file, _line, ...) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, _file, _line, __VA_ARGS__); oofem_exit(1); }
#define __OOFEM_ERROR(_file, _line, ...) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, _file, _line, __VA_ARGS__); oofem_exit(1); }
#define __OOFEM_WARNING(_file, _line, ...) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, _file, _line, __VA_ARGS__); }




#define __OOFEM_FATAL1(_file, _line, _1) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, _file, _line, _1); oofem_exit(1); }
#define __OOFEM_FATAL2(_file, _line, _1, _2) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, _file, _line, _1, _2); oofem_exit(1); }
#define __OOFEM_FATAL3(_file, _line, _1, _2, _3) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, _file, _line, _1, _2, _3); oofem_exit(1); }


#define __OOFEM_ERROR1(_file, _line, _1) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, _file, _line, _1); oofem_exit(1); }
#define __OOFEM_ERROR2(_file, _line, _1, _2) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, _file, _line, _1, _2); oofem_exit(1); }
#define __OOFEM_ERROR3(_file, _line, _1, _2, _3) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, _file, _line, _1, _2, _3); oofem_exit(1); }
#define __OOFEM_ERROR4(_file, _line, _1, _2, _3, _4) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, _file, _line, _1, _2, _3, _4); oofem_exit(1); }
#define __OOFEM_ERROR5(_file, _line, _1, _2, _3, _4, _5) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, _file, _line, _1, _2, _3, _4, _5); oofem_exit(1); }
#define __OOFEM_ERROR6(_file, _line, _1, _2, _3, _4, _5, _6) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, _file, _line, _1, _2, _3, _4, _5, _6); oofem_exit(1); }

#define __OOFEM_WARNING1(_file, _line, _1) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, _file, _line, _1); }
#define __OOFEM_WARNING2(_file, _line, _1, _2) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, _file, _line, _1, _2); }
#define __OOFEM_WARNING3(_file, _line, _1, _2, _3) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, _file, _line, _1, _2, _3); }
#define __OOFEM_WARNING4(_file, _line, _1, _2, _3, _4) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, _file, _line, _1, _2, _3, _4); }
#define __OOFEM_WARNING5(_file, _line, _1, _2, _3, _4, _5) { oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, _file, _line, _1, _2, _3, _4, _5); }
//@}
} // end namespace oofem
#endif // error_h
