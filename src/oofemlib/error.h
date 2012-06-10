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

//   ************************************
//   *** Error macros                 ***
//   ************************************

#ifndef error_h
#define error_h

#include "oofemcfg.h"
#include "logger.h"

namespace oofem {
/** Cause oofem program termination by calling exit. */
void oofem_exit(int code);

/**@
 * Macros calling the error/warning service, which is defined for all classes
 * derived form FEMComponent and some others classes.
 * It inserts automatically file and line informations into call.
 * This macro can be used only within classes that implement error service.
 * The corresponding error service is assumed to have fist argument the filename, second argument line number,
 * followed by format string and optional arguments (as in printf family of functions). It typically formats the
 * final message and uses oofem default loggers to report message.
 *
 * These macros add to given message also file and line information.
 * They can be implemented in a very elegant way using macro with variable
 * number of arguments (__VA_ARGS__). But since the many compilers do not support
 * macros with variable number of arguments we have two choices:
 * to implement this as a function with variable number of arguments, but then we can not
 * add file and line info via __FILE__ and __LINE__ macros. Or if file and line
 * info is preferred, then instead of single macro with variable number of arguments we can have
 * series of "classical" macros with increasing number of parameters.
 * The latter approach is used here.
 *
 */
//@{
#define OOFEM_CLASS_ERROR1(_1) error(__FILE__, __LINE__, _1);
#define OOFEM_CLASS_ERROR2(_1, _2) error(__FILE__, __LINE__, _1, _2);
#define OOFEM_CLASS_ERROR3(_1, _2, _3) error(__FILE__, __LINE__, _1, _2, _3);
#define OOFEM_CLASS_ERROR4(_1, _2, _3, _4) error(__FILE__, __LINE__, _1, _2, _3, _4);
#define OOFEM_CLASS_ERROR5(_1, _2, _3, _4, _5) error(__FILE__, __LINE__, _1, _2, _3, _4, _5);
#define OOFEM_CLASS_ERROR6(_1, _2, _3, _4, _5, _6) error(__FILE__, __LINE__, _1, _2, _3, _4, _5, _6);

#define _error1(_1) OOFEM_CLASS_ERROR1(_1)
#define _error2(_1, _2) OOFEM_CLASS_ERROR2(_1, _2)
#define _error3(_1, _2, _3) OOFEM_CLASS_ERROR3(_1, _2, _3)
#define _error4(_1, _2, _3, _4) OOFEM_CLASS_ERROR4(_1, _2, _3, _4)
#define _error5(_1, _2, _3, _4, _5) OOFEM_CLASS_ERROR5(_1, _2, _3, _4, _5)
#define _error6(_1, _2, _3, _4, _5, _6) OOFEM_CLASS_ERROR6(_1, _2, _3, _4, _5, _6)

#define OOFEM_CLASS_WARNING1(_1) warning(__FILE__, __LINE__, _1);
#define OOFEM_CLASS_WARNING2(_1, _2) warning(__FILE__, __LINE__, _1, _2);
#define OOFEM_CLASS_WARNING3(_1, _2, _3) warning(__FILE__, __LINE__, _1, _2, _3);
#define OOFEM_CLASS_WARNING4(_1, _2, _3, _4) warning(__FILE__, __LINE__, _1, _2, _3, _4);

#define _warning1(_1) OOFEM_CLASS_WARNING1(_1)
#define _warning2(_1, _2) OOFEM_CLASS_WARNING2(_1, _2)
#define _warning3(_1, _2, _3) OOFEM_CLASS_WARNING3(_1, _2, _3)
#define _warning4(_1, _2, _3, _4) OOFEM_CLASS_WARNING4(_1, _2, _3, _4)

#define OOFEM_CLASS_ERROR(_1) OOFEM_CLASS_ERROR1(_1)
#define _error(_1) _error1(_1)
#define OOFEM_CLASS_WARNING(_1)  OOFEM_CLASS_WARNING1(_1)
#define _warning(_1) _warning1(_1)
//@}

/*
 * // log-family macros that use OOFEM loggers and exit
 *#define OOFEM_FATAL(...) {oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_FATAL, __FILE__, __LINE__, __VA_ARGS__); oofem_exit(1);}
 *#define OOFEM_ERROR(...) {oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_ERROR, __FILE__, __LINE__, __VA_ARGS__); oofem_exit(1);}
 *#define OOFEM_WARNING(...) {oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_WARNING, __FILE__, __LINE__, __VA_ARGS__);}
 *
 *
 * //log-family macros that allow to pass file and line info. They use OOFEM loggers and exit
 *#define __OOFEM_FATAL(_file,_line,...) {oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_FATAL,_file,_line,__VA_ARGS__);oofem_exit(1);}
 *#define __OOFEM_ERROR(_file,_line,...) {oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_ERROR,_file,_line,__VA_ARGS__);oofem_exit(1);}
 *#define __OOFEM_WARNING(_file,_line,...) {oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_WARNING,_file,_line,__VA_ARGS__);}
 *
 */

/**@
 * Log-family macros that use OOFEM loggers and exit (for errors and fatals).
 */
//@{
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

#define OOFEM_FATAL(_1) OOFEM_FATAL1(_1)
#define OOFEM_ERROR(_1) OOFEM_ERROR1(_1)
#define OOFEM_WARNING(_1) OOFEM_WARNING1(_1)
//@}

/**@
 * Log-family macros that allow to pass file and line info. They use OOFEM loggers and exit (fatals and errors).
 */
//@{
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

#define __OOFEM_FATAL(_file, _line, _1)  __OOFEM_FATAL1(_file, _line, _1)
#define __OOFEM_ERROR(_file, _line, _1) __OOFEM_ERROR1(_file, _line, _1)
#define __OOFEM_WARNING(_file, _line, _1) __OOFEM_WARNING1(_file, _line, _1)
//@}
} // end namespace oofem
#endif // error_h
