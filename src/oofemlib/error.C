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

#include "logger.h"

#ifndef __MAKEDEPEND
 #include <cstdio>
 #include <cstdlib>
#endif

namespace oofem {
/**
 * Global variable containing warning level, should be in interval (0,3).
 * Zero level suppress all warning messages, level 3 causes to report all
 * warning messages. Default level is set to maximum level.
 */
//int oofem_warningLevel = 3;


void oofem_exit(int code)
{
    oofem_errLogger.printStatistics();
    fprintf(stderr, "oofem exit code %d\n", code);
    exit(code);
}

/*
 * #ifndef HAVE_MACRO_VA_ARGS
 *
 *#ifndef __MAKEDEPEND
 *#include <stdarg.h>
 *#endif
 *
 *
 *
 *#define __PROCESS_LOG \
 * char buff[MAX_ERROR_MSG_LENGTH]; \
 * va_list args; \
 * va_start(args, format); \
 * vsnprintf(buff, MAX_ERROR_MSG_LENGTH, format, args); \
 * va_end(args);
 *
 * void _error(const char *format, ...)
 * {
 * __PROCESS_LOG;
 * oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_ERROR, NULL, 0,buff);
 * exit (1);
 * }
 * void _warning(const char *format, ...)
 * {
 * __PROCESS_LOG;
 * oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_WARNING, NULL ,0,buff);
 * }
 *
 *
 * void OOFEM_FATAL(const char *format, ...)
 * {
 * __PROCESS_LOG;
 * oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_FATAL, NULL ,0,buff);
 * exit (1);
 * }
 *
 * void OOFEM_ERROR(const char *format, ...)
 * {
 * __PROCESS_LOG;
 * oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_ERROR, NULL,0,buff);
 * exit (1);
 *
 * }
 *
 * void OOFEM_WARNING(const char *format, ...)
 * {
 * __PROCESS_LOG;
 * oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_WARNING, NULL,0,buff);
 * }
 *
 * void  __OOFEM_FATAL(const char* _file,int _line,const char *format, ...)
 * {
 * __PROCESS_LOG;
 * oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_FATAL, _file,_line,buff);
 * exit (1);
 * }
 *
 * void __OOFEM_ERROR(const char* _file,int _line,const char *format, ...)
 * {
 * __PROCESS_LOG;
 * oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_ERROR, _file,_line,buff);
 * exit (1);
 * }
 * void __OOFEM_WARNING(const char* _file,int _line,const char *format, ...)
 * {
 * __PROCESS_LOG;
 * oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_WARNING, _file,_line,buff);
 * }
 *
 *#endif
 */
} // end namespace oofem
