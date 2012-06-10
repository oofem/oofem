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

#ifndef logger_h
#define logger_h

#ifdef __PARALLEL_MODE
 #ifndef __MAKEDEPEND
  #include <mpi.h>
 #endif
#endif

#include "oofemcfg.h"
#ifndef __MAKEDEPEND
 #include <cstdio>
#endif

namespace oofem {
class Logger
{
public:
    /// Type defining basic log levels.
    enum logLevelType {
        LOG_LEVEL_FORCED=-1,
        LOG_LEVEL_FATAL=0, LOG_LEVEL_ERROR=0,
        LOG_LEVEL_WARNING = 1,
        LOG_LEVEL_RELEVANT = 2,
        LOG_LEVEL_INFO = 3,
        LOG_LEVEL_ALL = 4, LOG_LEVEL_DEBUG = 4
    };
protected:
    /// Stream used for logging.
    FILE *mylogStream;
    /// flag indicating whether to close mylogStream.
    bool closeFlag;
    /// Current log level, messages with higher level are not reported.
    logLevelType logLevel;
    /// Counter of all warning and error messages.
    int numberOfWrn, numberOfErr;
public:
    Logger(logLevelType level, FILE *stream);
    ~Logger();
    /// Redirects log output to given file name (with path).
    void appendlogTo(char *fname);

    /// Writes the normal log message.
    void writeLogMsg(logLevelType level, const char *format, ...);
    /// Writes extended log message with file and line info.
    void writeELogMsg(logLevelType level, const char *_file, int _line, const char *format, ...);
    /// Flushes the log stream.
    void flush() { fflush(mylogStream); }

    /// Sets log level to given one. Only log messages with level less or equal given threshold will be printed.
    void setLogLevel(logLevelType level) { logLevel = level; }
    /// Sets log level to given one. Only log messages with level less or equal given threshold will be printed.
    void setLogLevel(int level);
    /// Prints number of errors and warning logged.
    void printStatistics();

protected:
    const char *giveLevelName(logLevelType l) const;
};

extern Logger oofem_logger;
extern Logger oofem_errLogger;


/**@
 * General log-family macros (those for error and warning reporting).
 * These macros add to given message also file and line information.
 * They can be implemented in a very elegant way using macro with variable
 * number of arguments (__VA_ARGS__). But since the many compilers do not support
 * macros with variable number of arguments we have two choices:
 * to implement this as a function with variable number of arguments, but then we can not
 * add file and line info via __FILE__ and __LINE__ macros. Or if file and line
 * info is preferred, then instead of single macro with variable number of arguments we can have
 * series of "classical" macros with increasing number of parameters.
 * The latter approach is used here.
 */
//@{
#define LOG_FATAL1(logger, _1) logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, _1)
#define LOG_FATAL2(logger, _1, _2) logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, _1, _2)
#define LOG_FATAL3(logger, _1, _2, _3) logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, _1, _2, _3)
#define LOG_ERROR1(logger, _1) logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1)
#define LOG_ERROR2(logger, _1, _2) logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1, _2)
#define LOG_ERROR3(logger, _1, _2, _3) logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1, _2, _3)
#define LOG_WARNING1(logger, _1) logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1)
#define LOG_WARNING2(logger, _1, _2) logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1, _2)
#define LOG_WARNING3(logger, _1, _2, _3) logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1, _2, _3)

#define LOG_FATAL(logger, _1) LOG_FATAL1(logger, _1)
#define LOG_ERROR(logger, _1) LOG_ERROR1(logger, _1)
#define LOG_WARNING(logger, _1) LOG_WARNING1(logger, _1)
//@}


/**@
 * log-family macros that use default OOFEM loggers
 */
//@{
#define OOFEM_LOG_FATAL1(_1) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, _1)
#define OOFEM_LOG_FATAL2(_1, _2) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, _1, _2)
#define OOFEM_LOG_FATAL3(_1, _2, _3) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, _1, _2, _3)
#define OOFEM_LOG_ERROR1(_1) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1)
#define OOFEM_LOG_ERROR2(_1, _2) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1, _2)
#define OOFEM_LOG_ERROR3(_1, _2, _3) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, _1, _2, _3)
#define OOFEM_LOG_WARNING1(_1) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1)
#define OOFEM_LOG_WARNING2(_1, _2) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1, _2)
#define OOFEM_LOG_WARNING3(_1, _2, _3) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, _1, _2, _3)

#define OOFEM_LOG_FATAL(_1) OOFEM_LOG_FATAL1(_1)
#define OOFEM_LOG_ERROR(_1) OOFEM_LOG_ERROR1(_1)
#define OOFEM_LOG_WARNING(_1) OOFEM_LOG_WARNING1(_1)
//@}

#ifdef HAVE_MACRO_VA_ARGS
/// used internally
/**@
 * Log reporting macros (those that do add file and line info)
 */
//@{
 #define __LOG_E_MESSAGE(logger, level, _file, _line, ...) logger.writeELogMsg(level, _file, _line, __VA_ARGS__)
 #define LOG_FORCED_MSG(logger, ...) logger.writeLogMsg(Logger :: LOG_LEVEL_FORCED, __VA_ARGS__)
/*
 * // General log-family macros
 * // Alternative definition using VARARGS
 *
 *#define LOG_FATAL(logger, ...) logger.writeELogMsg(Logger::LOG_LEVEL_FATAL, __FILE__, __LINE__, __VA_ARGS__)
 *#define LOG_ERROR(logger, ...) logger.writeELogMsg(Logger::LOG_LEVEL_ERROR, __FILE__, __LINE__, __VA_ARGS__)
 *#define LOG_WARNING(logger, ...) logger.writeELogMsg(Logger::LOG_LEVEL_WARNING, __FILE__, __LINE__, __VA_ARGS__)
 */
 #define LOG_RELEVANT(logger, ...) logger.writeLogMsg(Logger :: LOG_LEVEL_RELEVANT, __VA_ARGS__)
 #define LOG_INFO(logger, ...) logger.writeLogMsg(Logger :: LOG_LEVEL_INFO, __VA_ARGS__)
 #define LOG_DEBUG(logger, ...) logger.writeLogMsg(Logger :: LOG_LEVEL_DEBUG, __VA_ARGS__)

/* log-family macros that use default OOFEM loggers */
/*
 * #define OOFEM_LOG_FATAL(...) oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_FATAL, __FILE__, __LINE__, __VA_ARGS__)
 * #define OOFEM_LOG_ERROR(...) oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_ERROR, __FILE__, __LINE__, __VA_ARGS__)
 * #define OOFEM_LOG_WARNING(...) oofem_errLogger.writeELogMsg(Logger::LOG_LEVEL_WARNING, __FILE__, __LINE__, __VA_ARGS__)
 */
 #define OOFEM_LOG_RELEVANT(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_RELEVANT, __VA_ARGS__)
 #define OOFEM_LOG_INFO(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_INFO, __VA_ARGS__)
 #define OOFEM_LOG_DEBUG(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_DEBUG, __VA_ARGS__)
//@}
#else

void LOG_FORCED_MSG(Logger &logger, const char *format, ...);
void LOG_RELEVANT(Logger &logger, const char *format, ...);
void LOG_INFO(Logger &logger, const char *format, ...);
void LOG_DEBUG(Logger &logger, const char *format, ...);

void OOFEM_LOG_RELEVANT(const char *format, ...);
void OOFEM_LOG_INFO(const char *format, ...);
void OOFEM_LOG_DEBUG(const char *format, ...);

#endif
} // end namespace oofem
#endif // logger_h
