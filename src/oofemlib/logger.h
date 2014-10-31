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

#ifndef logger_h
#define logger_h

#include "oofemcfg.h"

#include <cstdio>
#include <string>

// MSVC doesn't properly implement C99. (might need to wrap __func__ behind a macro to support all platforms correctly(?))
#ifdef _MSC_VER
# define __func__ __FUNCTION__
#endif

namespace oofem {
/**
 * Logger class used by OOFEM to print information during analysis.
 * Prints warnings and errors into a separate stream from normal output.
 */
class OOFEM_EXPORT Logger
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
    FILE *logStream, *errStream;
    /// flag indicating whether to close mylogStream.
    bool closeFlag, errCloseFlag;
    /// Current log level, messages with higher level are not reported.
    logLevelType logLevel;
    /// Counter of all warning and error messages.
    int numberOfWrn, numberOfErr;
public:
    Logger(logLevelType level);
    ~Logger();
    /// Redirects log output to given file name (with path).
    void appendLogTo(const std :: string &fname);
    /// Redirects error output to given file name (with path).
    void appendErrorTo(const std :: string &fname);
    /// Redirects log output to given stream.
    void appendLogTo(FILE* stream);
    /// Redirects error output to given stream.
    void appendErrorTo(FILE* stream);


    /// Writes the normal log message.
    void writeLogMsg(logLevelType level, const char *format, ...);
    /// Writes extended log message with file and line info.
    void writeELogMsg(logLevelType level, const char *_func, const char *_file, int _line, const char *format, ...);
    /// Flushes the log stream.
    void flush() { fflush(logStream); fflush(errStream); }

    /// Sets log level to given one. Only log messages with level less or equal given threshold will be printed.
    void setLogLevel(logLevelType level) { logLevel = level; }
    /// Sets log level to given one. Only log messages with level less or equal given threshold will be printed.
    void setLogLevel(int level);
    /// Prints number of errors and warning logged.
    void printStatistics();

protected:
    const char *giveLevelName(logLevelType l) const;
};

extern OOFEM_EXPORT Logger oofem_logger;

/**
 * Log reporting macros
 */
//@{
#define OOFEM_LOG_FATAL(...) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __func__, __FILE__, __LINE__, __VA_ARGS__)
#define OOFEM_LOG_ERROR(...) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __func__, __FILE__, __LINE__, __VA_ARGS__)
#define OOFEM_LOG_WARNING(...) oofem_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING,  __func__, __FILE__, __LINE__, __VA_ARGS__)

#define OOFEM_LOG_FORCED(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_FORCED, __VA_ARGS__)
#define OOFEM_LOG_RELEVANT(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_RELEVANT, __VA_ARGS__)
#define OOFEM_LOG_INFO(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_INFO, __VA_ARGS__)
#define OOFEM_LOG_DEBUG(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_DEBUG, __VA_ARGS__)
//@}
} // end namespace oofem
#endif // logger_h
