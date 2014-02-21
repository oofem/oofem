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

#ifdef __PARALLEL_MODE
 #include <mpi.h>
#endif

#include <cstdio>
#include <string>

namespace oofem {
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
    FILE *mylogStream;
    /// flag indicating whether to close mylogStream.
    bool closeFlag;
    /// Current log level, messages with higher level are not reported.
    logLevelType logLevel;
    /// Counter of all warning and error messages.
    int numberOfWrn, numberOfErr;
public:
    Logger(logLevelType level, FILE * stream);
    ~Logger();
    /// Redirects log output to given file name (with path).
    void appendlogTo(const std :: string &fname);

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

extern OOFEM_EXPORT Logger oofem_logger;
extern OOFEM_EXPORT Logger oofem_errLogger;

/**
 * Log reporting macros
 */
//@{
#ifdef HAVE_MACRO_VA_ARGS
#define OOFEM_LOG_FATAL(...) oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __FILE__, __LINE__, __VA_ARGS__)
#define OOFEM_LOG_ERROR(...) oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __FILE__, __LINE__, __VA_ARGS__)
#define OOFEM_LOG_WARNING(...) oofem_errLogger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __FILE__, __LINE__, __VA_ARGS__)

#define OOFEM_LOG_FORCED(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_FORCED, __VA_ARGS__)
#define OOFEM_LOG_RELEVANT(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_RELEVANT, __VA_ARGS__)
#define OOFEM_LOG_INFO(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_INFO, __VA_ARGS__)
#define OOFEM_LOG_DEBUG(...) oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_DEBUG, __VA_ARGS__)
#else
void OOFEM_LOG_FATAL(const char *format, ...);
void OOFEM_LOG_ERROR(const char *format, ...);
void OOFEM_LOG_WARNING(const char *format, ...);
void OOFEM_LOG_FORCED(const char *format, ...);
void OOFEM_LOG_RELEVANT(const char *format, ...);
void OOFEM_LOG_INFO(const char *format, ...);
void OOFEM_LOG_DEBUG(const char *format, ...);
#endif

//@}
} // end namespace oofem
#endif // logger_h
