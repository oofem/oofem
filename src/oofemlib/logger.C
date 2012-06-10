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
#include "error.h"

#ifndef __MAKEDEPEND
 #include <cstdarg>
#endif

#ifdef __GNUC__
#include <cxxabi.h>
#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
// Taken from https://idlebox.net/2008/0901-stacktrace-demangled/ which indicated free usage.
/** Print a demangled stack backtrace of the caller function to FILE* out. */
static inline void print_stacktrace(FILE *out = stderr, unsigned int max_frames = 63)
{
    fprintf(out, "stack trace:\n");

    // storage array for stack trace address data
    void* addrlist[max_frames+1];

    // retrieve current stack addresses
    int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

    if (addrlen == 0) {
        fprintf(out, "  <empty, possibly corrupt>\n");
        return;
    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char** symbollist = backtrace_symbols(addrlist, addrlen);

    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
    char* funcname = (char*)malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < addrlen; i++)
    {
        char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

        // find parentheses and +address offset surrounding the mangled name:
        // ./module(function+0x15c) [0x8048a6d]
        for (char *p = symbollist[i]; *p; ++p)
        {
            if (*p == '(')
                begin_name = p;
            else if (*p == '+')
                begin_offset = p;
            else if (*p == ')' && begin_offset) {
                end_offset = p;
                break;
            }
        }

        if (begin_name && begin_offset && end_offset
            && begin_name < begin_offset)
        {
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';

            // mangled name is now in [begin_name, begin_offset) and caller
            // offset in [begin_offset, end_offset). now apply
            // __cxa_demangle():

            int status;
            char* ret = abi::__cxa_demangle(begin_name,
                                            funcname, &funcnamesize, &status);
            if (status == 0) {
                funcname = ret; // use possibly realloc()-ed string
                fprintf(out, "  %s : %s+%s\n",
                        symbollist[i], funcname, begin_offset);
            }
            else {
                // demangling failed. Output function name as a C function with
                // no arguments.
                fprintf(out, "  %s : %s()+%s\n",
                        symbollist[i], begin_name, begin_offset);
            }
        }
        else
        {
            // couldn't parse the line? print the whole line.
            fprintf(out, "  %s\n", symbollist[i]);
        }
    }

    free(funcname);
    free(symbollist);
}
#else
static inline void print_stacktrace(FILE *out = stderr, unsigned int max_frames = 63)
{
    fprintf(out, "No backtrace available\n");
}
#endif

namespace oofem {
#define LOG_ERR_HEADER "_______________________________________________________"
#define LOG_ERR_TAIL   "_______________________________________________________\a\n"

Logger :: Logger(logLevelType level, FILE *stream)
{
    this->logLevel = level;
    if ( stream ) {
        this->mylogStream = stream;
    } else {
        this->mylogStream = stdout;
    }

    this->closeFlag = false;
    numberOfWrn = numberOfErr = 0;
}

Logger :: ~Logger()
{
    if ( this->closeFlag ) {
        fclose(this->mylogStream);
    }
}

void
Logger :: appendlogTo(char *fname)
{
    FILE *stream = NULL;
    if ( this->closeFlag ) {
        stream = freopen(fname, "w", mylogStream);
    } else {
        stream = fopen(fname, "w");
    }

    if ( stream == NULL ) {
        OOFEM_WARNING2("Logger::appendlogTo : file opening error (%s)", fname);
    } else {
        mylogStream = stream;
    }

    this->closeFlag = true;
}

void
Logger :: writeLogMsg(logLevelType level, const char *format, ...)
{
    int rank = 0;

#ifdef __PARALLEL_MODE
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if ( rank == 0 )
    {
        va_list args;

        if ( level <= this->logLevel ) {
            va_start(args, format);
            vfprintf(mylogStream, format, args);
            va_end(args);
        }
    }

    if ( ( level == LOG_LEVEL_FATAL ) || ( level == LOG_LEVEL_ERROR ) ) {
        numberOfErr++;
    } else if ( level == LOG_LEVEL_WARNING ) {
        numberOfWrn++;
    }
}

void
Logger :: writeELogMsg(logLevelType level, const char *_file, int _line, const char *format, ...)
{
    va_list args;

    if  ( level <= this->logLevel ) {
        if ( _file ) {
            fprintf(mylogStream, "%s\n%s: (%s:%d)\n", LOG_ERR_HEADER, giveLevelName(level), _file, _line);
        } else {
            fprintf( mylogStream, "%s\n%s:\n", LOG_ERR_HEADER, giveLevelName(level) );
        }

        va_start(args, format);
        vfprintf(mylogStream, format, args);
        va_end(args);
        fprintf(mylogStream, "\n%s", LOG_ERR_TAIL);
    }

    if ( ( level == LOG_LEVEL_FATAL ) || ( level == LOG_LEVEL_ERROR ) ) {
        numberOfErr++;
        print_stacktrace(mylogStream, 10);
    } else if ( level == LOG_LEVEL_WARNING ) {
        numberOfWrn++;
    }
}

const char *
Logger :: giveLevelName(logLevelType l) const
{
    switch ( l ) {
    //case LOG_LEVEL_FATAL:
    case LOG_LEVEL_ERROR:
        return "Error";

    case LOG_LEVEL_WARNING:
        return "Warning";

    default:
        return "Info";
    }
}

void
Logger :: setLogLevel(int level)
{
    if ( ( level >= ( int ) LOG_LEVEL_FATAL ) && ( level <= ( int ) LOG_LEVEL_DEBUG ) ) {
        this->logLevel = ( logLevelType ) level;
    }
}


void
Logger :: printStatistics()
{
    int rank = 0;

#ifdef __PARALLEL_MODE
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    int totalNumberOfErr = numberOfErr, totalNumberOfWrn = numberOfWrn;
#ifdef __PARALLEL_MODE
    MPI_Reduce(& numberOfErr, & totalNumberOfErr, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(& numberOfWrn, & totalNumberOfWrn, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    if ( rank == 0 )
    {
        // force output
        fprintf(mylogStream, "Total %d error(s) and %d warning(s) reported\n", totalNumberOfErr, totalNumberOfWrn);
    }
}


#ifndef HAVE_MACRO_VA_ARGS


 #ifdef _MSC_VER

  #define __PROCESS_LOG \
    char buff [ MAX_ERROR_MSG_LENGTH ]; \
    va_list args; \
    va_start(args, format); \
    _vsnprintf(buff, MAX_ERROR_MSG_LENGTH, format, args); \
    va_end(args);

 #else

  #define __PROCESS_LOG \
    char buff [ MAX_ERROR_MSG_LENGTH ]; \
    va_list args; \
    va_start(args, format); \
    vsnprintf(buff, MAX_ERROR_MSG_LENGTH, format, args); \
    va_end(args);

 #endif

void LOG_FORCED_MSG(Logger &logger, const char *format, ...)
{
    __PROCESS_LOG;
    logger.writeLogMsg(Logger :: LOG_LEVEL_FORCED, buff);
}

void LOG_RELEVANT(Logger &logger, const char *format, ...)
{
    __PROCESS_LOG;
    logger.writeLogMsg(Logger :: LOG_LEVEL_RELEVANT, buff);
}


void LOG_INFO(Logger &logger, const char *format, ...)
{
    __PROCESS_LOG;
    logger.writeLogMsg(Logger :: LOG_LEVEL_INFO, buff);
}

void LOG_DEBUG(Logger &logger, const char *format, ...)
{
    __PROCESS_LOG;
    logger.writeLogMsg(Logger :: LOG_LEVEL_DEBUG, buff);
}

void OOFEM_LOG_RELEVANT(const char *format, ...)
{
    __PROCESS_LOG;
    oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_RELEVANT, buff);
}


void OOFEM_LOG_INFO(const char *format, ...)
{
    __PROCESS_LOG;
    oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_INFO, buff);
}

void OOFEM_LOG_DEBUG(const char *format, ...)
{
    __PROCESS_LOG;
    oofem_logger.writeLogMsg(Logger :: LOG_LEVEL_DEBUG, buff);
}

#endif
} // end namespace oofem
