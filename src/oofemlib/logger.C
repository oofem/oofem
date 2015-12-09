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

#include "logger.h"
#include "error.h"

#include <cstdarg>
#ifdef __PARALLEL_MODE
 #include <mpi.h>
#endif


#if defined ( __GNUC__ ) && defined ( HAVE_EXECINFO_H )
 #include <cxxabi.h>
 #include <execinfo.h>
 #include <cstdio>
 #include <cstdlib>
// Taken from https://idlebox.net/2008/0901-stacktrace-demangled/ which indicated free usage.
/** Print a demangled stack backtrace of the caller function to FILE* out. */
static inline void print_stacktrace(FILE *out = stderr, unsigned int max_frames = 63)
{
    int addrlen = 0;
    fprintf(out, "stack trace:\n");

    // storage array for stack trace address data
    void *addrlist [ max_frames + 1 ];

    // retrieve current stack addresses
    addrlen = backtrace( addrlist, sizeof( addrlist ) / sizeof( void * ) );
    if ( addrlen == 0 ) {
        fprintf(out, "  <empty, possibly corrupt>\n");
        return;
    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char **symbollist;
    symbollist = backtrace_symbols(addrlist, addrlen);
    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
    char *funcname = ( char * ) malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for ( int i = 2; i < addrlen; i++ ) {
        char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

        // find parentheses and +address offset surrounding the mangled name:
        // ./module(function+0x15c) [0x8048a6d]
        for ( char *p = symbollist [ i ]; * p; ++p ) {
            if ( * p == '(' ) {
                begin_name = p;
            } else if ( * p == '+' ) {
                begin_offset = p;
            } else if ( * p == ')' && begin_offset ) {
                end_offset = p;
                break;
            }
        }

        if ( begin_name && begin_offset && end_offset &&
             begin_name < begin_offset ) {
            * begin_name++ = '\0';
            * begin_offset++ = '\0';
            * end_offset = '\0';

            // mangled name is now in [begin_name, begin_offset) and caller
            // offset in [begin_offset, end_offset). now apply
            // __cxa_demangle():

            int status;
            char *ret = abi :: __cxa_demangle(begin_name,
                                              funcname, & funcnamesize, & status);
            if ( status == 0 ) {
                funcname = ret; // use possibly realloc()-ed string
                fprintf(out, "  %s : %s+%s\n",
                        symbollist [ i ], funcname, begin_offset);
            } else {
                // demangling failed. Output function name as a C function with
                // no arguments.
                fprintf(out, "  %s : %s()+%s\n",
                        symbollist [ i ], begin_name, begin_offset);
            }
        } else {
            // couldn't parse the line? print the whole line.
            fprintf(out, "  %s\n", symbollist [ i ]);
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

// Default log output
Logger oofem_logger(Logger :: LOG_LEVEL_INFO);

Logger :: Logger(logLevelType level)
{
    this->logLevel = level;
    this->logStream = stdout;
    this->errStream = stderr;

    this->closeFlag = false;
    this->errCloseFlag = false;
    numberOfWrn = numberOfErr = 0;
}

Logger :: ~Logger()
{
    if ( this->closeFlag ) {
        fclose(this->logStream);
    }
    if ( this->errCloseFlag ) {
        fclose(this->errStream);
    }
}

void
Logger :: appendLogTo(const std :: string &fname)
{
    FILE *stream = NULL;
    if ( this->closeFlag ) {
        stream = freopen(fname.c_str(), "a", this->logStream);
    } else {
        stream = fopen(fname.c_str(), "a");
    }

    if ( stream == NULL ) {
        OOFEM_WARNING( "file opening error (%s)", fname.c_str() );
    } else {
        this->logStream = stream;
    }

    this->closeFlag = true;
}

void
Logger :: appendErrorTo(const std :: string &fname)
{
    FILE *stream = NULL;
    if ( this->errCloseFlag ) {
        stream = freopen(fname.c_str(), "a", this->errStream);
    } else {
        stream = fopen(fname.c_str(), "a");
    }

    if ( stream == NULL ) {
        OOFEM_WARNING( "file opening error (%s)", fname.c_str() );
    } else {
        this->errStream = stream;
    }

    this->errCloseFlag = true;
}

void
Logger :: appendLogTo(FILE *stream)
{
    if ( this->closeFlag ) {
      fclose (this->logStream);
    }

    if ( stream == NULL ) {
        OOFEM_ERROR( "Logger::appendLogTo : null stream given" );
    } else {
        this->logStream = stream;
    }

    this->closeFlag = false;
}

void
Logger :: appendErrorTo(FILE *stream)
{
    if ( this->errCloseFlag ) {
        fclose (this->errStream);
    }

    if ( stream == NULL ) {
        OOFEM_ERROR( "Logger::appendLogTo : null stream given" );
    } else {
        this->errStream = stream;
    }

    this->errCloseFlag = false;
}



void
Logger :: writeLogMsg(logLevelType level, const char *format, ...)
{
    int rank = 0;

#ifdef __PARALLEL_MODE
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
#endif
    (void)rank;//prevent a warning about unused variable
    FILE *stream = this->logStream;
    if ( level == LOG_LEVEL_FATAL || level == LOG_LEVEL_ERROR ) {
        numberOfErr++;
        stream = this->errStream;
    } else if ( level == LOG_LEVEL_WARNING ) {
        numberOfWrn++;
        stream = this->errStream;
    }


    //    if ( rank == 0 ) {
    if (1) {
        va_list args;

        if ( level <= this->logLevel ) {
            va_start(args, format);
            vfprintf(stream, format, args);
            va_end(args);
        }
    }
}

void
Logger :: writeELogMsg(logLevelType level, const char *_func, const char *_file, int _line, const char *format, ...)
{
    va_list args;

    FILE *stream = this->logStream;
    if ( level == LOG_LEVEL_FATAL || level == LOG_LEVEL_ERROR ) {
        numberOfErr++;
        stream = this->errStream;
    } else if ( level == LOG_LEVEL_WARNING ) {
        numberOfWrn++;
        stream = this->errStream;
    }

    if  ( level <= this->logLevel ) {
        if ( _file ) {
            fprintf(stream, "%s\n%s: (%s:%d)\n", LOG_ERR_HEADER, giveLevelName(level), _file, _line);
        } else {
            fprintf(stream, "%s\n%s:\n", LOG_ERR_HEADER, giveLevelName(level) );
        }
        if ( _func ) {
            fprintf(stream, "In %s:\n", _func );
        }

        va_start(args, format);
        vfprintf(stream, format, args);
        va_end(args);
        fprintf(stream, "\n%s", LOG_ERR_TAIL);
    }

    if ( level == LOG_LEVEL_FATAL || level == LOG_LEVEL_ERROR ) {
        print_stacktrace(this->errStream, 10);
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
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
#endif

    int totalNumberOfErr = numberOfErr, totalNumberOfWrn = numberOfWrn;
#ifdef __PARALLEL_MODE
    MPI_Reduce(& numberOfErr, & totalNumberOfErr, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(& numberOfWrn, & totalNumberOfWrn, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    if ( rank == 0 ) {
        // force output
        fprintf(logStream, "Total %d error(s) and %d warning(s) reported\n", totalNumberOfErr, totalNumberOfWrn);
    }
}

} // end namespace oofem
