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
#include "util.h"

#include <cstdarg>
#ifdef __PARALLEL_MODE
 #include <mpi.h>
#endif


namespace oofem {
#define LOG_ERR_HEADER "_______________________________________________________"
#define LOG_ERR_TAIL   "_______________________________________________________\a\n"

// Default log output
Logger oofem_logger(Logger :: LOG_LEVEL_INFO);

Logger :: Logger(logLevelType level) :
    logStream(stdout),
    errStream(stderr),
    closeFlag(false),
    errCloseFlag(false),
    logLevel(level),
    numberOfWrn(0),
    numberOfErr(0)
#ifdef __PARALLEL_MODE
    ,comm(MPI_COMM_SELF)
#endif
{}

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
    MPI_Comm_rank(this->comm, & rank);
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
            fprintf(stream, "%s\n%s:", LOG_ERR_HEADER, giveLevelName(level));
        } else {
            fprintf(stream, "%s\n%s:", LOG_ERR_HEADER, giveLevelName(level));
        }

        va_start(args, format);
        vfprintf(stream, format, args);
        va_end(args);
        fprintf(stream, "\n");    

        if ( _func ) {
            fprintf(stream, "In %s ", _func );
        }

        if ( _file ) {
            fprintf(stream, "(%s:%d)", _file, _line);
        } 

        fprintf(stream, "\n%s", LOG_ERR_TAIL);
    }

    if ( level == LOG_LEVEL_FATAL || level == LOG_LEVEL_ERROR ) {
#ifndef CEMPY
        print_stacktrace(this->errStream, 10);
#endif
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

#ifdef __PARALLEL_MODE
void
Logger :: setComm(MPI_Comm comm)
{
    this->comm = comm;
}
#endif

void
Logger :: printStatistics()
{
    int rank = 0;

#ifdef __PARALLEL_MODE
    MPI_Comm_rank(this->comm, & rank);
#endif

    int totalNumberOfErr = numberOfErr, totalNumberOfWrn = numberOfWrn;
#ifdef __PARALLEL_MODE
    MPI_Reduce(& numberOfErr, & totalNumberOfErr, 1, MPI_INT, MPI_SUM, 0, this->comm);
    MPI_Reduce(& numberOfWrn, & totalNumberOfWrn, 1, MPI_INT, MPI_SUM, 0, this->comm);
#endif
    if ( rank == 0 ) {
        // force output
        fprintf(logStream, "Total %d error(s) and %d warning(s) reported\n", totalNumberOfErr, totalNumberOfWrn);
    }
}

} // end namespace oofem
