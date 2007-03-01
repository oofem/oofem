/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/logger.C,v 1.1.2.3 2004/04/16 13:04:05 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2003   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. 
*/

#include "logger.h"
#include "error.h"
#include "cltypes.h"

#ifndef __MAKEDEPEND
#include <stdarg.h>
#endif

#define LOG_ERR_HEADER "_______________________________________________________"
#define LOG_ERR_TAIL   "_______________________________________________________\a\n"

Logger :: Logger (logLevelType level, FILE* stream)
{
  this->logLevel = level;
  if (stream) this->mylogStream = stream;
  else this->mylogStream = stdout;
  this->closeFlag = false;
  numberOfWrn = numberOfErr = 0;
}

Logger::~Logger () 
{
  if (this->closeFlag) fclose (this->mylogStream);
}

void
Logger :: appendlogTo (char* fname)
{
  FILE *stream=NULL;
  if (this->closeFlag) 
    freopen (fname, "w", stream);
  else 
    stream = fopen (fname, "w");
  
  if (stream == NULL) {
    OOFEM_WARNING2 ("Logger::appendlogTo : file opening error (%s)", fname);
  } else mylogStream = stream;
  this->closeFlag = true;
}

void
Logger :: writeLogMsg (logLevelType level, char *format, ...)
{
  va_list args;
  
  if (level <= this->logLevel) {
    va_start(args, format);
    vfprintf(mylogStream, format, args);
    va_end(args);
  }

  if ((level == LOG_LEVEL_FATAL) || (level == LOG_LEVEL_ERROR))
    numberOfErr++;
  else if (level == LOG_LEVEL_WARNING) numberOfWrn++;
}


void
Logger :: writeELogMsg (logLevelType level, const char* _file, int _line, char *format, ...)
{
  va_list args;
  
  if  (level <= this->logLevel) {
    if (_file)
      fprintf (mylogStream, "%s\n%s: (%s:%d)\n", LOG_ERR_HEADER, giveLevelName(level), _file, _line);
    else
      fprintf (mylogStream, "%s\n%s:\n", LOG_ERR_HEADER, giveLevelName(level));

    va_start(args, format);
    vfprintf(mylogStream, format, args);
    va_end(args);
    fprintf (mylogStream, "\n%s", LOG_ERR_TAIL);
  }

  if ((level == LOG_LEVEL_FATAL) || (level == LOG_LEVEL_ERROR))
    numberOfErr++;
  else if (level == LOG_LEVEL_WARNING) numberOfWrn++;
}

const char*
Logger :: giveLevelName (logLevelType l) const 
{
  switch (l) {
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
Logger::setLogLevel (int level) 
{
  if ((level >= (int)LOG_LEVEL_FATAL) && (level <= (int)LOG_LEVEL_DEBUG))
    this->logLevel = (logLevelType) level;
}


void 
Logger :: printStatistics () 
{
  // force output
  fprintf(mylogStream, "Total %d error(s) and %d warning(s) reported\n", numberOfErr, numberOfWrn);
}


#ifndef HAVE_MACRO_VA_ARGS


#ifdef _MSC_VER 

#define __PROCESS_LOG \
char buff[MAX_ERROR_MSG_LENGTH];\
va_list args;\
va_start(args, format);\
_vsnprintf(buff, MAX_ERROR_MSG_LENGTH, format, args);\
va_end(args);

#else

#define __PROCESS_LOG \
char buff[MAX_ERROR_MSG_LENGTH];\
va_list args;\
va_start(args, format);\
vsnprintf(buff, MAX_ERROR_MSG_LENGTH, format, args);\
va_end(args);

#endif

void LOG_FORCED_MSG(Logger& logger, const char *format, ...)
{
  __PROCESS_LOG;
  logger.writeLogMsg(Logger::LOG_LEVEL_FORCED, buff);
}

void LOG_RELEVANT(Logger& logger, const char *format, ...)
{
  __PROCESS_LOG;
  logger.writeLogMsg(Logger::LOG_LEVEL_RELEVANT, buff);
}


void LOG_INFO(Logger& logger, const char *format, ...)
{
  __PROCESS_LOG;
  logger.writeLogMsg(Logger::LOG_LEVEL_INFO, buff);
}

void LOG_DEBUG(Logger& logger, const char *format, ...)
{
  __PROCESS_LOG;
  logger.writeLogMsg(Logger::LOG_LEVEL_DEBUG, buff);
}

void OOFEM_LOG_RELEVANT(const char *format, ...)
{
  __PROCESS_LOG;
  oofem_logger.writeLogMsg(Logger::LOG_LEVEL_RELEVANT, buff);
}


void OOFEM_LOG_INFO(const char *format, ...)
{
  __PROCESS_LOG;
  oofem_logger.writeLogMsg(Logger::LOG_LEVEL_INFO, buff);
}

void OOFEM_LOG_DEBUG(const char *format, ...)
{
  __PROCESS_LOG;
  oofem_logger.writeLogMsg(Logger::LOG_LEVEL_DEBUG, buff);
}


#endif
