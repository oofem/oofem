/* $Header: /home/cvs/bp/oofem/oofemlib/src/clock.h,v 1.6 2003/05/19 13:03:57 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

/*
 * //   ********************
 * //   *** CLOCK MODULE ***
 * //   ********************
 */

#ifndef clock_h
#define clock_h

#include "compiler.h"

#ifndef _MSC_VER // If no Microsoft C

 #ifndef __MAKEDEPEND
  #ifdef TIME_WITH_SYS_TIME
   # include <sys/time.h>
   # include <time.h>
  #else
   # ifdef HAVE_SYS_TIME_H
    #  include <sys/time.h>
   # else
    #  include <time.h>
   # endif
  #endif
 #endif

// for getrusage - user time reporting
 #ifndef __MAKEDEPEND
  #include <sys/resource.h>
 #endif

/*
 * oofem timeval structure used to measure user time.
 * struct timeval
 * {
 *  __time_t tv_sec;            // Seconds.
 *  __time_t tv_usec;           // Microseconds.
 * }
 */
 #define OOFEM_USEC_LIM 1000000

#else // _MSC_VER active
 #ifndef __MAKEDEPEND
  #include <time.h>
 #endif

//oofem timeval structure used to measure user time.
struct timeval
{
    unsigned long tv_sec;          // Seconds.
    unsigned long tv_usec;         // Microseconds.
};

 #define OOFEM_USEC_LIM 1

#endif // end of _MSC_VER

namespace oofem {
typedef timeval oofem_timeval;

/**
 * Returns current time in seconds.
 * @return current time in oofem_timeval structure.
 */
void getTime(oofem_timeval &answer);

/**
 * Returns current time in seconds in time_t.
 */
time_t getTime();
/**
 * Function returns user time stored in oofem_timeval structure.
 */
void getUtime(oofem_timeval &answer);

/**
 * Function returning the elapsed user-time from given point, identified by oofem_timeval variable
 */
void getRelativeUtime(oofem_timeval &answer, oofem_timeval &from);
/**
 * Function returning the elapsed user-time startig at "from" and ending at "to", identified by oofem_timeval variables
 */
void getRelativeUtime(oofem_timeval &answer, oofem_timeval &from, oofem_timeval &to);
/**
 * Function converts total seconds into hours, minutes and remaining seconds
 */
void convertTS2HMS(int &nhrs, int &nmin, int &nsec, long int tsec);
void convertTS2HMS(int &nhrs, int &nmin, int &nsec, double tsec);
} // end namespace oofem
#endif // clock_h

