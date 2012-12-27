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

#include "timer.h"
#include "compiler.h"
#include <cstdio>

#ifndef _MSC_VER
// for getrusage - user time reporting
#include <sys/resource.h>
#endif

namespace oofem {

#ifndef _MSC_VER

void Timer :: getUtime(oofem_timeval &answer)
{
    struct rusage rsg;
    getrusage(RUSAGE_SELF, & rsg);
    answer.tv_sec  = rsg.ru_utime.tv_sec;
    answer.tv_usec = rsg.ru_utime.tv_usec;
}

void Timer :: getTime(oofem_timeval &answer)
{
    gettimeofday(&answer, NULL);
}

#else // #ifndef _MSC_VER

void Timer :: getUtime(oofem_timeval &answer)
{
    clock_t utime = clock();
    answer.tv_sec = utime / CLOCKS_PER_SEC;
    answer.tv_usec = 0;
}

void Timer :: getTime(oofem_timeval &answer)
{
    time_t t;
    t = time(NULL);
    answer.tv_sec = ( unsigned long ) t;
    answer.tv_usec = 0;
}

#endif


Timer :: Timer()
{
    initTimer();
}

void Timer :: startTimer()
{
    this->initTimer();
    this->getTime(start_wtime);
    this->getUtime(start_utime);
    running = true;
}

void Timer :: stopTimer()
{
    this->pauseTimer();
    running = false;
}

void Timer :: pauseTimer()
{
    this->getTime(end_wtime);
    this->getUtime(end_utime);
    running = false;
    this->updateElapsedTime();
}

void Timer :: resumeTimer()
{
    this->getTime(start_wtime);
    this->getUtime(start_utime);
    running = true;
}

void Timer :: initTimer()
{
    elapsedWTime.tv_sec = elapsedWTime.tv_usec = elapsedUTime.tv_sec = elapsedUTime.tv_usec = 0;
    running = false;
}

double Timer :: getUtime()
{
    this->updateElapsedTime();
    return ( double ) elapsedUTime.tv_sec + ( double ) elapsedUTime.tv_usec * 1.e-6;
}

double Timer :: getWtime()
{
    updateElapsedTime();
    return ( double ) elapsedWTime.tv_sec + ( double ) elapsedWTime.tv_usec * 1.e-6;
}

void Timer :: convert2HMS(int& nhrs, int& nmin, int& nsec, long int tsec)
{
    long int _nsec = tsec;
    if ( _nsec > 60 ) {
        nmin = _nsec / 60;
        _nsec %= 60;
    }
    
    if ( nmin > 60 ) {
        nhrs = nmin / 60;
        nmin %= 60;
    }
    
    nsec = _nsec;
}

void Timer :: convert2HMS(int& nhrs, int& nmin, int& nsec, double tsec)
{
    Timer  ::  convert2HMS(nhrs, nmin, nsec, (long int)tsec);
}

void Timer :: toString(char* buff)
{
    sprintf( buff, "ut: %f.3s, wt: %f.3s", getUtime(), getWtime() );
}

void Timer :: updateElapsedTime()
{
    if ( running ) {
        pauseTimer();
        resumeTimer();
    }

#ifndef _MSC_VER
    oofem_timeval etime;
    timersub(& end_wtime, & start_wtime, & etime);
    timeradd(& etime, & elapsedWTime, & elapsedWTime);

    timersub(& end_utime, & start_utime, & etime);
    timeradd(& etime, & elapsedUTime, & elapsedUTime);
#endif

    start_utime = end_utime;
    start_wtime = end_wtime;
}

double EngngModelTimer :: getUtime(EngngModelTimer :: EngngModelTimerType t)
{
    return timers [ t ].getUtime();
}

double EngngModelTimer :: getWtime(EngngModelTimer :: EngngModelTimerType t)
{
    return timers [ t ].getWtime();
}

void EngngModelTimer :: convert2HMS(int& nhrs, int& nmin, int& nsec, long int tsec) const
{
    Timer  ::  convert2HMS(nhrs, nmin, nsec, tsec);
}

void EngngModelTimer :: convert2HMS(int& nhrs, int& nmin, int& nsec, double tsec) const
{
    Timer  ::  convert2HMS(nhrs, nmin, nsec, tsec);
}

void EngngModelTimer :: toString(EngngModelTimer :: EngngModelTimerType t, char* buff)
{
    return timers [ t ].toString(buff);
}
    
}

