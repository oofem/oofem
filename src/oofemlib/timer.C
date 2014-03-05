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

#include "timer.h"

#include <cstdio>

#ifndef _WIN32 //_MSC_VER and __MINGW32__ included
//for getrusage - user time reporting
 #include <sys/resource.h>
#else
 #include <ctime>
#endif

namespace oofem {
void Timer :: getUtime(std :: chrono :: duration< double > &answer)
{
#ifdef _WIN32 //_MSC_VER and __MINGW32__ included
    clock_t utime = clock();
    answer = std :: chrono :: seconds(utime / CLOCKS_PER_SEC);
#else
    struct rusage rsg;
    getrusage(RUSAGE_SELF, & rsg);
    answer = std :: chrono :: seconds(rsg.ru_utime.tv_sec) + std :: chrono :: microseconds(rsg.ru_utime.tv_usec);
#endif
}

void Timer :: getTime(std :: chrono :: time_point< std :: chrono :: high_resolution_clock > &answer)
{
    answer = std :: chrono :: high_resolution_clock :: now();
}

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
    elapsedWTime = elapsedWTime.zero();
    elapsedUTime = elapsedUTime.zero();
    running = false;
}

double Timer :: getUtime()
{
    this->updateElapsedTime();
    return elapsedUTime.count();
}

double Timer :: getWtime()
{
    updateElapsedTime();
    return elapsedWTime.count();
}

void Timer :: convert2HMS(int &nhrs, int &nmin, int &nsec, double tsec)
{
    long int _nsec = ( long int ) tsec;
    nhrs = 0;
    nmin = 0;
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

void Timer :: toString(char *buff)
{
    std :: sprintf( buff, "ut: %f.3s, wt: %f.3s", getUtime(), getWtime() );
}

void Timer :: updateElapsedTime()
{
    if ( running ) {
        pauseTimer();
        resumeTimer();
    }

    elapsedWTime += end_wtime - start_wtime;
    elapsedUTime += end_utime - start_utime;

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

void EngngModelTimer :: convert2HMS(int &nhrs, int &nmin, int &nsec, double tsec)
{
    Timer :: convert2HMS(nhrs, nmin, nsec, tsec);
}

void EngngModelTimer :: toString(EngngModelTimer :: EngngModelTimerType t, char *buff)
{
    return timers [ t ].toString(buff);
}
}
