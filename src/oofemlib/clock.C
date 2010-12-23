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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

#include "clock.h"
#include "compiler.h"
#ifndef _MSC_VER

namespace oofem {
void getUtime(oofem_timeval &answer)
{
    struct rusage rsg;
    getrusage(RUSAGE_SELF, & rsg);
    answer.tv_sec  = rsg.ru_utime.tv_sec;
    answer.tv_usec = rsg.ru_utime.tv_usec;
}

void getRelativeUtime(oofem_timeval &answer, oofem_timeval &from, oofem_timeval &to)
{
    if ( to.tv_usec < from.tv_usec ) {
        answer.tv_usec = ( OOFEM_USEC_LIM - from.tv_usec ) + to.tv_usec;
        answer.tv_sec  = to.tv_sec - from.tv_sec - 1;
    } else {
        answer.tv_usec = to.tv_usec - from.tv_usec;
        answer.tv_sec  = to.tv_sec - from.tv_sec;
    }
}


void getRelativeUtime(oofem_timeval &answer, oofem_timeval &from)
{
    oofem_timeval to;
    getUtime(to);

    getRelativeUtime(answer, from, to);
}

void getTime(oofem_timeval &answer)
{
    gettimeofday(& answer, NULL);
}


time_t getTime()
{
    time_t t;
    t = time(NULL);
    return t;
}

#else // #ifndef _MSC_VER

namespace oofem {
void getUtime(oofem_timeval &answer)
{
    clock_t utime = clock();
    answer.tv_sec = utime / CLOCKS_PER_SEC;
    answer.tv_usec = 0;
}

void getRelativeUtime(oofem_timeval &answer, oofem_timeval &from, oofem_timeval &to)
{
    clock_t utime = clock();

    answer.tv_sec = to.tv_sec - from.tv_sec;
    answer.tv_usec = 0;
}

void getRelativeUtime(oofem_timeval &answer, oofem_timeval &from)
{
    oofem_timeval utime;
    getUtime(utime);
    getRelativeUtime(answer, from, utime);
}

void getTime(oofem_timeval &answer)
{
    time_t t;
    t = time(NULL);
    answer.tv_sec = ( unsigned long ) t;
    answer.tv_usec = 0;
}

time_t getTime()
{
    time_t t;
    t = time(NULL);
    return t;
}

#endif

void convertTS2HMS(int &nhrs, int &nmin, int &nsec, long int tsec)
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

void convertTS2HMS(int &nhrs, int &nmin, int &nsec, double tsec)
{
    long int _nsec = ( long int ) tsec;
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
} // end namespace oofem
