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

#ifndef timer_h
#define timer_h

#include "clock.h"

namespace oofem {
/**
 * Class implementing single timer, providing wall clock and user time capabilities.
 */
class Timer
{
    /// Wall clock time structures.
    oofem_timeval start_wtime, end_wtime;
    /// User time struct.
    oofem_timeval start_utime, end_utime;
    /// Accumulated wtime and utime (in seconds) from start.
    oofem_timeval elapsedWTime, elapsedUTime;
    /// Flag indicating whether timer is running.
    bool running;

public:
    Timer() { initTimer(); }
    void startTimer() {
        this->initTimer();
        oofem :: getTime(start_wtime);
        oofem :: getUtime(start_utime);
        running = true;
    }
    void stopTimer() {
        this->pauseTimer();
        running = false;
    }
    void pauseTimer() {
        oofem :: getTime(end_wtime);
        oofem :: getUtime(end_utime);
        running = false;
        this->updateElapsedTime();
    }
    void resumeTimer() {
        oofem :: getTime(start_wtime);
        oofem :: getUtime(start_utime);
        running = true;
    }
    void initTimer() {
        elapsedWTime.tv_sec = elapsedWTime.tv_usec = elapsedUTime.tv_sec = elapsedUTime.tv_usec = 0;
        running = false;
    }
    bool isRunning() { return running; }

    /// Returns total user time elapsed in seconds
    double getUtime() {
        this->updateElapsedTime();
        return ( double ) elapsedUTime.tv_sec + ( double ) elapsedUTime.tv_usec / OOFEM_USEC_LIM;
    }
    /// Returns total elapsed wall clock time in seconds
    double getWtime() {
        updateElapsedTime();
        return ( double ) elapsedWTime.tv_sec + ( double ) elapsedWTime.tv_usec / OOFEM_USEC_LIM;
    }

    // Converts total seconds into hours, minutes, and seconds.
    void convert2HMS(int &nhrs, int &nmin, int &nsec, long int tsec) const { oofem :: convertTS2HMS(nhrs, nmin, nsec, tsec); }
    // Converts total seconds into hours, minutes, and seconds.
    void convert2HMS(int &nhrs, int &nmin, int &nsec, double tsec) const { oofem :: convertTS2HMS(nhrs, nmin, nsec, tsec); }

    /// Prints receiver state into a string.
    void toString(char *buff)  { sprintf( buff, "ut: %f.3s, wt: %f.3s", getUtime(), getWtime() ); }

    void updateElapsedTime()  {
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
};

/**
 * Timer class, assumed to be an attribute of engineering model, serving stop-watch facility for engineering model.
 * It can handle several timers independently, each corresponding to different solution stage, etc.
 * Each timer is capable to track elapsed wall clock time as well as user time.
 */
class EngngModelTimer
{
public:
    /**
     * Enumeration to distinguish different type of timers.
     *
     * EMTT_NetComputationalStepTimer timer (and particularly its wall clock time) should measure only computation itself,
     * no communication, therefore it should be measure of workload (in terms of wall clock time) on particular processors.
     * It also typically not include time needed to solve the system of equations, since this has to be done in parallel,
     * so solution takes the same time on all processors and include unwanted synchronization.
     */
    enum EngngModelTimerType {
        EMTT_AnalysisTimer,
        EMTT_SolutionStepTimer,
        EMTT_NetComputationalStepTimer,
        EMTT_LoadBalancingTimer,
        EMTT_DataTransferTimer,
        EMTT_LastTimer
    };

protected:
    /// Array of Timer classes.
    Timer timers [ EMTT_LastTimer ];

public:
    EngngModelTimer() { }
    virtual ~EngngModelTimer() { }

    /**@name Profiling routines. */
    //@{
    void startTimer(EngngModelTimerType t) { timers [ t ].startTimer(); }
    void stopTimer(EngngModelTimerType t) { timers [ t ].stopTimer(); }
    void pauseTimer(EngngModelTimerType t) { timers [ t ].pauseTimer(); }
    void resumeTimer(EngngModelTimerType t) { timers [ t ].resumeTimer(); }
    void initTimer(EngngModelTimerType t) { timers [ t ].initTimer(); }
    //@}

    /**@name Reporting routines. */
    //@{
    /// Returns total user time elapsed.
    double getUtime(EngngModelTimerType t)  { return timers [ t ].getUtime(); }
    /// Returns elapsed wall clock time.
    double getWtime(EngngModelTimerType t)   { return timers [ t ].getWtime(); }
    /// Returns pointer to timer determined by EngngModelTimerType.
    const Timer *getTimer(EngngModelTimerType t)  { return timers + t; }
    /// Converts total seconds into hours, mins, and seconds.
    void convert2HMS(int &nhrs, int &nmin, int &nsec, long int tsec) const { oofem :: convertTS2HMS(nhrs, nmin, nsec, tsec); }
    void convert2HMS(int &nhrs, int &nmin, int &nsec, double tsec) const { oofem :: convertTS2HMS(nhrs, nmin, nsec, tsec); }
    /// Printing & formatting.
    void toString(EngngModelTimerType t, char *buff) { return timers [ t ].toString(buff); }
    //@}
};
} // end namespace oofem
#endif // timer_h
